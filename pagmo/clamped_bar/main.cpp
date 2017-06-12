// eva
# include "frame.hpp"
# include "thermo.hpp"
# include "utils.hpp"
# include "io.hpp"
# include "pagmo_utils.hpp"
# include <iostream>

# include <boost/math/constants/constants.hpp>
# include <boost/container/flat_map.hpp>
# include "boost_flat_container_serialization.hpp"
# include <chrono>


/// Builds the map containing containing the problem configuration
/// by reaading the specified file, throwing the proper exceptions 
/// when  a field is missing or has an invalid value.
/// @param[in] filename name of the file containing the configuration
/// @returns a map containing the options red
boost::program_options::variables_map
read_config_file(const std::string & filename);


int main(int argc, char * argv [])
{
    
    // Namespaces
    using namespace utils;
    using eva::real;
    
    // Aliases
    using struct_kind_t = eva::frame_kind<2>;
    using thermo_kind_t = eva::thermo_kind<struct_kind_t>;
    using structure_t =  eva::thermo_structure<struct_kind_t>;;
    using joint_t     = typename eva::joint_of  <structure_t>::type;
    using element_t   = typename eva::element_of<structure_t>::type;
    // using kind_t      = typename eva::kind_of<structure_t>::type;
    
    // Read cmd line options
    auto cmd_line_opts  = eva::utils::handle_cmd_line_options(argc, argv);
    
    // Read config file
    auto config_file = cmd_line_opts["config"].as<std::string>();
    auto config = read_config_file(config_file);
    
    auto l = config["length"].as<real>();//100.e-3;
    auto h = config["height"].as<real>();//25.0e-3;
    auto m = config["rows"].as<int>();//30;
    auto n = config["cols"].as<int>();//50;

    auto E = config["E"].as<real>();
    auto k = config["k"].as<real>();
    
    // Create grid representing the beam
    auto grid = make_grid<joint_t>(l, h, m, n, true); // <--Jitter

    // Fix first and last columns joints (set bcs)
    for (auto i = 0u; i < m; ++i)
    {
        //                      x , y , theta
        grid[i][ 0 ].bcs << 0., 0., 0.;
        grid[i][n-1].bcs << 0., 0., 0.;
    }
    
    // Apply mechanical load on joint [0][3] (lowest row)
    grid[0][3].load << 0, 16 * 1.e4;

    // Apply mechanical load on joint [m-1][n-3] (highest row)
    // grid[m-1][n-3].load << 0, -16 * 1.e4;

    // Apply mechanical load on joint [m-1][n-3] (highest row)
    // grid[0][n/2].load << 0, -160 * 1.e4;

    // Apply load on the entire lowest row
    // for (auto i = 0u; i < n; ++i) grid[0][i].load << 0, -16 * 1.e4;

    // Apply homogenous Dirchlet BC on the outer elements
    for (auto i = 0u; i < m; ++i)
    {
        grid[i][ 0 ].T_bc = 0.;
        grid[i][n-1].T_bc = 0.;
    }
    
    // Apply thermal load on joint [2][2]
    grid[m-1][2].flux_bc = 0.5 * 1.3;

    
    // Build structure and tag nodes indicating on 
    // which boundary they lie:
    // -] 0 -> internal node
    // -] 1 -> lies on a boundary row
    // -] 2 -> lies on a boundary column
    // -] 3 -> 1+2 hold
    auto structure    = structure_t();
    auto vertex_idx   = 0u;
    auto boundary_map = boost::container::flat_map<size_t, int>();
    
    for (auto i = 0u; i < m; ++i)
        for (auto j = 0u; j < n; ++j)
        {
            add_vertex(grid[i][j], structure);

            // Init value as 0 by default
            auto & value = boundary_map[vertex_idx];
            value = 0;
            
            if (i == 0 || i == m-1) value += 1;
            if (j == 0 || j == n-1) value += 2;

            // Also keep track of nodes with loads applied ?
            ++vertex_idx;
        }
    // for (auto vidx : boost::make_iterator_range(vertices(structure)))
    //     std::cout << vidx << " -> " << boundary_map[vidx] << std::endl;
    

    // Select feasible connections
    auto topology    = make_topology(grid, neumann_rule());
    auto cardinality = topology.nonZeros();
    std::cout << "Num elements = " << cardinality << std::endl;

    // Write to file structure, boundary map & topology
    eva::write_graphviz  (   structure,   "base_frame.dot");
    eva::utils::serialize(boundary_map, "boundary_map.bin");
    eva::utils::serialize(    topology,     "topology.bin");

    std::cout << "Delta X = " << l / (n-1) << "\n";
    std::cout << "Delta Y = " << h / (m-1) << "\n";
    
    
    // Add elements to the structure
    for (int kk = 0; kk < topology.outerSize(); ++kk)
        for (topology_t::InnerIterator it(topology, kk); it; ++it)
            if (it.value())
            {
                auto r = 4.e-3;
                if (boundary_map[it.row()] > 0 && boundary_map[it.col()] > 0)
                    r += 2.e-3;
                
                constexpr double pi = boost::math::constants::pi<double>();

                auto element = element_t();
                element.E = E;
                element.A = pi * r*r;
                element.I = pi/2 * r*r*r*r;
                element.k = k;
                
                add_edge(it.row(), it.col(), element, structure);
            }
    
    // Take start time
    auto beg_solve = std::chrono::high_resolution_clock::now();
    
    // Solve structural problem & compute compliance
    using struct_solver_t = eva::sparse_solver_params<>;
    auto struct_results  = solve(structure, struct_kind_t(), struct_solver_t());
    auto compliance  = eva::compute_compliance(structure, struct_results);
    std::cout << "Compliance = " << compliance << std::endl;

    // Compute volume
    auto volume = compute_mass(structure, 1.);
    std::cout << "Volume = " << volume << std::endl;
    
    // Solve thermal problem & compute max temperature
    using thermo_solver_t = eva::dense_solver_params<thermo_kind_t::default_dense_solver_t>;
    auto thermo_results = solve(structure, thermo_kind_t(), thermo_solver_t());
    auto max_t = eva::get_max_temperature(structure, thermo_results);
    std::cout << "Max temperature = " << max_t << std::endl;

    // Take time & print it
    auto end_solve = std::chrono::high_resolution_clock::now();
    auto solving_time = std::chrono::duration<double>(end_solve - beg_solve);    
    std::cout << "Solving time = " << solving_time.count() << std::flush << std::endl;;

    // // Compute volume
    // auto volume = eva::compute_mass(structure, 1.);
    // std::cout << "Volume = " << volume << std::endl;    

    // Display
    // eva::display(structure);

    // // // Save results
    // // write_vtu(structure, results, "test.vtu");

    // Assemble & display displaced structure
    auto dsp_structure = structure;
    for (auto v : boost::make_iterator_range(vertices(dsp_structure)))
    {
        dsp_structure[v].coords += struct_results[v].displacement;
    }
    eva::display(dsp_structure);

    
    return 0;
}


boost::program_options::variables_map
read_config_file(const std::string & filename)
{
	namespace po = boost::program_options;
    using eva::real;

	// Setup available options
	po::options_description opts("Problem Config");
	opts.add_options()
        
        ("length,l", po::value<real>()->required(), "Bar length [m]")

        ("height,h", po::value<real>()->required(), "Bar height [m]")

        ("rows,m", po::value<int>()->required(), "Number of rows in the discretization grid")

        ("cols,n", po::value<int>()->required(), "Number of cols in the discretization grid")
        
		("E", po::value<real>()->default_value(69.e9), "Elasticity module")

        ("k", po::value<real>()->default_value(237), "Thermal conductivity")
        
		;
    
	// Read config
    std::ifstream ifile(filename);    
	po::variables_map vmap;
	po::store(po::parse_config_file(ifile, opts, true), vmap);

	// Otherwise notify errors
	po::notify(vmap);
	return vmap;
}



