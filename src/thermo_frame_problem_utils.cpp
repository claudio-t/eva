# include "thermo_frame_problem_utils.hpp"

namespace utils { namespace thermo_frame_problem {


size_t
compute_dof_nr(const size_t m , const size_t n, const std::string & rule)
{
    // Compute #dof associated to vertices
    // (common for both rules)
    auto n_v = n*m;
    auto v_dof = 2*n_v - 2*(n+m-2) - 4;

    
    // Compute #dof associated to edges
    auto n_e = 0u;
    if (rule == "neumann" || rule == "Neumann")
    {
        auto n_e = n*(m-1) + m*(n-1);
    }
    else if (rule == "moor" || rule == "Moor")
    {
        auto n_e = n*(m-1) + m*(n-1) + 2*(m-1)*(n-1);
    }
    else
    {
        throw std::invalid_argument("rule must be either 'Neumann' or 'Moor'");
    }
    auto e_dof = 1 * n_e;

    return v_dof + e_dof;
}

size_t
compute_dof_nr(
    const topology_t & topology,
    const boundary_map_t & boundary_map)
{    
    // Compute problem dimension (nr of dof)
    auto nr_links = topology.nonZeros(); // nr of allowed elements
    auto nr_inner = 0u; // nr of inner nodes
    auto nr_outer = 0u; // nr of outer nodes
    for (const auto & el : boundary_map) el.second == 0 ? ++nr_inner : ++nr_outer;
    
    auto prb_dim = nr_links + 2*nr_inner + nr_outer - 4;

    return prb_dim;
}


joint_grid_t
make_grid(eva::real length, eva::real height, size_t m, size_t n, bool jitter)
{
    assert(joint_t().coords.size() == 2);
        
    if (m < 2) throw std::invalid_argument("Parameter m must be >= 2");
    if (n < 2) throw std::invalid_argument("Parameter n must be >= 2");
    
    // Prealloc joint rows
    auto joints = joint_grid_t(m);
    joints.reserve(m);

    // Compute grid x and y spacings
    auto delta_x = length / (n - 1);
    auto delta_y = height / (m - 1);

    // // Add some randomess
    std::mt19937 rng;
    rng.seed(std::random_device()());
    auto dx = std::uniform_real_distribution<float>(-delta_x / 2, delta_x / 2);
    auto dy = std::uniform_real_distribution<float>(-delta_y / 2, delta_y / 2);

    for (auto i = 0u; i < m; ++i)
    {        
        // Reserve space
        auto n_nodes = n;
        joints[i].resize(n_nodes);
        
        for (auto j = 0u; j < n_nodes; ++j)
        {
            // Init coordinates
            auto x = j*delta_x;
            auto y = i*delta_y;

            // Add random jitter
            if (jitter)
            {                
                auto is_boundary_row = i == 0 || i == m-1;
                auto is_boundary_col = j == 0 || j == n-1;

                // Add random jitter
                if (!is_boundary_row) y += 0.8*dy(rng);
                if (!is_boundary_col) x += 0.8*dx(rng);
            }
            
            
            // Set joint coords
            joints[i][j].coords << x, y;
        }
    }

    return joints;
}


void add_bcs(joint_grid_t & grid, const po::variables_map & config)
{    
    const auto m = grid.size();
    const auto n = grid.front().size();
    
    // Fix first and last columns joints (set bcs)
    for (auto i = 0u; i < m; ++i)
    {
        //                      x , y , theta
        grid[i][ 0 ].bcs << 0., 0., 0.;
        grid[i][n-1].bcs << 0., 0., 0.;
    }
    
    // Apply homogenous Dirchlet BC on the outer elements
    for (auto i = 0u; i < m; ++i)
    {
        grid[i][ 0 ].T_bc = 0.;
        grid[i][n-1].T_bc = 0.;
    }
}



void add_loads(joint_grid_t & grid, const po::variables_map & config)
{ 
    const auto m = grid.size();
    const auto n = grid.front().size();
    
    // Apply mechanical load on joint [lowest_row][3] (lowest row)
    grid[0][3].load << 0, 16 * 1.e4;

    // Apply thermal load on joint [highest_row][2]
    grid[m-1][2].flux_bc = 0.5 * 1.3;
}



structure_t
make_frame(
    const joint_grid_t & grid,
    const topology_t & topology,
    const po::variables_map & config)
{
    auto E = config["E"].as<real>(); assert(E >= 0.);
    auto k = config["k"].as<real>(); assert(k >= 0.);
    
    auto structure = structure_t();

    const auto m = grid.size();
    const auto n = grid.front().size();
    
    for (auto i = 0u; i < m; ++i)
        for (auto j = 0u; j < n; ++j)
            add_vertex(grid[i][j], structure);

    // Add elements to the structure
    for (int kk = 0; kk < topology.outerSize(); ++kk)
        for (topology_t::InnerIterator it(topology, kk); it; ++it)
            if (it.value())
            {
                auto r = 2.e-3;  
                constexpr double pi = boost::math::constants::pi<double>();

                auto element = element_t();
                element.E = E;
                element.A = pi * r*r;
                element.I = pi/2 * r*r*r*r;
                element.k = k;
                
                add_edge(it.row(), it.col(), element, structure);
            }
    
    return structure;
}


boundary_map_t
make_boundary_map(const joint_grid_t & grid)
{
   // Build structure and tag nodes indicating on 
    // which boundary they lie:
    // -] 0 -> internal node
    // -] 1 -> lies on a boundary row
    // -] 2 -> lies on a boundary column
    // -] 3 -> 1+2 hold
    auto vertex_idx   = 0u;
    auto boundary_map = boundary_map_t();
    
    const auto m = grid.size();
    const auto n = grid.front().size();  
    
    for (auto i = 0u; i < m; ++i)
        for (auto j = 0u; j < n; ++j)
        {
            // Init value as 0 by default
            auto & value = boundary_map[vertex_idx];
            value = 0;
            
            if (i == 0 || i == m-1) value += 1;
            if (j == 0 || j == n-1) value += 2;

            // Also keep track of nodes with loads applied ?
            ++vertex_idx;
        }
    return boundary_map;
}


po::variables_map
read_config_file(const std::string & filename)
{
    using eva::real;

	// Setup available options
	po::options_description structure_opts("Structure generation");
	structure_opts.add_options()        

        ("length,l", po::value<real>()->required(), "Bar length [m]")

        ("height,h", po::value<real>()->required(), "Bar height [m]")

        ("rows,m", po::value<int>()->required(), "Number of rows in the discretization grid")

        ("cols,n", po::value<int>()->required(), "Number of cols in the discretization grid")

        ;
    
	po::options_description problem_opts("Problem configuration");
	problem_opts.add_options()
        
        ("x-jitter", po::value<real>()->required(), "Maximum jitter along the x-axis")

        ("y-jitter", po::value<real>()->required(), "Maximum jitter along the y-axis")

        ("min-radius", po::value<real>()->required(), "Minimum radius allowed")

        ("max-radius", po::value<real>()->required(), "Maximum radius allowed")
        
        ("min-volume", po::value<real>()->default_value(-1.), "Minimum volume allowed")

        ("max-volume", po::value<real>()->default_value(-1.), "Maximum volume allowed")

        ("tolerance", po::value<real>()->default_value(0.), "Volume constraint tolerance")
        
		("E", po::value<real>()->default_value(69.e9), "Elasticity module")

        ("k", po::value<real>()->default_value(237), "Thermal conductivity")
        ;

    po::options_description optimization_opts("Optimization options");
    optimization_opts.add_options()

        ("evolutions", po::value<size_t>()->default_value(5), "Number of archipelago evolutions")

        ("islands", po::value<size_t>()->default_value(6), "Number of archipelago islands")
        
        ("individuals", po::value<size_t>()->default_value(1000), "Number of individuals per island")

        ("generations", po::value<size_t>()->default_value(6), "Number of NSGA2 generations")
        
		;
    
    po::options_description all_opts("Optimization Problem Configuration");
    all_opts.add(structure_opts).add(problem_opts).add(optimization_opts);
    
	// Read config
    std::ifstream ifile(filename);
    if (ifile.bad()) throw std::ios_base::failure("Unable to open config file");
    
	po::variables_map vmap;
	po::store(po::parse_config_file(ifile, all_opts, true), vmap);

	// Otherwise notify errors
    try
    {
        po::notify(vmap);
    }
    catch (std::exception & ex)
    {
        std::cerr
            << "Error when parsing config file.\n Options are:\n"
            << all_opts << std::endl;
        throw ex;
    }
	return vmap;
}



problem_t
make_problem(
    const structure_t & structure,
    const boundary_map_t & boundary_map,
    const topology_t & topology,
    const po::variables_map & config)
{    
   // Compute problem dimension (nr of dof)
    auto nr_links = topology.nonZeros(); // nr of allowed elements
    auto nr_inner = 0u; // nr of inner nodes
    auto nr_outer = 0u; // nr of outer nodes
    for (const auto & el : boundary_map) el.second == 0 ? ++nr_inner : ++nr_outer;
    
    auto prb_dim = nr_links + 2*nr_inner + nr_outer - 4;    

    
    // Build problem bounds
    auto x_jitter   = config[  "x-jitter"].as<real>();
    auto y_jitter   = config[  "y-jitter"].as<real>();
    auto min_radius = config["min-radius"].as<real>();
    auto max_radius = config["max-radius"].as<real>();
    // auto tolerance  = config[ "tolerance"].as<real>();
    
    // auto min_volume = min_radius * nr_outer;
    // if (config.count("min-volume"))
    //     min_volume = config["min-volume"].as<real>();

    // auto max_volume = max_radius * nr_links;
    // if (config.count("max-volume"))
    //     max_volume = config["max-volume"].as<real>();
    
    assert(  x_jitter >= 0.);
    assert(  y_jitter >= 0.);
    assert(min_radius >= 0.);
    assert(max_radius >= 0.);
    // assert(min_volume >= 0.);
    // assert(max_volume >= 0.);
    // assert( tolerance >= 0.);
    
    auto bounds = problem_t::bounds {
        x_jitter, y_jitter, min_radius, max_radius//, min_volume, max_volume
    };

    // // Build problem physical properties
    // auto E = config["E"].as<real>(); assert(E >= 0.);
    // auto k = config["k"].as<real>(); assert(k >= 0.);
    // auto phys_props = problem_t::physical_properties {E, k};
    
    // Build problem
    return  problem_t(
        structure, boundary_map, prb_dim, bounds//, phys_props
        );
}


int get_champion_island_idx(const pagmo::archipelago & archi)
{
    int idx = 0;
    auto min_glob = archi.get_island(0)->get_population().champion().f[0];
    auto problem  = archi.get_island(0)->get_population().problem().clone();
    
	for (size_t i = 1;i<archi.get_size(); ++i)
    {
		auto champ_loc = archi.get_island(i)->get_population().champion();
		if (champ_loc.f[0] < min_glob && problem->feasibility_x(champ_loc.x))
        {
			min_glob = champ_loc.f[0];
			idx = i;
		}
	}
	return idx;
}


pagmo::population::champion_type
get_champion(const pagmo::archipelago & archi)
{
    auto champ_isl_idx = get_champion_island_idx(archi);
    return archi.get_island(champ_isl_idx)->get_population().champion();
}



void
display_displaced(const structure_t & frame)
{
    // Solve structural problem
    auto struct_results = solve(frame, struct_kind_t());

    auto dsp_frame = frame;
    for (auto v : boost::make_iterator_range(vertices(dsp_frame)))
    {
        dsp_frame[v].coords += struct_results[v].displacement;
    }
    eva::display(dsp_frame);
}


void
export_to_vtu(const structure_t & frame, const std::string & filename)
{
    // Solve structural problem
    auto struct_results = solve(frame, struct_kind_t());
    
    // Solve thermal problem
    using solver_t = eva::dense_solver_params<thermo_kind_t::default_dense_solver_t>;
    auto thermo_results = solve(frame, thermo_kind_t(), solver_t());    
    
    // Save results to vtu file
    auto vtk_grid = eva::to_vtk_unstructured_grid(frame);

    eva::vtk_add_joint_results(struct_results, vtk_grid);
    eva::vtk_add_joint_results(thermo_results, vtk_grid);
    
    eva::write_vtu(vtk_grid, filename);
}



}}//end namespaces
