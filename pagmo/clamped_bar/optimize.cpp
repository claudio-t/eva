/////////////////////////////////////////////////////
//eva
# include "clamped_bar_problem.hpp"
# include "utils.hpp"
# include "pagmo_utils.hpp"
# include "io.hpp"
// std
# include <cassert>
# include <chrono>
// pagmo
# include <pagmo/src/population.h>
# include <pagmo/src/algorithms.h>
# include <pagmo/src/archipelago.h>
# include <pagmo/src/topologies.h>
# include <pagmo/src/algorithm/gsl_fr.h>
# include <pagmo/src/algorithm/nlopt_aug_lag.h>
/////////////////////////////////////////////////////

// Namespaces
namespace po = boost::program_options;
using eva::real;
    
// Aliases
using problem_t =  pagmo::problem::clamped_bar;

using struct_kind_t = problem_t::struct_kind_t;
using thermo_kind_t = problem_t::thermo_kind_t;
using structure_t   = problem_t::structure_t;
using joint_t       = problem_t::joint_t;
using element_t     = problem_t::element_t;

using boundary_map_t        = problem_t::boundary_map_t;
using topology_t            = problem_t::topology_t;
using bounds_t              = problem_t::bounds;
using physical_properties_t = problem_t::physical_properties;


/// Builds the map containing containing the problem configuration
/// by reaading the specified file, throwing the proper exceptions 
/// when  a field is missing or has an invalid value.
///
/// @param[in] Filename name of the file containing the configuration
///
/// @returns A map containing the options red
po::variables_map
read_config_file(const std::string & filename);

/// Builds the problem object.
///
/// @param[in] structure The Base frame containing the structure nodes
/// @param[in] boundary_map The map which indicates boundary nodes
/// @param[in] topology The network representing the allowed set of
///                     elements
/// @param[in] config The map representing the config file content
///
/// @returns An instance of the problem
problem_t
make_problem(
    const structure_t & structure,
    const boundary_map_t & boundary_map,
    const topology_t & topology,
    const po::variables_map & config);


int main(int argc, char * argv [])
{
    try
    {
        
        // Read cmd line options
        auto cmd_line_opts  = eva::utils::handle_cmd_line_options(argc, argv);
    
        // Read config file
        auto config_file = cmd_line_opts["config"].as<std::string>();
        auto config = read_config_file(config_file);

        // Read frame file
        auto structure_file = config["structure"].as<std::string>();
        auto structure = eva::read_from_graphviz<structure_t>("base_frame.dot");    

        // Read boundary map file
        auto boundary_map_file = config["boundary-map"].as<std::string>();
        auto boundary_map = boundary_map_t();
        eva::utils::deserialize(boundary_map, boundary_map_file);
        
        // Read boundary map file
        auto topology_file = config["topology"].as<std::string>();
        auto topology = topology_t();
        eva::utils::deserialize(topology, topology_file);
        
        // Build problem
        auto problem = make_problem(structure, boundary_map, topology, config);

        // Setup NSGA2 algorithm
        auto gens  = config["generations"].as<size_t>(); assert(gens > 0);
        auto cr    = 0.9;
        auto eta_c = 10;
        auto m     = 0.05;
        auto eta_m = 50;

        auto algorithm = pagmo::algorithm::nsga2(gens, cr,eta_c, m, eta_m).clone();
        algorithm->set_screen_output(true);

                                             
        // Create the archipelago
        auto n_isls = config[  "islands"  ].as<size_t>(); assert(n_isls > 0);
        auto n_inds = config["individuals"].as<size_t>(); assert(n_inds > 0);
        
        auto archi_topo = pagmo::topology::ring();
        auto archi = pagmo::archipelago(*algorithm, problem, n_isls, n_inds, archi_topo);

        // Evolve archipelago
        auto n_evo = config["evolutions"].as<size_t>(); assert(n_evo > 0);
        for (auto i = 0u; i < n_evo; ++i)
        {
            std::cout << "Starting generation " << i << "\n";
            
            // GLOBALLY optimize
            auto beg_evo = std::chrono::high_resolution_clock::now();
            archi.evolve(1);
            archi.join();
            auto end_evo = std::chrono::high_resolution_clock::now();

            // Print elapsed time
            auto solving_time = std::chrono::duration<double>(end_evo - beg_evo);    
            std::cout << "Elapsed time = " << solving_time.count() << std::endl;
            
            // Retrieve champion //
            auto champ_isl_idx = utils::get_champion_island_idx(archi);
            auto champion = archi.get_island(champ_isl_idx)->get_population().champion();

            // Print fitness
            std::cout << "Fitness = " << champion.f << std::endl;
        }
        
        // eva::utils::serialize(problem, "prb.bin");    
        // auto prb2 = problem;
        // eva::utils::deserialize(prb2, "prb.bin");
    }
    catch (std::exception & ex)
    {
        std::cerr << "An error occurred: " << ex.what() << std::endl;
        return 1;
    }
    catch (...)
    {
        std::cerr << "Unknown error occurred" << std::endl;
        return 1;
    }
    

    return 0;
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
    auto tolerance  = config[ "tolerance"].as<real>();
    
    auto min_volume = min_radius * nr_outer;
    if (config.count("min-volume"))
        min_volume = config["min-volume"].as<real>();

    auto max_volume = max_radius * nr_links;
    if (config.count("max-volume"))
        max_volume = config["max-volume"].as<real>();
    
    assert(  x_jitter >= 0.);
    assert(  y_jitter >= 0.);
    assert(min_radius >= 0.);
    assert(max_radius >= 0.);
    assert(min_volume >= 0.);
    assert(max_volume >= 0.);
    assert( tolerance >= 0.);
    
    auto bounds = problem_t::bounds {
        x_jitter, y_jitter, min_radius, max_radius, min_volume, max_volume
    };

    // Build problem physical properties
    auto E = config["E"].as<real>();
    auto k = config["k"].as<real>();

    assert(E >= 0.);
    assert(k >= 0.);
    
    auto phys_props = problem_t::physical_properties {E, k};
    
    // Build problem
    return  problem_t(
        structure, topology, boundary_map, prb_dim, bounds, phys_props
        );
}


po::variables_map
read_config_file(const std::string & filename)
{
	namespace po = boost::program_options;
    using eva::real;

	// Setup available options
	po::options_description structure_opts("Structure configuration");
	structure_opts.add_options()        

        ("structure", po::value<std::string>()->required(),
         "Name of the file which contains the base structure")
        
        ("boundary-map", po::value<std::string>()->required(),
         "Name of the file which contains the serialized boundary map")

        ("topology", po::value<std::string>()->required(),
         "Name of the file which contains the serialized topology")
        ;
    
	po::options_description problem_opts("Problem options");
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
