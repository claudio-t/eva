# include "clamped_bar_problem.hpp"
# include "utils.hpp"
# include "pagmo_utils.hpp"
# include "io.hpp"
# include "boost_flat_container_serialization.hpp"
# include <cassert>


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
    using problem_t =  pagmo::problem::clamped_bar;
    using struct_kind_t = problem_t::struct_kind_t;
    using thermo_kind_t = problem_t::thermo_kind_t;
    using structure_t   = problem_t::thermo_kind_t;
    
    using joint_t     = problem_t::joint_t;
    using element_t   = problem_t::element_t;
    
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
    auto boundary_map = problem_t::boundary_map_t();
    eva::utils::deserialize(boundary_map, boundary_map_file);
        
    // Read boundary map file
    auto topology_file = config["topology"].as<std::string>();
    auto topology = problem_t::topology_t();
    eva::utils::deserialize(topology, topology_file);

     // Compute problem dimension (nr of dof)
    auto nr_links = topology.nonZeros(); // nr of allowed elements
    auto nr_inner = 0u; // nr of inner nodes
    auto nr_outer = 0u; // nr of outer nodes
    for (const auto & el : boundary_map) el.second == 0 ? ++nr_inner : ++nr_outer;
    
    auto prb_dim = nr_links + 2*nr_inner + nr_outer - 4;   

    
    // Build problem bounds
    auto x_jitter = config["x-jitter"].as<real>();
    auto y_jitter = config["y-jitter"].as<real>();
    auto min_radius = config["min-radius"].as<real>();
    auto max_radius = config["max-radius"].as<real>();
    
    auto min_volume = real(0.);
    if (config.count("min-volume"))
        min_volume = config["min-volume"].as<real>();
    else
        min_volume = min_radius * nr_outer;

    auto max_volume = real(0.);
    if (config.count("max-volume"))
        max_volume = config["max-volume"].as<real>();
    else
        max_volume = max_radius * nr_links;
    
    assert(x_jitter >= 0.);
    assert(y_jitter >= 0.);
    assert(min_radius >= 0.);
    assert(max_radius >= 0.);
    assert(min_volume >= 0.);
    assert(max_volume >= 0.);
    
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
    auto prb = problem_t(
        structure, topology, boundary_map, prb_dim, bounds, phys_props
        );

    eva::utils::serialize(prb, "prb.bin");
    
    auto prb2 = prb;
    eva::utils::deserialize(prb2, "prb.bin");

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

        ("structure", po::value<std::string>()->required(),
         "Name of the file which contains the base structure")
        
        ("boundary-map", po::value<std::string>()->required(),
         "Name of the file which contains the serialized boundary map")

        ("topology", po::value<std::string>()->required(),
         "Name of the file which contains the serialized topology")

        ("x-jitter", po::value<real>()->required(), "Maximum jitter along the x-axis")

        ("y-jitter", po::value<real>()->required(), "Maximum jitter along the y-axis")

        ("min-radius", po::value<real>()->required(), "Minimum radius allowed")

        ("max-radius", po::value<real>()->required(), "Maximum radius allowed")
        
        ("min-volume", po::value<real>()->default_value(-1.), "Minimum volume allowed")

        ("max-volume", po::value<real>()->default_value(-1.), "Maximum volume allowed")
        
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
