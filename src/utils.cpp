# include "utils.hpp"
# include <iostream>

namespace eva { namespace utils {

boost::program_options::variables_map
handle_cmd_line_options(int argc, char * argv[])
{
	namespace po = boost::program_options;
    
	// Setup available options
	po::options_description opts("Program options");
	opts.add_options()
		("help,h", "Produce this help message")

        ("config,c", po::value<std::string>()->default_value("config.conf"),
         "File containing the problem configuration")

        ("structure,s", po::value<std::string>()->default_value("structure.dot"),
         "File containing the problem structure")
		;

	// Store args
	po::variables_map vmap;
	po::store(po::parse_command_line(argc, argv, opts), vmap);

	// Check for help
	if (vmap.count("help"))
	{
		std::cout << opts << std::endl;
		throw help_exception();
	}

	// Otherwise notify errors
	po::notify(vmap);
	return vmap;
}

}}// end namespaces
