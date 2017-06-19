//eva
# include "thermo_frame_problem.hpp"
# include "thermo_frame_problem_utils.hpp"
# include "utils.hpp"
# include "io.hpp"
# include <iomanip>

# include <pagmo/src/population.h>

using namespace utils::thermo_frame_problem;

boost::program_options::variables_map
handle_cmd_line_options(int argc, char * argv[]);

int main(int argc, char * argv[])
{
    try
    {
        // Read cmd line options
        auto cmd_line_opts  = handle_cmd_line_options(argc, argv);

        // Retrieve problem components
        auto problem_file = cmd_line_opts["problem"].as<std::string>();
        auto problem = problem_t();
        eva::utils::deserialize(problem, problem_file);

        // Retrieve individuals
        auto individuals_file = cmd_line_opts["population"].as<std::string>();
        auto individuals = std::vector<pagmo::population::individual_type>();
        eva::utils::deserialize(individuals, individuals_file);

        if (individuals.size() == 0)
        {
            std::cout << "Warning: empty population (check input file)" << std::endl;
            return 0;
        }

        // Get individual indices to be displayed and/or printed
        auto display_idxs = cmd_line_opts["display"]     .as< std::vector<size_t> >();
        auto export_idxs  = cmd_line_opts["export-vtu"]  .as< std::vector<size_t> >();
        auto print_idxs   = cmd_line_opts["show-fitness"].as< std::vector<size_t> >();
        
        if (cmd_line_opts.count("show-all"))
        {
            print_idxs = std::vector<std::size_t>(individuals.size());
            std::iota(begin(print_idxs), end(print_idxs), 0u);
        }
        
        // Display individuals
        for (auto idx : display_idxs)
        {
            auto frame = problem.encode_genes((individuals[idx].cur_x));
            eva::display(frame);
            display_displaced(frame);
        }

        // Export individuals
        for (auto idx : export_idxs)
        {
            auto vtu_name = "frame_" + std::to_string(idx) + ".vtu";
            std::cout
                << "Exporting individual " << idx
                << " to " << vtu_name << std::endl;
            
            auto frame = problem.encode_genes((individuals[idx].cur_x));
            export_to_vtu(frame, vtu_name);
        }
        
        
        // Print fitnesses
        auto print_pop = pagmo::population(problem);
        for (auto & idx : print_idxs) print_pop.push_back(individuals[idx].cur_x);

        auto ii = 0u;
        auto num_digits = (print_pop.size() > 0) ? int(log10(double(print_pop.size() )))+1 : 1;
        for (auto & ind : print_pop)
            std::cout
                << "Idx: " << std::setw(num_digits) <<print_idxs[ii++] << "  "
                << "Fitness:" << ind.cur_f << std::endl;

        auto show_champion = cmd_line_opts.count("show-champion");
        auto export_csv    = cmd_line_opts.count("export-csv");
        // auto export_fronts = cmd_line_opts.count("export-pareto-fronts");
        
        if (show_champion || export_csv)
        {
          auto pop = pagmo::population(problem);
          for (auto & ind : individuals) pop.push_back(ind.cur_x);

          if (show_champion)
              std::cout << "Champion:\n" << pop.champion() << std::endl;

          // if (export_csv)
          // {
          //     auto csv_name = cmd_line_opts["export-csv"].as<std::string>();
          //     std::ofstream ofs(csv_name);
          //     auto idx = 0u;
          //     ofs << "index,compliance,volume,max_temperature\n";
          //     for (const auto & ind : pop)\
          //     {
          //         ofs
          //             << idx++ << ", "
          //             << ind.cur_f[0] << ","
          //             << ind.cur_f[1] << ","
          //             << ind.cur_f[2] << "\n";
          //     }
          // }

          if (export_csv)
          {
              auto csv_name = cmd_line_opts["export-csv"].as<std::string>();
              std::ofstream ofs(csv_name);
              auto idx = 0u;
              ofs << "index,front,compliance,volume,max_temperature\n";

              // Compute fronts
              auto fronts = pop.compute_pareto_fronts();

              // Write them
              for (auto front_idx = 0u; front_idx < fronts.size(); ++front_idx)
                  for (auto idx : fronts[front_idx])
                  {
                      const auto & ind = pop.get_individual(idx);
                      ofs
                      << idx++ << ","
                      << front_idx << ","
                      << ind.cur_f[0] << ","
                      << ind.cur_f[1] << ","
                      << ind.cur_f[2] << "\n";
                  }
          }
        }


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



boost::program_options::variables_map
handle_cmd_line_options(int argc, char * argv[])
{
	namespace po = boost::program_options;
    
	// Setup available options
	po::options_description opts("Program options");
	opts.add_options()
		("help,h", "Produce this help message")

        ("problem", po::value<std::string>()->default_value("problem.bin"),
         "File containing the serialized problem")

        ("population", po::value<std::string>()->default_value("population.bin"),
         "File containing the serialized population")
        
        ("display", po::value< std::vector<size_t> >()->multitoken()
         ->default_value(std::vector<size_t>())
         /*->implicit_value(std::vector<size_t>())*/,
         "List of individuals indices to display")

        ("export-vtu", po::value< std::vector<size_t> >()->multitoken()
         ->default_value(std::vector<size_t>())
         /*->implicit_value(std::vector<size_t>())*/,
         "List of individuals indices to export")

        ("export-csv", po::value<std::string>()->zero_tokens()->implicit_value("population.csv"),
         "Compute Pareto fronts and export them into a csv file")

        
        ("show-fitness", po::value< std::vector<size_t> >()->multitoken()
         ->default_value(std::vector<size_t>()),
         "List of individual indices to print fitnesses")

        ("show-all", po::value<bool>()->zero_tokens(), "Prints fitness of all individuals")
         
        ("show-champion", po::value<bool>()->zero_tokens(), "Prints the population champion")
		;

	// Store args
	po::variables_map vmap;
	po::store(po::parse_command_line(argc, argv, opts), vmap);

	// Check for help
	if (vmap.count("help"))
	{
		std::cout << opts << std::endl;
		throw eva::utils::help_exception();
	}

	// Otherwise notify errors
	po::notify(vmap);
	return vmap;
}
