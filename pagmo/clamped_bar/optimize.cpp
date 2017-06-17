/////////////////////////////////////////////////////
//eva
# include "thermo_frame_problem.hpp"
# include "thermo_frame_problem_utils.hpp"
# include "utils.hpp"
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
using namespace utils::thermo_frame_problem;

// These will be removed once a proper config file entry will be available
void custom_add_bcs  (joint_grid_t & grid, const po::variables_map & config);
void custom_add_loads(joint_grid_t & grid, const po::variables_map & config);


int main(int argc, char * argv [])
{
    try
    {
        // Read cmd line options
        auto cmd_line_opts  = eva::utils::handle_cmd_line_options(argc, argv);
    
        // Read config file
        auto config_file = cmd_line_opts["config"].as<std::string>();
        auto config = read_config_file(config_file);

        // Make regular grid
        auto l = config["length"].as<real>();//100.e-3;
        auto h = config["height"].as<real>();//25.0e-3;
        auto m = config["rows"].as<int>();//30;
        auto n = config["cols"].as<int>();//50;
        auto grid = make_grid(l, h, m, n);

        // Add bcs & loads
        custom_add_bcs  (grid, config);
        custom_add_loads(grid, config);


        // Make topology
        auto topology = make_topology(grid, neumann_rule());
        
        // Make base frame
        auto structure = make_frame(grid, topology, config);
        eva::display(structure);
        
        // Read boundary map file
        auto boundary_map = make_boundary_map(grid);
        
        // Build problem
        auto problem = make_problem(structure, boundary_map, topology, config);

        // Setup NSGA2 algorithm
        auto gens  = config["generations"].as<size_t>(); assert(gens > 0);
        auto cr    = 0.9;
        auto eta_c = 10;
        auto mut   = 0.05;
        auto eta_m = 50;

        auto algorithm = pagmo::algorithm::nsga2(gens, cr,eta_c, mut, eta_m).clone();
        algorithm->set_screen_output(true);

                                             
        // Create the archipelago
        auto n_isls = config[  "islands"  ].as<size_t>(); assert(n_isls > 0);
        auto n_inds = config["individuals"].as<size_t>(); assert(n_inds > 0);
        
        auto archi_topo = pagmo::topology::ring();
        auto archi = pagmo::archipelago(*algorithm, problem, n_isls, n_inds, archi_topo);

        // Evolve archipelago
        auto n_evo   = config["evolutions"].as<size_t>(); assert(n_evo > 0);
        auto da_best = pagmo::population::champion_type();
        
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
            da_best = get_champion(archi);
            
            // Print fitness
            std::cout << "Best Fitness = " << da_best.f << std::endl;
        }

        // Save population & problem
        auto pop = std::vector<pagmo::population::individual_type>();
        for (auto i = 0u; i < archi.get_size(); ++i)
        {
            auto ind_idx = 0u;
                for (const auto & ind : archi.get_island(i)->get_population())
                    pop.emplace_back(ind);
        }
        eva::utils::serialize(pop, "population.bin");
        eva::utils::serialize(problem, "problem.bin");
        export_to_vtu(problem.encode_genes(da_best.x), "da_best.vtu");
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


void custom_add_bcs(joint_grid_t & grid, const po::variables_map & config)
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


void custom_add_loads(joint_grid_t & grid, const po::variables_map & config)
{ 
    const auto m = grid.size();
    const auto n = grid.front().size();
    
    // Apply mechanical load on joint [lowest_row][3] (lowest row)
    grid[0][3].load << 0, 16 * 1.e4;

    // Apply thermal load on joint [highest_row][2]
    grid[m-1][2].flux_bc = 0.5 * 1.3;
}

