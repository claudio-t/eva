# include <iostream>
# include "truss.hpp"
# include "frame.hpp"
# include "io.hpp"

int main(int argc, char * argv [])
{
    auto filename  = std::string();
    
    // Check for input file
    if (argc < 2) {
        // std::cout << "Input file missing!\n";
        // return 0;
        filename = std::string("/home/claudio/eva/input/frame2d.dot");
    }
    else
        filename = std::string(argv[1]);

    // Read structure from file
    auto structure = eva::read_from_graphviz<eva::frame2d>(filename);

    // Solve structure
    auto results = solve(structure, eva::sparse_solver_params<>());
    
    for (auto idx = 0u; idx < num_vertices(structure); ++idx)
    {
        std::cout
            << "Displacement = "
            << results[idx].displacement.transpose()
            << std::endl;

        std::cout
            << "Reaction = "
            << results[idx].reaction.transpose()
            << std::endl;

        std::cout
            << "Reaction Torque = "
            << results[idx].react_torque.transpose()
            << std::endl;
        
        std::cout
            << "BCs = "
            << structure[idx].bcs.transpose()
            << std::endl;

        std::cout
            << "Load = "
            << structure[idx].load.transpose()
            << std::endl;

        std::cout
            << "Reaction Torque = "
            << structure[idx].torque.transpose()
            << std::endl;
        
        std::cout
            << "Coords = "
            << structure[idx].coords.transpose()
            << std::endl;
        
        std::cout << std::endl;
    }
    
    // // Compute&print compliance
    // std::cout << "-----------------------------------------------" 
    //           << std::endl << "Compliance:\n";
    // auto compliance = 0.;
    // for (const auto& r : res) compliance += r.displacement.transpose() * r.force;
    // std::cout << compliance << std::endl;
    
    // Display undeformed structure
    display(structure);
    
    // Save to vtu for paraview
    write_vtu(structure, results, "test.vtu");
    
    return 0;
}

