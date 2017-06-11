# include "thermo.hpp"
# include "truss.hpp"
# include "frame.hpp"
# include "io.hpp"
# include <iostream>

# include <boost/archive/binary_iarchive.hpp>
# include <boost/archive/binary_oarchive.hpp>


template <typename T>
bool serialize(const T& data, const std::string & filename)
{
    std::ofstream ofs(filename.c_str(), std::ios::out);
    if (!ofs.is_open())
        return false;
    {
        boost::archive::binary_oarchive oa(ofs);
        oa << data;
    }
    ofs.close();
    return true;
}

template <typename T>
bool deserialize(T& data, const std::string & filename)
{
    std::ifstream ifs(filename.c_str(), std::ios::in);
    if (!ifs.is_open())
        return false;
    {
        boost::archive::binary_iarchive ia(ifs);
        ia >> data;
    }
    ifs.close();
    return true;
}

int main(int argc, char * argv [])
{
    auto filename  = std::string();
    
    // Check for input file
    if (argc < 2) {
        // std::cout << "Input file missing!\n";
        // return 0;
        filename = std::string("/home/claudio/eva/input/thermo_simple.dot");
    }
    else
        filename = std::string(argv[1]);

    // Setup aliases
    using struct_kind_t = eva::frame_kind<2>;
    using thermo_kind_t = eva::thermo_kind<struct_kind_t>;
    using thermo_structure_t = eva::thermo_structure<struct_kind_t>;

    // Read thermal+structural data & build structure from file
    auto structure = eva::read_from_graphviz<thermo_structure_t>(filename);    
    
    // Solve structural problem
    auto struct_results = solve(structure, struct_kind_t());

    // Solve thermal problem
    using solver_t = eva::dense_solver_params<thermo_kind_t::default_dense_solver_t>;
    auto thermo_results = solve(structure, thermo_kind_t(), solver_t());
    
    for (auto idx = 0u; idx < num_vertices(structure); ++idx)
    {
        // std::cout
        //     << "Temperature = "
        //     << results[idx].T
        //     << std::endl;

        // std::cout
        //     << "Flux = "
        //     << results[idx].flux
        //     << std::endl;
        // std::cout
        //     << "Displacement = "
        //     << struct_results[idx].displacement.transpose()
        //     << std::endl;

        // std::cout
        //     << "Reaction = "
        //     << struct_results[idx].reaction.transpose()
        //     << std::endl;

        // std::cout
        //     << "BCs = "
        //     << structure[idx].bcs.transpose()
        //     << std::endl;

        
        // std::cout
        //     << "Load = "
        //     << structure[idx].load.transpose()
        //      << "[N]"
        //     << std::endl;
        std::cout
            << "Idx = " << idx << "\n";
        
        std::cout
            << "Dirichlet BC = "
            << structure[idx].T_bc
            << "[K]"
            << std::endl;

        std::cout
            << "Heat Flux BC= "
            << structure[idx].flux_bc
             << "[W]"
            << std::endl;
        
        std::cout
            << "Coords = "
            << structure[idx].coords.transpose()
            << std::endl;

    std::cout
            << "Temperature = "
            << thermo_results[idx].T
            << "[K]"
            << std::endl;
        
        
        std::cout << std::endl;
    }

    for (auto e : boost::make_iterator_range(edges(structure)))
        std::cout << "(" << source(e, structure) << ", " << target(e, structure) << ") --> "
          << compute_heat_flow(e, structure, thermo_results) << "\n";
    

    // // Display
    // eva::display(structure);

    // // Save results to vtu file
    // auto vtk_grid = eva::to_vtk_unstructured_grid(structure);
    
    // eva::vtk_add_joint_properties  (structure, vtk_grid);
    // eva::vtk_add_element_properties(structure, vtk_grid);
    
    // eva::vtk_add_joint_results(thermo_results, vtk_grid);
    // eva::vtk_add_joint_results(struct_results, vtk_grid);
    
    // eva::write_vtu(vtk_grid, "test.vtu");

    // // Write ref file
    // boost::write_graphviz(
    //     std::cout, structure,
    //     eva::make_joint_properties_writer(structure),
    //     eva::make_element_properties_writer(structure));
    
    // eva::write_graphviz(structure, "ref.dot");

    // // Serialize
    // serialize(structure, "serialized_structure.bin");

    // // Read from binary
    // auto des_struct = thermo_structure_t();
    // deserialize(des_struct, "serialized_structure.bin");

    // boost::write_graphviz(
    //     std::cout, structure,
    //     eva::make_joint_properties_writer(structure),
    //     eva::make_element_properties_writer(structure));
    
    // eva::write_graphviz(des_struct, "deserialized.dot");
    return 0;
}
