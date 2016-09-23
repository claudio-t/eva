# include <iostream>
# include "truss.hpp"
# include "frame.hpp"
# include "io.hpp"
# include "utils.hpp"

int main(int argc, char * argv [])
{
    // Check for input file
    if (argc < 2) {
        std::cout << "Input file missing!\n";
        return 0;
    }
    
    // Read structure from file
    auto filename  = std::string(argv[1]);
    auto structure = eva::read_from_graphviz<eva::frame2d>(filename);

    // Solve structure
    auto results = solve(structure, eva::sparse_solver_params<>());

    // Print properties
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
            << "BCs = "
            << structure[idx].bcs.transpose()
            << std::endl;

        std::cout
            << "Load = "
            << structure[idx].load.transpose()
            << std::endl;
        
        std::cout
            << "Coords = "
            << structure[idx].coords.transpose()
            << std::endl;
        
        std::cout << std::endl;
    }

    // Compute compliance
    auto compliance = 0.;
    for (size_t vidx = 0u; vidx < num_vertices(structure); ++vidx)
        compliance += structure[vidx].load.transpose() * results[vidx].displacement;
    // std::cout << "Compliance = "<< compliance << std::endl;
    
    // Display structure
    display(structure);
    
    // Print properties (test)
    // auto ugrid = to_vtk_unstructured_grid(structure);
    // auto * cellData = ugrid->GetCellData();
    // for (int i = 0; i < cellData->GetNumberOfArrays(); ++i)
    // {
    //     vtkDataArray* data = cellData->GetArray(i);
    //     cout << "name " << data->GetName() << endl;
    //     for (int j = 0; j < data->GetNumberOfTuples(); ++j)
    //     {
    //         double value = data->GetTuple1(i);
    //         cout << "  value " << j << "th is " << value << endl;
    //     }
    // 

    // // Retrieve element member forces
    // auto element        = decltype(structure)::edge_descriptor(); 
    // auto element_exists = false;
    // std::tie(element, element_exists) = edge(1, 2, structure);
    
    // if (element_exists) {
    //     auto f_element = get_internal_forces(element, structure, u, f);
    //     std::cout << "NA ->  " << std::get<0>(f_element) << std::endl;
    //     std::cout << "NB ->  " << std::get<1>(f_element) << std::endl;
    // }
    // else {
    //     std::cout << "Warning: element " << element << " does not exist!\n";
    // }
    
    ///////////////////////////////////////////////////////////////

    // using v_prop = decltype(structure)::vertex_bundled;
    // std::ofstream os;
    // os.open("output/frame2d.dot");    
    // write_graphviz(os, structure, make_joint_properties_writer(structure));
    
    
    // Init properties to write
    // auto dyna_props = boost::dynamic_properties();
    // dyna_props.property("node_id", get(boost::vertex_index, structure)); 
    // // Node props
    // using vertex_prop_t = decltype(structure)::vertex_bundled;
    // dyna_props.property("coords", get(&vertex_prop_t::coords, structure)); 
    // // dyna_props.property("torque", get(&vertex_prop_t::torque, structure)); 
    // dyna_props.property(  "load",   get(&vertex_prop_t::load, structure)); 
    // dyna_props.property(   "bcs",    get(&vertex_prop_t::bcs, structure));
    // // dyna_props.property( "dummy",  get(&vertex_prop_t::dummy, structure)); 
    // // Edge props
    // using edge_prop_t = decltype(structure)::edge_bundled;
    // dyna_props.property("E", get(&edge_prop_t::E, structure)); 
    // dyna_props.property("A", get(&edge_prop_t::A, structure)); 
    // // dyna_props.property("I", get(&edge_prop_t::I, structure)); 
    // // Write to sd::out
    // write_graphviz_dp(std::cout, structure, dyna_props);
    
    return 0;
}

