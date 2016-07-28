# include <iostream>
# include "truss.hpp"
# include "frame.hpp"
# include "io.hpp"


int main(int argc, char * argv [])
{
    // Read from file
    if (argc < 2) {
        std::cout << "Input file missing!\n";
        return 0;
    }
    auto filename = std::string(argv[1]);
    
    auto structure = eva::read_from_graphviz<eva::truss3d>(filename);
    
    auto u = std::vector<eva::real> ();
    auto f = std::vector<eva::real> ();
    std::tie(u,f) = solve(structure, eva::dense_solver_params());
    
    std::cout << "-----------------------------------------------" 
              << std::endl << "Displacements:\n";
    for (const auto& el : u) std::cout << el << std::endl;
    //~ 
    std::cout << "-----------------------------------------------" 
              << std::endl << "Reactions:\n";
    for (const auto& el : f) std::cout << el << std::endl;
    
    std::cout << "-----------------------------------------------" 
              << std::endl << "Compliance:\n";
    auto compliance = 0.;
    for(auto i = 0u; i < f.size(); ++i) compliance += u[i]*f[i];
    std::cout << compliance << std::endl;

    auto element        = decltype(structure)::edge_descriptor(); 
    auto element_exists = false;
    std::tie(element, element_exists) = edge(1, 2, structure);
    
    if (element_exists) {
        auto f_element = get_internal_forces(element, structure, u, f);
        std::cout << "NA ->  " << std::get<0>(f_element) << std::endl;
        std::cout << "NB ->  " << std::get<1>(f_element) << std::endl;
    }
    else {
        std::cout << "Warning: element " << element << " does not exists!\n";
    }
    
    ///////////////////////////////////////////////////////////////
    // // Init properties to write
    // auto dyna_props = boost::dynamic_properties();
    // dyna_props.property("node_id", get(boost::vertex_index, structure)); 
    // // Node props
    // using vertex_prop_t = decltype(structure)::vertex_bundled;
    // dyna_props.property("coords", get(&vertex_prop_t::coords, structure)); 
    // dyna_props.property("torque", get(&vertex_prop_t::torque, structure)); 
    // dyna_props.property(  "load",   get(&vertex_prop_t::load, structure)); 
    // dyna_props.property(   "bcs",    get(&vertex_prop_t::bcs, structure)); 
    // // Edge props
    // using edge_prop_t = decltype(structure)::edge_bundled;
    // dyna_props.property("E", get(&edge_prop_t::E, structure)); 
    // dyna_props.property("A", get(&edge_prop_t::A, structure)); 
    // dyna_props.property("I", get(&edge_prop_t::I, structure)); 
    // // Write to sd::out
    // write_graphviz_dp(std::cout, structure, dyna_props);
    
    return 0;
}

