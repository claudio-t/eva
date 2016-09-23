# include <iostream>
# include <fstream>
# include <boost/graph/adjacency_list.hpp>
# include <boost/graph/graph_traits.hpp>
# include <boost/range/iterator_range.hpp>
# include <boost/graph/graphviz.hpp>

// Vertex properties
typedef boost::property< boost::vertex_name_t, std::string > vertex_p;
// Edge properties
typedef boost::property< boost::edge_weight_t, double > edge_p;
// Graph properties
typedef boost::property< boost::graph_name_t, std::string > graph_p;
// Adjacency_list-based type
typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS, 
                                vertex_p, edge_p > graph_t;

int main(int argc, char * argv[]) {

    // Init an empty graph and the dynamic_property_maps.
    auto graph = graph_t (0);
    auto dynap = boost::dynamic_properties();
        
    // Setup properties to be read
    {
        auto name   = get(boost::vertex_name, graph);
        auto weight = get(boost::edge_weight, graph);
        dynap.property("node_id",   name);
        dynap.property( "weight", weight);
        
        
        // Read
        std::ifstream ifile ("graph.dot");
        if (!ifile.is_open()) std::cerr << "Cannot open file" << std::endl;
        
        bool status = boost::read_graphviz(ifile, graph, dynap, "node_id");
    }
    
    auto name   = get(boost::vertex_name, graph);
    auto weight = get(boost::edge_weight, graph);
    std::cout << typeid(name).name() << std::endl;
    for(auto&& v : make_iterator_range(vertices(graph))) {
        std::cout << name[v] << std::endl;
    }
    for(auto&& e : make_iterator_range(edges(graph))) {
        std::cout << weight[e] << std::endl;
    }
    
    return 0;
}


