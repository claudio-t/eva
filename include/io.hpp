# ifndef __EVA_IO__
# define __EVA_IO__

/**
 * @file io.hpp
 * @brief Contains utilities to read and write structures.
 */
 
// std
# include <string>
# include <fstream>
 
// boost
# include <boost/graph/graphviz.hpp>
# include <boost/spirit/include/qi.hpp>
# include <boost/spirit/include/phoenix.hpp>


namespace eva {
    
//####################################### DECLARATIONS #############################################

// -- Read -- //

/// Reads and builds a structure from a proper input file
template <typename Structure>
Structure
read_from_graphviz(const std::string& filename);

/// Adds to the dynamic property map the properties required for reading a given structure
template <typename Structure> 
void 
setup_joint_properties(Structure& structure, boost::dynamic_properties& dyna_props);

/// Functor that selects the joint properties that have to be red from the input file.
/// Has to be specialized for each joint type
template <typename JointProps> struct joint_properties_selector;

/// Specializes joint_properties_selector for reading a truss_joint 
template <int N> struct joint_properties_selector<truss_joint<N>>;

/// Specializes joint_properties_selector for reading a frame_joint<2> 
template <> struct joint_properties_selector<frame_joint<2>>;


template <typename Structure> 
void 
setup_element_properties(Structure& structure, boost::dynamic_properties& dyna_props);

template <typename ElementProps> struct element_properties_selector;

template <int N> struct element_properties_selector<truss_element<N>>;

//####################################### DEFINITIONS ##############################################

// -- Read -- //
template <typename Structure>
Structure
read_from_graphviz(const std::string& filename) 
{    
    // Build empty structure //
    auto structure = Structure();
    
    // Setup properties to be read
    auto dyna_props = boost::dynamic_properties(boost::ignore_other_properties);
    
    setup_joint_properties  (structure, dyna_props);
    setup_element_properties(structure, dyna_props);
    
    // Read from file
    std::ifstream ifile(filename.c_str());
    boost::read_graphviz(ifile, structure, dyna_props);

    return structure;
}


template <int N>
struct joint_properties_selector<truss_joint<N>> 
{    
    template <typename Structure> 
    void 
    operator()(Structure& s, boost::dynamic_properties& dps) 
    {
        // Type checking
        using vert_prop_type = typename Structure::vertex_bundled;
        static_assert(std::is_base_of<truss_joint<N>, vert_prop_type>::value, 
                      "Unfeasible structure type");
        
        // Read coordinates, loads and boundary conditions
        dps.property("coords", get(&vert_prop_type::coords, s));
        dps.property(  "load", get(&vert_prop_type::load,   s));
        dps.property(   "bcs", get(&vert_prop_type::bcs,    s));
    }
};


template <int N>
struct element_properties_selector<truss_element<N>> 
{    
    template <typename Structure> 
    void 
    operator()(Structure& s, boost::dynamic_properties& dps)
    {
        // Type checking
        using edge_prop_type = typename Structure::edge_bundled;
        static_assert(std::is_base_of<truss_element<N>, edge_prop_type>::value, 
                      "Unfeasible structure type");
        
        // Read Young's modulus and cross sectional area
        dps.property("E", get(&edge_prop_type::E, s));
        dps.property("A", get(&edge_prop_type::A, s));
    }
};


template <>
struct joint_properties_selector<frame_joint<2>> 
{    
    template <typename Structure> 
    void 
    operator()(Structure& s, boost::dynamic_properties& dps) 
    {
        // Type checking
        using vert_prop_type = typename Structure::vertex_bundled;
        
        static_assert(std::is_base_of<frame_joint<2>, vert_prop_type>::value, 
                      "Unfeasible structure type");
        
        // Read coordinates, torques, loads and boundary conditions
        dps.property("coords", get(&vert_prop_type::coords, s));
        dps.property("torque", get(&vert_prop_type::torque, s));
        dps.property(  "load", get(&vert_prop_type::load,   s));
        dps.property(   "bcs", get(&vert_prop_type::bcs,    s));

    }
};

template <>
struct element_properties_selector<frame_element<2>> 
{    
    template <typename Structure> 
    void 
    operator()(Structure& s, boost::dynamic_properties& dps)
    {
        // Type checking
        using edge_prop_type = typename Structure::edge_bundled;
        static_assert(std::is_base_of<frame_element<2>, edge_prop_type>::value, 
                      "Unfeasible structure type");
        
        // Read Young's modulus and cross sectional area
        dps.property("E", get(&edge_prop_type::E, s));
        dps.property("A", get(&edge_prop_type::A, s));
        dps.property("I", get(&edge_prop_type::I, s));
    }
};


template <typename Structure> 
void
setup_joint_properties(Structure& s, boost::dynamic_properties& dps) 
{
    joint_properties_selector<typename Structure::vertex_bundled>()(s, dps);
}


template <typename Structure> 
void
setup_element_properties(Structure& s, boost::dynamic_properties& dps) 
{
    element_properties_selector<typename Structure::edge_bundled>()(s, dps);
}


} //end namespace eva

namespace boost { namespace detail {
    
template <int N, typename T>
inline void 
assign_value_action(const double value, eva::fixed_vector<N,T>& vec, const size_t idx) 
{
    vec(idx) = value;
}



template<int N, typename T>
struct lexical_cast_do_cast<eva::fixed_vector<N,T>, std::string>
{         
    /// Alias for parser input iterator type
    using iterator_type = std::string::const_iterator;
    
    static inline
    eva::fixed_vector<N,T>
    lexical_cast_impl(const std::string& str) 
    {
        namespace qi    = boost::spirit::qi;
        namespace phx   = boost::phoenix;
        //~ namespace ascii = boost::spirit::ascii;
        
        qi::rule<iterator_type, eva::fixed_vector<N,T>(),
                 qi::locals<int>, qi::blank_type> rule = 
            qi::eps [ qi::_a = 0 ]
            >> '['
            >> qi::repeat(N-1) [
                qi::double_[phx::bind(&assign_value_action<N,T>, qi::_1, qi::_val, qi::_a++)]
                >> ',' 
            ]
            >> qi::double_[phx::bind(&assign_value_action<N,T>, qi::_1, qi::_val, qi::_a)]
            >> ']'
        ;
                  
        auto ret = eva::fixed_vector<N,T>();
        
        auto iter(str.begin()), iend(str.end());
        
        bool r = phrase_parse(iter, iend, rule, qi::blank, ret);
        
        //~ if (r && iter == iend) { 
            //~ std::cout << "Parsing succeded\n";
        //~ } else {
        if (!r || iter != iend) {
            std::cout << "-------------------------\n";
            std::cout << "Parsing failed\n";
            std::cout << str << std::endl;
            //~ std::cout << r << "\n" << (iter==iend) << "\n";
            std::cout << "-------------------------\n";
        }
        return ret;
    }
    /* OLD, UGLY, SLOW
    static inline
    eva::fixed_vector<N,T>
    lexical_cast_impl(const std::string& str) 
    {
        auto line = str;                         // Temporary
        std::istringstream iss (std::move(line));// Move construct istringstream
        
        auto vec = eva::fixed_vector<N>();
        for (auto i = 0; i < N; ++i) {
            std::string str;
            iss >> str;
            // FIXME: what for size_t, int etc..??
            if (std::is_same<T, double>::value) vec(i) = std::stod(str);
            if (std::is_same<T,  float>::value) vec(i) = std::stof(str);
            if (std::is_same<T,    int>::value) vec(i) = std::stoi(str);
        }
        return vec;
    }*/
};
} } // end details and boost

# endif //__EVA_IO__

