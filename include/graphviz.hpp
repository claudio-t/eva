# ifndef __EVA_GRAPHVIZ__
# define __EVA_GRAPHVIZ__

/**
 * @file graphviz.hpp
 * @brief Contains utilities to read and write structures in graphviz
 * format (.dot files)
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

//------------------------------------------- Read -----------------------------------------------//
/// Reads and builds a structure from a proper input file
template <typename Structure>
Structure
read_from_graphviz(const std::string& filename);

/// Adds to the dynamic property map the properties required for
/// reading a given structure
template <typename Structure> 
void 
setup_joint_properties(Structure& structure, boost::dynamic_properties& dyna_props);

/// Functor that selects the joint properties that have to be red
/// from the input file. Has to be specialized for each joint type
template <typename JointProps> struct joint_properties_selector;

/// Specialization for reading a truss joint 
template <int N> struct joint_properties_selector<truss_joint<N>>;

/// Specialization for reading a 2D frame joint 
template <> struct joint_properties_selector<frame_joint<2>>;

/// Adds to the dynamic property map the properties required for
/// reading a given structure
template <typename Structure> 
void 
setup_element_properties(Structure& structure, boost::dynamic_properties& dyna_props);

/// Functor that selects the element properties that have to be red
/// from the input file. Has to be specialized for each element type
template <typename ElementProps> struct element_properties_selector;

/// Specialization for reading a truss_element
template <int N> struct element_properties_selector<truss_element<N>>;

/// Specialization for reading a 2D frame element
template <> struct element_properties_selector<frame_element<2>>;


//------------------------------------------ Write -----------------------------------------------//
template <typename Structure>
auto make_joint_properties_writer(const Structure& s);

/// Functor that selects the joint properties that have to be written
/// to the output file. Has to be specialized for each joint type
template <typename Structure, typename JointProperties>
class joint_properties_writer;

/// Specialization for writing a truss joint
template <typename Structure, int N>
class joint_properties_writer<Structure, truss_joint<N>>;

/// Specialization for writing a 2D frame joint
template <typename Structure>
class joint_properties_writer<Structure, frame_joint<2>>;

/// Utility function for pretty printing a fixed vector
template <int N> std::string
to_string(
    const fixed_vector<N>& v,
    const std::string lx_delim = "[",
    const std::string rx_delim = "]");



//####################################### DEFINITIONS ##############################################

//------------------------------------------- Read -----------------------------------------------//
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


//------------------------------------------ Write -----------------------------------------------//
template <typename S>
auto
make_joint_properties_writer(const S& s)
{
    return joint_properties_writer<S, typename S::vertex_bundled>(s);
}

template <typename Structure, int N>
class joint_properties_writer<Structure, truss_joint<N> >
{
public:
    using props = typename Structure::vertex_bundled;

    template <typename P>
    using property_map_t = boost::vec_adj_list_vertex_property_map<
        Structure, const Structure *,
        P, const P&,
        P eva::truss_joint<N>::*
        >;

    /// Constructor that initializes member data
    joint_properties_writer(const Structure& s)
        : coords_list_(get(&props::coords, s)),
          loads_list_ (get(&props::load,   s)),
          bcs_list_   (get(&props::bcs,    s))
        {}
    
    template <class Vertex>
    void operator()(std::ostream& os, const Vertex& v) const
    {
        os << '['
           << "coords = \"" << to_string(coords_list_[v]) << "\", " 
           << "load = \""   << to_string(loads_list_[v])  << "\", "
           << "bcs = \""    << to_string(bcs_list_[v])    << "\", "
           // << "pos = \""    << to_string(coords_list_[v], "", "") << "!\"" 
           << ']';
    }
    
private:
    const property_map_t<fixed_vector<N> > coords_list_;
    const property_map_t<fixed_vector<N> >  loads_list_;
    const property_map_t<fixed_vector<N> >    bcs_list_;
};

template <typename Structure>
class joint_properties_writer<Structure, frame_joint<2> >
{
public:
    using props = typename Structure::vertex_bundled;

    template <typename P>
    using property_map_t = boost::vec_adj_list_vertex_property_map<
        Structure, const Structure *,
        P, const P&,
        P eva::frame_joint<2>::*
        >;

    /// Constructor that initializes member data
    joint_properties_writer(const Structure& s)
        : coords_list_ (get(&props::coords, s)),
          loads_list_  (get(&props::load,   s)),
          bcs_list_    (get(&props::bcs,    s)),
          torques_list_(get(&props::torque, s))
        {}

    /// Actual functor implementation
    template <class Vertex>
    void operator()(std::ostream& os, const Vertex& v) const
    {
        os << '['
           << "coords = \"" << to_string(coords_list_[v]) << "\", " 
           << "load = \""   << to_string(loads_list_[v])  << "\", "
           << "bcs = \""    << to_string(bcs_list_[v])    << "\", "
           // << "torque = \"" << torques_list[v] << "\"" 
           << ']';
    }
    
private:
    const property_map_t<fixed_vector<2> >  coords_list_;
    const property_map_t<fixed_vector<2> >   loads_list_;
    const property_map_t<fixed_vector<3> >     bcs_list_;
    const property_map_t<fixed_vector<1> > torques_list_;
    
};

template <int N>
std::string
to_string(const fixed_vector<N>& v,
          const std::string lx_delim,
          const std::string rx_delim)
{
    std::ostringstream oss;
    oss << lx_delim;
    for(auto i = 0u; i < N-1; ++i) oss << v[i] << ", ";
    oss << v[N-1];
    oss << rx_delim;

    return oss.str();
}

} //end namespace eva

namespace boost { namespace detail {

template<int N, typename T>
struct lexical_cast_do_cast<eva::fixed_vector<N,T>, std::string>
{     
    static void  
    do_assign(const double value, eva::fixed_vector<N,T>& vec, const size_t idx) 
    {
        vec(idx) = value;
    }
    
    static inline
    eva::fixed_vector<N,T>
    lexical_cast_impl(const std::string& str) 
    {
        namespace qi    = boost::spirit::qi;
        namespace phx   = boost::phoenix;
        // namespace ascii = boost::spirit::ascii;
        
        // Alias for parser input iterator type
        using iterator_type = std::string::const_iterator;
        
        // Alias for a parser that parses a T
        auto t_ = qi::real_parser<T>();
        static_assert(std::is_floating_point<T>::value, 
                       "Integral types are not supported yet.");
        
        qi::rule<iterator_type, eva::fixed_vector<N,T>(),
                 qi::locals<size_t>, qi::blank_type> rule = 
            qi::eps [ qi::_a = 0 ]
            > '['
            > qi::repeat(N-1) [
                t_[phx::bind(&lexical_cast_do_cast::do_assign, qi::_1, qi::_val, qi::_a++)]
                >> ',' 
                ]
            > t_[phx::bind(&lexical_cast_do_cast::do_assign, qi::_1, qi::_val, qi::_a)]
            > ']'
            ;
        // qi::rule<
        //     iterator_type, eva::fixed_vector<N,T>(),
        //     qi::locals<size_t>, qi::blank_type
        //     > rule = 
        //     qi::eps [ qi::_a = 0 ]
        //     > '['
        //     > qi::repeat(N-1) [
        //         t_[phx::bind(&assign_value_action<N,T>, qi::_1, qi::_val, qi::_a++)]
        //         > ',' 
        //     ]
        //     > t_[phx::bind(&assign_value_action<N,T>, qi::_1, qi::_val, qi::_a)]
        //     > ']'
        // ;
        
        rule.name("vector rule");
        qi::on_error<qi::fail> (
            rule,
            std::cerr
                << phx::val("Error! Expecting ")
                << qi::_4                               // what failed?
                << phx::val(" here: \"")
                << phx::construct<std::string>(qi::_3, qi::_2)  // iterators to error-pos, end
                << phx::val("\"")
                << std::endl
        );
                  
        auto ret = eva::fixed_vector<N,T>();
        
        auto iter(str.begin()), iend(str.end());
        
        bool r = phrase_parse(iter, iend, rule, qi::blank, ret);
        
        if (!r || iter != iend) {
            std::cout << "-------------------------\n";
            std::cout << "Parsing failed\n";
            std::cout << str << std::endl;
            std::cout << "-------------------------\n";
        }
        return ret;
    }
};


template<>
struct lexical_cast_do_cast<eva::fixed_vector<1, double>, std::string>
{ 
    static inline
    eva::fixed_vector<1, double>
    lexical_cast_impl(const std::string& str) 
    {
        auto ret = eva::fixed_vector<1, double>();
        ret << std::stod(str);
        return ret;
    }
};
template<>
struct lexical_cast_do_cast<eva::fixed_vector<1, float>, std::string>
{ 
    static inline
    eva::fixed_vector<1, float>
    lexical_cast_impl(const std::string& str) 
    {
        auto ret = eva::fixed_vector<1, float>();
        ret << std::stof(str);
        return ret;
    }
};
template<>
struct lexical_cast_do_cast<eva::fixed_vector<1, int>, std::string>
{ 
    static inline
    eva::fixed_vector<1, int>
    lexical_cast_impl(const std::string& str) 
    {
        auto ret = eva::fixed_vector<1, int>();
        ret << std::stoi(str);
        return ret;
    }
};

} } // end details and boost

# endif //__EVA_GRAPHVIZ__
