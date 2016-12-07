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
Structure read_from_graphviz(const std::string& filename);

/// Adds to the dynamic property map the properties required for
/// reading a given structure
template <typename S> 
void 
setup_joint_properties(S& structure, boost::dynamic_properties& dyna_props);

/// Functor that selects the joint properties that have to be red
/// from the input file. Has to be specialized for each joint type
template <typename JointProps>
struct properties_selector;

/// Specialization for reading a truss joint 
template <int N>
struct properties_selector< sa::truss_joint<N>>;

/// Specialization for reading a 2D frame joint 
template <>
struct properties_selector< sa::frame_joint<2>>;

/// Adds to the dynamic property map the properties required for
/// reading a given structure
template <typename S> 
void 
setup_element_properties(S& structure, boost::dynamic_properties& dyna_props);

/// Functor that selects the element properties that have to be red
/// from the input file. Has to be specialized for each element type
template <typename ElementProps>
struct properties_selector;

/// Specialization for reading a truss_element
template <int N>
struct properties_selector< sa::truss_element<N> >;

/// Specialization for reading a 2D frame element
template <>
struct properties_selector< sa::frame_element<2> >;


//------------------------------------------ Write -----------------------------------------------//
template <typename S>
auto make_joint_properties_writer(const S& structure);

/// Functor that selects the joint properties that have to be written
/// to the output file. Has to be specialized for each joint type
template <typename S, typename JointProperties>
class properties_writer;

/// Specialization for writing a truss joint
template <typename S, int N>
class properties_writer< S, sa::truss_joint<N> >;

/// Specialization for writing a 2D frame joint
template <typename S>
class properties_writer< S, sa::frame_joint<2> >;

template <typename S>
auto make_element_properties_writer(const S& structure);

/// Utility function for pretty printing a fixed vector
template <int N> std::string
to_string(
    const fixed_vector<N>& v,
    const std::string lx_delim = "[",
    const std::string rx_delim = "]");



//####################################### DEFINITIONS ##############################################

//------------------------------------------- Read -----------------------------------------------//
template <typename S>
S read_from_graphviz(const std::string& filename) 
{    
    // Build empty structure //
    auto structure = S();
    
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
struct properties_selector< sa::truss_joint<N> > 
{    
    template <typename S> 
    void 
    operator()(S& structure, boost::dynamic_properties& dps) 
    {
        using vert_prop_type = typename joint_of<S>::type;
        
        // Read coordinates, loads and boundary conditions
        dps.property("coords", get(&vert_prop_type::coords, structure));
        dps.property(  "load", get(&vert_prop_type::load,   structure));
        dps.property(   "bcs", get(&vert_prop_type::bcs,    structure));
    }
};


template <int N>
struct properties_selector< sa::truss_element<N> > 
{    
    template <typename S> 
    void 
    operator()(S& structure, boost::dynamic_properties& dps)
    {
        using edge_prop_type = typename element_of<S>::type;
        
        // Read Young's modulus and cross sectional area
        dps.property("E", get(&edge_prop_type::E, structure));
        dps.property("A", get(&edge_prop_type::A, structure));
    }
};


template <>
struct properties_selector< sa::frame_joint<2> > 
{    
    template <typename S> 
    void 
    operator()(S& structure, boost::dynamic_properties& dps) 
    {
        using vert_prop_type = typename joint_of<S>::type;
        
        // Read coordinates, torques, loads and boundary conditions
        dps.property("coords", get(&vert_prop_type::coords, structure));
        dps.property("torque", get(&vert_prop_type::torque, structure));
        dps.property(  "load", get(&vert_prop_type::load,   structure));
        dps.property(   "bcs", get(&vert_prop_type::bcs,    structure));

    }
};

template <>
struct properties_selector< sa::frame_element<2> > 
{    
    template <typename S> 
    void 
    operator()(S& structure, boost::dynamic_properties& dps)
    {
        using edge_prop_type = typename element_of<S>::type;
        
        // Read Young's modulus and cross sectional area
        dps.property("E", get(&edge_prop_type::E, structure));
        dps.property("A", get(&edge_prop_type::A, structure));
        dps.property("I", get(&edge_prop_type::I, structure));
    }
};


template <typename S> 
void
setup_joint_properties(S& s, boost::dynamic_properties& dps) 
{
    properties_selector<typename joint_of<S>::type>()(s, dps);
}


template <typename S> 
void
setup_element_properties(S& s, boost::dynamic_properties& dps) 
{
    properties_selector<typename element_of<S>::type>()(s, dps);
}


//------------------------------------------ Write -----------------------------------------------//
template <typename S>
auto
make_joint_properties_writer(const S& s)
{
    return properties_writer<S, typename joint_of<S>::type>(s);
}

template <typename S, int N>
class properties_writer< S, sa::truss_joint<N> >
{
public:
    using props = typename joint_of<S>::type;

    template <typename P>
    using property_map_t = boost::vec_adj_list_vertex_property_map<
        S, const S *,
        P, const P&,
        P eva::sa::truss_joint<N>::*
        >;

    /// Constructor that initializes member data
    properties_writer(const S& structure)
        : coords_list_(get(&props::coords, structure))
        , loads_list_ (get(&props::load,   structure))
        , bcs_list_   (get(&props::bcs,    structure))
        {}

    template <class Vertex>
    void write_impl(std::ostream& os, const Vertex& v) const
    {
        os << "coords = \"" << to_string(coords_list_[v]) << "\", " 
           << "load = \""   << to_string(loads_list_[v])  << "\", "
           << "bcs = \""    << to_string(bcs_list_[v])    << "\", ";
    }
    
    template <class Vertex>
    void operator()(std::ostream& os, const Vertex& v) const
    {
        os << '[';
        write_impl(os, v);
        os << ']';
    }
    
protected:
    const property_map_t< fixed_vector<N> > coords_list_;
    const property_map_t< fixed_vector<N> >  loads_list_;
    const property_map_t< fixed_vector<N> >    bcs_list_;
};

template <typename S>
class properties_writer< S, sa::frame_joint<2> >
    // : public properties_writer< S, sa::truss_joint<2> >
{
public:
    using props = typename joint_of<S>::type;
    using base_t = properties_writer< S, sa::truss_joint<2> >;

    template <typename P>
    using property_map_t = boost::vec_adj_list_vertex_property_map<
        S, const S *,
        P, const P&,
        P eva::sa::frame_joint<2>::*
        >;

    /// Constructor that initializes member data
    properties_writer(const S& structure)
        // : base_t(structure),
        : coords_list_ (get(&props::coords, structure)),
          loads_list_  (get(&props::load,   structure)),
          bcs_list_    (get(&props::bcs,    structure)),
          torques_list_(get(&props::torque, structure))
        {}

    /// Actual functor implementation
    template <class Vertex>
    void write_impl(std::ostream& os, const Vertex& v) const
    {
        os << "coords = \"" << to_string(coords_list_[v]) << "\", " 
           << "load = \""   << to_string(loads_list_[v])  << "\", "
           << "bcs = \""    << to_string(bcs_list_[v])    << "\",";
        // base_t::write_impl(os, v);
        os << "torque = \"" << torques_list_[v] << "\"";
        
        // return os;
    }

    template <class Vertex>
    void operator()(std::ostream& os, const Vertex& v) const
    {
        os << '[';
        write_impl(os, v);
        os << ']';
    }
    
protected:
    const property_map_t< fixed_vector<2> >  coords_list_;
    const property_map_t< fixed_vector<2> >   loads_list_;
    const property_map_t< fixed_vector<3> >     bcs_list_;
    const property_map_t< fixed_vector<1> > torques_list_;
    
};


//------------------------------------------ Write -----------------------------------------------//

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
struct lexical_converter_impl< eva::fixed_vector<N,T>, std::string >
{     
    static void  
    do_assign(const double value, eva::fixed_vector<N,T>& vec, const size_t idx) 
    {
        vec(idx) = value;
    }
    
    static inline bool
    try_convert(const std::string& str, eva::fixed_vector<N,T>& vec) 
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

        // The actual parser
        qi::rule<iterator_type, eva::fixed_vector<N,T>(),
                 qi::locals<size_t>, qi::blank_type> rule = 
            qi::eps [ qi::_a = 0 ]
            > '['
            > qi::repeat(N-1) [
                t_[phx::bind(&lexical_converter_impl::do_assign, qi::_1, qi::_val, qi::_a++)]
                >> ',' 
                ]
            > t_[phx::bind(&lexical_converter_impl::do_assign, qi::_1, qi::_val, qi::_a)]
            > ']'
            ;
        
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
        
        auto iter(str.begin()), iend(str.end());
        
        bool r = phrase_parse(iter, iend, rule, qi::blank, vec);
        
        if (!r || iter != iend)
        {
            std::cout << "-------------------------\n";
            std::cout << "Parsing failed\n";
            std::cout << str << std::endl;
            std::cout << "-------------------------\n";
            return false;
        }
        else
        {
            return true;
        }
    }
};


template<>
struct lexical_converter_impl< eva::fixed_vector<1, double>, std::string >
{ 
    static inline bool
    try_convert(const std::string& str, eva::fixed_vector<1, double>& vec) 
    {
        vec << std::stod(str);
        return true;
    }
};
template<>
struct lexical_converter_impl< eva::fixed_vector<1, float>, std::string >
{ 
    static inline bool
    try_convert(const std::string& str, eva::fixed_vector<1, float>& vec) 
    {   
        vec << std::stof(str);
        return true;
    }
};
template<>
struct lexical_converter_impl< eva::fixed_vector<1, int>, std::string >
{ 
    static inline bool
    try_convert(const std::string& str, eva::fixed_vector<1, int>& vec) 
    {
        vec << std::stoi(str);
        return true;
    }
};

} } // end details and boost




// template<typename Target, typename Source>
// struct lexical_converter_impl
// {
//     typedef lexical_cast_stream_traits<Source, Target>  stream_trait;

//     typedef detail::lexical_istream_limited_src<
//         BOOST_DEDUCED_TYPENAME stream_trait::char_type,
//         BOOST_DEDUCED_TYPENAME stream_trait::traits,
//         stream_trait::requires_stringbuf,
//         stream_trait::len_t::value + 1
//         > i_interpreter_type;

//     typedef detail::lexical_ostream_limited_src<
//         BOOST_DEDUCED_TYPENAME stream_trait::char_type,
//         BOOST_DEDUCED_TYPENAME stream_trait::traits
//         > o_interpreter_type;

//     static inline bool try_convert(const Source& arg, Target& result) {
//         i_interpreter_type i_interpreter;

//         // Disabling ADL, by directly specifying operators.
//         if (!(i_interpreter.operator <<(arg)))
//             return false;

//         o_interpreter_type out(i_interpreter.cbegin(), i_interpreter.cend());

//         // Disabling ADL, by directly specifying operators.
//         if(!(out.operator >>(result)))
//             return false;

//         return true;
//     }
// };





// namespace boost { namespace detail {

// template<int N, typename T>
// struct lexical_cast_do_cast< eva::fixed_vector<N,T>, std::string >
// {     
//     static void  
//     do_assign(const double value, eva::fixed_vector<N,T>& vec, const size_t idx) 
//     {
//         vec(idx) = value;
//     }
    
//     static inline
//     eva::fixed_vector<N,T>
//     lexical_cast_impl(const std::string& str) 
//     {
//         namespace qi    = boost::spirit::qi;
//         namespace phx   = boost::phoenix;
//         // namespace ascii = boost::spirit::ascii;
        
//         // Alias for parser input iterator type
//         using iterator_type = std::string::const_iterator;
        
//         // Alias for a parser that parses a T
//         auto t_ = qi::real_parser<T>();
//         static_assert(std::is_floating_point<T>::value, 
//                        "Integral types are not supported yet.");

//         // The actual parser
//         qi::rule<iterator_type, eva::fixed_vector<N,T>(),
//                  qi::locals<size_t>, qi::blank_type> rule = 
//             qi::eps [ qi::_a = 0 ]
//             > '['
//             > qi::repeat(N-1) [
//                 t_[phx::bind(&lexical_cast_do_cast::do_assign, qi::_1, qi::_val, qi::_a++)]
//                 >> ',' 
//                 ]
//             > t_[phx::bind(&lexical_cast_do_cast::do_assign, qi::_1, qi::_val, qi::_a)]
//             > ']'
//             ;
        
//         rule.name("vector rule");
//         qi::on_error<qi::fail> (
//             rule,
//             std::cerr
//                 << phx::val("Error! Expecting ")
//                 << qi::_4                               // what failed?
//                 << phx::val(" here: \"")
//                 << phx::construct<std::string>(qi::_3, qi::_2)  // iterators to error-pos, end
//                 << phx::val("\"")
//                 << std::endl
//         );
                  
//         auto ret = eva::fixed_vector<N,T>();
        
//         auto iter(str.begin()), iend(str.end());
        
//         bool r = phrase_parse(iter, iend, rule, qi::blank, ret);
        
//         if (!r || iter != iend) {
//             std::cout << "-------------------------\n";
//             std::cout << "Parsing failed\n";
//             std::cout << str << std::endl;
//             std::cout << "-------------------------\n";
//         }
//         return ret;
//     }
// };


// template<>
// struct lexical_cast_do_cast< eva::fixed_vector<1, double>, std::string >
// { 
//     static inline
//     eva::fixed_vector<1, double>
//     lexical_cast_impl(const std::string& str) 
//     {
//         auto ret = eva::fixed_vector<1, double>();
//         ret << std::stod(str);
//         return ret;
//     }
// };
// template<>
// struct lexical_cast_do_cast< eva::fixed_vector<1, float>, std::string >
// { 
//     static inline
//     eva::fixed_vector<1, float>
//     lexical_cast_impl(const std::string& str) 
//     {
//         auto ret = eva::fixed_vector<1, float>();
//         ret << std::stof(str);
//         return ret;
//     }
// };
// template<>
// struct lexical_cast_do_cast< eva::fixed_vector<1, int>, std::string >
// { 
//     static inline
//     eva::fixed_vector<1, int>
//     lexical_cast_impl(const std::string& str) 
//     {
//         auto ret = eva::fixed_vector<1, int>();
//         ret << std::stoi(str);
//         return ret;
//     }
// };

// } } // end details and boost

# endif //__EVA_GRAPHVIZ__
