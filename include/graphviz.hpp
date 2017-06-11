# ifndef __EVA_GRAPHVIZ__
# define __EVA_GRAPHVIZ__

/**
 * @file graphviz.hpp
 * @brief Contains utilities to read and write structures in graphviz
 * format (.dot files)
 */

// self
# include "core.hpp"

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

/// Adds to the dynamic property map the properties required for
/// reading a given structure
template <typename S> 
void 
setup_element_properties(S& structure, boost::dynamic_properties& dyna_props);

/// Functor that selects the joint properties that have to be red
/// from the input file. Has to be specialized for each joint and
/// element type
template <typename Props>
struct properties_selector;


#ifdef __EVA_TRUSS__
/// Specialization for reading a truss joint 
template <int N>
struct properties_selector< truss_joint<N>>;

/// Specialization for reading a truss_element
template <int N>
struct properties_selector< truss_element<N> >;
#endif//__EVA_TRUSS__


#ifdef __EVA_FRAME__
/// Specialization for reading a 2D frame joint 
template <>
struct properties_selector< frame_joint<2> >;

/// Specialization for reading a 2D frame element
template <>
struct properties_selector< frame_element<2> >;
#endif//__EVA_FRAME__

#ifdef __EVA_THERMO__
/// Specialization for reading a 2D frame joint 
template <typename BaseKind>
struct properties_selector< thermo_joint<BaseKind> >;

/// Specialization for reading a 2D frame element
template <typename BaseKind>
struct properties_selector< thermo_element<BaseKind> >;
#endif//__EVA_THERMO__

//------------------------------------------ Write -----------------------------------------------//
template <typename S>
auto write_graphviz(const S& structure, const std::string & filename);

template <typename S>
auto make_joint_properties_writer(const S& structure);

template <typename S>
auto make_element_properties_writer(const S& structure);


template <typename S>
auto make_element_properties_writer(const S& structure);

/// Utility function for pretty printing a fixed vector
template <int N> std::string
to_string(
    const fixed_vector<N>& v,
    const std::string lx_delim = "[",
    const std::string rx_delim = "]");


/// Functor that selects the joint properties that have to be written
/// to the output file. Has to be specialized for each joint type
template <typename S, typename JointProperties>
class properties_writer;

#ifdef __EVA_TRUSS__
/// Specialization for writing a truss joint
template <typename S, int N>
class properties_writer< S, truss_joint<N> >;

/// Specialization for writing a truss joint
template <typename S, int N>
class properties_writer< S, truss_element<N> >;
#endif//__EVA_TRUSS__


#ifdef __EVA_FRAME__
/// Specialization for writing a 2D frame joint
template <typename S>
class properties_writer< S, frame_joint<2> >;

/// Specialization for writing a truss joint
template <typename S, int N>
class properties_writer< S, truss_element<N> >;
#endif//__EVA_FRAME__

#ifdef __EVA_THERMO__

template <typename S, typename BaseKind>
class properties_writer< S, thermo_joint<BaseKind> >;


#endif//__EVA_THERMO__


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


#ifdef __EVA_TRUSS__

template <int N>
struct properties_selector< truss_joint<N> > 
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
struct properties_selector< truss_element<N> > 
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
#endif//__EVA_TRUSS__


#ifdef __EVA_FRAME__
template <>
struct properties_selector< frame_joint<2> > 
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
struct properties_selector< frame_element<2> > 
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
#endif//__EVA_FRAME__


#ifdef __EVA_THERMO__

template <typename BaseKind>
struct properties_selector< thermo_joint<BaseKind> > 
{    
    template <typename S> 
    void 
    operator()(S& structure, boost::dynamic_properties& dps) 
    {
        using vert_prop_type = typename joint_of<S>::type;

        // Fill base kind props
        using base_kind_type = typename kind_of<S>::type::base_kind_t;
        using base_vert_prop = typename base_kind_type::joint_type;
        properties_selector<base_vert_prop>()(structure, dps);
        
        // Read coordinates and boundary conditions
        // dps.property( "coords", get(&vert_prop_type::coords, structure));
        dps.property(   "T_bc", get(&vert_prop_type::T_bc,   structure));
        dps.property("flux_bc", get(&vert_prop_type::flux_bc,structure));
    }
};

template <typename BaseKind>
struct properties_selector< thermo_element<BaseKind> > 
{    
    template <typename S> 
    void 
    operator()(S& structure, boost::dynamic_properties& dps)
    {
        using edge_prop_type = typename element_of<S>::type;

        // Fill base kind props
        using base_edge_prop = typename edge_prop_type::base_type;   
        properties_selector<base_edge_prop>()(structure, dps);        
        
        // Read thermal conductity coefficient and cross sectional area
        dps.property("k", get(&edge_prop_type::k, structure));        
        // dps.property("A", get(&edge_prop_type::A, structure));
    }
};
#endif//__EVA_THERMO__

//------------------------------------------ Write -----------------------------------------------//
template <typename S>
auto write_graphviz(const S& structure, const std::string & filename)
{
    std::ofstream ofs(filename);
    
    boost::write_graphviz(
        ofs, structure,
        make_joint_properties_writer(structure),
        make_element_properties_writer(structure));
}

template <typename S>
auto
make_joint_properties_writer(const S& s)
{
    return properties_writer<S, typename joint_of<S>::type>(s);
}

template <typename S>
auto
make_element_properties_writer(const S& s)
{
    return properties_writer<S, typename element_of<S>::type>(s);
}

#ifdef __EVA_TRUSS__
template <typename S, int N>
class properties_writer< S, truss_joint<N> >
{
public:
    // using props_t = typename joint_of<S>::type;
    // template <typename P>
    // using property_map_t = boost::vec_adj_list_vertex_property_map<
    //     S, const S *,
    //     P, const P &,
    //     P eva::joint_of<S>::type::*
    //     // P eva::truss_joint<N>::*
    //     >;

    /// Constructor that initializes member data
    properties_writer(const S & structure)
        : structure_(structure){}

    template <class Vertex>
    void write_impl(std::ostream & os, const Vertex & v) const
    {
        os << "coords = \"" << to_string(structure_[v].coords) << "\", " 
           << "load = \""   << to_string(structure_[v].load)   << "\", "
           << "bcs = \""    << to_string(structure_[v].bcs)    << "\"";
    }
    
    template <class Vertex>
    void operator()(std::ostream & os, const Vertex & v) const
    {
        os << '[';
        write_impl(os, v);
        os << ']';
    }
    
protected:
    const S & structure_;
};

template <typename S, int N>
class properties_writer< S, truss_element<N> >
{
public:
    /// Constructor that initializes member data
    properties_writer(const S & structure)
        : structure_(structure){}

    template <class Edge>
    void write_impl(std::ostream & os, const Edge & e) const
    {
        os << "A = \"" << structure_[e].A << "\", " 
           << "E = \"" << structure_[e].E << "\"";
        
    }
    
    template <class Edge>
    void operator()(std::ostream & os, const Edge & e) const
    {
        os << '[';
        write_impl(os, e);
        os << ']';
    }
    
protected:
    const S & structure_;
};
#endif//__EVA_TRUSS__


#ifdef __EVA_FRAME__
template <typename S>
class properties_writer< S, frame_joint<2> >
    : public properties_writer< S, truss_joint<2> >
{
public:
    using props_t = typename joint_of<S>::type;
    using base_t  = properties_writer< S, truss_joint<2> >;

    using base_t::structure_;

    /// Constructor that initializes member data
    properties_writer(const S & structure)
        : base_t(structure) {}

    /// Actual functor implementation
    template <class Vertex>
    void write_impl(std::ostream & os, const Vertex & v) const
    {
        base_t::write_impl(os, v);
        os << ", " << "torque = \"" << structure_[v].torque << "\"";
    }

    template <class Vertex>
    void operator()(std::ostream & os, const Vertex & v) const
    {
        os << '[';
        write_impl(os, v);
        os << ']';
    }
    
protected:
};

template <typename S>
class properties_writer< S, frame_element<2> >
    : public properties_writer< S, truss_element<2> >
{
public:
    using props_t = typename joint_of<S>::type;
    using base_t  = properties_writer< S, truss_element<2> >;

    using base_t::structure_;

    /// Constructor that initializes member data
    properties_writer(const S & structure)
        : base_t(structure) {}

    /// Actual functor implementation
    template <class Edge>
    void write_impl(std::ostream & os, const Edge & e) const
    {
        base_t::write_impl(os, e);
        os << ", " << "I = \"" << structure_[e].I << "\"";
    }

    template <class Edge>
    void operator()(std::ostream & os, const Edge & e) const
    {
        os << '[';
        write_impl(os, e);
        os << ']';
    }
    
protected:
};
#endif//__EVA_FRAME__


#ifdef __EVA_THERMO__
template <typename S, typename BaseKind>
class properties_writer< S, thermo_joint<BaseKind> >
    : public properties_writer< S, typename thermo_joint<BaseKind>::base_type >
{
public:
    using base_t  = properties_writer< S, typename thermo_joint<BaseKind>::base_type >;
    using base_t::structure_;

    /// Constructor that initializes member data
    properties_writer(const S & structure)
        : base_t(structure) {}

    /// Actual functor implementation
    template <class Vertex>
    void write_impl(std::ostream & os, const Vertex & v) const
    {
        base_t::write_impl(os, v);
        os << ", " << "T_bc = \"" << structure_[v].   T_bc << "\", "
           << "flux_bc = \"" << structure_[v].flux_bc << "\"";
    }

    template <class Vertex>
    void operator()(std::ostream& os, const Vertex& v) const
    {
        os << '[';
        write_impl(os, v);
        os << ']';
    }
    
protected:
};


template <typename S, typename BaseKind>
class properties_writer< S, thermo_element<BaseKind> >
    : public properties_writer< S, typename thermo_element<BaseKind>::base_type >
{
public:
    using base_t  = properties_writer< S, typename thermo_element<BaseKind>::base_type >;
    using base_t::structure_;

    /// Constructor that initializes member data
    properties_writer(const S & structure)
        : base_t(structure) {}

    /// Actual functor implementation
    template <class Edge>
    void write_impl(std::ostream & os, const Edge & e) const
    {
        base_t::write_impl(os, e);
        os << ", " << "k = \"" << structure_[e].k << "\"";
    }

    template <class Edge>
    void operator()(std::ostream & os, const Edge & e) const
    {
        os << '[';
        write_impl(os, e);
        os << ']';
    }
    
protected:
};
#endif//__EVA_THERMO__

//------------------------------------------ Helpers -----------------------------------------------//

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
