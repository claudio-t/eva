# ifndef __EVA_PARSER__
# define __EVA_PARSER__

#define BOOST_SPIRIT_DEBUG

# include <eigen3/Eigen/Dense>
# include <boost/spirit/include/qi.hpp>
# include <boost/spirit/include/phoenix.hpp>
# include <boost/fusion/include/adapt_struct.hpp>

#include <boost/fusion/include/tuple.hpp>
#include <boost/fusion/adapted/boost_array.hpp>

# include <iostream>
# include <string>
# include <array>
# include <vector>

namespace eva { 
    
typedef Eigen::Matrix<double, 3, 1> vector3;
//~ typedef boost::array<double, 3> vector3;
//~ typedef std::vector<double> vector3;
struct test_element 
{    
    size_t  id;
    //~ double dummy;
    vector3 coords;
    //~ vector3 bcs;
    //~ vector3 load;
}; }
//~ namespace boost { namespace spirit { namespace traits {
    //~ template <typename T, size_t N>
        //~ struct is_container<boost::array<T, N>, void> : mpl::false_ { };
//~ } } }
// MUST BE AT GLOBAL SCOPE!!!!!!!!!!!
BOOST_FUSION_ADAPT_STRUCT(
    eva::test_element,
    (size_t, id)
    //~ (double, dummy)
    (eva::vector3, coords)
    //~ (eva::vector3, bcs)
    //~ (eva::vector3, load)
)


namespace eva {
    
namespace qi    = boost::spirit::qi;
namespace phx   = boost::phoenix;
namespace ascii = boost::spirit::ascii;

namespace details { void assign_value(const double value, vector3& vec, const size_t idx) 
{
    vec(idx) = value;
} }

// Parse a vector
//~ template <typename Iterator>
//~ bool parse_vector(Iterator first, Iterator last, std::vector<double>& v) {
    //~ 
    //~ using qi::double_;
    //~ using qi::phrase_parse;
    //~ using qi::_1;
    //~ using ascii::space;
    //~ using phoenix::push_back;
    //~ 
    //~ auto parser = '[' >> double_[push_back(phoenix::ref(v), _1)] % ',' >> ']';
    //~ bool r = phrase_parse(first, last, parser, space);
    //~ 
    //~ // Or more simply
    //~ auto parser = '[' >> double_ % ',' >> ']';
    //~ bool r = phrase_parse(first, last, parser, space, v);
            //~ 
    //~ if (first != last) // fail if we did not get a full match
        //~ return false;
    //~ 
    //~ return r;
//~ }
//~ template <typename It>
//~ qi::rule<It, vector3(), qi::locals<size_t>, ascii::space_type>
//~ make_array_rule(const char delim_lx, const char delim_rx)
//~ {
    //~ qi::rule<It, vector3(), qi::locals<size_t>, ascii::space_type> rule = 
        //~ qi::eps [ qi::_a = 0u ]
        //~ >> delim_lx
        //~ >> qi::repeat(2) [
                //~ qi::double_[phx::bind(&details::assign_value, qi::_1, qi::_val, qi::_a)]
                //~ >> ',' 
        //~ ]
        //~ >> qi::double_[phx::bind(&details::assign_value, qi::_1, qi::_val, qi::_a)]
        //~ >> delim_rx
    //~ ;
    //~ 
    //~ return rule;
//~ }
template <typename It>
qi::rule<It, vector3(), qi::locals<size_t>, ascii::space_type>
make_array_rule(const char delim_lx, const char delim_rx)
{
    qi::rule<It, vector3(), qi::locals<size_t>, ascii::space_type> rule = 
        qi::eps [ qi::_a = 0u ]
        >> delim_lx
        >> qi::repeat(2) [
                qi::double_[phx::bind(&details::assign_value, qi::_1, qi::_val, qi::_a++)]
                >> ',' 
        ]
        >> qi::double_[phx::bind(&details::assign_value, qi::_1, qi::_val, qi::_a)]
        >> delim_rx
    ;
    
    return rule;
}


template <typename Iterator>
struct element_parser : qi::grammar<Iterator, test_element(), ascii::space_type> 
{
    element_parser() : element_parser::base_type(element_rule) 
    {
        //~ using size_t_parser = qi::uint_parser<size_t>;
        //~ auto size_t_ = size_t_parser; 
        auto el_id_rule = qi::int_;
        
        auto coords_rule = make_array_rule<Iterator>('(', ')');
        //~ auto vector_rule = make_array_rule<Iterator>('[', ']');
        //~ auto coords_rule = '(' >> qi::double_ >> qi::double_ >> qi::double_ >> ')';
                
        element_rule %= el_id_rule              // ID rule
                     >> ':'
                     //~ >> ('(' >> qi::double_ >> ',' >> qi::double_ >> ',' >> qi::double_ >> ')')
                     >> coords_rule             // Coordinates
                     //~ >> ( ',' >> vector_rule)   // Optional bcs
                     //~ >> ( ',' >> vector_rule)   // Optional load
                     >> ';'
        ;
        
        //~ using qi::on_error;
        //~ using qi::fail;
        //~ using qi::debug;
        //~ using phx::construct;
        //~ using phx::val;
        //~ coords_rule.name("coords");
        //~ element_rule.name("element");
        //~ on_error<fail>
        //~ (
            //~ element_rule
          //~ , std::cout
                //~ << val("Error! Expecting ")
                //~ << qi::_4                               // what failed?
                //~ << val(" here: \"")
                //~ << construct<std::string>(qi::_3, qi::_2)   // iterators to error-pos, end
                //~ << val("\"")
                //~ << std::endl
        //~ );
        //~ debug(coords_rule);
        //~ debug(element_rule);
    }
    
private:

    //~ qi::rule<Iterator, size_t(), ascii::space_type>       el_id_rule;
    //~ qi::rule<Iterator, vector3(), qi::space_type>      coords_rule;
    //~ qi::rule<Iterator, vector3(), ascii::space_type>      vector_rule;
    qi::rule<Iterator, test_element(), ascii::space_type> element_rule;
    
};


} // end namespace eva


# endif //__EVA_PARSER__
