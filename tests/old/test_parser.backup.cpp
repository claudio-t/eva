# include <string>
# include <vector>
# include <array>

//~ # include <boost/array.hpp>
//~ # include <boost/spirit/include/qi.hpp>

# include <eigen3/Eigen/Dense>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
//~ #include <boost/bind.hpp>

//~ #define BOOST_RESULT_OF_USE_DECLTYPE
//~ #define BOOST_SPIRIT_USE_PHOENIX_V3

namespace qi    = boost::spirit::qi;
namespace phx   = boost::phoenix;
namespace ascii = boost::spirit::ascii;

typedef Eigen::Matrix<double, 3, 1> fixed_vector3;
//~ typedef std::array<double, 2> fixed_vector2;


void assign_value(const double value, fixed_vector3& vec, const size_t idx) 
{
    vec(idx) = value;
}


int main() {
    using It = std::string::const_iterator;

    std::string const src("(0.1, 1.2, 2.3)");
    auto f(src.begin()), l(src.end());

    //~ using data_t = std::array<double, 2>;
    //~ data_t dst {};
//~ 
        //~ qi::rule<It, data_t(), qi::locals<data_t::iterator>, qi::space_type> rule = 
            //~ qi::eps [ qi::_a = phx::begin(qi::_val) ]
            //~ >> qi::repeat(2) [
                    //~ qi::double_[ *qi::_a++ = qi::_1 ]
            //~ ];
    using data_t = fixed_vector3;
    data_t dst {};
    //~ auto assign_value2 = [](auto value, auto vec, auto idx) {vec(idx) = value;};
    
        qi::rule<It, data_t(), qi::locals<size_t>, qi::space_type> rule = 
            qi::eps [ qi::_a = 0u ]
            >> '('
            >> qi::repeat(2) [
                    qi::double_[phx::bind(&assign_value, qi::_1, qi::_val, qi::_a++)]
                    >> ',' 
            ]
            >> qi::double_[phx::bind(&assign_value, qi::_1, qi::_val, qi::_a)]
            >> ')'
        ;
            

    bool r = qi::phrase_parse(f, l, rule, qi::space, dst);

    if (r) {
        std::cout << "Parse succeeded\n";

        //~ for(auto&& el : dst) std::cout << el << " ";
        std:: cout << dst;
        std::cout << "\n";
    } else {
        std::cout << "Parse failed at '" << std::string(f,l) << "'\n";
    }
}

//~ namespace qi = boost::spirit::qi;
//~ namespace ascii = boost::spirit::ascii;
//~ 
//~ template <typename T, int N> using fixed_vector = Eigen::Matrix<T, N, 1>;
//~ 
//~ ///////////////////////////////////////////////////////////////////////////////
//~ // create a wrapper holding the fixed_vector and a current insertion point
//~ namespace client  { namespace detail {
    //~ 
//~ template <typename T>
//~ struct adapt_array;
//~ 
//~ template <typename T, int N>
//~ struct adapt_array<fixed_vector<T, N> >
//~ {
    //~ typedef fixed_vector<T, N> array_type;
    //~ 
    //~ typedef T value_type;
//~ 
    //~ adapt_array(array_type& arr)
      //~ : arr_(arr), current_(0) {}
//~ 
    //~ // expose a push_back function compatible with std containers
    //~ bool push_back(value_type const& val)
    //~ {
        //~ // if the array is full, we need to bail out
        //~ // returning false will fail the parse
        //~ if (current_ >= N) 
            //~ return false;
//~ 
        //~ arr_(current_++) = val;
        //~ return true;
    //~ }
//~ 
    //~ array_type& arr_;
    //~ std::size_t current_;
//~ };
//~ }
//~ 
//~ namespace result_of
//~ {
//~ template <typename T>
//~ struct adapt_array;
//~ 
//~ template <typename T, int N>
//~ struct adapt_array<fixed_vector<T, N> >
//~ {
    //~ typedef detail::adapt_array<fixed_vector<T, N> > type;
//~ };
//~ }
//~ 
//~ template <typename T, int N>
//~ inline detail::adapt_array<fixed_vector<T, N> > 
//~ adapt_array(fixed_vector<T, N>& arr) 
//~ {    
    //~ return detail::adapt_array<fixed_vector<T, N> >(arr);
//~ }
//~ 
//~ }
//~ 
//~ ///////////////////////////////////////////////////////////////////////////////
//~ // specialize Spirit's container specific customization points for our adaptor
//~ namespace boost { namespace spirit { namespace traits 
//~ {
    //~ template <typename T, int N>
     //~ struct is_container<client::detail::adapt_array<fixed_vector<T, N> > > 
       //~ : boost::mpl::true_
     //~ {};
//~ 
    //~ template <typename T, int N>
    //~ struct container_value<client::detail::adapt_array<fixed_vector<T, N> > >
    //~ {
        //~ typedef T type;     // value type of container
    //~ };
//~ 
    //~ template <typename T, int N>
    //~ struct push_back_container<
      //~ client::detail::adapt_array<fixed_vector<T, N> >, T>
    //~ {
        //~ static bool call(client::detail::adapt_array<fixed_vector<T, N> >& c
          //~ , T const& val)
        //~ {
            //~ return c.push_back(val);
        //~ }
    //~ };
//~ }}}
//~ 
//~ int main()
//~ {
    //~ typedef std::string::const_iterator iterator_type;
    //~ typedef fixed_vector<int, 2> array_type; 
    //~ typedef client::result_of::adapt_array<array_type>::type adapted_type; 
//~ 
    //~ array_type arr;
//~ 
    //~ std::string str = "1 2";
    //~ iterator_type iter = str.begin();
    //~ iterator_type end = str.end();
//~ 
    //~ qi::rule<iterator_type, adapted_type(), ascii::space_type> r = *qi::int_;
//~ 
    //~ adapted_type attr = client::adapt_array(arr);
    //~ bool result = qi::phrase_parse(iter, end, r, ascii::space, attr);
//~ 
    //~ if (result) 
        //~ std::cout << "Parsed: " << arr[0] << ", " << arr[1] << std::endl;
//~ 
    //~ return 0;
//~ }





//~ # include "parser.hpp"
//~ 
//~ int main(int argc, char * argv[]) {
    //~ 
    //~ std::cout << "/////////////////////////////////////////////////////////\n\n";
    //~ std::cout << "\t\tAn test parser for Spirit...\n\n";
    //~ std::cout << "/////////////////////////////////////////////////////////\n\n";
//~ 
    //~ std::cout << "Type [q or Q] to quit\n\n";
//~ 
    //~ using boost::spirit::ascii::space;
    //~ typedef std::string::const_iterator iterator_type;
//~ 
    //~ eva::element_parser<iterator_type> g; // Our grammar
    //~ std::string str;
    //~ while (getline(std::cin, str))
    //~ {
        //~ if (str.empty() || str[0] == 'q' || str[0] == 'Q')
            //~ break;
//~ 
        //~ eva::test_element el;
        //~ std::string::const_iterator iter = str.begin();
        //~ std::string::const_iterator end = str.end();
        //~ bool r = phrase_parse(iter, end, g, space, el);
//~ 
        //~ if (r && iter == end)
        //~ {
            //~ std::cout << boost::fusion::tuple_open('[');
            //~ std::cout << boost::fusion::tuple_close(']');
            //~ std::cout << boost::fusion::tuple_delimiter(", ");
//~ 
            //~ std::cout << "-------------------------\n";
            //~ std::cout << "Parsing succeeded\n";
            //~ std::cout << el.id <<  "\t" << el.coords[0] << " " << el.coords[1] << std::endl;
            //~ std::cout << "\n-------------------------\n";
        //~ }
        //~ else
        //~ {
            //~ std::cout << "-------------------------\n";
            //~ std::cout << "Parsing failed\n";
            //~ std::cout << "-------------------------\n";
        //~ }
    //~ }
//~ 
    //~ std::cout << "Bye... :-) \n\n";
    //~ return 0;
//~ }
