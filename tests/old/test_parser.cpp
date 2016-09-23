# include "parser.hpp"
# include <string>





int main(int argc, char * argv[]) {
    
    namespace qi    = boost::spirit::qi;
    namespace phx   = boost::phoenix;
    namespace ascii = boost::spirit::ascii;
    
    std::cout << "/////////////////////////////////////////////////////////\n\n";
    std::cout << "\t\tAn test parser for Spirit...\n\n";
    std::cout << "/////////////////////////////////////////////////////////\n\n";

    std::cout << "Type [q or Q] to quit\n\n";

    typedef std::string::const_iterator iterator_type;

    eva::element_parser<iterator_type> g; // Our grammar
    //~ auto g = eva::make_array_rule<iterator_type>('(', ')');
    std::string str;
    while (getline(std::cin, str)) {
        
        if (str.empty() || str[0] == 'q' || str[0] == 'Q')
            break;

        eva::test_element el;
        //~ eva::vector3 el;
        std::string::const_iterator iter = str.begin();
        std::string::const_iterator end = str.end();
        bool r = phrase_parse(iter, end, g, ascii::space, el);

        if (r && iter == end) {
            //~ std::cout << boost::fusion::tuple_open('[');
            //~ std::cout << boost::fusion::tuple_close(']');
            //~ std::cout << boost::fusion::tuple_delimiter(", ");

            std::cout << "-------------------------\n";
            std::cout << "Parsing succeeded\n";
            //~ std::cout << el;
            //~ std::cout << el.id <<  "\t" << el.coords[0] << " " << el.coords[1] << std::endl;
            std::cout << "\n-------------------------\n";
        }
        else {
            std::cout << "-------------------------\n";
            std::cout << "Parsing failed\n";
            std::cout << "-------------------------\n";
        }
    }

    std::cout << "Bye... :-) \n\n";
    
    return 0;
}

