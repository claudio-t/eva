# ifndef __UTILS__
# define __UTILS__

namespace utils {

/////////// BRACE-INITIALIZABLE //////////////
/// @brief Aggregator class that allows brace initialization
/// @tparam Ts arbitrary number of base classes
/// @attention Each base class must support brace initialization
template <class... Ts>
struct brace_initializable : public Ts... {
    
    /// Default constructor
    brace_initializable() = default;
    
    /// Brace-initialization supporting constructor
    brace_initializable(const Ts&... ts) : Ts(ts)... {}
    
};

////////// STATIC FOR LOOP //////////////
template<int First, int Last>
struct static_for {
    
    template<class F>
    void operator()(const F& f) const {
        f(First);
        static_for<First+1, Last>()(f);
    }
};

template<int Last>
struct static_for<Last, Last> {
    
    template<class F>
    void operator()(const F& f) const {}
};


////////// FAST ARRAY NORM //////////////


/////////// SEQUENTIALLY PRINTABLE ////////////
//~ template <class... Ts>
//~ struct sequentially_printable {
    //~ using 
//~ }

} // end namespace utils

# endif //__UTILS__
