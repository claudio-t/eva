# ifndef __EVA_NUM_ARRAY__
# define __EVA_NUM_ARRAY__

# include <array>
# include <iostream>

namespace eva {
    
//######################################## Declarations ############################################

//------------------------------------------ Classes ---------------------------------------------//
/// Numerical fixed size array type
template <typename T, size_t N, template<typename> class Initializer>
class num_array;

/// Struct used for initialize a num_array to zero
template <typename T> struct zero_init;

/// Struct used for initialize a num_array to quiet NaN
template <typename T> struct nan_init;

/// Struct used for initialize a num_array to infinity
template <typename T> struct inf_init;

//----------------------------------------- Functions --------------------------------------------//
/// Creates an std::array filled with a given value (can be evaluated at compile-time)
template <size_t N, typename T> constexpr std::array<T, N> make_array(T val);

// FIXME!!! INCOMPLETE AND TOTALLY UNEFFICIENT!!!
/// Makes a num_array istream-able
template <typename T, size_t N, template<typename> class I>
std::istream& operator>>(std::istream& is, num_array<T, N, I>& arr);

/// Makes a num_array ostream-able
template <typename T, size_t N, template<typename> class I>
std::ostream& operator<<(std::ostream& os, const num_array<T, N, I>& arr);


//######################################## Definitions #############################################

//------------------------------------------ Classes ---------------------------------------------//

template <typename T, size_t N, template<typename> class I = zero_init >
class num_array {
public:
    // Member Types //
    using value_type      = T;                  ///< Value type
    using iterator        = T*;                 ///< Iterator type
    using const_iterator  = const T*;           ///< Const iterator type
    using reference       = value_type&;        ///< Iterator type
    using const_reference = const value_type&;  ///< Const iterator type
    
    // Constructors //
    /// Default constructor that initializes all values to the default provided.
    constexpr num_array() : vals_(make_array<N>(I<T>::value)) {}
    
    /// Constructor that allows the usage of list initialization
    template <typename... Ts> constexpr num_array(Ts... ts) : vals_{ts...} {}
    
    // Iterators //
    /// Returns an iterator pointing to the beginning of the num_array
    iterator begin() { return vals_.begin(); }
    /// Returns an iterator pointing to the end of the num_array (:= invalid position)
    iterator end()   { return vals_.end(); }
    /// Returns a constant iterator pointing to the beginning of the num_array
    constexpr const_iterator cbegin() const { return vals_.begin(); }
    /// Returns a constant iterator pointing to the end of the num_array (:= invalid position)
    constexpr const_iterator cend()   const { return vals_.end(); }
    
    // Data access //
    /// Provides random access to the stored values
                     reference operator[](size_t idx)       { return vals_[idx]; }
    /// Provides random (constant/RO) access to the stored values
    constexpr const value_type operator[](size_t idx) const { return vals_[idx]; }
    
    /// Returns the array size
    constexpr size_t size() const { return vals_.size(); }
    
private:
    
    /// The actual array content
    std::array<T,N> vals_;
    
    // Type check: T must represent a number //
    static_assert(std::is_integral<T>() || std::is_floating_point<T>(),
                  "The value type T in num_array must represent a number");
};
    
template <typename T> struct zero_init {
    constexpr static T value = 0; 
};

template <typename T> struct nan_init {
    constexpr static T value = std::numeric_limits<T>::quiet_NaN(); 
};

template <typename T> struct inf_init {
    constexpr static T value = std::numeric_limits<T>::infinity(); 
};



//----------------------------------------- Functions --------------------------------------------//

namespace details {
template <typename T, size_t... Idxs>
constexpr std::array<T, sizeof...(Idxs)>
make_array_helper(T val, std::index_sequence<Idxs...>) {
// clang disable missing braces
# pragma clang diagnostic push
# pragma clang diagnostic ignored "-Wmissing-braces"
// gcc disable unused value
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wunused-value"
    return { (Idxs, val)... };
# pragma GCC diagnostic pop
# pragma clang diagnostic pop
} }

template <size_t N, typename T>
constexpr std::array<T, N> make_array(T val) {
    return details::make_array_helper(val, std::make_index_sequence<N>());
}
// Does not work with gcc < 5.0
//~ template <size_t N, typename T>
//~ constexpr std::array<T, N> make_array(T val) {
    //~ 
    //~ auto arr = std::array<T,N>();
    //~ for (const auto& el : arr) el = val;
    //~ return arr;
//~ }

  

// FIXME!!! INCOMPLETE AND TOTALLY UNEFFICIENT!!!
template <typename T, size_t N, template<typename> class I>
std::istream& operator>>(std::istream& is, num_array<T, N, I>& arr) {
    
    auto line = std::string ();                 // Temporary
    std::getline(is, line);                     // Read line from istream
    std::istringstream iss (std::move(line));   // Move construct istringstream
    
    for (auto& el : arr) {
        std::string str;
        iss >> str;
        // FIXME: what for size_t, int etc..??
        if (std::is_same<T, double>::value) el = std::stod(str);
        if (std::is_same<T,  float>::value) el = std::stof(str);
        if (std::is_same<T,    int>::value) el = std::stoi(str);
    }
    return is;
}

template <typename T, size_t N, template<typename> class I>
std::ostream& operator<<(std::ostream& os, const num_array<T, N, I>& arr) {
    // Copy directly to ostream
    std::copy(arr.cbegin(), arr.cend(), std::ostream_iterator<T>(os, " "));
    return os;
}

} // end namespace eva


# endif //__EVA_NUM_ARRAY__
