# ifndef __EVA_SOLVER_PARAMETERS__
# define __EVA_SOLVER_PARAMETERS__

/// Generic interface
template <typename... Params>
struct solver_parameters;


/// Specialization for handling
template <typename T, typename... Ts>
struct solver_parameters<T, Ts...> : solver_parameters<Ts...>
{
    using std::conditional;
    using std::is_base_of;
    
    using solver_t  = typename conditional<is_base_of<structure_tag, T>::value>::type;
    using algebra_t = typename conditional<is_base_of<algebra_tag, T>::value>::type;
    
};


# endif //__EVA_SOLVER_PARAMETERS__
