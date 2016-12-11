# ifndef __EVA_THERMO__
# define __EVA_THERMO__
/**
 * @file thermo.hpp
 * @brief Contains methods and classes to solve thermal problems on structures.
 */
// eva
# include "core.hpp"

namespace eva {


//-------------------------------------- Thermo tag types ----------------------------------------//
/// Thermo structure tag type
template <typename StructureKind>
struct thermo_kind;

//----------------------------------------- Properties -------------------------------------------//
/// Thermo structure joint (node) properties
template <typename StructureKind>
struct thermo_joint;

/// Thermo structure element (edge) properties
template <typename StructureKind>
struct thermo_element;

/// Thermo structure type
template <typename StructureKind>
using thermo_structure = generic_structure<
    thermo_joint,
    thermo_element,
    thermo_kind
    >;


//------------------------------------- Problem Assembling ---------------------------------------//
template <>
struct system_submatrices_assembler<dense_algebra_t, thermo_kind>;

template <>
struct system_submatrices_assembler<sparse_algebra_t, thermo_kind>;

template <>
struct known_terms_assembler<thermo_kind>

//------------------------------------- Results Assembling ---------------------------------------//

struct result<thermo_kind>;



//####################################### DEFINITIONS ##############################################

//-------------------------------------- Thermo tag types ----------------------------------------//
struct thermo_kind
{
    constexpr static int ndof = 1;
};


//----------------------------------------- Properties -------------------------------------------//
struct thermo_joint
{
    real bcs;  ///< Thermal Dirichlet BC (fixed temperature)    
    real flux; ///< Thermal flux (Neumann condition)
};


struct thermo_element
{
    real k; // Thermal conductivity [W/(m*K)]
};


//------------------------------------- Problem Assembling ---------------------------------------//
template <>
struct system_submatrices_assembler<dense_algebra_t, thermo_kind>
{
    template <typename Structure>
    std::array<dense_matrix, 3>
    operator()(const Structure& s,
               const std::vector<index_t>& dofmap,
               const size_t n_f, const size_t n_b)
    {
        // Allocate & init submatrices
        dense_matrix K_ff = dense_matrix::Zero(n_f, n_f);
        dense_matrix K_fb = dense_matrix::Zero(n_f, n_b);
        dense_matrix K_bb = dense_matrix::Zero(n_b, n_b);

        // Loop over all vertices
        for (auto&& v : make_iterator_range(vertices(s)))
        {
            // Get current vertex mapped value
            auto avv = dofmap[av];

            // Init diagonal value
            auto dval = real(0);
            
            // For each of the adjacent nodes
            for (auto&& av : make_iterator_range(adjacent_vertices(v, s)))
            {
                // Compute element length
                auto l = (s[av].coords - s[v]).norm();

                // Compute element thermal resistance
                const auto& el = edge(v, av, s);
                auto R = l / (el.A * el.k);

                // Get adj vertex mapped value
                auto avv = dofmap[av];

                // Assign value
                switch (vv < n_f)
                {
                    case 0 : K_bb(vv-n_f, avv-n_f) = 1/R; break;
                        
                    case 1 : K_fb(vv, avv-n_f)     = 1/R; break;
                        
                    case 2 : /* Does nothing (no K_bf) */ break;
                        
                    case 3 : K_ff(vv, avv)         = 1/R; break;
                }                

                // Update diagonal value
                dval += 1/R;            
            }

            // Assign diagonal value
            switch (vv < n_f)
            {
                case 0 : K_bb(vv-n_f, vv-n_f) = -dval; break;
                        
                case 1 : K_fb(vv, vv-n_f)     = -dval; break;
                        
                case 2 : /* Does nothing (no K_bf) */  break;
                    
                case 3 : K_ff(vv, vv)         = -dval; break;
            } 

        }
# pragma clang diagnostic push
# pragma clang diagnostic ignored "-Wmissing-braces"
        return { std::move(K_ff), std::move(K_fb), std::move(K_bb) };
# pragma clang diagnostic pop  
    }    
};

template <>
struct system_submatrices_assembler<sparse_algebra_t, thermo_kind>
{
    template <typename Structure>
    std::array<sparse_matrix, 3>
    operator()(const Structure& s,
               const std::vector<index_t>& dofmap,
               const size_t n_f, const size_t n_b)
    {
        auto dms = system_submatrices_assembler<
            dense_algebra_t, thermo_kind>()(s, dofmap, n_f, n_b);
        
        return {sparse_matrix(dms[0].sparseView()),
                sparse_matrix(dms[1].sparseView()),
                sparse_matrix(dms[2].sparseView())}
    }
}
    


template <>
struct known_terms_assembler<thermo_kind>
{
    template <typename S>
    std::array<dense_vector, 2>
    operator()(const S& s, const size_t n_f, const size_t n_b) 
    {
        // Init return var
# pragma clang diagnostic push
# pragma clang diagnostic ignored "-Wmissing-braces"
        auto ret = std::array<dense_vector, 2> {dense_vector(n_f), dense_vector(n_b)};
# pragma clang diagnostic pop
        
        auto& f_f = std::get<0>(ret);   // Fixed fluxes vector (Free DOF only)
        auto& u_b = std::get<1>(ret);   // Temperatures vector (BC DOF only)
        
        size_t pos_f = 0u;     // f_f incremental iterator
        size_t pos_u = n_b;    // u_b decremental iterator
        
        for (auto v : boost::make_iterator_range(vertices(s)))
        {
            // Get bcs & loads
            const auto bc_d = s[v].bc_d;
            const auto flux = s[v].flux;

            // If free DOF => write to f & post-increase position
            if (std::isnan(bc_d))
                f_f(pos_f++) = load(i);
            // If BC DOF => write to u & pre-decrease position
            else 
                u_b(--pos_u) = bcs(i);
        }
        
        return ret;
    }
};


//------------------------------------- Results Assembling ---------------------------------------//

struct result<thermo_kind> 
{
    real T;   ///< Temperature
    real flux ///< Thermal flux

    result(const fixed_vector<ndof>& t,
           const fixed_vector<ndof>& f,
           const thermo_joint&)
        : T    (t[0])
        , flux (f[0])
    { }
};


} // end namespaces

# endif//__EVA_THERMO__

