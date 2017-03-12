# ifndef __EVA_THERMO__
# define __EVA_THERMO__
/**
 * @file thermo.hpp
 * @brief Contains methods and classes to solve thermal problems on structures.
 */
// eva
# include "core.hpp"
#include <iostream>
#include <iomanip>

namespace eva {


//-------------------------------------- Thermo tag types ----------------------------------------//
/// Thermo structure tag type
template <typename BaseKind>
struct thermo_kind;

//----------------------------------------- Properties -------------------------------------------//
/// Thermo structure joint (node) properties
template <typename BaseKind>
struct thermo_joint;

/// Thermo structure element (edge) properties
// template <typename StructureKind>
template <typename BaseKind>
struct thermo_element;

/// Thermo structure type
template <typename BaseKind>
using thermo_structure = generic_structure<
    thermo_joint<BaseKind>,
    thermo_element<BaseKind>,
    thermo_kind<BaseKind>
    >;


//------------------------------------- Problem Assembling ---------------------------------------//
template <typename BaseKind>
struct system_submatrices_assembler< dense_algebra_t, thermo_kind<BaseKind> >;

template <typename BaseKind>
struct system_submatrices_assembler< sparse_algebra_t, thermo_kind<BaseKind> >;

template <typename BaseKind>
struct known_terms_assembler< thermo_kind<BaseKind> >;

//------------------------------------- Results Assembling ---------------------------------------//

template <typename BaseKind>
struct result< thermo_kind<BaseKind> >;

//---------------------------------------- Post-process ------------------------------------------//

template <typename Structure, typename Result>
real compute_heat_flow(
    const typename Structure::edge_descriptor,
    const Structure & structure,
    const std::vector<Result> & results);


//####################################### DEFINITIONS ##############################################

//-------------------------------------- Thermo tag types ----------------------------------------//
template <typename BaseKind>
struct thermo_kind : public BaseKind
{
    using base_kind_t = BaseKind;
    constexpr static int ndof = 1;
    
    using default_dense_solver_t  = Eigen::PartialPivLU<dense_matrix>;
    using default_sparse_solver_t = Eigen::BiCGSTAB<sparse_matrix>;
};


//----------------------------------------- Properties -------------------------------------------//
// template <typename BaseKind>
template <typename BaseKind>
struct thermo_joint : public BaseKind::joint_type
{
    using base_type = typename BaseKind::joint_type;
    
    real T_bc;  ///< Thermal Dirichlet BC (fixed temperature)    
    real flux_bc; ///< Thermal flux (Neumann condition)

    thermo_joint()
        :
        // coords(fixed_vector<N>::Zero()),
        T_bc(std::numeric_limits<real>::quiet_NaN()),
        flux_bc(0) {}
          
};

template <typename BaseKind>
struct thermo_element : public BaseKind::element_type
{
    using base_type = typename BaseKind::element_type;
    real k; ///< Thermal conductivity [W/(m*K)]
};


//------------------------------------- Problem Assembling ---------------------------------------//
template <typename BaseKind>
struct system_submatrices_assembler< dense_algebra_t, thermo_kind<BaseKind> >
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

        for (auto e : make_iterator_range(edges(s)))
        {
            auto src = source(e, s);
            auto trg = target(e, s);

            // Compute element thermal resistance
            real l = (s[src].coords - s[trg].coords).norm();
            real R = l / (s[e].A * s[e].k);

            // Add contributes to the stiffness matrix
            auto src_m = dofmap[src];
            auto trg_m = dofmap[trg];

            if (src_m < n_f)
            {
                // Update Diagonal
                K_ff(src_m, src_m) += 1/R;

                // Update extra-diagonal
                if (trg_m < n_f)
                    K_ff(src_m, trg_m) -= 1/R;
                else
                    K_fb(src_m, trg_m - n_f) -=1/R;
                
            }
            else
            {
                K_bb(src_m - n_f, src_m - n_f) = 1/R;
                
                // if (trg_m >= n_f)
                    // K_fb(src_m - n_f, trg_m - n_f) -= 1/R;
                // Do nothing on K_bf
            }
            
            if (trg_m < n_f)
            {
                K_ff(trg_m, trg_m) += 1/R;

                if (src_m < n_f)
                    K_ff(trg_m, src_m) -= 1/R;
                else
                    K_fb(trg_m, src_m - n_f) -=1/R;
            }
            else
            {
                K_bb(trg_m - n_f, trg_m - n_f) = 1/R;
                
                // if (src_m < n_f)
                //     K_fb(trg_m - n_f, src_m) += 1/R;
            }    
        }
        

        // // Loop over all vertices
        // for (auto&& v : boost::make_iterator_range(vertices(s)))
        // {
        //     // Get current vertex mapped value
        //     auto vv = dofmap[v];

        //     // Init diagonal value
        //     auto dval = real(0);
            
        //     // For each of the adjacent nodes
        //     for (auto&& av : make_iterator_range(adjacent_vertices(v, s)))
        //     {
        //         // Compute element length
        //         real l = (s[av].coords - s[v].coords).norm();

        //         // Compute element thermal resistance
        //         const auto el = edge(v, av, s).first;
        //         real R = l / (s[el].A * s[el].k);

        //         // Get adj vertex mapped value
        //         auto avv = dofmap[av];

        //         // Assign value: if edge (v, av) is first encountered
        //         // flux is outflow (negative), otherwise is inflow
        //         // (positive).
        //         // Edge is first encountered if v < av or if av has
        //         // dirichlet conditions.
        //         bool outflow = (v < av) || (avv >= n_f);
        //         switch((vv < n_f) + 2*(avv < n_f))
        //         {
        //             case 0 : K_bb(vv-n_f, avv-n_f) = (outflow ? -1/R : 1/R); break;
                        
        //             case 1 : K_fb(vv, avv-n_f)     = (outflow ? -1/R : 1/R); break;
                        
        //             case 2 : /* Does nothing (no K_bf) */ break;
                        
        //             case 3 : K_ff(vv, avv)         = (outflow ? -1/R : 1/R); break;
        //         }                

        //         // Update diagonal value
        //         dval += 1. / (outflow ? R : -R);//(v < av ? 1. : -1.) * 1/R;            
        //     }

        //     // Assign diagonal value
        //     if (vv < n_f)  K_ff(vv, vv)         = dval;
        //     else           K_bb(vv-n_f, vv-n_f) = dval;
            
        // }
        std::cout << "DOFS:\n"
                  << "n_f = " << n_f << "\n"
                  << "n_b = " << n_b << "\n";

        std::cout << "Map = ";
        for (const auto & el : dofmap) std::cout << el << " ";
        std::cout << "\n\n";
        
        std::cout << std::setprecision(16);
        
        std::cout << "K_ff:\n" << K_ff << "\n\n"
                  << "K_fb:\n" << K_fb << "\n\n"
                  << "K_bb:\n" << K_bb << "\n\n";
        
            
# pragma clang diagnostic push
# pragma clang diagnostic ignored "-Wmissing-braces"
        return { std::move(K_ff), std::move(K_fb), std::move(K_bb) };
# pragma clang diagnostic pop  
    }    
};

template <typename BaseKind>
struct system_submatrices_assembler< sparse_algebra_t, thermo_kind<BaseKind> >
{
    template <typename Structure>
    std::array<sparse_matrix, 3>
    operator()(const Structure& s,
               const std::vector<index_t>& dofmap,
               const size_t n_f, const size_t n_b)
    {
        auto dms = system_submatrices_assembler<
            dense_algebra_t, thermo_kind<BaseKind> >()(s, dofmap, n_f, n_b);
        
        return {sparse_matrix(dms[0].sparseView()),
                sparse_matrix(dms[1].sparseView()),
                sparse_matrix(dms[2].sparseView())};
    }
};
    


template <typename BaseKind>
struct known_terms_assembler< thermo_kind<BaseKind> >
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
        
        auto & f_f = std::get<0>(ret);   // Fixed fluxes vector (Free DOF only)
        auto & u_b = std::get<1>(ret);   // Temperatures vector (BC DOF only)
        
        size_t pos_f = 0u;     // f_f incremental iterator
        size_t pos_u = n_b;    // u_b decremental iterator
        
        for (auto v : boost::make_iterator_range(vertices(s)))
        {
            // Get bcs & loads
            const auto T_bc = s[v].T_bc;
            const auto flux_bc = s[v].flux_bc;

            // If free DOF => write to f & post-increase position
            if (std::isnan(T_bc))
                f_f(pos_f++) = flux_bc;
            // If BC DOF => write to u & pre-decrease position
            else 
                u_b(--pos_u) = T_bc;
        }

        // std::cout << "f_f:\n" << f_f << "\n\n"
        //           << "u_b:\n" << u_b << "\n\n";
        
        return ret;
    }
};



//------------------------------------- Results Assembling ---------------------------------------//

template <typename BaseKind>
struct result< thermo_kind<BaseKind> > 
{
    real T;    ///< Temperature
    real flux; ///< Thermal flux

    result(const fixed_vector<1>& t,
           const fixed_vector<1>& f,
           const thermo_joint<BaseKind>& data)
        : T (t[0]), flux (f[0])
    {}
};


//---------------------------------------- Post-process ------------------------------------------//
template <typename Structure, typename Result>
real compute_heat_flow(
    const typename Structure::edge_descriptor e,
    const Structure & s,
    const std::vector<Result> & results)
{
    auto src = source(e, s);
    auto trg = target(e, s);
    
    // Retrieve nodal temperatures
    auto T_src = results[src].T;
    auto T_trg = results[trg].T;

    // Compute element thermal resistance
    real l = (s[src].coords - s[trg].coords).norm();
    real R = l / (s[e].A * s[e].k);

    // std::cout << std::setprecision(16);
    // std::cout << "(" << src << ", " << trg << ") --> "
    //           << (T_src - T_trg) / R << "\n";
    // std::cout << "L = " << l << "\t" << "R = " << R << "\n";
    // std::cout << "T_src = " << T_src << "\t" << "T_trg = " << T_trg << "\n\n";

    // Compute flux
    return (T_src - T_trg) / R;
}


} // end namespaces

# endif//__EVA_THERMO__

