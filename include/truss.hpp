# ifndef __EVA_TRUSS__
# define __EVA_TRUSS__

/**
 * @file truss.hpp
 * @brief Contains methods and classes to solve truss structures.
 */
// eva
# include "core_sa.hpp"


namespace eva {
    
//####################################### DECLARATIONS #############################################

//----------------------------------------- Properties -------------------------------------------//
/// Truss tag type
template <int N> struct truss_kind;

/// Truss joint (node) properties
template <int N> struct truss_joint;

/// Truss element (edge) properties
template <int N> struct truss_element;


//--------------------------------------- Truss  types -------------------------------------------//
/// 2-dimensional truss structure type
using truss2d = generic_structure< truss_joint<2>, truss_element<2>, truss_kind<2> >;
    
/// 3-dimensional truss structure type
using truss3d = generic_structure< truss_joint<3>, truss_element<3>, truss_kind<3> >;


//-------------------------------------- Post-processing -----------------------------------------//
// template <int N> 
// struct internal_forces_getter< truss_kind<N> >;


//------------------------------------- Problem Assembling ---------------------------------------//

/// Specializes system_submatrices_assembler functor for both 2D and
/// 3D trusses.
/// This specialization simply calls the stiffness_submatrices_assembler
/// functor, which -in turn- automatically selects the proper
/// specialization of the element_matrix_assembler functor.
template <typename A, int N>
struct system_submatrices_assembler<A, truss_kind<N> >;

/// Specializes known_terms_assembler functor for both 2D and 3D trusses
template <int N>
struct known_terms_assembler< truss_kind<N> >;

/// Specializes element_matrix_assembler functor for 2D trusses
template <>
struct element_matrix_assembler< truss_kind<2> >;

/// Specializes element_matrix_assembler functor for 3D trusses
template <>
struct element_matrix_assembler< truss_kind<3> >;


//------------------------------------- Results Assembling ---------------------------------------//
template <int N>
struct result< truss_kind<N> >;


template <int N>
struct result< truss_kind<N> >
{
    using kind_type = truss_kind<N>;
    
    constexpr static int ndof = truss_kind<N>::ndof;
    constexpr static int sdim = truss_kind<N>::sdim;
    
    fixed_vector<sdim> displacement;
    fixed_vector<sdim> reaction;

    result(const fixed_vector<ndof>& u,
           const fixed_vector<ndof>& f,
           const truss_joint<N>& p)
        : displacement(u)
        , reaction    (f)
    {
        // Keep result only if there is no load applied
        if (p.load != fixed_vector<sdim>::Zero())
            reaction = fixed_vector<sdim>::Zero();
    }
};

//####################################### DEFINITIONS ##############################################
  
//----------------------------------------- Properties -------------------------------------------//
template <int N> struct truss_kind 
{     
    constexpr static int ndof = N; 
    constexpr static int sdim = N;

    using joint_type   = truss_joint<N>;
    using element_type = truss_element<N>;
    using result_type  = result< truss_kind<N> >;
    
    using default_dense_solver_t  = Eigen::LDLT<dense_matrix>;
    using default_sparse_solver_t = Eigen::ConjugateGradient<sparse_matrix>;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {/* Do nothing! */}    
};


template <int N> struct truss_joint 
{
    fixed_vector<N> coords; ///< Joint coordinates
    fixed_vector<N> load;   ///< Applied load [N]
    fixed_vector<N> bcs;    ///< Boundary conditions
    
    /// Default constructor. Initializes
    /// coords to zero, load to zero and bcs to nan
    truss_joint()
        : coords(fixed_vector<N>::Zero()) 
        , load  (fixed_vector<N>::Zero())
        , bcs   (fixed_vector<N>::Constant(std::numeric_limits<real>::quiet_NaN()))
        {}

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & coords & load & bcs;
    }        
};


template <int N> struct truss_element 
{    
    real E; ///< Young's modulus      [GPa]
    
    real A; ///< Cross sectional area [m^2]

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & E & A;
    }
};



//------------------------------------- Problem Assembling ---------------------------------------//


template <typename A, int N>
struct system_submatrices_assembler<A, truss_kind<N> >
{
    template <typename S>
    auto
    operator()(const S& s,
               const std::vector<index_t>& dofmap,
               const size_t n_f, const size_t n_b)
    {   
        return stiffness_submatrices_assembler<
            truss_kind<N>, A
            >()(s, dofmap, n_f, n_b);
    }    
};
    

template <int N>
struct known_terms_assembler< truss_kind<N> > 
{
    template <typename S>
    std::array<dense_vector, 2>
    operator()(const S& s, const size_t n_f, const size_t n_b) 
    {
        // Aux vars
        const static int ndof = truss_kind<N>::ndof;   // #DOF per node (2D => 3, 3D => 6) 

        // Init return var
# pragma clang diagnostic push
# pragma clang diagnostic ignored "-Wmissing-braces"
        auto ret = std::array<dense_vector, 2> {dense_vector(n_f), dense_vector(n_b)};
# pragma clang diagnostic pop
        
        auto& f_f = std::get<0>(ret);   // Force vector (Free DOF only)
        auto& u_b = std::get<1>(ret);   // Displacement vector (BC DOF only)
        
        size_t pos_f = 0u;     // f_f incremental iterator
        size_t pos_u = n_b;    // u_b decremental iterator
        
        for (auto v : boost::make_iterator_range(vertices(s)))
        {
            // Get bcs & loads
            const auto& bcs  = s[v].bcs;
            const auto& load = s[v].load;
            
            // Loop on spatial components
            for (size_t i = 0u; i < ndof; ++i)
            {
                // If free DOF => write to f & post-increase position
                if (std::isnan(bcs(i)))
                    f_f(pos_f++) = load(i);
                // If BC DOF => write to u & pre-decrease position
                else 
                    u_b(--pos_u) = bcs(i);
            }
        }
        return ret;
    }
};


template <>
struct element_matrix_assembler< truss_kind<2> > 
{
    template <typename S>
    fixed_matrix<4,4>
    operator()(const typename S::edge_descriptor& e, const S& s) 
    { 
        // Get node indexes
        const auto na = source(e, s);
        const auto nb = target(e, s);
        
        // Get needed properties
        const auto EA = s[e].E * s[e].A;  // Young modulus times cross sectional area
        const auto a  = s[na].coords;     // Node A coordinates
        const auto b  = s[nb].coords;     // Node B coordinates
        
        // Compute truss length, sine and cosine
        real L  = (a-b).norm();
        real cc = (b[0]-a[0]) / L;
        real ss = (b[1]-a[1]) / L;
        
        // Cache needed values
        real c2 = EA / L * cc*cc;
        real cs = EA / L * cc*ss;
        real s2 = EA / L * ss*ss;
        
        // Build element matrix in GLOBAL coordinates
        fixed_matrix<4,4> k;
        
        k << c2, cs,-c2,-cs,
             cs, s2,-cs,-s2,
            -c2,-cs, c2, cs,
            -cs,-s2, cs, s2;
            
        return k;
    }
};

template <> 
struct element_matrix_assembler< truss_kind<3> > 
{
    template <typename S>
    fixed_matrix<6,6>
    operator()(const typename S::edge_descriptor& e, const S& s) 
    { 
        // Get node indexes
        const auto na = source(e, s);
        const auto nb = target(e, s);
        
        // Get needed properties
        const auto EA = s[e].E * s[e].A;  // Young modulus times cross sectional area
        const auto a  = s[na].coords;     // Node A coordinates
        const auto b  = s[nb].coords;     // Node B coordinates
        
        real L  = (a-b).norm();
        real cc = sqrt(EA)/L; 
        
        // Compute EA / L * cosine(theta_x/y/z) 
        real cx = cc*(b[0]-a[0]);
        real cy = cc*(b[1]-a[1]);
        real cz = cc*(b[2]-a[2]);
                
        // Build matrix
        fixed_matrix<6,6> k;
        
        k << cx*cx,  cx*cy,  cx*cz, -cx*cx, -cx*cy, -cx*cz,
             cx*cy,  cy*cy,  cy*cz, -cx*cy, -cy*cy, -cy*cz,
             cx*cz,  cy*cz,  cz*cz, -cx*cz, -cy*cz, -cz*cz,
            -cx*cx, -cx*cy, -cx*cz,  cx*cx,  cx*cy,  cx*cz,
            -cx*cy, -cy*cy, -cy*cz,  cx*cy,  cy*cy,  cy*cz,
            -cx*cz, -cy*cz, -cz*cz,  cx*cz,  cy*cz,  cz*cz;
            
        return k;
    }
};


//-------------------------------------- Post-processing -----------------------------------------//
// template <int N> 
// struct internal_forces_getter< truss_kind<N> > 
// {
//     template <typename S>
//     std::array<fixed_vector<kind_of<S>::type::ndof>,2>
//     operator()(const typename S::edge_descriptor& e, const S& s, 
//                const std::vector<real>& u, const std::vector<real>& f) 
//     {
//         const static size_t ndof = kind_of<S>::type::ndof; 
        
//         // Get local element matrix
//         auto K_el = assemble_element_matrix(e, s);
        
//         // Get bar endings idxs and applied loads
//         auto a = source(e, s);
//         auto b = target(e, s);
        
//         const auto& load_a = s[a].load;
//         const auto& load_b = s[b].load;
        
//         // NOTE: maybe use a fixed size eigen vector
//         dense_vector u_el(2*ndof); // Nodal displacements
//         //~ dense_vector q_el(2*ndof); // Nodal applied loads
//         auto q_el  = fixed_vector<2*ndof>(); 
//              q_el << s[a].load, s[b].load;
    
//         for (size_t i = 0u; i < ndof; ++i)
//         {
//             // Displacements
//             u_el(i)      = u[ndof*a+i];
//             u_el(i+ndof) = u[ndof*b+i];
//             // Loads
//             q_el(i)      = load_a(i);
//             q_el(i+ndof) = load_b(i);
//         }
//         // std::cout << std::endl;
//         auto forces = K_el*u_el - q_el;
        
//         // std::cout << "K_el:\n" << K_el << std::endl;
//         // std::cout << "u_el:\n" << u_el << std::endl;
//         // std::cout << "q_el:\n" << q_el << std::endl;
        
//         // Fill result
//         auto res = std::array<fixed_vector<ndof>,2>{};
//         auto& res_a = std::get<0>(res);
//         auto& res_b = std::get<1>(res);
        
//         for (size_t i = 0u; i < ndof; ++i)
//         {
//             res_a(i) = forces(i);
//             res_b(i) = forces(i+ndof);
//         }
//         return res;
//     }
// };

    
} //end namespace eva

# endif //__EVA_TRUSS__
