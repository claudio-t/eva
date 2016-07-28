# ifndef __EVA_FRAME__
# define __EVA_FRAME__

/**
 * @file frame.hpp
 * @brief Contains methods and classes to solve frame structures.
 */
// eva
# include "structure.hpp"

namespace eva {
    
//####################################### DECLARATIONS #############################################

//--------------------------------------- Properties ---------------------------------------------//
/// Frame joint (node) properties
template <int N> struct frame_joint;

/// Frame element (edge) properties
template <int N> struct frame_element;

/// Frame tag type
template <int N> struct frame_type 
{         
    constexpr static int ndof = 3*(N-1); 
    constexpr static int sdim = N; 
};

//--------------------------------- 2/3-D truss typedefs -----------------------------------------//
/// 2-dimensional frame structure type
using frame2d = generic_structure<frame_joint<2>, frame_element<2>, frame_type<2>>;

/// 3-dimensional frame structure type
using frame3d = generic_structure<frame_joint<3>, frame_element<3>, frame_type<3>>;

//----------------------------------- Functors and Functions -------------------------------------//
// -- Problem Assembling -- //

/// Specializes element_matrix_assembler functor for a 2D frame
template <> 
struct element_matrix_assembler<frame_type<2>>;

/// Specializes element_matrix_assembler functor for a 3D frame
template <> 
struct element_matrix_assembler<frame_type<3>>;

/// Specializes known_terms_assembler functor for both 2D and 3D frames
template <int N>
struct known_terms_assembler<frame_type<N>>;



//####################################### DEFINITIONS ##############################################

//--------------------------------------- Properties ---------------------------------------------//
template <int N> struct frame_joint 
{    
    fixed_vector<N> coords;     ///< Joint coordinates
    
    fixed_vector<N> load;       ///< Nodal load [N]
    
    fixed_vector<2*N-3> torque; ///< Nodal torque [Nm]
    
    fixed_vector<3*(N-1)> bcs;  ///< Boundary conditions [[displ.s], [rot.s]]
    
    /// Default constructor. Initializes coords to zero,
    /// load to zero, torque to zero and bcs to nan
    frame_joint()
        : coords(fixed_vector<N>::Zero()) 
        , load  (fixed_vector<N>::Zero())
        , torque(fixed_vector<2*N-3>::Zero())
        , bcs   (fixed_vector<3*(N-1)>::Constant(std::numeric_limits<real>::quiet_NaN()))
        {}
};


template <> struct frame_element<2> 
{    
    real E; ///< Young's modulus       [Pa]
    
    real A; ///< Cross sectional area  [m^2]
    
    real I; ///< Second moment of area [m^4]
};


template <> struct frame_element<3> 
{    
    real E;  ///< Young's modulus       [GPa]
    
    real A;  ///< Cross sectional area  [m^2]
    
    real G;  ///< Shear modulus [GPa]
    
    real Ix; ///< Second moment of area wrt (local) x-axis [m^4]
    
    real Iy; ///< Second moment of area wrt (local) y-axis [m^4]
    
    real Iz; ///< Second moment of area wrt (local) z-axis [m^4]
};


//---------------------------------------- Functions ---------------------------------------------//
    
template <> 
struct element_matrix_assembler<frame_type<2>> 
{
    template <typename S>
    fixed_matrix<6,6>
    operator()(const typename S::edge_descriptor& e, const S& structure) 
    { 
        
        // Get node indexes
        auto na = source(e, structure);
        auto nb = target(e, structure);
        
        // Get needed properties
        auto EA = structure[e].E * structure[e].A;  // Young modulus times cross sectional area
        auto EI = structure[e].E * structure[e].I;  // Young modulus times momentum of inertia
        auto a  = structure[na].coords;             // Node A coordinates
        auto b  = structure[nb].coords;             // Node B coordinates
        
        // Compute truss length, sine and cosine
        real L = (a-b).norm();
        real c = (b[0]-a[0]) / L;
        real s = (b[1]-a[1]) / L;
        
        // Cache needed values
        real cc = c*c;
        real cs = c*s;
        real ss = s*s;
        
        real L2 = L*L;
        real L3 = L2*L;
                    
        // Build element matrix in GLOBAL coordinates
        // DOF: [ [F_x, F_y, M]_{node A}, [F_x, F_y, M]_{node B} ] 
        auto k = fixed_matrix<6,6> ();
        
        k <<  EA/L*cc + 12.*EI/L3*ss,  EA/L*cs - 12.*EI/L3*cs, -6.*EI/L2*s,
             -EA/L*cc - 12.*EI/L3*ss, -EA/L*cs + 12.*EI/L3*cs, -6.*EI/L2*s,
            
              EA/L*cs - 12.*EI/L3*cs,  EA/L*ss + 12.*EI/L3*cc,  6.*EI/L2*c,
             -EA/L*cs + 12.*EI/L3*cs, -EA/L*ss - 12.*EI/L3*cc,  6.*EI/L2*c,
                  
                  -6.*EI/L2*s      ,        6.*EI/L2*c     ,    4.*EI/L   ,
                   6.*EI/L2*s      ,       -6.*EI/L2*c     ,    2.*EI/L   ,
                  
             -EA/L*cc - 12.*EI/L3*ss, -EA/L*cs + 12.*EI/L3*cs,  6.*EI/L2*s,
              EA/L*cc + 12.*EI/L3*ss,  EA/L*cs - 12.*EI/L3*cs,  6.*EI/L2*s,
                                     
             -EA/L*cs + 12.*EI/L3*cs, -EA/L*ss - 12.*EI/L3*cc, -6.*EI/L2*c,
              EA/L*cs - 12.*EI/L3*cs,  EA/L*ss + 12.*EI/L3*cc, -6.*EI/L2*c,
                   
                  -6.*EI/L2*s      ,        6.*EI/L2*c     ,    2.*EI/L   ,
                   6.*EI/L2*s      ,       -6.*EI/L2*c     ,    4.*EI/L   ;
        
        return k;
    }
};


template <> 
struct element_matrix_assembler<frame_type<3>> 
{
    template <typename S>
    fixed_matrix<6,6>
    operator()(const typename S::edge_descriptor& e, const S& structure) 
    { 
        // Get node indexes
        auto na = source(e, structure);
        auto nb = target(e, structure);
        
        // Get needed properties
        auto EA = structure[e].E * structure[e].A;  // Young modulus times cross sectional area
        auto EI = structure[e].E * structure[e].I;  // Young modulus times momentum of inertia
        auto a  = structure[na].coords;     // Node A coordinates
        auto b  = structure[nb].coords;     // Node B coordinates
        
        // Compute truss length, sine and cosine
        real L = (a-b).norm();
        real c = (b[0]-a[0]) / L;
        real s = (b[1]-a[1]) / L;
        
        // Cache needed values
        real cc = c*c;
        real cs = c*s;
        real ss = s*s;
        
        real L2 = L*L;
        real L3 = L2*L;   
        
        // Build element matrix in GLOBAL coordinates
        // DOF: [ [F_x, F_y, M]_{node A}, [F_x, F_y, M]_{node B} ] 
        auto k = fixed_matrix<6,6> ();
        
        k << EA/L*cc + 12.*EI/L3*ss,  EA/L*cs - 12.*EI/L3*cs, -6.*EI/L2*s,
            -EA/L*cc - 12.*EI/L3*ss, -EA/L*cs + 12.*EI/L3*cs, -6.*EI/L2*s,
            
             EA/L*cs - 12.*EI/L3*cs,  EA/L*ss + 12.*EI/L3*cc,  6.*EI/L2*c,
            -EA/L*cs + 12.*EI/L3*cs, -EA/L*ss - 12.*EI/L3*cc,  6.*EI/L2*c,
                  
                 -6.*EI/L2*s      ,        6.*EI/L2*c     ,    4.*EI/L   ,
                  6.*EI/L2*s      ,       -6.*EI/L2*c     ,    2.*EI/L   ,
                  
            -EA/L*cc - 12.*EI/L3*ss, -EA/L*cs + 12.*EI/L3*cs,  6.*EI/L2*s,
             EA/L*cc + 12.*EI/L3*ss,  EA/L*cs - 12.*EI/L3*cs,  6.*EI/L2*s,
                                     
            -EA/L*cs + 12.*EI/L3*cs, -EA/L*ss - 12.*EI/L3*cc, -6.*EI/L2*c,
             EA/L*cs - 12.*EI/L3*cs,  EA/L*ss + 12.*EI/L3*cc, -6.*EI/L2*c,
                   
                 -6.*EI/L2*s      ,        6.*EI/L2*c     ,    2.*EI/L   ,
                  6.*EI/L2*s      ,       -6.*EI/L2*c     ,    4.*EI/L   ;
        
        return k;
    }
};


template <int N>
struct known_terms_assembler<frame_type<N>> 
{
    template <typename S>
    std::array<dense_vector, 2>
    operator()(const S& s, const size_t n_f, const size_t n_b) 
    {
        
        // Aux vars
        //~ const static size_t sdim = S::graph_bundled::sdim;  // Spatial dim
        const static size_t ndof = S::graph_bundled::ndof;  // #DOF per node (2D => 3, 3D => 6) 
        
        // Init ret var
        auto ret = std::array<dense_vector, 2> {dense_vector(n_f),
                                                dense_vector(n_b)};
                                                
        auto& f_f = std::get<0>(ret);   // Force vector (Free DOF only)
        auto& u_b = std::get<1>(ret);   // Displacement vector (BC DOF only)
        
        size_t pos_f = 0u;     // f_f incremental iterator
        size_t pos_u = n_b;    // u_b decremental iterator
        
        for (auto v : make_iterator_range(vertices(s))) {
            // Get bcs, loads and torques
            const auto& bcs    = s[v].bcs;
            const auto& load   = s[v].load;
            const auto& torque = s[v].torque;
            
            // Source term storage for DOF associated with vertex v:
            // concatenate load and torque
            auto rhs_v = fixed_vector<ndof>();
            rhs_v << load, torque;
            
            // Compute equivalent load configuration for 
            // a distributed load a/o a concentrated load
            // ...
            // ...
            
            // Fill known term vectors
            for (size_t i = 0u; i < ndof; ++i) {
                // If free DOF => write to f & post-increase position
                if (std::isnan(bcs(i)))
                    f_f(pos_f++) = rhs_v(i);
                // If BC DOF => write to u & pre-decrease position
                else 
                    u_b(--pos_u) = bcs(i);
            }
        }
        
        return ret;
    }
};
    
} //end namespace eva

# endif //__EVA_FRAME__
