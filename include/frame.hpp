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
/// Truss joint (node) properties
template <size_t N> struct frame_joint;

/// Truss element (edge) properties
template <size_t N> struct frame_element;

/// Truss tag type
template <size_t N> struct frame_type {     
    
    constexpr static size_t ndof = 3u*(N-1u); 
    constexpr static size_t sdim = N; 
};

//--------------------------------- 2/3-D truss typedefs -----------------------------------------//
/// 2-dimensional truss structure type
typedef generic_structure<frame_joint<2>, frame_element<2>, frame_type<2>> frame2d; 

/// 3-dimensional truss structure type
typedef generic_structure<frame_joint<3>, frame_element<3>, frame_type<3>> frame3d; 



//####################################### DEFINITIONS ##############################################

//--------------------------------------- Properties ---------------------------------------------//
template <size_t N> struct frame_joint {
    
    fixed_array1d<N, real, zero_init> coords;       ///< Joint coordinates
    
    fixed_array1d<N, real, zero_init> load;         ///< Nodal load [N]
    
    fixed_array1d<2u*N-3u, real, zero_init> torque; ///< Nodal torque [Nm]
    
    fixed_array1d<3u*(N-1u), real, nan_init> bcs;   ///< Boundary conditions [[displ.s], [rot.s]]
};


template <> struct frame_element<2> {
    
    real E; ///< Young's modulus       [GPa]
    
    real A; ///< Cross sectional area  [m^2]
    
    real I; ///< Second moment of area [m^4]
};


template <> struct frame_element<3> {
    
    real E;  ///< Young's modulus       [GPa]
    
    real A;  ///< Cross sectional area  [m^2]
    
    real G;  ///< Shear modulus [GPa]
    
    real Ix; ///< Second moment of area wrt (local) x-axis [m^4]
    
    real Iy; ///< Second moment of area wrt (local) y-axis [m^4]
    
    real Iz; ///< Second moment of area wrt (local) z-axis [m^4]
};


//---------------------------------------- Functions ---------------------------------------------//
    
template <> 
struct element_matrix_assembler<frame_type<2>> {
    template <typename S>
    fixed_array2d<6,6>
    operator()(const typename S::edge_descriptor& e, const S& structure) { 
        
        // Get node indexes
        auto na = source(e, structure);
        auto nb = target(e, structure);
        
        // Get needed properties
        auto EA = structure[e].E * structure[e].A;  // Young modulus times cross sectional area
        auto EI = structure[e].E * structure[e].I;  // Young modulus times momentum of inertia
        auto a  = structure[na].coords;             // Node A coordinates
        auto b  = structure[nb].coords;             // Node B coordinates
        
        // Compute truss length, sine and cosine
        real L = sqrt((b[0]-a[0])*(b[0]-a[0]) + (b[1]-a[1])*(b[1]-a[1]));
        real c = (b[0]-a[0]) / L;
        real s = (b[1]-a[1]) / L;
        
        // Cache needed values
        real cc = c*c;
        real cs = c*s;
        real ss = s*s;
        
        real L2 = L*L;
        real L3 = L2*L;
            
        // Convenient alias
        typedef  fixed_array2d<6,6>::value_type arr1d;    
        
        // Build element matrix in GLOBAL coordinates
        // DOF: [ [F_x, F_y, M]_{node A}, [F_x, F_y, M]_{node B} ] 
        auto k = fixed_array2d<6,6> {

            arr1d{ EA/L*cc + 12.*EI/L3*ss,  EA/L*cs - 12.*EI/L3*cs, -6.*EI/L2*s,
                  -EA/L*cc - 12.*EI/L3*ss, -EA/L*cs + 12.*EI/L3*cs, -6.*EI/L2*s },
            
            arr1d{ EA/L*cs - 12.*EI/L3*cs,  EA/L*ss + 12.*EI/L3*cc,  6.*EI/L2*c,
                  -EA/L*cs + 12.*EI/L3*cs, -EA/L*ss - 12.*EI/L3*cc,  6.*EI/L2*c },
                  
            arr1d{      -6.*EI/L2*s      ,         6.*EI/L2*c     ,    4.*EI/L ,
                         6.*EI/L2*s      ,        -6.*EI/L2*c     ,    2.*EI/L  },
                  
            arr1d{-EA/L*cc - 12.*EI/L3*ss, -EA/L*cs + 12.*EI/L3*cs,  6.*EI/L2*s,
                   EA/L*cc + 12.*EI/L3*ss,  EA/L*cs - 12.*EI/L3*cs,  6.*EI/L2*s },
                                     
            arr1d{-EA/L*cs + 12.*EI/L3*cs, -EA/L*ss - 12.*EI/L3*cc, -6.*EI/L2*c,
                   EA/L*cs - 12.*EI/L3*cs,  EA/L*ss + 12.*EI/L3*cc, -6.*EI/L2*c },
                   
            arr1d{      -6.*EI/L2*s      ,         6.*EI/L2*c     ,    2.*EI/L ,
                         6.*EI/L2*s      ,        -6.*EI/L2*c     ,    4.*EI/L  },
                                     
        };
        //~ std::cout << std::endl
        //~ 
                  //~ << "------------------------ K_" << na << nb << "------------------------"
                  //~ << std::endl;
        //~ 
        //~ for (auto& row : k) {
            //~ for (const auto el : row) std::cout << el << "\t";
            //~ std::cout << std::endl;
        //~ }
//~ 
        //~ std::cout << std::endl
                  //~ << std::endl;
        
        return k;
    }
};


//~ template <> 
//~ struct element_matrix_assembler<frame_type<3>> {
    //~ template <typename S>
    //~ fixed_array2d<6,6>
    //~ operator()(const typename S::edge_descriptor& e, const S& structure) { 
        //~ 
        //~ // Get node indexes
        //~ auto na = source(e, structure);
        //~ auto nb = target(e, structure);
        //~ 
        //~ // Get needed properties
        //~ auto EA = structure[e].E * structure[e].A;  // Young modulus times cross sectional area
        //~ auto EI = structure[e].E * structure[e].I;  // Young modulus times momentum of inertia
        //~ auto a  = structure[na].coords;     // Node A coordinates
        //~ auto b  = structure[nb].coords;     // Node B coordinates
        //~ 
        //~ // Compute truss length, sine and cosine
        //~ real L = sqrt((b[0]-a[0])*(b[0]-a[0]) + (b[1]-a[1])*(b[1]-a[1]));
        //~ real c = (b[0]-a[0]) / L;
        //~ real s = (b[1]-a[1]) / L;
        //~ 
        //~ // Cache needed values
        //~ real cc = c*c;
        //~ real cs = c*s;
        //~ real ss = s*s;
        //~ 
        //~ real L2 = L*L;
        //~ real L3 = L2*L;
            //~ 
        //~ // Convenient alias
        //~ typedef  fixed_array2d<6,6>::value_type arr1d;    
        //~ 
        //~ // Build element matrix in GLOBAL coordinates
        //~ // DOF: [ [F_x, F_y, M]_{node A}, [F_x, F_y, M]_{node B} ] 
        //~ auto k = fixed_array2d<6,6> {
//~ 
            //~ arr1d{ EA/L*cc + 12.*EI/L3*ss,  EA/L*cs - 12.*EI/L3*cs, -6.*EI/L2*s,
                  //~ -EA/L*cc - 12.*EI/L3*ss, -EA/L*cs + 12.*EI/L3*cs, -6.*EI/L2*s },
            //~ 
            //~ arr1d{ EA/L*cs - 12.*EI/L3*cs,  EA/L*ss + 12.*EI/L3*cc,  6.*EI/L2*c,
                  //~ -EA/L*cs + 12.*EI/L3*cs, -EA/L*ss - 12.*EI/L3*cc,  6.*EI/L2*c },
                  //~ 
            //~ arr1d{      -6.*EI/L2*s      ,         6.*EI/L2*c     ,    4.*EI/L ,
                         //~ 6.*EI/L2*s      ,        -6.*EI/L2*c     ,    2.*EI/L  },
                  //~ 
            //~ arr1d{-EA/L*cc - 12.*EI/L3*ss, -EA/L*cs + 12.*EI/L3*cs,  6.*EI/L2*s,
                   //~ EA/L*cc + 12.*EI/L3*ss,  EA/L*cs - 12.*EI/L3*cs,  6.*EI/L2*s },
                                     //~ 
            //~ arr1d{-EA/L*cs + 12.*EI/L3*cs, -EA/L*ss - 12.*EI/L3*cc, -6.*EI/L2*c,
                   //~ EA/L*cs - 12.*EI/L3*cs,  EA/L*ss + 12.*EI/L3*cc, -6.*EI/L2*c },
                   //~ 
            //~ arr1d{      -6.*EI/L2*s      ,         6.*EI/L2*c     ,    2.*EI/L ,
                         //~ 6.*EI/L2*s      ,        -6.*EI/L2*c     ,    4.*EI/L  },
                                     //~ 
        //~ };
        //~ return k;
    //~ }
//~ };


template <size_t N>
struct known_terms_assembler<frame_type<N>> {
    template <typename S>
    std::array<dense_vector, 2>
    operator()(const S& s, const size_t n_f, const size_t n_b) {
        
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
            
            // Source term storage for DOF associated with vertex v
            auto rhs_v = fixed_vector<ndof>();
            
            // Store LOAD components
            const auto l_size = load.size();
            for (size_t i = 0u; i < l_size; ++i) 
                rhs_v(i) = load(i);
            
            // Store TORQUE components
            const auto t_size = torque.size();
            for (size_t i = 0u; i < t_size; ++i) 
                rhs_v(i+l_size) = torque(i);
            
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