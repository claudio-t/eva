# ifndef __EVA_FRAME__
# define __EVA_FRAME__

/**
 * @file frame.hpp
 * @brief Contains methods and classes to solve frame structures.
 */
// eva
# include "core_sa.hpp"
# include "truss.hpp"

namespace eva {
    
//######################################## DECLARATIONS ############################################

//----------------------------------------- Properties -------------------------------------------//
/// Frame tag type
template <int N> struct frame_kind;

/// Frame joint (node) properties
template <int N> struct frame_joint;

/// Frame element (edge) properties
template <int N> struct frame_element;

//-------------------------------------- Frame tag types -----------------------------------------//
/// 2-dimensional frame structure tag type
using frame2d = generic_structure< frame_joint<2>, frame_element<2>, frame_kind<2> >;

/// 3-dimensional frame structure tag type
using frame3d = generic_structure< frame_joint<3>, frame_element<3>, frame_kind<3> >;


//------------------------------------- Problem Assembling ---------------------------------------//

/// Specializes system_submatrices_assembler functor for both 2D and
/// 3D trusses.
/// This specialization simply calls the stiffness_submatrices_assembler
/// functor, which -in turn- automatically selects the proper
/// specialization of the element_matrix_assembler functor.
template <typename A, int N>
struct system_submatrices_assembler<A, frame_kind<N> >;

/// Specializes known_terms_assembler functor for both 2D and 3D frames
template <int N>
struct known_terms_assembler< frame_kind<N> >;

/// Specializes element_matrix_assembler functor for a 2D frame
template <> 
struct element_matrix_assembler< frame_kind<2> >;

/// Specializes element_matrix_assembler functor for a 3D frame
template <> 
struct element_matrix_assembler< frame_kind<3> >;


//------------------------------------- Results Assembling ---------------------------------------//
template <int N>
struct result< frame_kind<N> >;

//------------------------------------- Results Assembling ---------------------------------------//
template <int N>
struct result< frame_kind<N> >
{
    using kind_type = frame_kind<N>;
    
    constexpr static int ndof = frame_kind<N>::ndof;
    constexpr static int sdim = frame_kind<N>::sdim;
    constexpr static int rdim = frame_kind<N>::rdim;
    
    fixed_vector<sdim> displacement;
    fixed_vector<rdim> rotation;
    fixed_vector<sdim> reaction;
    fixed_vector<rdim> react_torque;
    
    // Initialize members
    result(const fixed_vector<ndof>& u,
           const fixed_vector<ndof>& f,
           const frame_joint<N>& p)
        : displacement(u.head(sdim))
        , rotation    (u.tail(rdim))
        , reaction    (f.head(sdim))
        , react_torque(f.tail(rdim))
    {        
        // Keep reaction only if there is no load applied
        if (p.load != fixed_vector<sdim>::Zero())
            reaction = fixed_vector<sdim>::Zero();
        
        // Keep torque reaction only if there is no load applied
        if (p.torque != fixed_vector<rdim>::Zero())
            reaction = fixed_vector<sdim>::Zero();
    }
};


//######################################## DEFINITIONS #############################################

//---------------------------------------- Properties --------------------------------------------//
template <int N> struct frame_kind
{         
    constexpr static int ndof = 3*(N-1); // #(tot DOF)
    constexpr static int sdim = N;       // #(spatial coordinates)
    constexpr static int rdim = 2*N-3;   // #(rotations/torque
                                         // #dimensions)

    using joint_type   = frame_joint<N>;
    using element_type = frame_element<N>;
    using result_type  = result< frame_kind<N> >;

    
    using default_dense_solver_t  = Eigen::LDLT<dense_matrix>;
    using default_sparse_solver_t = Eigen::ConjugateGradient<sparse_matrix>;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {/* Do nothing with constexpr static stuff! */}
};

template <int N>
std::ostream & operator<<(std::ostream & os, const frame_kind<N> & k)
{
    return os
        << "ndof: " << k.ndof << "\n"
        << "sdim: " << k.sdim << "\n"
        << "rdim: " << k.rdim << "\n"
        ;
}


template <int N> struct frame_joint 
{

    // fixed_vector<N>     coords; ///< Joint coordinates [m]
    // fixed_vector<N>     load;   ///< Nodal load   [N]
    // fixed_vector<2*N-3> torque; ///< Nodal torque [Nm]
    // fixed_vector<3*N-3> bcs;    ///< Boundary conditions [[displ.s], [rot.s]]
    
    // /// Default constructor. Initializes coords to zero,
    // /// load to zero, torque to zero and bcs to nan
    // frame_joint()
    //     : coords(fixed_vector<N>::Zero()) 
    //     , load  (fixed_vector<N>::Zero())
    //     , torque(fixed_vector<2*N-3>::Zero())
    //     , bcs   (fixed_vector<3*N-3>::Constant(std::numeric_limits<real>::quiet_NaN()))
    //     {}   
    
    using k_ = frame_kind<N>;
    // constexpr static int ndof = frame_kind<N>::ndof;
    // constexpr static int sdim = frame_kind<N>::sdim;
    // constexpr static int rdim = frame_kind<N>::rdim;
    
    fixed_vector<k_::sdim> coords; ///< Joint coordinates [m]
    fixed_vector<k_::sdim> load;   ///< Nodal load   [N]
    fixed_vector<k_::rdim> torque; ///< Nodal torque [Nm]
    fixed_vector<k_::ndof> bcs;    ///< Boundary conditions [[displ.s], [rot.s]]
    
    /// Default constructor. Initializes coords to zero,
    /// load to zero, torque to zero and bcs to nan
    frame_joint()
        : coords(fixed_vector<k_::sdim>::Zero()) 
        , load  (fixed_vector<k_::sdim>::Zero())
        , torque(fixed_vector<k_::rdim>::Zero())
        , bcs   (fixed_vector<k_::ndof>::Constant(std::numeric_limits<real>::quiet_NaN()))
        {}

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        // ar & ndof & sdim & rdim;
        ar & coords & load & torque & bcs;
    }
};

template <int N>
std::ostream & operator<<(std::ostream & os, const frame_joint<N> & j)
{
    return os
        << "coords: " << j.coords.transpose() << "\n"
        << "  load: " <<   j.load.transpose() << "\n"
        << "torque: " << j.torque.transpose() << "\n"
        << "   bcs: " <<    j.bcs.transpose() << "\n"
        ;
}



template <> struct frame_element<2> 
{    
    real E; ///< Young's modulus       [Pa]
    real A; ///< Cross sectional area  [m^2]
    real I; ///< Second moment of area [m^4]

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & E & A & I;
    }
};

inline std::ostream & operator<<(std::ostream & os, const frame_element<2> & e)
{
    return os
        << "E: " << e.E << "\n"
        << "A: " << e.A << "\n"
        << "I: " << e.I << "\n"
        ;
}


template <> struct frame_element<3> 
{    
    real E; ///< Young's modulus       [Pa]    
    real A; ///< Cross sectional area  [m^2]
    real G; ///< Shear modulus         [Pa]
    
    real Ix; ///< Second moment of area wrt (local) x-axis [m^4]
    real Iy; ///< Second moment of area wrt (local) y-axis [m^4]
    real Iz; ///< Second moment of area wrt (local) z-axis [m^4]

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & E & A & G & Ix & Iy & Iz;
    }
};


inline std::ostream & operator<<(std::ostream & os, const frame_element<3> & e)
{
    return os
        << "E: " << e.E << "\n"
        << "A: " << e.A << "\n"
        << "G: " << e.G << "\n"
        << "Ix: " << e.Ix << "\n"
        << "Iy: " << e.Iy << "\n"
        << "Iz: " << e.Iz << "\n"
        ;
}


//-------------------------------------- Problem Assembling --------------------------------------//

template <typename A, int N>
struct system_submatrices_assembler<A, frame_kind<N> >
{
    template <typename S>
    auto
    operator()(const S& s,
               const std::vector<index_t>& dofmap,
               const size_t n_f, const size_t n_b)
    {
        return stiffness_submatrices_assembler<frame_kind<N>, A>()(s, dofmap, n_f, n_b);
    }    
};


template <int N>
struct known_terms_assembler< frame_kind<N> > 
{
    template <typename S>
    std::array<dense_vector, 2>
    operator()(const S& s, const size_t n_f, const size_t n_b) 
    {
        
        // Aux vars
        const static size_t ndof = frame_kind<N>::ndof;
        
        // Init ret var
# pragma clang diagnostic push
# pragma clang diagnostic ignored "-Wmissing-braces"
        auto ret = std::array<dense_vector, 2> {dense_vector(n_f),
                                                dense_vector(n_b)};
# pragma clang diagnostic pop
        
        auto& f_f = std::get<0>(ret);   // Force vector (Free DOF only)
        auto& u_b = std::get<1>(ret);   // Displacement vector (BC DOF only)
        
        size_t pos_f = 0u;     // f_f incremental iterator
        size_t pos_u = n_b;    // u_b decremental iterator
        
        for (auto v : boost::make_iterator_range(vertices(s)))
        {
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
            for (size_t i = 0u; i < ndof; ++i)
            {
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


template <> 
struct element_matrix_assembler< frame_kind<2> > 
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
struct element_matrix_assembler< frame_kind<3> > 
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
        // DOF: [ [F_x, F_y, F_z, M_x, M_y, M_z]_{node A},
        //        [F_x, F_y, F_z, M_x, M_y, M_z]_{node B} ] 
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

} //end namespace eva

# endif //__EVA_FRAME__
