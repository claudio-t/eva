# ifndef __EVA_TRUSS__
# define __EVA_TRUSS__

/**
 * @file truss.hpp
 * @brief Contains methods and classes to solve truss structures.
 */
// eva
# include "structure.hpp"


namespace eva {
    
//####################################### DECLARATIONS #############################################

//--------------------------------------- Properties ---------------------------------------------//
/// Truss joint (node) properties
template <int N> struct truss_joint;

/// Truss element (edge) properties
template <int N> struct truss_element;

/// Truss tag type
template <int N> struct truss_type 
{     
    constexpr static int ndof = N; 
    constexpr static int sdim = N; 
};


//--------------------------------- 2/3-D Truss Typedefs -----------------------------------------//
/// 2-dimensional truss structure type
typedef generic_structure<truss_joint<2>, truss_element<2>, truss_type<2>> truss2d; 

/// 3-dimensional truss structure type
typedef generic_structure<truss_joint<3>, truss_element<3>, truss_type<3>> truss3d; 



//----------------------------------- Functors and Functions -------------------------------------//

// -- Post-processing --//
template <int N> struct internal_forces_getter<truss_type<N>>;


// -- Problem Assembling -- //
/// Builds the truss element matrix in GLOBAL coordinates
template <int N> 
fixed_matrix<2*N, 2*N> 
truss_element_matrix(const fixed_vector<N>& a, const fixed_vector<N>& b, real EA);

template<> 
fixed_matrix<4,4>
truss_element_matrix<2>(const fixed_vector<2>& a, const fixed_vector<2>& b, real EA);

template<> 
fixed_matrix<6,6>
truss_element_matrix<3>(const fixed_vector<3>& a, const fixed_vector<3>& b, real EA);

template <int N> struct known_terms_assembler<truss_type<N>>;
    


//####################################### DEFINITIONS ##############################################
  
//--------------------------------------- Properties ---------------------------------------------//
template <int N> struct truss_joint 
{
    fixed_vector<N> coords; ///< Joint coordinates
    
    fixed_vector<N> load;   ///< Applied load [N]
    
    fixed_vector<N> bcs;     ///< Boundary conditions
    
    truss_joint()
    : coords(fixed_vector<N>::Zero()) 
    , load  (fixed_vector<N>::Zero())
    , bcs   (fixed_vector<N>::Constant(std::numeric_limits<real>::quiet_NaN()))
    {}
};


template <int N> struct truss_element 
{    
    real E; ///< Young's modulus      [GPa]
    
    real A; ///< Cross sectional area [m^2]
};



//----------------------------------- Functors and Functions -------------------------------------//

// -- Post-processing --//
template <int N> 
struct internal_forces_getter<truss_type<N>> 
{
    template <typename S>
    std::array<fixed_vector<S::graph_bundled::ndof>,2>
    operator()(const typename S::edge_descriptor& e, const S& s, 
               const std::vector<real>& u, const std::vector<real>& f) 
    {
        const static size_t ndof = S::graph_bundled::ndof; 
        
        // Get local element matrix 
        //~ std::cout <<  s[e].E * s[e].A << std::endl;
        auto K_el = assemble_element_matrix(e, s);
        
        // Get bar endings idxs and applied loads
        auto a = source(e, s);
        auto b = target(e, s);
        
        const auto& load_a = s[a].load;
        const auto& load_b = s[b].load;
        
        // NOTE: maybe use a fixed size eigen vector
        dense_vector u_el(2*ndof); // Nodal displacements
        //~ dense_vector q_el(2*ndof); // Nodal applied loads
        auto q_el  = fixed_vector<2*ndof>(); 
             q_el << s[a].load, s[b].load;
    
        for (size_t i = 0u; i < ndof; ++i) {
            // Displacements
            u_el(i)      = u[ndof*a+i];
            u_el(i+ndof) = u[ndof*b+i];
            // Loads
            q_el(i)      = load_a(i);
            q_el(i+ndof) = load_b(i);
        }
        std::cout << std::endl;
        auto forces = K_el*u_el - q_el;
        
        std::cout << "K_el:\n" << K_el << std::endl;
        std::cout << "u_el:\n" << u_el << std::endl;
        std::cout << "q_el:\n" << q_el << std::endl;
        
        // Fill result
        auto res = std::array<fixed_vector<ndof>,2>{};
        auto& res_a = std::get<0>(res);
        auto& res_b = std::get<1>(res);
        
        for (size_t i = 0u; i < ndof; ++i) {
            res_a(i) = forces(i);
            res_b(i) = forces(i+ndof);
        }
        
        return res;
    }
};


// -- Problem Assembling -- //
template <>
struct element_matrix_assembler<truss_type<2>> 
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
        real L  = sqrt((b[0]-a[0])*(b[0]-a[0]) + (b[1]-a[1])*(b[1]-a[1]));
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
struct element_matrix_assembler<truss_type<3>> 
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
        
        real L  = sqrt((b[0]-a[0])*(b[0]-a[0]) + (b[1]-a[1])*(b[1]-a[1]) + (b[2]-a[2])*(b[2]-a[2]));
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


template <int N>
struct known_terms_assembler<truss_type<N>> 
{
    template <typename S>
    std::array<dense_vector, 2>
    operator()(const S& s, const size_t n_f, const size_t n_b) 
    {
        const static int ndof = S::graph_bundled::ndof;   // #DOF per node (2D => 3, 3D => 6) 
        auto ret = std::array<dense_vector, 2> {dense_vector(n_f),
                                                dense_vector(n_b)};
                                                
        auto& f_f = std::get<0>(ret);   // Force vector (Free DOF only)
        auto& u_b = std::get<1>(ret);   // Displacement vector (BC DOF only)
        
        size_t pos_f = 0u;     // f_f incremental iterator
        size_t pos_u = n_b;    // u_b decremental iterator
        
        for (auto v : make_iterator_range(vertices(s))) {
            // Get bcs & loads
            const auto& bcs  = s[v].bcs;
            const auto& load = s[v].load;
            
            // Loop on spatial components
            for (size_t i = 0u; i < ndof; ++i) {
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

//######################################## Unused/old ##############################################
//~ template <> fixed_array1d<4, index_t> 
//~ build_local_to_global_dofmap<2>(index_t na, index_t nb, const std::vector<index_t>& dm) {
    //~ 
    //~ return  fixed_array1d<4, index_t> {dm[2*na], dm[2*na+1],
                                       //~ dm[2*nb], dm[2*nb+1]};
//~ }
//~ 
//~ template <> fixed_array1d<6, index_t> 
//~ build_local_to_global_dofmap<3>(index_t na, index_t nb, const std::vector<index_t>& dm) {
    //~ 
    //~ return  fixed_array1d<6,index_t> {dm[3*na], dm[3*na+1], dm[3*na+2],
                                      //~ dm[3*nb], dm[3*nb+1], dm[3*nb+2]};
//~ }

//~ template<> fixed_array2d<4,4>
//~ truss_element_matrix<2>(const fixed_array1d<2>& a,
                        //~ const fixed_array1d<2>& b,
                        //~ real EA) {
//~ 
    //~ // Compute truss length, sine and cosine
    //~ real L = sqrt((b[0]-a[0])*(b[0]-a[0]) + (b[1]-a[1])*(b[1]-a[1]));
    //~ real c = (b[0]-a[0]) / L;
    //~ real s = (b[1]-a[1]) / L;
    //~ 
    //~ // Cache needed values
    //~ real cc = EA / L * c*c;
    //~ real cs = EA / L * c*s;
    //~ real ss = EA / L * s*s;
        //~ 
    //~ // Convenient alias
    //~ typedef  fixed_array2d<4,4>::value_type arr1d;    
    //~ 
    //~ // Build element matrix in GLOBAL coordinates
    //~ auto k = fixed_array2d<4,4> {arr1d{ cc, cs,-cc,-cs },
                                 //~ arr1d{ cs, ss,-cs,-ss },
                                 //~ arr1d{-cc,-cs, cc, cs },
                                 //~ arr1d{-cs,-ss, cs, ss }};
    //~ return k;
//~ }


//~ template<> fixed_array2d<6,6>
//~ truss_element_matrix<3>(const fixed_array1d<3>& a,
                        //~ const fixed_array1d<3>& b,
                        //~ real EA) {
     //~ 
    //~ real L  = sqrt((b[0]-a[0])*(b[0]-a[0]) + (b[1]-a[1])*(b[1]-a[1]) + (b[2]-a[2])*(b[2]-a[2]));
    //~ real cc = sqrt(EA)/L; 
    //~ 
    //~ // Compute EA / L * cosine(theta_x/y/z) 
    //~ real cx = cc*(b[0]-a[0]);
    //~ real cy = cc*(b[1]-a[1]);
    //~ real cz = cc*(b[2]-a[2]);
    //~ 
    //~ // Convenient alias
    //~ typedef  fixed_array2d<6,6>::value_type arr1d;
    //~ 
    //~ // Build matrix
    //~ auto k = fixed_array2d<6,6> {arr1d{ cx*cx,  cx*cy,  cx*cz, -cx*cx, -cx*cy, -cx*cz },
                                 //~ arr1d{ cx*cy,  cy*cy,  cy*cz, -cx*cy, -cy*cy, -cy*cz },
                                 //~ arr1d{ cx*cz,  cy*cz,  cz*cz, -cx*cz, -cy*cz, -cz*cz },
                                 //~ arr1d{-cx*cx, -cx*cy, -cx*cz,  cx*cx,  cx*cy,  cx*cz },
                                 //~ arr1d{-cx*cy, -cy*cy, -cy*cz,  cx*cy,  cy*cy,  cy*cz },
                                 //~ arr1d{-cx*cz, -cy*cz, -cz*cz,  cx*cz,  cy*cz,  cz*cz }};
    //~ return k;                          
//~ }

//~ template <size_t N>
//~ std::ostream& operator<<(std::ostream& os, const eva::truss_type<N>& tt) {
    //~ return os << "I'm a truss";
//~ }
//~ 
//~ 
//~ template <size_t N>
//~ std::ostream& operator<<(std::ostream& os, const eva::truss_joint<N>& tj) {
    //~ std::copy(begin(tj.coords), end(tj.coords), std::ostream_iterator<eva::real>(os, " "));
    //~ return os;
//~ }
//~ 
//~ template <size_t N>
//~ std::ostream& operator<<(std::ostream& os, const eva::truss_element<N>& te) {
    //~ return os << te.E;
//~ }
//~ 
//~ 
//~ template <size_t N>
//~ std::istream& operator>>(std::istream& is, eva::truss_type<N>& tt) {
    //~ return is;
//~ }
//~ 
//~ template <size_t N>
//~ std::istream& operator>>(std::istream& is, eva::truss_joint<N>& te) {
    //~ return is;
//~ }
//~ 
//~ template <size_t N>
//~ std::istream& operator>>(std::istream& is, eva::truss_element<N>& tj) {
    //~ return is;
//~ }
    
} //end namespace eva

# endif //__EVA_TRUSS__
