# ifndef __EVA_CORE__
# define __EVA_CORE__

/// @file  core.hpp
/// @brief Contains common classes, functions and functors used to deal
///        with any of the implemented kind problems.


// boost
# include <boost/graph/adjacency_list.hpp>
# include <boost/graph/graph_traits.hpp>
# include <boost/range/iterator_range.hpp>
// eigen
# include <eigen3/Eigen/Dense>
# include <eigen3/Eigen/Sparse>

// std
# include <array>
// /* remove me! */ # include <chrono> /* remove me! */
/* remove me! */ # include <iostream> /* remove me! */



namespace eva { //namespace core {

    
//######################################## Declarations ############################################


//------------------------------------- Aliases and Typedefs -------------------------------------//

/// Alias for type representing a real number
using real = double;

/// Alias for a type representing an index
using index_t = std::size_t;

/// Alias for a fixed size vector
template <int N, typename T = real>
using fixed_vector = Eigen::Matrix<T, N, 1>;

/// Alias for a fixed size matrix
template <int M, int N, typename T = real>
using fixed_matrix = Eigen::Matrix<T, M, N>;

/// Alias for a dense matrix used for computation
using dense_matrix = Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic>;

/// Alias for a dense vector used for computation
using dense_vector = Eigen::Matrix<real, Eigen::Dynamic, 1>;

/// Alias for a sparse matrix used for computation
using sparse_matrix = Eigen::SparseMatrix<real>; // Does not work with CRS!

/// Alias for a sparse vector used for computation (CCS format)
using sparse_vector = Eigen::SparseVector<real>;

/// Alias for a generic structure
template < typename JointProps, typename ElementProps, typename Kind >
using generic_structure = boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS,
                                                 JointProps, ElementProps, Kind, boost::vecS >;


//-------------------------------------- Problem Solving -----------------------------------------//
/// Specifies the algebra (dense/sparse) and solver (see Eigen docs) types that have to be used
template <typename Algebra, typename Solver>
struct solver_params;

/// Dense algebra tag type
struct dense_algebra_t  {};

/// Sparse algebra tag type
struct sparse_algebra_t {};

/// Default settings for a dense solver
template <typename S = Eigen::LDLT<dense_matrix> >
using dense_solver_params = solver_params<dense_algebra_t, S>;
// struct dense_solver_params;

/// Specifies the algebra and solver types that the sparse solver has to use
template <typename S = Eigen::ConjugateGradient<sparse_matrix> >
using sparse_solver_params = solver_params<sparse_algebra_t, S>;


//---------------------------------------- DOF handling ------------------------------------------//
/// Builds a DOF map where the positions of DOF associated 
/// to BCs are located at the end
template <typename S>
std::tuple<std::vector<index_t>, size_t, size_t>
build_global_dofmap(const S& s);

/// Builds the (element) local to global DOF map
template <int N> 
fixed_vector<2*N, index_t> 
build_local_to_global_dofmap(const index_t idx_a, const index_t idx_b,
                             const std::vector<index_t>& dm);



/// Takes two vectors -one containing a quantity on free DOF, the other on BC DOF-
/// and merge them together reordering DOF accordingly to the supplied DOF map
dense_vector 
merge_and_reorder(const dense_vector& v_f, const dense_vector& v_b, 
                  const std::vector<index_t>& dofmap);



//------------------------------------- Results Assembling ---------------------------------------//
/// Represents the results of the solving procedure applied to a
/// certain structure. Has to be specialized for every kind of structure.
template <typename StructureKind>
struct result;


//---------------------------------------- Type helpers -----------------------------------------//
/// Trait used to recognize the kind of a given structure
template <typename S>
struct kind_of;

/// Trait used to recognize the joint (vertex properties) type of a given structure
template <typename S>
struct joint_of;

/// Trait used to recognize the element (edge properties) type of a given structure
template <typename S>
struct element_of;

/// Trait used to recognize the result type of a given structure
template <typename S>
struct result_of;


//######################################## Definitions #############################################      
//-------------------------------------- Problem Solving -----------------------------------------//
template <typename A, typename S>
struct solver_params
{
    using algebra_t = A;
    using solver_t  = S;
};



 
//---------------------------------------- DOF handling ------------------------------------------//
template <typename S>
std::tuple<std::vector<index_t>, size_t, size_t>
build_global_dofmap(const S& s) 
{
    constexpr size_t dim = kind_of<S>::type::ndof;
    
    // Initialize the DOF map
    size_t n_verts = num_vertices(s);
    auto  dofmap   = std::vector<index_t> (dim*n_verts);
    
    size_t n_bcs  = 0u;
    size_t bc_pos = dofmap.size()-1u;
    size_t ff_pos = 0u;
    
    for (size_t v = 0u; v < n_verts; ++v)
    {
        const auto& bcs = s[v].bcs;
        
        for (size_t i = 0u; i < dim; ++i)
        {
            // If coordinate = NaN => Free DOF else BC DOF
            if (std::isnan(bcs[i])) {
                // Put DOF at the beginning and increment position
                dofmap[dim*v+i] = ff_pos++;
            }
            else {
                // Put DOF at the end and increment position
                dofmap[dim*v+i] = bc_pos--;
                ++n_bcs;
            }
        }
    }
    // Return (dofmap, #(free DOF), #(BC DOF))
    return std::make_tuple(dofmap, dofmap.size()-n_bcs, n_bcs);
}


template <int N> 
fixed_vector<2*N, index_t> 
build_local_to_global_dofmap(const index_t idx_a, const index_t idx_b,
                             const std::vector<index_t>& dm) 
{    
    auto ret = fixed_vector<2*N, index_t> ();
    
    for (auto i = 0; i < N; ++i)
    {
        ret(i)   = dm[N*idx_a + i];    // Node A -> positions from 0 to  N-1
        ret(i+N) = dm[N*idx_b + i];    // Node B -> positions from N to 2N-1
    }
    
    return ret;
}


inline
dense_vector
merge_and_reorder(const dense_vector& v_f, 
                  const dense_vector& v_b, 
                  const std::vector<index_t>& dofmap) 
{
    // DOF sizes
    size_t n_f = v_f.size(); // Free DOF
    size_t n_b = v_b.size(); // BC DOF
    size_t n_t = n_f + n_b;  // Tot nr of DOF

    // Pre-allocate return variable
    auto rv = dense_vector(n_t);

    // For each DOF
    for (size_t i = 0u; i < n_t; ++i)
    {
        // Retrieve mapped position
        auto ii = dofmap[i];
        // Fill position by picking the value from the proper vector
         rv[i] = (ii < n_f) ? v_f(ii) : v_b(ii-n_f);
    }
    return rv;
}


// dense_vector
// reorder(const dense_vector& v, const std::vector<index_t>& dofmap)
// {
//     // Pre-allocate results
//     auto n_t = v.size();
//     auto rv  = dense_vector(n_t);
    
//     // Reorder
//     for (size_t i = 0u; i < n_t; ++i)
//         rv[i] = v[dofmap[i]];
    
//     return rv;
// }


//---------------------------------------- Type helpers -----------------------------------------//
template <typename S>
struct kind_of
{
    using type = typename S::graph_bundled;
};

template <typename S>
struct joint_of
{
    using type = typename S::vertex_bundled;
};

template <typename S>
struct element_of
{
    using type = typename S::edge_bundled;
};
    
template <typename S>
struct result_of
{
    using type = result<typename kind_of<S>::type>;
};

} // end namespace eva


# endif //__EVA_CORE__
