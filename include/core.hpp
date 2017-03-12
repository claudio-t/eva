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

namespace eva {

    
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
/// Generic solver functor: has to be specialized for every kind of structure
template <typename StructureKind, typename Params>
struct solver;

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

/// Solves a given problem automatically deducing the structure type
template <
    typename Structure,
    typename Kind,
    typename Params//sparse_solver_params<>,
    >
auto
solve(const Structure& structure, const Kind kind = Kind(), const Params p = Params());


//------------------------------------- Problem Assembling ---------------------------------------//
/// Assembles the system stiffness matrix in global coordinates.
template <typename A, typename S, typename Kind = void>
auto
assemble_system_submatrices(const S& s, const std::vector<index_t>& dofmap,
                               const size_t n_f, const size_t n_b);

/// Functor that assembles the system stiffness matrix. 
/// Has to be specialized for every kind of structure.
template <typename StructureKind, typename AlgebraType> 
struct system_submatrices_assembler;

/// Builds the known terms, i.e. the force vector portion related to the free DOF and
/// the displacement vector portion associated to BC DOF
template <typename S>
std::array<dense_vector, 2>
assemble_known_terms(const S& s, const size_t n_f, const size_t n_b);

/// Functor that assembles the known terms vectors (source term and
/// BCs related DOF). Has to be specialized for every kind of structure.
template <typename StructureKind> 
struct known_terms_assembler;


//---------------------------------------- DOF handling ------------------------------------------//
/// Builds a DOF map where the positions of DOF associated 
/// to BCs are located at the end
template <typename Kind = void, typename S>
std::tuple<std::vector<index_t>, size_t, size_t>
build_global_dofmap(const S& s);

/// Actual implementation of build_global_dofmap.
template <typename Kind, size_t Dim>
struct global_dofmap_builder;

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


/// Assembles the output of the solving procedure for a generic
/// given structure.
template <typename StructureKind, typename Structure>
std::vector< result<StructureKind> >
assemble_results(const dense_vector& u,
                 const dense_vector& f,
                 const Structure& s);

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

} // end namespace eva

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
                                 //  IMPLEMENTATION  //
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
# include "../src/core.tcc"


# endif //__EVA_CORE__
