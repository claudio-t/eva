# ifndef __EVA_STRUCTURAL_ANALYSIS__
# define __EVA_STRUCTURAL_ANALYSIS__

/// @file  core_sa.hpp
/// @brief Contains common functions and functors used to solve
///        structural analysis.

# include "core.hpp"

namespace eva { namespace sa {

//-------------------------------------- Problem Solving -----------------------------------------//
/// Solves a given problem automatically deducing the structure type
template <
    typename Params = sparse_solver_params<>,
    typename Kind = void, typename Structure
    >
auto
solve(const Structure& structure, Params p = Params());

/// Generic solver functor: has to be specialized for every kind of structure
template <typename StructureKind, typename Params>
struct solver;


//------------------------------------- Problem Assembling ---------------------------------------//
/// Assembles the system stiffness matrix in global coordinates.
template <typename A, typename S>
auto
assemble_stiffness_submatrices(const S& s, const std::vector<index_t>& dofmap,
                               const size_t n_f, const size_t n_b);

/// Functor that assembles the system stiffness matrix. 
/// Has to be specialized for every kind of structure.
template <typename AlgebraType> 
struct stiffness_submatrices_assembler;

/// Builds the known terms, i.e. the force vector portion related to the free DOF and
/// the displacement vector portion associated to BC DOF
template <typename S>
std::array<dense_vector, 2>
assemble_known_terms(const S& s, const size_t n_f, const size_t n_b);

/// Functor that assembles the element (local) stiffness matrix. 
/// Has to be specialized for every kind of structure.
template <typename StructureKind> 
struct element_matrix_assembler;

/// Assembles the element (local) stiffness matrix
template <typename Structure> 
auto
assemble_element_matrix(const typename Structure::edge_descriptor& e,
                        const Structure& s);

/// Functor that assembles the known terms vectors (source term and
/// BCs related DOF). Has to be specialized for every kind of structure.
template <typename StructureKind> 
struct known_terms_assembler;


/// Counts number of non-zero entries of the stiffness matrix
/// submatrices K_ff, K_fb, K_bb
template <typename S>
std::array<std::vector<size_t>, 3>
count_nnz_entries(const S& s, const std::vector<index_t>& dofmap,
                  const size_t nf, const size_t nb);


//------------------------------------- Results Assembling ---------------------------------------//

/// Assembles the output of the solving procedure for a generic
/// given structure.
template <typename StructureKind, typename Structure>
std::vector<result<StructureKind>>
assemble_results(const dense_vector& u,
                 const dense_vector& f,
                 const Structure& s);

//-------------------------------------- Post-processing -----------------------------------------//
// Get member internal forces
template <typename Structure> 
auto
get_internal_forces(const typename Structure::edge_descriptor& e, const Structure& s, 
                    const std::vector<real>& u, const std::vector<real>& f);

/// Functor that computes member internal forces.
/// Has to be specialized for every structure type.
template <typename StructureKind> 
struct internal_forces_getter;

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
                                 //  IMPLEMENTATION  //
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
# include "../src/core_sa.tcc"


}}// end namespace eva and sa




# endif //__EVA_STRUCTURAL_ANALYSIS__
