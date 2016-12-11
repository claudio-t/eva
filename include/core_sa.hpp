# ifndef __EVA_STRUCTURAL_ANALYSIS__
# define __EVA_STRUCTURAL_ANALYSIS__

/// @file  core_sa.hpp
/// @brief Contains common functions and functors used to solve
///        structural analysis.

# include "core.hpp"

namespace eva { namespace sa {

// In core_sa ci va:
// specializzazione per stiffness_submatrices_assembler e
// known_terms_assembler da usare con frame e truss.
// (Quindi anche le varie assemble element matrix..)


//------------------------------------- Problem Assembling ---------------------------------------//

/// Functor that assembles the system stiffness matrix.
template <typename AlgebraType> 
struct stiffness_submatrices_assembler;


/// Assembles the element (local) stiffness matrix. It is a
/// conveniency for calling the proper element_matrix_assembler
/// functor without explicitly specifying template parameters.
template <typename Structure> 
auto
assemble_element_matrix(const typename Structure::edge_descriptor& e,
                        const Structure& s);

/// Functor that assembles the element (local) stiffness matrix. 
/// It has been specialized for truss and frame kind types.
template <typename StructureKind> 
struct element_matrix_assembler;


/// Counts number of non-zero entries of the stiffness matrix
/// submatrices K_ff, K_fb, K_bb
template <typename S>
std::array<std::vector<size_t>, 3>
count_nnz_entries(const S& s, const std::vector<index_t>& dofmap,
                  const size_t nf, const size_t nb);




//-------------------------------------- Post-processing -----------------------------------------//
// // Get member internal forces
// template <typename Structure> 
// auto
// get_internal_forces(const typename Structure::edge_descriptor& e, const Structure& s, 
//                     const std::vector<real>& u, const std::vector<real>& f);

// /// Functor that computes member internal forces.
// /// Has to be specialized for every structure type.
// template <typename StructureKind> 
// struct internal_forces_getter;

}} //end namespace eva and sa

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
                                 //  IMPLEMENTATION  //
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
# include "../src/core_sa.tcc"



# endif //__EVA_STRUCTURAL_ANALYSIS__