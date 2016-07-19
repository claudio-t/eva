# ifndef __EVA_STRUCTURE__
# define __EVA_STRUCTURE__

/**
 * @file structure.hpp
 * @brief Contains structure definition
 */

// std
# include <array>
/* remove me! */ # include <chrono> /* remove me! */
// boost
# include <boost/graph/adjacency_list.hpp>
# include <boost/graph/graph_traits.hpp>
# include <boost/range/iterator_range.hpp>
// eigen
# include <eigen3/Eigen/Dense>
# include <eigen3/Eigen/Sparse>

// eva
//~ # include "num_array.hpp"


namespace eva {
    
//######################################## Declarations ############################################
    
//------------------------------------- Aliases and Typedefs -------------------------------------//

/// Alias for type representing a real number
typedef double real; //using double = real;

/// Alias for a type representing an index
typedef size_t index_t;

/// Alias for a fixed size vector
template <int N, typename T = real>
using fixed_vector = Eigen::Matrix<T, N, 1>;

/// Alias for a fixed size matrix
template <int M, int N, typename T = real>
using fixed_matrix = Eigen::Matrix<T, M, N>;

/// Alias for a dense matrix used for computation
typedef Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> dense_matrix;

/// Alias for a dense vector used for computation
typedef Eigen::Matrix<real, Eigen::Dynamic, 1> dense_vector;

/// Alias for a sparse matrix used for computation
typedef Eigen::SparseMatrix<real> sparse_matrix; // Does not work with CRS!

/// Alias for a sparse vector used for computation (CCS format)
typedef Eigen::SparseVector<real> sparse_vector;

/// Alias for a generic structure
template < typename NodeProps, typename EdgeProps, typename Kind >
using generic_structure = boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS,
                                                 NodeProps, EdgeProps, Kind, boost::vecS >;

//----------------------------------- Functors and Functions -------------------------------------//

// -- Post-processing --//
template <typename StructureType> struct internal_forces_getter;

template <typename Structure> 
decltype(auto)
get_internal_forces(const typename Structure::edge_descriptor& e, const Structure& s, 
                    const std::vector<real>& u, const std::vector<real>& f);

// -- Problem Solving -- //
/// Generic solver functor: has to be specialized for every structure type
template <typename StructureType> struct solver;

/// Convenient wrapper for the solver functor
template <typename Structure> decltype(auto) solve(const Structure& s);


// -- Problem Assembling -- //
/// Builds the system stiffness matrix in global coordinates
template <typename S>
std::array<dense_matrix, 3>
assemble_stiffness_matrix(const S& s, const std::vector<index_t>& dofmap,
                          const size_t n_f, const size_t n_b);

/// Builds the known terms, i.e. the force vector portion related to the free DOF and
/// the displacement vector portion associated to BC DOF
template <typename S>
std::array<dense_vector, 2>
assemble_known_terms(const S& s, const size_t n_f, const size_t n_b);

/// Assembles the element (local) stiffness matrix. 
/// Has to be specialized for every structure type
template <typename StructureType> struct element_matrix_assembler;

/// Convenient wrapper for the element_matrix_assembler functor
template <typename Structure> 
decltype(auto)
assemble_element_matrix(const typename Structure::edge_descriptor& e, const Structure& s);

/// Assembles the known terms vectors (source term and BCs related DOF).
/// Has to be specialized for every structure type
template <typename StructureType> struct known_terms_assembler;


// -- DOF handling -- //
/// Builds a DOF map where the positions of DOF associated 
/// to BCs are located at the end
template <typename S>
std::tuple<std::vector<index_t>, size_t, size_t>
build_global_dofmap(const S& s);

/// Builds the (element) local to global DOF map
template <int N> 
fixed_vector<2*N, index_t> 
build_local_to_global_dofmap(index_t na, index_t nb, const std::vector<index_t>& dm);

/// Count number of non-zero entries of the stiffness matrix submatrices K_ff, K_fb, K_bb
template <typename S>
std::array<std::vector<size_t>, 3>
count_nnz_entries(const S& s, const std::vector<index_t>& dofmap,
                  const size_t nf, const size_t nb);

/// Takes two vectors -one containing a quantity on free DOF, the other on BC DOF-
/// and merge them together reordering DOF accordingly to the supplied DOF map
std::vector<real> 
merge_and_reorder(const dense_vector& v_f, const dense_vector& v_b, 
                  const std::vector<index_t>& dofmap);
                  
std::vector<real> 
merge_and_reorder(const sparse_vector& v_f, const sparse_vector& v_b, 
                  const std::vector<index_t>& dofmap);


//######################################## Definitions #############################################

//----------------------------------- Functors and Functions -------------------------------------//
      
// -- Post-processing --//
template <typename StructureType> 
struct internal_forces_getter 
{
    template <typename S>
    void 
    operator()(const typename S::edge_descriptor& e, const S& s,
               const std::vector<real>& u, const std::vector<real>& f) 
    { 
        static_assert(!std::is_same<S,S>::value, 
                      "Method not implemented");
    }
};

template <typename S> 
decltype(auto)
get_internal_forces(const typename S::edge_descriptor& e, const S& s, 
                    const std::vector<real>& u, const std::vector<real>& f) 
{
    return internal_forces_getter<typename S::graph_bundled>()(e, s, u, f);
}

      
// -- Problem Solving -- //
template <typename StructureType>
struct solver 
{    
    template <class S>
    decltype(auto) operator()(const S& s) 
    {
        // Assemble system: K*u = f
        // 
        // | K_ff  K_fb |*| u_f |  =  | f_f | 
        // | K_bf  K_bb | | u_b |     | f_b |
        //
        // K := Stiffness matrix
        // u .= displacement vector
        // f := load/reaction vector
        //
        // $_f := Free DOF
        // $_b := BC   DOF
        
        // Build (global) DOF map with BCs related nodes located to the back
        auto   dofmap = std::vector<index_t>();
        size_t n_f, n_b;
        std::tie(dofmap, n_f, n_b) = build_global_dofmap(s);
        
        // Timing
        auto start1 = std::chrono::high_resolution_clock::now();
        //
        
        // Assemble system stiffness matrices (in global coordinates)
        auto matrices = assemble_stiffness_matrix(s, dofmap, n_f, n_b);
        const auto& K_ff = std::get<0>(matrices);
        const auto& K_fb = std::get<1>(matrices);
        const auto& K_bb = std::get<2>(matrices);
        
        //~ std::cout << std::endl
                  //~ << "------------------------ K_ff ------------------------"
                  //~ << std::endl
                  //~ << K_ff
                  //~ << std::endl
                  //~ << std::endl;
        //~ const sparse_matrix K_ff = std::get<0>(matrices).sparseView();
        //~ const sparse_matrix K_fb = std::get<1>(matrices).sparseView();
        //~ const sparse_matrix K_bb = std::get<2>(matrices).sparseView();
        
        // Assemble known terms (in global coordinates)
        auto known_terms = assemble_known_terms(s, n_f, n_b);
        const auto& f_f = std::get<0>(known_terms); // displacements on BC nodes
        const auto& u_b = std::get<1>(known_terms); // applied loads
        
        // Timing
        auto end1 = std::chrono::high_resolution_clock::now();
        auto elapsed_time1 = std::chrono::duration<double,std::micro>(end1-start1);
        std::cout << "Assembling time = " << elapsed_time1.count() << std::endl;
        //
        
        // Timing
        auto start2 = std::chrono::high_resolution_clock::now();
        //
        // Solve condensed system:
        // K_ff * u_f = f_f - K_fb * u_b
        Eigen::LDLT<dense_matrix> eigen_solver(K_ff);        
        //~ Eigen::SimplicialLLT<sparse_matrix> eigen_solver(K_ff);
        //~ Eigen::ConjugateGradient<sparse_matrix> eigen_solver(K_ff);
        auto u_f = eigen_solver.solve(f_f - K_fb * u_b);
        
        // Timing
        auto end2 = std::chrono::high_resolution_clock::now();
        auto elapsed_time2 = std::chrono::duration<double,std::micro>(end2-start2);
        std::cout << "Solving time = " << elapsed_time2.count() << std::endl;
        //
        
        // Compute reactions on BC DOF:
        // f_b = K_bf * u_f + K_bb * u_b
        auto f_b = K_fb.transpose() * u_f + K_bb * u_b;
        
        // Asemble displacement and force vectors 
        // using the original DOF numbering
        auto uu = merge_and_reorder(u_f, u_b, dofmap);
        auto ff = merge_and_reorder(f_f, f_b, dofmap);
        
        return std::make_pair(std::move(uu), std::move(ff));
    }
};


template <typename Structure> 
decltype(auto)
solve(const Structure& s) 
{
    return solver<typename Structure::graph_bundled>()(s);
}



// -- Problem Assembling -- //
template <typename StructureType> 
struct element_matrix_assembler 
{
    template <typename S>
    void operator()(const typename S::edge_descriptor& e, const S& s) 
    { 
        static_assert(!std::is_same<S,S>::value, 
                      "Assemble method not implemented");
    }
};

template <typename S> 
decltype(auto)
assemble_element_matrix(const typename S::edge_descriptor& e, const S& s) 
{
    return element_matrix_assembler<typename S::graph_bundled>()(e, s);
}


template <typename StructureType>
struct known_terms_assembler 
{
    template <typename S>
    void
    operator()(const S& s, const size_t n_f, const size_t n_b) 
    {
        static_assert(!std::is_same<S,S>::value, 
                      "Assemble method not implemented");
    }
};

template <typename S>
std::array<dense_vector, 2>
assemble_known_terms(const S& s, const size_t n_f, const size_t n_b) 
{
    return known_terms_assembler<typename S::graph_bundled>()(s, n_f, n_b);
}


template <typename S>
std::array<dense_matrix, 3>
assemble_stiffness_matrix(const S& s, const std::vector<index_t>& dofmap,
                          const size_t n_f, const size_t n_b) 
{
    const static size_t dim = S::graph_bundled::ndof;
    
    // Init SYSTEM SUB-matrices
    auto matrices = std::array<dense_matrix, 3> {dense_matrix::Zero(n_f, n_f),
                                                 dense_matrix::Zero(n_f, n_b),
                                                 dense_matrix::Zero(n_b, n_b)};
    
    auto& K_ff = std::get<0>(matrices);
    auto& K_fb = std::get<1>(matrices);
    auto& K_bb = std::get<2>(matrices);
    
    // Loop over all edges
    for (auto&& e : make_iterator_range(edges(s))) {
        
        // Get node indexes
        auto na = source(e, s);
        auto nb = target(e, s);
        
        // Build local to global DOF map
        auto loc_to_glob = build_local_to_global_dofmap<dim>(na, nb, dofmap);
                
        // Get ELEMENT stiffness matrix in GLOBAL coordinates
        auto K_e = assemble_element_matrix(e, s);
        
        // Compose SYSTEM submatrices
        for (size_t i = 0; i < 2*dim; ++i)
            for (size_t j = 0; j < 2*dim; ++j) {
                
                auto ii = loc_to_glob(i);
                auto jj = loc_to_glob(j);
                
                switch((ii < n_f) + 2*(jj < n_f)) {
                    
                    case 0 : K_bb(ii-n_f, jj-n_f) += K_e(i, j); break;
                    
                    case 1 : K_fb(ii, jj-n_f) += K_e(i, j);     break;
                    
                    case 2 : /* Does nothing (no K_bf)  */      break;
                    
                    case 3 : K_ff(ii, jj) += K_e(i, j);         break;
                }
            }
    }
    // FIXME: Check for NRVO! (should work)
    return matrices;
}


template <typename S>
std::array<sparse_matrix, 3>
assemble_stiffness_matrix_sp(const S& s, const std::vector<index_t>& dofmap,
                             const size_t n_f, const size_t n_b) 
{    
    const static size_t dim = S::graph_bundled::ndof; 
        
    // Compute number of non-zero entries for each matrix
    auto nnzs = count_nnz_entries(s, dofmap, n_f, n_b);
    const auto& nnz_ff = std::get<0>(nnzs);
    const auto& nnz_fb = std::get<1>(nnzs);
    const auto& nnz_bb = std::get<2>(nnzs);
    
    // Init submatrices
    auto matrices = std::array<sparse_matrix, 3> { sparse_matrix(n_f, n_f), 
                                                   sparse_matrix(n_f, n_b), 
                                                   sparse_matrix(n_b, n_b) }; 
    // Reserve exactly the required space
    for (size_t idx = 0u; idx < dim; ++idx) matrices[idx].reserve(nnzs[idx]);
        
    auto& K_ff = std::get<0>(matrices); 
    auto& K_fb = std::get<1>(matrices);
    auto& K_bb = std::get<2>(matrices);
    
    // Loop over all edges
    for (auto&& e : make_iterator_range(edges(s))) {
        
        // Get node indexes
        auto na = source(e, s);
        auto nb = target(e, s);
        
        // Build local to global DOF map
        auto loc_to_glob = build_local_to_global_dofmap<dim>(na, nb, dofmap);
                
        // Get ELEMENT stiffness matrix in GLOBAL coordinates
        auto K_e = assemble_element_matrix(e, s);
    
        // Compose SYSTEM submatrices
        for (size_t i = 0; i < 2*dim; ++i)
            for (size_t j = 0; j < 2*dim; ++j) {
                
                auto ii = loc_to_glob(i);
                auto jj = loc_to_glob(j);
                
                switch((ii < n_f) + 2*(jj < n_f)) {
                    
                    case 0: K_bb.coeffRef(ii-n_f, jj-n_f) += K_e(i, j); break;
                    
                    case 1: K_fb.coeffRef(ii, jj-n_f) += K_e(i, j);     break;
                    
                    case 2: /* Does nothing (no K_bf)  */               break;
        
                    case 3: K_ff.coeffRef(ii, jj) += K_e(i, j);         break;
                }
            }
    }
    // Compress, then return matrices (K_ff, K_fb, K_bb)
    for (size_t idx = 0u; idx < dim; ++idx) matrices[idx].makeCompressed();
    
    return matrices;
}
template <typename S>
std::array<sparse_matrix, 3>
assemble_stiffness_matrix_sp2(const S& s, const std::vector<index_t>& dofmap,
                              const size_t n_f, const size_t n_b) 
{
    const static size_t dim = S::graph_bundled::ndof; 
    
    // Compute number of non-zero entries for each matrix
    auto nnzs = count_nnz_entries(s, dofmap, n_f, n_b);
    
    std::array<size_t, 3> nnz_tot {std::accumulate(begin(nnzs[0]), end(nnzs[0]), 0u),
                                   std::accumulate(begin(nnzs[1]), end(nnzs[1]), 0u), 
                                   std::accumulate(begin(nnzs[2]), end(nnzs[2]), 0u)};
    
    // Allocate Triplet lists reserving the required space
    auto triplet_lists = std::array<std::vector<Eigen::Triplet<real>>, 3> ();
    for (size_t idx = 0u; idx < 3u; ++idx) 
        triplet_lists[idx].reserve(nnz_tot[idx]);
        
    auto& K_ff_list = std::get<0>(triplet_lists); 
    auto& K_fb_list = std::get<1>(triplet_lists);
    auto& K_bb_list = std::get<2>(triplet_lists);
    
    
    // Loop over all edges
    for (auto&& e : make_iterator_range(edges(s))) {
        
        // Get node indexes
        auto na = source(e, s);
        auto nb = target(e, s);
        
        // Build local to global DOF map
        auto loc_to_glob = build_local_to_global_dofmap<dim>(na, nb, dofmap);
                
        // Get ELEMENT stiffness matrix in GLOBAL coordinates
        auto K_e = assemble_element_matrix(e, s);
    
        // Compose SYSTEM submatrices
        for (size_t i = 0; i < 2*dim; ++i)
            for (size_t j = 0; j < 2*dim; ++j) {
                
                auto ii = loc_to_glob(i);
                auto jj = loc_to_glob(j);
                
                switch((ii < n_f) + 2*(jj < n_f)) {
                    
                    case 0 /*K_bb*/: K_bb_list.emplace_back(ii-n_f, jj-n_f, K_e(i, j)); break;
                    
                    case 1 /*K_fb*/: K_fb_list.emplace_back(ii, jj-n_f, K_e(i, j));     break;
                    
                    case 2 /*K_bf*/: /* Does nothing (no K_bf)  */                      break;
        
                    case 3 /*K_ff*/: K_ff_list.emplace_back(ii, jj, K_e(i, j));         break;
                }
            }
    }
    
    // Allocate empty submatrices reserving the required space
    auto matrices = std::array<sparse_matrix, 3> { sparse_matrix(n_f, n_f), 
                                                   sparse_matrix(n_f, n_b), 
                                                   sparse_matrix(n_b, n_b) }; 

    
    for (size_t idx = 0u; idx < dim; ++idx)
        matrices[idx].setFromTriplets(triplet_lists[idx].begin(), triplet_lists[idx].end());
    
    // Compress, then return matrices (K_ff, K_fb, K_bb)
    for (size_t idx = 0u; idx < dim; ++idx) matrices[idx].makeCompressed();
    
    return matrices;
}


 
// -- DOF Handling -- //
template <typename S>
std::tuple<std::vector<index_t>, size_t, size_t>
build_global_dofmap(const S& s) 
{
    const static size_t dim = S::graph_bundled::ndof;
    
    // Initialize the DOF map
    size_t n_verts = num_vertices(s);
    auto  dofmap   = std::vector<index_t> (dim*n_verts);
    
    size_t n_bcs  = 0u;                 // #(BC DOF)
    size_t bc_pos = dofmap.size()-1u;   //
    size_t ff_pos = 0u;
    
    for (size_t v = 0u; v < n_verts; ++v) {
        // Get BC field
        const auto& bcs = s[v].bcs;
        
        for (size_t i = 0u; i < dim; ++i) {
            // If coordinate = NaN => Free DOF else BC DOF
            if (std::isnan(bcs[i]))
                // Put DOF at the beginning
                dofmap[dim*v+i] = ff_pos++;
            else {
                // Put DOF at the end and increment counter
                dofmap[dim*v+i] = bc_pos--;
                ++n_bcs;
            }
        }
    }
    // Return (dofmap, #(free DOF), #(BC DOF))
    return std::make_tuple(dofmap, dofmap.size()-n_bcs, n_bcs);
}


template <typename S>
std::array<std::vector<size_t>, 3>
count_nnz_entries(const S& s, const std::vector<index_t>& dofmap,
                  const size_t n_f, const size_t n_b) 
{    
    // Aux vars
    const static size_t dim = S::graph_bundled::ndof;
    
    // Init return variables
    auto ret = std::array<std::vector<size_t>, 3> {std::vector<size_t>(n_f),
                                                   std::vector<size_t>(n_f), 
                                                   std::vector<size_t>(n_b)};
    auto& nnz_ff = std::get<0>(ret);    // #(non-zero entries) of K_ff rows
    auto& nnz_fb = std::get<1>(ret);    // #(non-zero entries) of K_fb rows
    auto& nnz_bb = std::get<2>(ret);    // #(non-zero entries) of K_bb rows
    
    // For each vertex
    for (auto&& v : make_iterator_range(vertices(s))) {
        // #(non-zero entries) associated to vertex DOF 
        size_t nr_f = 0u;
        size_t nr_b = 0u;
        
        // Check if DOF is free or restrained (BC)
        auto dofs = s[v].bcs;
        // Loop on current vertex DOF
        for (size_t i = 0u; i < dim; ++i) {
            // If nan =>  BC DOF  => increase bc non-zero entries count 
            if(std::isnan(dofs[i])) ++nr_b;
            //  Else  => free DOF => increase free non-zero entries count
            else                   ++nr_f;
        }
        
        // Loop on neighbour nodes and do the same thing for each one of them
        for (auto&& nn : make_iterator_range(adjacent_vertices(v, s))) {
            // Check if DOF is free or restrained (BC) 
            auto bcs = s[nn].bcs;
            // Loop on bcs
            for (size_t j = 0u; j < dim; ++j) {
                // If nan =>  BC DOF  => increase bc non-zero entries count 
                if(std::isnan(bcs[j])) ++nr_b;
                //  Else  => free DOF => increase free non-zero entries count
                else                   ++nr_f;
            }
        }
        
        // Now for each DOF belonging to the current vertex
        for (size_t dof = 0u; dof < dim; ++dof) {
            auto ii = dofmap[dof];
            if (ii < n_f) {
                // Write #(nnz entries) of K_ff and K_fb
                nnz_ff[ii] = nr_f;
                nnz_fb[ii] = nr_b;
                
            } else {
                // Write #(nnz entries) of K_ff and K_fb
                nnz_bb[ii-n_f] = nr_b;
                //nnz_bf[ii-nf] = nr_f;
            }
        }
    }
    // #(non-zero entries) of (K_ff, K_fb, K_bb)
    return ret;
}



template <int N> 
fixed_vector<2*N, index_t> 
build_local_to_global_dofmap(index_t na, index_t nb, const std::vector<index_t>& dm) 
{    
    auto ret = fixed_vector<2*N, index_t> ();
    
    for (auto i = 0; i < N; ++i) {
        ret(i)   = dm[N*na + i];    // Node A -> positions from 0 to  N-1
        ret(i+N) = dm[N*nb + i];    // Node B -> positions from N to 2N-1
    }
    
    return ret;
}


std::vector<real> 
merge_and_reorder(const dense_vector& v_f, 
                  const dense_vector& v_b, 
                  const std::vector<index_t>& dofmap) 
{    
    // DOF sizes
    size_t n_f = v_f.size();  // Free DOF
    size_t n_b = v_b.size();  // BC DOF
    size_t n_t = n_f  + n_b;  // Tot nr of DOF
    
    assert(dofmap.size() == n_t);
    
    auto rv = std::vector<real>(); rv.reserve(n_t);
    
    for (size_t i = 0u; i < n_t; ++i) {
        // Pick element to add from v_f if ii < n_f and from v_b otherwise
        auto ii = dofmap[i]; 
        rv.emplace_back( (ii < n_f) ? v_f(ii) : v_b(ii-n_f) );
    }
    return rv;
}


std::vector<real> 
merge_and_reorder(const sparse_vector& v_f, 
                  const sparse_vector& v_b, 
                  const std::vector<index_t>& dofmap) 
{    
    // DOF sizes
    size_t n_f = v_f.size();  // Free DOF
    size_t n_b = v_b.size();  // BC DOF
    size_t n_t = n_f  + n_b;  // Tot nr of DOF
    
    assert(dofmap.size() == n_t);
    
    // Auxiliary vector that just concatenates v_f and v_b
    auto auxv = std::vector<real>(n_t);
    
    // Fill auxiliary vector with free DOF values
    for (auto it = Eigen::SparseVector<double>::InnerIterator(v_f); it; ++it)
        auxv[it.index()] = it.value();
    
    // Fill auxiliary vector with BC DOF values
    for (auto it = Eigen::SparseVector<double>::InnerIterator(v_b); it; ++it)
        auxv[it.index() + n_f] = it.value();
    
    // Build return vector by reordering back DOF
    auto rv = std::vector<real>(); rv.reserve(n_t);
    for (size_t i = 0u; i < n_t; ++i) 
        rv.emplace_back(auxv[dofmap[i]]);
    
    return rv;
}

} // end namespace eva


# endif //__EVA_STRUCTURE__