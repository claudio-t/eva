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


namespace eva {
    
//######################################## Declarations ############################################
    
//------------------------------------- Aliases and Typedefs -------------------------------------//

/// Alias for type representing a real number
using real = double;

/// Alias for a type representing an index
using index_t = size_t;

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
template < typename NodeProps, typename EdgeProps, typename Kind >
using generic_structure = boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS,
                                                 NodeProps, EdgeProps, Kind, boost::vecS >;

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

/// Solves a given problem automatically deducing the structure type
template <
    typename Params = sparse_solver_params<>,
    typename Kind = void, typename Structure
    >
auto
solve(const Structure& s, Params p = Params());

/// Generic solver functor: has to be specialized for every kind of structure
template <typename StructureKind, typename Params>
struct solver;


//------------------------------------- Problem Assembling ---------------------------------------//
/// Assembles the system stiffness matrix in global coordinates.
template <typename A, typename S>
auto
assemble_stiffness_matrix(const S& s, const std::vector<index_t>& dofmap,
                          const size_t n_f, const size_t n_b);

/// Functor that assembles the system stiffness matrix. 
/// Has to be specialized for every structure type
template <typename AlgebraType> 
struct stiffness_matrix_assembler;

/// Builds the known terms, i.e. the force vector portion related to the free DOF and
/// the displacement vector portion associated to BC DOF
template <typename S>
std::array<dense_vector, 2>
assemble_known_terms(const S& s, const size_t n_f, const size_t n_b);

/// Functor that assembles the element (local) stiffness matrix. 
/// Has to be specialized for every kind of structure
template <typename StructureKind> 
struct element_matrix_assembler;

/// Assembles the element (local) stiffness matrix
template <typename Structure> 
auto
assemble_element_matrix(const typename Structure::edge_descriptor& e, const Structure& s);

/// Functor that assembles the known terms vectors (source term and BCs related DOF).
/// Has to be specialized for every kind of structure
template <typename StructureKind> 
struct known_terms_assembler;


//---------------------------------------- DOF handling ------------------------------------------//
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
dense_vector 
merge_and_reorder(const dense_vector& v_f, const dense_vector& v_b, 
                  const std::vector<index_t>& dofmap);


//------------------------------------- Results Assembling ---------------------------------------//
/// Represents the results of the solvin procedure applied to a
/// certain structure. Has to be specialized for every kind of structure.
template <typename StructureKind>
struct result;

/// Assembles the output of the solving procedure for a generic
/// given structure.
template <typename StructureKind>
std::vector<result<StructureKind>>
assemble_results(
    const dense_vector& u,
    const dense_vector& f,
    const std::vector<index_t>& dofmap);


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


//######################################## Definitions #############################################      
//-------------------------------------- Problem Solving -----------------------------------------//
template <typename Params, typename Kind, typename Structure>
auto
solve(const Structure& s, Params p) 
{
    // If Kind is void use Structure kind, else directly use Kind
    using kind_t = typename std::conditional< 
        std::is_void<Kind>::value, 
        typename Structure::graph_bundled, Kind
        >::type;
                    
    return solver<kind_t, Params>()(s);
}

template <typename A, typename S>
struct solver_params
{
    using algebra_t = A;
    using solver_t  = S;
};


template <typename StructureKind, typename Params>
struct solver 
{    
    using algebra_t = typename Params::algebra_t;
    using  solver_t = typename  Params::solver_t;
    
    template <class S>
    auto
    operator()(const S& s) 
    {
        // Assemble system: K*u = f
        // 
        // | K_ff  K_fb |*| u_f |  =  | f_f | 
        // | K_bf  K_bb | | u_b |     | f_b |
        //
        // K := Stiffness matrix
        // u := displacement vector
        // f := load/reaction vector
        //
        // $_f := Free DOF
        // $_b := BC   DOF
        
        // Build (global) DOF map with BCs related nodes located to the back
        auto dofmap = std::vector<index_t>();
        size_t n_f, n_b;
        std::tie(dofmap, n_f, n_b) = build_global_dofmap(s);
        
        // Timing
        auto start1 = std::chrono::high_resolution_clock::now();
        //
        
        // Assemble system stiffness matrices (in global coordinates)
        auto matrices = assemble_stiffness_matrix<algebra_t>(s, dofmap, n_f, n_b);
        const auto& K_ff = std::get<0>(matrices);
        const auto& K_fb = std::get<1>(matrices);
        const auto& K_bb = std::get<2>(matrices);
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
        solver_t eigen_solver(K_ff);
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
        // using the original DOF ordering
        auto uu = merge_and_reorder(u_f, u_b, dofmap);
        auto ff = merge_and_reorder(f_f, f_b, dofmap);

        // auto uu = dense_vector(n_f + n_b);
        // uu << std::move(u_f), std::move(u_b);
        // reorder(uu, dofmap, n_f);

        return assemble_results<StructureKind>(uu, ff, dofmap);
        // return std::make_pair(std::move(uu), std::move(ff));
    }
};


//------------------------------------- Problem Assembling ---------------------------------------//
template <typename S> 
auto
assemble_element_matrix(const typename S::edge_descriptor& e, const S& s) 
{
    return element_matrix_assembler<typename S::graph_bundled>()(e, s);
}

template <typename StructureKind> 
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
std::array<dense_vector, 2>
assemble_known_terms(const S& s, const size_t n_f, const size_t n_b) 
{
    return known_terms_assembler<typename S::graph_bundled>()(s, n_f, n_b);
}

template <typename StructureKind>
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


template <typename A, typename S>
auto
assemble_stiffness_matrix(const S& s, const std::vector<index_t>& dofmap,
                          const size_t n_f, const size_t n_b) 
{
    return stiffness_matrix_assembler<A>()(s, dofmap, n_f, n_b);
}

template <typename AlgebraKind>
struct stiffness_matrix_assembler
{
    template <typename S>
    void
    operator()(const S& s, const std::vector<index_t>& dofmap,
                             const size_t n_f, const size_t n_b)
    {
        static_assert(
            std::is_same<S,S>::value, 
            "Algebra type must be either dense_algebra_t or sparse_algebra_t"
        );
    }
};

template <>
struct stiffness_matrix_assembler<dense_algebra_t>
{
    template <typename S>
    std::array<dense_matrix, 3>
    operator()(const S& s, const std::vector<index_t>& dofmap,
                             const size_t n_f, const size_t n_b)
    {
        const static size_t dim = S::graph_bundled::ndof;
        
        // Init SYSTEM SUB-matrices
        dense_matrix K_ff = dense_matrix::Zero(n_f, n_f);
        dense_matrix K_fb = dense_matrix::Zero(n_f, n_b);
        dense_matrix K_bb = dense_matrix::Zero(n_b, n_b);
        
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
        return {std::move(K_ff), std::move(K_fb), std::move(K_bb)};
    }
};

template <>
struct stiffness_matrix_assembler<sparse_algebra_t>
{
    template <typename S>
    std::array<sparse_matrix, 3>
    operator()(const S& s, const std::vector<index_t>& dofmap,
                             const size_t n_f, const size_t n_b)
    {
        return triplets_impl(s, dofmap, n_f, n_b);
        //~ return matrices_impl(s, dofmap, n_f, n_b);
    }
    
    template <typename S>
    std::array<sparse_matrix, 3>
    triplets_impl(const S& s, const std::vector<index_t>& dofmap,
                 const size_t n_f, const size_t n_b)
    {
        const static size_t dim = S::graph_bundled::ndof; 
        
        // Compute number of non-zero entries for each matrix
        auto nnzs = count_nnz_entries(s, dofmap, n_f, n_b);

        auto nnz_tot = std::array<size_t, 3> ();
        for (size_t i = 0u; i < 3u; ++i)
            std::accumulate(begin(nnzs[i]), end(nnzs[i]), 0u);

        // Allocate Triplet lists reserving the required space
        using triplet_list_t = std::vector<Eigen::Triplet<real>>;
        auto triplet_lists = std::array<triplet_list_t, 3> ();
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
        auto matrices = std::array<sparse_matrix, 3> {
            sparse_matrix(n_f, n_f),
            sparse_matrix(n_f, n_b),
            sparse_matrix(n_b, n_b)};
        
        for (size_t idx = 0u; idx < dim; ++idx)
            matrices[idx].setFromTriplets(triplet_lists[idx].begin(), triplet_lists[idx].end());
        
        // Compress, then return matrices (K_ff, K_fb, K_bb)
        for (size_t idx = 0u; idx < dim; ++idx) matrices[idx].makeCompressed();
        
        return matrices;
    }
    
    template <typename S>
    std::array<sparse_matrix, 3>
    matrices_impl(const S& s, const std::vector<index_t>& dofmap,
                  const size_t n_f, const size_t n_b)
    {    
        const static size_t dim = S::graph_bundled::ndof; 
            
        // Compute number of non-zero entries for each matrix
        auto nnzs = count_nnz_entries(s, dofmap, n_f, n_b);
        
        // Init submatrices
        auto matrices = std::array<sparse_matrix, 3> {
            sparse_matrix(n_f, n_f), 
            sparse_matrix(n_f, n_b), 
            sparse_matrix(n_b, n_b) }; 
        // Reserve exactly the required space
        for (size_t idx = 0u; idx < dim; ++idx) 
            matrices[idx].reserve(nnzs[idx]);
            
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
};

 
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
            if (std::isnan(bcs[i])) {
                // Put DOF at the beginning
                dofmap[dim*v+i] = ff_pos++;
            }
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
    auto ret = std::array<std::vector<size_t>, 3> {
        std::vector<size_t>(n_f),
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
            if (std::isnan(dofs[i])) ++nr_b;
            //  Else  => free DOF => increase free non-zero entries count
            else                    ++nr_f;
        }
        
        // Loop on neighbour nodes and do the same thing for each one of them
        for (auto&& nn : make_iterator_range(adjacent_vertices(v, s))) {
            // Check if DOF is free or restrained (BC) 
            auto bcs = s[nn].bcs;
            // Loop on bcs
            for (size_t j = 0u; j < dim; ++j) {
                // If nan =>  BC DOF  => increase bc non-zero entries count 
                if (std::isnan(bcs[j])) ++nr_b;
                //  Else  => free DOF => increase free non-zero entries count
                else                    ++nr_f;
            }
        }
        
        // Now for each DOF belonging to the current vertex
        for (size_t dof = 0u; dof < dim; ++dof) {
            auto ii = dofmap[dof];
            if (ii < n_f) {
                // Write #(nnz entries) of K_ff and K_fb
                nnz_ff[ii] = nr_f;
                nnz_fb[ii] = nr_b;
                
            }
            else {
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



dense_vector
merge_and_reorder(const dense_vector& v_f, 
                  const dense_vector& v_b, 
                  const std::vector<index_t>& dofmap) 
{
    // DOF sizes
    size_t n_f = v_f.size();  // Free DOF
    size_t n_b = v_b.size();  // BC DOF
    size_t n_t = n_f + n_b;  // Tot nr of DOF

    // Pre-allocate return variable
    auto rv = dense_vector(n_t);

    // For each DOF
    for (size_t i = 0u; i < n_t; ++i) {
        // Retrieve mapped position
        auto ii = dofmap[i];
        // Fill position by picking the value from the proper vector
         rv[i] = (ii < n_f) ? v_f(ii) : v_b(ii-n_f);
    }
    return rv;
}


dense_vector
reorder(const dense_vector& v, const std::vector<index_t>& dofmap)
{
    // Pre-allocate results
    auto n_t = v.size();
    auto rv  = dense_vector(n_t);
    
    // Reorder
    for (size_t i = 0u; i < n_t; ++i)
        rv[i] = v[dofmap[i]];
    
    return rv;
}



//------------------------------------- Results Assembling ---------------------------------------//
template <typename StructureKind>
std::vector<result<StructureKind>>
assemble_results(
    const dense_vector& u,
    const dense_vector& f,
    const std::vector<index_t>& dofmap)
{
    // Aux vars
    constexpr static int n_loc_dof = StructureKind::ndof;
    auto n_nodes = dofmap.size() / n_loc_dof;
    
    // Pre-allocate results
    typedef result<StructureKind> result_t;
    auto results = std::vector<result_t> ();
    results.reserve(n_nodes);

    // For each joint/node
    for(size_t node_it = 0u; node_it < n_nodes*n_loc_dof; node_it += n_loc_dof) {
        // Assemble the corresponding results
        results.emplace_back(result_t(u.segment<n_loc_dof>(node_it),
                                      f.segment<n_loc_dof>(node_it)));
    }
    
    return results;
}



//-------------------------------------- Post-processing -----------------------------------------//
template <typename StructureKind> 
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
auto
get_internal_forces(const typename S::edge_descriptor& e, const S& s, 
                    const std::vector<real>& u, const std::vector<real>& f) 
{
    return internal_forces_getter<typename S::graph_bundled>()(e, s, u, f);
}



} // end namespace eva


# endif //__EVA_STRUCTURE__
