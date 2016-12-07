/// @file  core_sa.tcc
/// @brief Contains definitions of template functions and class methods
///        declared in file core_sa.hpp .

//-------------------------------------- Problem Solving -----------------------------------------//
template <typename Params, typename Kind, typename S>
auto
solve(const S& s, Params p) 
{
    // If Kind is void use Structure kind, else directly use Kind
    using kind_t = typename std::conditional< 
        std::is_void<Kind>::value, 
        typename kind_of<S>::type, Kind
        >::type;

    auto r = solver<kind_t, Params>()(s);
    const auto& uu = std::get<0>(r);
    const auto& ff = std::get<1>(r);
    
    return assemble_results<kind_t>(uu, ff, s);
    // return solver<kind_t, Params>()(s);
}

template <typename StructureKind, typename Params>
struct solver 
{    
    using algebra_t = typename Params::algebra_t;
    using solver_t  = typename Params::solver_t;
    
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
        
        // Build (global) DOF map with BCs related nodes located at the back
        auto dofmap = std::vector<index_t>();
        size_t n_f, n_b;
        std::tie(dofmap, n_f, n_b) = build_global_dofmap(s);
        
        // Assemble system stiffness matrices (in global coordinates)
        // NOTE: use .sparseView()?
        auto matrices = assemble_stiffness_submatrices<algebra_t>(s, dofmap, n_f, n_b);
        const auto& K_ff = std::get<0>(matrices);
        const auto& K_fb = std::get<1>(matrices);
        const auto& K_bb = std::get<2>(matrices);
        
        // Assemble known terms (in global coordinates)
        auto known_terms = assemble_known_terms(s, n_f, n_b);
        const dense_vector& f_f = std::get<0>(known_terms); // displacements on BC nodes
        const dense_vector& u_b = std::get<1>(known_terms); // applied loads

        // Solve condensed system:
        // K_ff * u_f = f_f - K_fb * u_b
        solver_t eigen_solver(K_ff);
        dense_vector rhs = f_f - K_fb * u_b;
        dense_vector u_f = eigen_solver.solve(rhs);
        
        // Compute reactions on BC DOF:
        // f_b = K_bf * u_f + K_bb * u_b
        dense_vector f_b = K_fb.transpose() * u_f + K_bb * u_b;
        
        // Asemble displacement and force vectors 
        // using the original DOF ordering
        dense_vector uu = merge_and_reorder(u_f, u_b, dofmap);
        dense_vector ff = merge_and_reorder(f_f, f_b, dofmap);

        return std::make_pair(std::move(uu), std::move(ff));
        // return assemble_results<StructureKind>(uu, ff, s);
    }
};


//------------------------------------- Problem Assembling ---------------------------------------//
template <typename S> 
auto
assemble_element_matrix(const typename S::edge_descriptor& e, const S& s) 
{
    return element_matrix_assembler<typename kind_of<S>::type>()(e, s);
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
    return known_terms_assembler<typename kind_of<S>::type>()(s, n_f, n_b);
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
assemble_stiffness_submatrices(const S& s, const std::vector<index_t>& dofmap,
                          const size_t n_f, const size_t n_b) 
{
    return stiffness_submatrices_assembler<A>()(s, dofmap, n_f, n_b);
}

template <typename AlgebraKind>
struct stiffness_submatrices_assembler
{
    template <typename S>
    void
    operator()(const S& s,
               const std::vector<index_t>& dofmap,
               const size_t n_f, const size_t n_b)
    {
        static_assert(
            std::is_same<S,S>::value, 
            "Algebra type must be either dense_algebra_t or sparse_algebra_t"
            );
    }
};

template <>
struct stiffness_submatrices_assembler<dense_algebra_t>
{
    template <typename S>
    std::array<dense_matrix, 3>
    operator()(const S& s,
               const std::vector<index_t>& dofmap,
               const size_t n_f, const size_t n_b)
    {
        constexpr size_t dim = kind_of<S>::type::ndof;
        
        // Init SYSTEM SUB-matrices
        dense_matrix K_ff = dense_matrix::Zero(n_f, n_f);
        dense_matrix K_fb = dense_matrix::Zero(n_f, n_b);
        dense_matrix K_bb = dense_matrix::Zero(n_b, n_b);
        
        // Loop over all edges
        for (auto&& e : make_iterator_range(edges(s)))
        {
            // Get node indexes
            auto na = source(e, s);
            auto nb = target(e, s);
            
            // Build local to global DOF map
            auto loc_to_glob = build_local_to_global_dofmap<dim>(na, nb, dofmap);
                    
            // Get ELEMENT stiffness matrix in GLOBAL coordinates
            auto K_e = assemble_element_matrix(e, s);
            
            // Compose SYSTEM submatrices
            for (size_t i = 0; i < 2*dim; ++i)
                for (size_t j = 0; j < 2*dim; ++j)
                {    
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
# pragma clang diagnostic push
# pragma clang diagnostic ignored "-Wmissing-braces"
        return { std::move(K_ff), std::move(K_fb), std::move(K_bb) };
# pragma clang diagnostic pop
    }
};

template <>
struct stiffness_submatrices_assembler<sparse_algebra_t>
{
    template <typename S>
    std::array<sparse_matrix, 3>
    operator()(const S& s, const std::vector<index_t>& dofmap,
                             const size_t n_f, const size_t n_b)
    {
        return triplets_impl(s, dofmap, n_f, n_b);
        // return matrices_impl(s, dofmap, n_f, n_b);
    }
    
    template <typename S>
    std::array<sparse_matrix, 3>
    triplets_impl(const S& s, const std::vector<index_t>& dofmap,
                 const size_t n_f, const size_t n_b)
    {
        constexpr size_t dim = kind_of<S>::type::ndof; 
        
        // Compute number of non-zero entries for each matrix
        auto nnzs = count_nnz_entries(s, dofmap, n_f, n_b);

        auto nnz_tot = std::array<size_t, 3> ();
        for (size_t i = 0u; i < 3u; ++i)
        {
            nnz_tot[i] = std::accumulate(begin(nnzs[i]), end(nnzs[i]), 0u);
        }

        // Allocate Triplet lists reserving the required space
        using triplet_list_t = std::vector<Eigen::Triplet<real>>;
        auto K_ff_list = triplet_list_t(nnz_tot[0]); 
        auto K_fb_list = triplet_list_t(nnz_tot[1]);
        auto K_bb_list = triplet_list_t(nnz_tot[2]);
        
        // Loop over all edges
        for (auto&& e : make_iterator_range(edges(s)))
        {    
            // Get node indexes
            auto na = source(e, s);
            auto nb = target(e, s);
            
            // Build local to global DOF map
            auto loc_to_glob = build_local_to_global_dofmap<dim>(na, nb, dofmap);
                    
            // Get ELEMENT stiffness matrix in GLOBAL coordinates
            auto K_e = assemble_element_matrix(e, s);
        
            // Compose SYSTEM submatrices
            for (size_t i = 0; i < 2*dim; ++i)
                for (size_t j = 0; j < 2*dim; ++j)
                {                    
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
# pragma clang diagnostic push
# pragma clang diagnostic ignored "-Wmissing-braces"
        // Gather lists
        auto triplet_lists = std::array<triplet_list_t, 3> {K_ff_list, K_fb_list, K_bb_list};

        // Allocate empty submatrices reserving the required space
        auto matrices = std::array<sparse_matrix, 3> {
            sparse_matrix(n_f, n_f),
            sparse_matrix(n_f, n_b),
            sparse_matrix(n_b, n_b)};
# pragma clang diagnostic pop
        
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
        constexpr size_t dim = kind_of<S>::type::ndof; 
            
        // Compute number of non-zero entries for each matrix
        auto nnzs = count_nnz_entries(s, dofmap, n_f, n_b);
        
        // Init submatrices
# pragma clang diagnostic push
# pragma clang diagnostic ignored "-Wmissing-braces"
        auto matrices = std::array<sparse_matrix, 3> {
            sparse_matrix(n_f, n_f), 
            sparse_matrix(n_f, n_b), 
            sparse_matrix(n_b, n_b) }; 
# pragma clang diagnostic pop
        
        // Reserve exactly the required space
        for (size_t idx = 0u; idx < dim; ++idx) 
            matrices[idx].reserve(nnzs[idx]);
            
        auto& K_ff = std::get<0>(matrices); 
        auto& K_fb = std::get<1>(matrices);
        auto& K_bb = std::get<2>(matrices);
        
        // Loop over all edges
        for (auto&& e : make_iterator_range(edges(s)))
        {    
            // Get node indexes
            auto na = source(e, s);
            auto nb = target(e, s);
            
            // Build local to global DOF map
            auto loc_to_glob = build_local_to_global_dofmap<dim>(na, nb, dofmap);
                    
            // Get ELEMENT stiffness matrix in GLOBAL coordinates
            auto K_e = assemble_element_matrix(e, s);
        
            // Compose SYSTEM submatrices
            for (size_t i = 0; i < 2*dim; ++i)
                for (size_t j = 0; j < 2*dim; ++j)
                {    
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


template <typename S>
std::array<std::vector<size_t>, 3>
count_nnz_entries(const S& s, const std::vector<index_t>& dofmap,
                  const size_t n_f, const size_t n_b) 
{    
    // Aux vars
    constexpr size_t dim = kind_of<S>::type::ndof;
    
    // Init return variables
// # pragma clang diagnostic push
// # pragma clang diagnostic ignored "-Wmissing-braces"
//     auto ret = std::array<std::vector<size_t>, 3> {
//         std::vector<size_t>(n_f),
//         std::vector<size_t>(n_f), 
//         std::vector<size_t>(n_b)};
// # pragma clang diagnostic pop

    auto nnz_ff = std::vector<size_t>(n_f); // #(non-zero entries) of K_ff rows
    auto nnz_fb = std::vector<size_t>(n_f); // #(non-zero entries) of K_fb rows
    auto nnz_bb = std::vector<size_t>(n_b); // #(non-zero entries) of K_bb rows
    
    // For each vertex
    for (auto&& v : boost::make_iterator_range(vertices(s)))
    {
        // #(non-zero entries) associated to vertex DOF 
        size_t nr_f = 0u;
        size_t nr_b = 0u;
        
        // Check if DOF is free or restrained (BC)
        auto dofs = s[v].bcs;
        // Loop on current vertex DOF
        for (size_t i = 0u; i < dim; ++i)
        {
            // If nan =>  BC DOF  => increase bc non-zero entries count 
            if (std::isnan(dofs[i])) ++nr_b;
            //  Else  => free DOF => increase free non-zero entries count
            else                    ++nr_f;
        }
        
        // Loop on neighbour nodes and do the same thing for each one of them
        for (auto&& nn : make_iterator_range(adjacent_vertices(v, s)))
        {
            // Check if DOF is free or restrained (BC) 
            auto bcs = s[nn].bcs;
            // Loop on bcs
            for (size_t j = 0u; j < dim; ++j)
            {
                // If nan =>  BC DOF  => increase bc non-zero entries count 
                if (std::isnan(bcs[j])) ++nr_b;
                //  Else  => free DOF => increase free non-zero entries count
                else                    ++nr_f;
            }
        }
        
        // Now for each DOF belonging to the current vertex
        for (size_t dof = 0u; dof < dim; ++dof)
        {
            auto ii = dofmap[dof];
            if (ii < n_f) {
                // Write #(nnz entries) of K_ff and K_fb
                nnz_ff[ii] = nr_f;
                nnz_fb[ii] = nr_b;
                
            }
            else {
                // Write #(nnz entries) of K_bb and K_bf
                nnz_bb[ii-n_f] = nr_b;
                //nnz_bf[ii-nf] = nr_f;
            }
        }
    }
# pragma clang diagnostic push
# pragma clang diagnostic ignored "-Wmissing-braces"
    // #(non-zero entries) of (K_ff, K_fb, K_bb)
    return {std::move(nnz_ff), std::move(nnz_fb), std::move(nnz_bb)};
# pragma clang diagnostic pop
}


//------------------------------------- Results Assembling ---------------------------------------//
template <typename StructureKind, typename Structure>
std::vector< result<StructureKind> >
assemble_results(
    const dense_vector& u,
    const dense_vector& f,
    const Structure& s)
{
    // Aux vars
    constexpr static int n_loc_dof = StructureKind::ndof;
    const auto n_nodes = num_vertices(s);
    
    // Pre-allocate results
    using result_t = result<StructureKind>;
    auto results = std::vector<result_t> ();
    results.reserve(n_nodes);

    // For each joint (node)
    // for(size_t node_it = 0u; node_it < n_nodes*n_loc_dof; node_it += n_loc_dof)
    for (size_t node_it = 0u; node_it < n_nodes; ++node_it)
    {
        // Get current node starting index in u and f
        size_t node_start_idx = node_it * n_loc_dof;
        
        // Assemble the corresponding results
        results.emplace_back(result_t(u.segment<n_loc_dof>(node_start_idx),
                                      f.segment<n_loc_dof>(node_start_idx),
                                      s[node_it]));
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
    return internal_forces_getter<typename kind_of<S>::type>()(e, s, u, f);
}
