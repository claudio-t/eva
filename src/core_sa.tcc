/// @file  core_sa.tcc
/// @brief Contains definitions of template functions and class methods
///        declared in file core_sa.hpp .


namespace eva {

//------------------------------------- Problem Assembling ---------------------------------------//


template <typename S, typename Kind> 
auto
assemble_element_matrix(
    const typename S::edge_descriptor& e,
    const S& s, const Kind kind) 
{

    // using kind_t = typename std::conditional<
    //     std::is_same<typename kind_of<S>::type, Kind>::value,
        
    //     >::type;
    
    return element_matrix_assembler<Kind>()(e, s);
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




template <typename Kind>
struct stiffness_submatrices_assembler<Kind, dense_algebra_t>
{
    template <typename Structure>
    std::array<dense_matrix, 3>
    operator()(const Structure& s,
               const std::vector<index_t>& dofmap,
               const size_t n_f, const size_t n_b)
    {
        constexpr size_t dim = Kind::ndof;
        
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
            auto K_e = element_matrix_assembler<Kind>()(e, s);
            
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

template <typename Kind>
struct stiffness_submatrices_assembler<Kind, sparse_algebra_t>
{
    template <typename Structure>
    std::array<sparse_matrix, 3>
    operator()(const Structure& s,
               const std::vector<index_t>& dofmap,
               const size_t n_f, const size_t n_b)
    {
        return triplets_impl(s, dofmap, n_f, n_b);
        // return matrices_impl(s, dofmap, n_f, n_b);
    }
    
    template <typename Structure>
    std::array<sparse_matrix, 3>
    triplets_impl(const Structure& s,
                  const std::vector<index_t>& dofmap,
                  const size_t n_f, const size_t n_b)
    {
        constexpr size_t dim = Kind::ndof;
        
        // Compute number of non-zero entries for each matrix
        auto nnzs = count_nnz_entries<dim>(s, dofmap, n_f, n_b);

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
            auto K_e = element_matrix_assembler<Kind>()(e, s);
        
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
    
    template <typename Structure>
    std::array<sparse_matrix, 3>
    matrices_impl(const Structure& s,
                  const std::vector<index_t>& dofmap,
                  const size_t n_f, const size_t n_b)
    {    
        constexpr size_t dim = Kind::ndof; 
            
        // Compute number of non-zero entries for each matrix
        auto nnzs = count_nnz_entries<dim>(s, dofmap, n_f, n_b);
        
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
            auto K_e = element_matrix_assembler<Kind>()(e, s);
        
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


template <size_t Dim, typename S>
std::array<std::vector<size_t>, 3>
count_nnz_entries(const S& s, const std::vector<index_t>& dofmap,
                  const size_t n_f, const size_t n_b) 
{    
    // Aux vars
    constexpr size_t dim = Dim;//kind_of<S>::type::ndof;
    
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

//------------------------------------- Post-processing -----------------------------------------//
template <typename Structure, typename Result>
real compute_compliance(
    const Structure & structure,
    const std::vector<Result> & results)
{
    real compliance = 0.;
    
    for (const auto v : boost::make_iterator_range(vertices(structure)))
        compliance += structure[v].load.transpose() * results[v].displacement;
    
    return compliance;
}


template <typename Structure>
real compute_mass(const Structure & s, const real rho)
{
    real mass = 0.0;

    for(auto e : boost::make_iterator_range(edges(s)))
    {
        // Compute element length
        auto src = source(e, s);
        auto trg = target(e, s);
        real l = (s[src].coords - s[trg].coords).norm();

        // Compute element mass & add it
        mass += rho * l * s[e].A;
    }
    return mass;
}


} //end namespace eva








// template <typename StructureKind> 
// struct internal_forces_getter 
// {
//     template <typename S>
//     void 
//     operator()(const typename S::edge_descriptor& e, const S& s,
//                const std::vector<real>& u, const std::vector<real>& f) 
//     { 
//         static_assert(!std::is_same<S,S>::value, 
//                       "Method not implemented");
//     }
// };

// template <typename S> 
// auto
// get_internal_forces(const typename S::edge_descriptor& e, const S& s, 
//                     const std::vector<real>& u, const std::vector<real>& f) 
// {
//     return internal_forces_getter<typename kind_of<S>::type>()(e, s, u, f);
// }
