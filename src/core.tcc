/// @file  core.tcc
/// @brief Contains definitions of template functions and class methods
///        declared in file core.hpp .
// #include <iostream>

namespace eva {


//######################################## Definitions #############################################      
//-------------------------------------- Problem Solving -----------------------------------------//
template <
    typename Structure,
    typename Kind   = typename kind_of<Structure>::type,
    typename Params = dense_solver_params<typename Kind::default_dense_solver_t>
    >
auto solve(const Structure& s, const Kind, const Params) 
{
    auto r = solver<Kind, Params>()(s);
    const auto& uu = std::get<0>(r);
    const auto& ff = std::get<1>(r);
    
    return assemble_results<Kind>(uu, ff, s);
    // return solver<kind_t, Params>()(s);
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

        // Define kind type
        using kind_t = StructureKind;
        
        // Build (global) DOF map with BCs related nodes located at the back
        auto dofmap = std::vector<index_t>();
        size_t n_f, n_b;
        std::tie(dofmap, n_f, n_b) = build_global_dofmap<kind_t>(s);
        
        // Assemble system stiffness matrices (in global coordinates)
        // NOTE: use .sparseView()?
        auto matrices = assemble_system_submatrices<algebra_t, S, kind_t>(s, dofmap, n_f, n_b);
        const auto & K_ff = std::get<0>(matrices);
        const auto & K_fb = std::get<1>(matrices);
        const auto & K_bb = std::get<2>(matrices);
        
        
        // Assemble known terms (in global coordinates)
        auto known_terms = known_terms_assembler<kind_t>()(s, n_f, n_b);
        const dense_vector & f_f = std::get<0>(known_terms); // displacements on BC nodes
        const dense_vector & u_b = std::get<1>(known_terms); // applied loads

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
template <typename A, typename S, typename Kind>
auto
assemble_system_submatrices(const S& s, const std::vector<index_t>& dofmap,
                          const size_t n_f, const size_t n_b) 
{ 
    using kind_t = typename std::conditional<
        std::is_same<Kind, void>::value,
        typename kind_of<S>::type,
        Kind
        >::type;
    
    
    return system_submatrices_assembler<A, kind_t>()(s, dofmap, n_f, n_b);
}

template <typename AlgebraKind, typename StructureKind>
struct system_submatrices_assembler
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


// template <typename S, typename Kind>
// std::array<dense_vector, 2>
// assemble_known_terms(const S& s, const size_t n_f, const size_t n_b) 
// {
//     return known_terms_assembler<typename kind_of<S>::type>()(s, n_f, n_b);
// }

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


 
//---------------------------------------- DOF handling ------------------------------------------//
template <typename Kind, typename S>
std::tuple<std::vector<index_t>, size_t, size_t>
build_global_dofmap(const S& s) 
{
    using kind_t = typename std::conditional<
        std::is_same<Kind, void>::value,
        typename kind_of<S>::type,
        Kind
        >::type;
    
    return global_dofmap_builder<kind_t, kind_t::ndof>()(s);
}

template <typename Kind, size_t Dim>
struct global_dofmap_builder
{
    template <typename S>
    std::tuple<std::vector<index_t>, size_t, size_t>
    operator()(const S& s)
    {
        constexpr size_t dim = Kind::ndof;
    
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
                if (std::isnan(bcs[i]))
                {
                    // Put DOF at the beginning and increment position
                    dofmap[dim*v+i] = ff_pos++;
                }
                else
                {
                    // Put DOF at the end and increment position
                    dofmap[dim*v+i] = bc_pos--;
                    ++n_bcs;
                }
            }
        }
        // Return (dofmap, #(free DOF), #(BC DOF))
        return std::make_tuple(dofmap, dofmap.size()-n_bcs, n_bcs);   
    }
};  

template <typename Kind>
struct global_dofmap_builder<Kind, 1>
{
    template <typename S>
    std::tuple<std::vector<index_t>, size_t, size_t>
    operator()(const S& s)
    {
        // Initialize the DOF map
        size_t n_verts = num_vertices(s);
        auto  dofmap   = std::vector<index_t> (n_verts);
    
        size_t n_bcs  = 0u;
        size_t bc_pos = dofmap.size()-1u;
        size_t ff_pos = 0u;
    
        for (size_t v = 0u; v < n_verts; ++v)
        {
            const auto & bc_d = s[v].T_bc;
        
            // If T = NaN => Free DOF else Dirichlet BC DOF
            if (std::isnan(bc_d))
            {
                // Put DOF at the beginning
                dofmap[v] = ff_pos++;
            }
            else
            {
                // Put DOF at the end and increment counter
                dofmap[v] = bc_pos--;
                ++n_bcs;
            }
        }
        // Return (dofmap, #(free DOF), #(BC DOF))
        return std::make_tuple(dofmap, dofmap.size()-n_bcs, n_bcs);   
    }
};  



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

}// end namespace eva
