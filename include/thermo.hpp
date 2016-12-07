/**
 * @file thermo.hpp
 * @brief Contains methods and classes to solve thermal problems on structures.
 */
// eva
# include "structure.hpp"

namespace eva { namespace thermo {


//----------------------------------------- Properties -------------------------------------------//
/// Thermo structure tag type
template <typename StructureKind>
struct kind;

/// Thermo structure joint (node) properties
template <typename StructureKind>
struct joint;

/// Thermo structure element (edge) properties
template <typename StructureKind>
struct element;

//-------------------------------------- Thermo tag types -----------------------------------------//
/// Thermo structure type
template <typename StructureKind>
using structure = generic_structure<
    joint  <StructureKind>,
    element<StructureKind>,
    kind   <StructureKind>
    >;

template <typename StructureKind>        
struct kind : public StructureKind
{
    using structure_kind = StructureKind;
};


template <typename StructureKind>        
struct joint : public sa::joint<StructureKind>
{
    real th_bc;   ///< Thermal Dricihlet BCs (fixed temperature)    
    real th_flux; ///< Thermal flux
};

template <typename StructureKind>        
struct element : public sa::joint<StructureKind>
{
    real k; // Thermal conductivity
};




namespace eva {
//------------------------------------- Results Assembling ---------------------------------------//

template <typename StructureKind>
struct result< thermo::kind<StructureKind> >
    : public results<StructureKind>
{
    real T; ///< Temperature

    result(const real temp,
           const fixed_vector<ndof>& u,
           const fixed_vector<ndof>& f,
           const sa::truss_joint<N>& p)
        : results<StructureKind>(u, f, p)
        , T(temp)
    { }
};

} // end namespace eva
        
template <typename Params, typename StructureKind, typename S>
auto
// std::vector< result<StructureKind> >
solve(const S& s, Params p) 
{
    // Solve structural problem
    auto r = solver<StructureKind::structure_kind, Params>()(s);
    const auto& u = std::get<0>(r);
    const auto& f = std::get<1>(r);

    // Solve thermal problem
    auto t = thermal_solver<Params>()(s);

    // Gather results results
    constexpr static int n_loc_dof = StructureKind::ndof;
    const auto n_nodes = num_vertices(s);
    
    // Pre-allocate results
    using result_t = result< thermo_result<StructureKind> >;
    auto results = std::vector<result_t> ();
    results.reserve(n_nodes);

    // For each joint (node)
    for (size_t node_it = 0u; node_it < n_nodes; ++node_it)
    {
        // Get current node starting index in u and f
        size_t node_start_idx = node_it * n_loc_dof;
        
        // Assemble the corresponding results
        results.emplace_back(result_t(u.segment<n_loc_dof>(node_start_idx),
                                      f.segment<n_loc_dof>(node_start_idx),
                                      t(node_start_idx),
                                      s[node_it]));
    }    
    return results;
}



template <typename Params>
struct solver<thermo_structure, Params> 
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
        auto matrices = assemble_system_matrices<algebra_t>(s, dofmap, n_f, n_b);
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
        
        return assemble_results<StructureKind>(uu, ff, s);
    }
};



template <typename S>
std::tuple<std::vector<index_t>, size_t, size_t>
build_global_dofmap(const S& s) 
{    
    // Initialize the DOF map
    size_t n_verts = num_vertices(s);
    auto  dofmap   = std::vector<index_t> (dim*n_verts);
    
    size_t n_bcs  = 0u;   
    size_t bc_pos = dofmap.size()-1u;
    size_t ff_pos = 0u;
    
    for (size_t v = 0u; v < n_verts; ++v)
    {
        // If T = NaN => Free DOF else BC DOF
        if (std::isnan(s[v].T_bc)) {
            // Put DOF at the beginning
            dofmap[dim*v+i] = ff_pos++;
        }
        else
        {
            // Put DOF at the end and increment counter
            dofmap[v+i] = bc_pos--;
            ++n_bcs;
        }
        
    }
    // Return (dofmap, #(free DOF), #(BC DOF))
    return std::make_tuple(dofmap, dofmap.size()-n_bcs, n_bcs);
}



template <typename S>
auto assemble_thermo_matrices(const S& s,
                              const size_t n_f,
                              const size_t n_b)
{
}



template <typename StructureKind>
struct known_terms_assembler< structure<StructureKind> >
{
    template <typename S>
    std::array<dense_vector, 2>
    operator()(const S& s, const size_t n_f, const size_t n_b) 
    {
        // Init return var
# pragma clang diagnostic push
# pragma clang diagnostic ignored "-Wmissing-braces"
        auto ret = std::array<dense_vector, 2> {dense_vector(n_f), dense_vector(n_b)};
# pragma clang diagnostic pop
        
        auto& f_f = std::get<0>(ret);   // Force vector (Free DOF only)
        auto& u_b = std::get<1>(ret);   // Displacement vector (BC DOF only)
        
        size_t pos_f = 0u;     // f_f incremental iterator
        size_t pos_u = n_b;    // u_b decremental iterator
        
        for (auto v : boost::make_iterator_range(vertices(s)))
        {
            // Get bcs & loads
            const auto th_bc   = s[v].th_bc;
            const auto th_flux = s[v].th_flux;

            // If free DOF => write to f & post-increase position
            if (std::isnan(th_bc))
                f_f(pos_f++) = load(i);
            // If BC DOF => write to u & pre-decrease position
            else 
                u_b(--pos_u) = bcs(i);
        }
        
        return ret;
    }
};


template <>
struct stiffness_matrix_assembler<dense_algebra_t>
{
    template <typename S>
    std::array<dense_matrix, 3>
    operator()(const S& s,
               const std::vector<index_t>& dofmap,
               const size_t n_f, const size_t n_b)
    {        
        // Allocate & init submatrices
        dense_matrix K_ff = dense_matrix::Zero(n_f, n_f);
        dense_matrix K_fb = dense_matrix::Zero(n_f, n_b);
        dense_matrix K_bb = dense_matrix::Zero(n_b, n_b);

        // Loop over all vertices
        for (auto&& v : make_iterator_range(vertices(s)))
        {
            // Get current vertex mapped value
            auto avv = dofmap[av];

            // Init diagonal value
            auto dval = real(0);
            
            // For each of the adjacent nodes
            for (auto&& av : make_iterator_range(adjacent_vertices(v, s)))
            {
                // Compute element length
                auto l = (s[av].coords - s[v]).norm();

                // Compute element thermal resistance
                const auto& el = edge(v, av, s);
                auto R = l / (el.A * el.k);

                // Get adj vertex mapped value
                auto avv = dofmap[av];

                // Assign value
                switch (vv < n_f)
                {
                    case 0 : K_bb(vv-n_f, avv-n_f) = 1/R; break;
                        
                    case 1 : K_fb(vv, avv-n_f)     = 1/R; break;
                        
                    case 2 : /* Does nothing (no K_bf) */ break;
                        
                    case 3 : K_ff(vv, avv)         = 1/R; break;
                }                

                // Update diagonal value
                dval += 1/R;            
            }

            // Assign diagonal value
            switch (vv < n_f)
            {
                case 0 : K_bb(vv-n_f, vv-n_f) = -dval; break;
                        
                case 1 : K_fb(vv, vv-n_f)     = -dval; break;
                        
                case 2 : /* Does nothing (no K_bf) */  break;
                    
                case 3 : K_ff(vv, vv)         = -dval; break;
            } 

        }
# pragma clang diagnostic push
# pragma clang diagnostic ignored "-Wmissing-braces"
        return { std::move(K_ff), std::move(K_fb), std::move(K_bb) };
# pragma clang diagnostic pop
    }
};



}} // end namespaces

