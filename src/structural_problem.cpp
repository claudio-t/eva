/**
 * @file structural_problem.cpp
 * @brief Contains structural_problem class member definitions.
 */
# include "structural_problem.hpp"


namespace pagmo { namespace problem {

structural_problem::structural_problem(
    unsigned int dim,
    const std::vector<joint_t>& joints,
    const gene_mask_t& gene_mask,
    const int n_max
    )
    : base(int(dim), // Global dimension
           int(dim), // Integer dimension
           int(1),   // Fitness dimension
           int(3),   // Global constraints dimension
           int(1),   // Inequality constraints dimension
           0.)  // Constraints tolerance
    , joints_   (joints)
    , gene_mask_(gene_mask)
    , n_max_    (n_max)
{
    // Check that dim and gene_mask dimensions are coherent
    auto ldim = 0u;
    for (const auto& el: gene_mask) ldim += el.size();
        
    if (ldim != dim)
        throw std::runtime_error("Gene mask size does not match input dimension!");

    
    // Save indices of vertices with either an applied load or a BC
    auto no_bcs  = [](const decltype(joint_t::bcs)& vbcs){
        for (auto i = 0u; i < vbcs.size(); ++i)
            if (!std::isnan(vbcs[i])) return false;
        
        return true;
    };
    auto no_load = decltype(joint_t::load)::Zero();


    // for (const auto& joint : joints_)
    auto n_verts = joints.size();
    for (auto i = 0u; i < n_verts; ++i)
        if (!no_bcs(joints[i].bcs) || joints[i].load != no_load)
            mandatory_vertices_.emplace_back(i);
    
    // Problem bounds
    set_lb(0.); // Lower bound (forall vars)
    set_ub(1.); // Upper bound (forall vars)
}

base_ptr structural_problem::clone() const
{
    return base_ptr(new structural_problem(*this));
}

std::string structural_problem::get_name() const
{
    return "Structural Problem";
}

std::string structural_problem::human_readable_extra() const
{
    std::ostringstream oss;
    // oss << "\n\tData member value: " << m_member;
    return oss.str();
}

void structural_problem::compute_constraints_impl(constraint_vector& constraints,
                                                  const decision_vector& genes) const
{
    // Check for structure feasibility:
    // the graph representing the structure must have a single
    // connected component containing all of the nodes with either an
    // applied load or a boundary condition.

    // The graph should be a single connected component
    const auto structure = encode_genes(genes);
    constraints[0] = mandatory_vertices_.size() * double(!eva::utils::is_connected(structure));
    
    // Each mandatory vertex has at least one edge
    auto vertices_counter = true;
    for (const auto v_idx : mandatory_vertices_)
        if (out_degree(v_idx, structure) > 0)
            ++vertices_counter;
    
    constraints[1] = double(mandatory_vertices_.size() - vertices_counter);
    
    // Maximum number of links must be lesser-equal than n_max
    int sum = 0;
    for (const auto gene : genes) sum += int(gene);

    // For inequalty constraints if negative OK, if positive KO
    constraints[2] = double(sum - n_max_);
}

structural_problem::structure_t
structural_problem::encode_genes(const decision_vector& genes) const
{
    // Init empty structure and default element
    auto structure = structure_t();
    auto element   = element_t {/*E=*/1., /*A=*/1., /*I=*/1.};
    
    // auto mandatory_vertices = std::vector<eva::index_t>();
    for (const auto& joint : joints_)
        add_vertex(joint, structure);
 
    auto gene_it = 0u;
    for (auto node_it = 0u; node_it < gene_mask_.size(); ++node_it)
        for (auto link_it = 0u; link_it < gene_mask_[node_it].size(); ++link_it)
        {
            auto src = node_it;
            auto trg = gene_mask_[node_it][link_it];
            if (genes[gene_it++]) add_edge(src, trg, element, structure);
        }
    
    return structure;
}

void structural_problem::objfun_impl(fitness_vector& fitness,
                                     const decision_vector& genes) const
{
    // Build structure according to the current genes
    const auto structure = encode_genes(genes);
    
    // Setup solver and solve
    auto results = solve(structure, eva::dense_solver_params<>());
    
    // Compute compliance
    auto n_nodes    = num_vertices(structure);
    auto compliance = 0.;
    for (size_t vidx = 0u; vidx < n_nodes; ++vidx)
        compliance += structure[vidx].load.transpose() * results[vidx].displacement;

    // Set fitness equal to compliance
    fitness[0] = compliance;
    std::cout << "Compliance = "<< compliance << std::endl;
}

}} // end namespaces 