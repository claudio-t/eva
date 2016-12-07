/**
 * @file structural_problem.cpp
 * @brief Contains structural_problem class member definitions.
 */
# include "structural_problem.hpp"

# include "continuous_structural_problem.hpp"
# include <boost/math/constants/constants.hpp>


namespace pagmo { namespace problem {

structural_problem::continuous_structural_problem(
    unsigned int dim,
    const std::vector<joint_t>& joints,
    )
    : base(int(dim), // Global dimension
           int(0),   // Integer dimension
           int(1),   // Fitness dimension
           int(3),   // Global constraints dimension
           int(1),   // Inequality constraints dimension
           0.)       // Constraints tolerance
    , joints_(joints)
{
    /// Store joints with BCs/Loads
    store_mandatory_joints();

    // Set variables bounds
    set_bounds();
}


// void
// structural_problem::set_bounds()
// {
//     set_lb(lower_bounds); // Lower bound (forall vars)
//     set_ub(upper_bounds); // Upper bound (forall vars)
// }


void
structural_problem::store_mandatory_joints()
{
    // No BCs applied means no NaN component
    auto no_bcs  = [](const decltype(joint_t::bcs)& vbcs){
        for (auto i = 0u; i < vbcs.size(); ++i)
            if (!std::isnan(vbcs[i])) return false;
        
        return true;
    };

    // No load means zero vector
    auto no_load = decltype(joint_t::load)::Zero();

    // Check each joint for both 
    auto n_verts = joints.size();
    for (auto i = 0u; i < n_verts; ++i)
        if (!no_bcs(joints[i].bcs) || joints[i].load != no_load)
            mandatory_joints_.emplace_back(i);
}


std::string
structural_problem::get_name() const
{
    return "Structural Optimization Problem";
}

std::string
structural_problem::human_readable_extra() const
{
    std::ostringstream oss;
    // oss << "\n\tData member value: " << m_member;
    return oss.str();
}

void
structural_problem::compute_constraints_impl(constraint_vector& constraints,
                                             const decision_vector& genes) const
{
    // The graph should be a single connected component
    const auto structure = encode_genes(genes);
    constraints[0] = mandatory_vertices_.size() *
        double(!eva::utils::is_connected(structure));
    
    // Each mandatory vertex has at least one edge
    auto vertices_counter = true;
    for (const auto v_idx : mandatory_vertices_)
        if (out_degree(v_idx, structure) > 0)
            ++vertices_counter;
    
    constraints[1] = double(mandatory_vertices_.size() - vertices_counter);

    // Mass constraint
    double volume = 0.;
    int gene_it   = 0;
    constexpr double pi = boost::math::constants::pi<double>();
    
    for (auto node_it = 0u; node_it < gene_mask_.size(); ++node_it)
        for (auto link_it = 0u; link_it < gene_mask_[node_it].size(); ++link_it)
        {   
            auto p1 = structure[node_it].coords;
            auto p2 = structure[gene_mask_[node_it][link_it]].coords;

            // Elemen volume = A * length
            const auto r = genes[gene_it++];
            const auto A = pi * r*r;
            const auto L = (p2 - p1).norm();
            
            volume += A * L;
        }
    
    // NB: for inequalty constraints negative means OK, positive KO
    constraints[2] = compute_mass(structure) - max_mass_;
}

structural_problem::structure_t
structural_problem::encode_genes(const decision_vector& genes) const
{
    // Init empty structure and default element
    auto structure = structure_t();

    // Add joints
    for (const auto& joint : joints_) add_vertex(joint, structure);
 
    auto gene_it = 0u;
    for (auto node_it = 0u; node_it < gene_mask_.size(); ++node_it)
        for (auto link_it = 0u; link_it < gene_mask_[node_it].size(); ++link_it)
        {
            auto src = node_it;
            auto trg = gene_mask_[node_it][link_it];

            constexpr double pi = boost::math::constants::pi<double>();
            const auto r = genes[gene_it++];
            const auto E = 2.e11; // Steel ASTM-A36
            const auto A = pi * r*r;
            const auto I = pi/2 * r*r*r*r;
            
            if (r > 0.) add_edge(src, trg, {/*E=*/E, /*A=*/A, /*I=*/I}, structure);
        }
    
    return structure;
}

void
structural_problem::objfun_impl(fitness_vector& fitness,
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
    // std::cout << "Compliance = "<< compliance*1.e10 << std::endl;
}

}} // end namespaces 
