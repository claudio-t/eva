/**
 * @file structural_problem.cpp
 * @brief Contains structural_problem class member definitions.
 */
# include "structural_problem.hpp"


namespace pagmo { namespace problem {

base_ptr structural_problem::clone() const
{
    std::cout << "Cloning...\n";
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

// void structural_problem::compute_constraints_impl(constraint_vector& constraints,
//                                                   const decision_vector& genes) const
// {
//     // Maximum number of links must be lesser-equal than n_max
//     double sum = 0.;
//     for (const auto gene : genes) sum += double(gene);

//     // For inequalty constraints if negative OK, if positive KO
//     constraints[0] = sum - n_max_;
//     std::cout << "Nr of links = " << sum << "\n";
// }

void structural_problem::objfun_impl(fitness_vector& fitness,
                                     const decision_vector& genes) const
{
    // std::cout << "Evaluating...\n";
    
    // Build structure
    auto structure = structure_t();
    auto element   = element_t {1., 1.};

    // Add vertices keeping track of bcs and loaded related vertex indices
    auto mandatory_vertices = std::vector<eva::index_t>();
    for (const auto& joint : joints_)
    {
        auto vidx = add_vertex(joint, structure);
        
        // Store mandatory vertex indices
        if (joint.bcs  != decltype(joint.bcs) ::Zero()
            ||
            joint.load != decltype(joint.load)::Zero())
            mandatory_vertices.emplace_back(vidx);
    }

    auto gene_it = 0u;
    for (auto node_it = 0u; node_it < gene_mask_.size(); ++node_it)
        for (auto link_it = 0u; link_it < gene_mask_[node_it].size(); ++link_it)
        {
            auto src = node_it;
            auto trg = gene_mask_[node_it][link_it];
            if (genes[gene_it++]) add_edge(src, trg, element, structure);
            // std::cout
            //     << " Adding edge: "
            //     << structure[src].coords.transpose()
            //     << " -- "
            //     << structure[trg].coords.transpose()
            //     << std::endl;
        }
    
    // Check for structure feasibility:
    // the graph representing the structure must have a single
    // connected component containing all of the nodes with either an
    // applied load or a boundary condition.

    auto feasible = true;
    // -] check that each mandatory vertices has at least one edge
    // for (auto&& el : std::vector<int>{0, 4, 8})
        // add_edge(el, el+4, element, structure);
        
    for (const auto v_idx : mandatory_vertices)
    {
        if (out_degree(v_idx, structure) == 0)
        {
            feasible = false;
            break;
        }    
    }
    // -] check that the graph representing the structure is connected
    feasible = feasible && eva::utils::is_connected(structure);

    // if (!feasibility_x(genes)) std::cout << "X NOT FEASIBLE!\n";
    double sum = 0.;
    for (const auto gene : genes) sum += double(gene);
    feasible = feasible && sum <= 10.;
    
    if (feasible)
    {
        // Setup solver and solve
        auto results = solve(structure, eva::dense_solver_params<>());
    
        // Compute compliance
        auto n_nodes    = num_vertices(structure);
        auto compliance = 0.;
        for (size_t vidx = 0u; vidx < n_nodes; ++vidx)
            compliance += structure[vidx].load.transpose() * results[vidx].displacement;    

        // Set fitness equal to compliance
        std::cerr << "Compliance = "<< compliance << std::endl;
        fitness[0] = compliance;
    } 
    else
    {
        // Set the predefined unfeasible value
        std::cout << "unfeasible structure!\n";
        fitness[0] = UNFEASIBLE_SCENARIO_FITNESS;
    }
}

}} // end namespaces 
