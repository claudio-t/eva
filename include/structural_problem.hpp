# ifndef __PAGMO_STRUCTURAL_PROBLEM__
# define __PAGMO_STRUCTURAL_PROBLEM__
/**
 * @file structural_problem.hpp
 * @brief Contains the structural_problem class declaration used in
 * the topological optimization procedure.
 */

// eva
# include "truss.hpp"
# include "utils.hpp"
// std
# include <string>
// pagmo
# include <pagmo/src/config.h>
# include <pagmo/src/serialization.h>
# include <pagmo/src/types.h>
# include <pagmo/src/problem/base.h>

// Value returned in case of an unfeasible structure scenario
# ifndef UNFEASIBLE_SCENARIO_FITNESS
# define UNFEASIBLE_SCENARIO_FITNESS std::numeric_limits<double>::infinity()
# endif

namespace pagmo{ namespace problem {

/// A random problem
/**
 * @author Dario Izzo (dario.izzo@gmail.com)
 */
class __PAGMO_VISIBLE structural_problem : public base
{
public:
    /* Member types */
    using structure_t = eva::truss2d;
    using joint_t     = typename eva::joint_of  <structure_t>::type;
    using element_t   = typename eva::element_of<structure_t>::type;
    using gene_mask_t = std::vector< std::vector<size_t> >;
    
    /* Constructor(s) */
    structural_problem(unsigned int dim,
                       const std::vector<joint_t>& joints,
                       const gene_mask_t& gene_mask,
                       const int n_max)
        : base(int(dim), // Global dimension
               int(dim), // Integer dimension
               int(1),   // Fitness dimension
               int(0),   // Global constraints dimension
               int(0),   // Inequality constraints dimension
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

        // Problem bounds
        set_lb(0.); // Lower bound (forall vars)
        set_ub(1.); // Upper bound (forall vars)

        std::cout << "Constructor called" << std::endl;        
    }
    
    /* Pagmo stuff */    
    base_ptr clone() const;
    std::string get_name() const;

protected:
    /* Pagmo stuff */
    void objfun_impl(fitness_vector&, const decision_vector&) const;
    // void compute_constraints_impl(fitness_vector&, const decision_vector&) const;
    std::string human_readable_extra() const;
    
private:
    /* Serialization related stuff (still Pagmo stuff))*/
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int)
    {
        ar & boost::serialization::base_object<base>(*this);
        ar & joints_;
        ar & gene_mask_;
        ar & n_max_;
        // ar & structure_;
    }
    
    /* Member Data */
    const std::vector<joint_t>& joints_;    ///< List of available structure joints
    const gene_mask_t&          gene_mask_; ///< Genes mask
    const int                   n_max_;     ///< Maximum number of links
    // structure_t structure_; ///<  structure
};

// base_ptr structural_problem::clone() const
// {
//     return base_ptr(new structural_problem(*this));
// }


// // void structural_problem::compute_constraints_impl(constraint_vector& constraints,
// //                                                   const decision_vector& decision) const
// // {
// //     // The graph representing the structure must have a single
// //     // connected component containing all of the nodes with either an
// //     // applied load or a boundary condition.

// //     // How to check:
// //     // -] nr of in/out edges of mandatory vertices must be >= 1
// //     // -] graph has to be connected
// //     if (!feasible) constraints[0] = 1.;   
// // }

// void structural_problem::objfun_impl(fitness_vector& fitness,
//                                      const decision_vector& genes) const
// {
//     // MOVE ME SOMWHERE ELSE!! //
    
//     // Build structure
//     auto structure = structure_t();
//     auto element   = element_t{1., 1.};

//     // Add vertices keeping track of bcs and loaded related ones
//     auto mandatory_vertices = std::vector<eva::index_t>();
//     for (const auto& joint : joints_)
//     {
//         auto vidx = add_vertex(joint, structure);
        
//         // Store mandatory vertex index
//         if (joint.bcs  != decltype(joint.bcs) ::Zero()
//             ||
//             joint.load != decltype(joint.load)::Zero())
//             mandatory_vertices.emplace_back(vidx);
//     }
    
//     for (auto node_it = 0u; node_it < gene_mask_.size(); ++node_it)
//         for (auto link_it = 0u; link_it < gene_mask_[node_it].size(); ++link_it)
//         {
//             auto src = node_it;
//             auto trg = gene_mask_[node_it][link_it];
//             add_edge(src, trg, element, structure);
//             // std::cout
//             //     << " Adding edge: "
//             //     << structure[src].coords.transpose()
//             //     << " -- "
//             //     << structure[trg].coords.transpose()
//             //     << std::endl;
//         }
    
//     // Check for structure feasibility:
//     auto feasible = true;
//     // -] check that each mandatory vertices has at least one edge
//     for (const auto v_idx : mandatory_vertices)
//     {
//         if (out_degree(v_idx, structure) == 0)
//         {
//             feasible = false;
//             break;
//         }    
//     }
//     // -] check that the graph representing the structure is connected
//     feasible = feasible && eva::utils::is_connected(structure);
    
    
//     if (feasible)
//     {
//         // Setup solver and solve
//         auto results = solve(structure, eva::sparse_solver_params<>());
    
//         // Compute compliance
//         auto n_nodes    = num_vertices(structure);
//         auto compliance = 0.;
//         for (size_t vidx = 0u; vidx < n_nodes; ++vidx)
//             compliance += structure[vidx].load.transpose() * results[vidx].displacement;    

//         // Set fitness equal to compliance
//         fitness[0] = compliance;
//     }
//     else
//         // Set the predefined unfeasible value
//         fitness[0] = UNFEASIBLE_SCENARIO_FITNESS;
// }


}} // namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::structural_problem)
# endif // __PAGMO_STRUCTURAL_PROBLEM__
