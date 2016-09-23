# ifndef __PAGMO_CONTINUOUS_STRUCTURAL_PROBLEM__
# define __PAGMO_CONTINUOUS_STRUCTURAL_PROBLEM__
/**
 * @file continuous structural_problem.hpp
 * @brief Contains the structural_problem class declaration used in
 * the topological optimization procedure.
 */

// eva
# include "truss.hpp"
# include "frame.hpp"
# include "utils.hpp"
# include "vtk.hpp"
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
class __PAGMO_VISIBLE continuous_structural_problem : public base
{
public:
    /* Member types */
    using structure_t = eva::frame2d;
    using joint_t     = typename eva::joint_of  <structure_t>::type;
    using element_t   = typename eva::element_of<structure_t>::type;
    using gene_mask_t = std::vector< std::vector<size_t> >;
    
    /* Constructor(s) */
    continuous_structural_problem(
        unsigned int dim,
        const std::vector<joint_t>& joints,
        const gene_mask_t& gene_mask,
        const double m_max
        );
    
    /* Pagmo stuff */    
    base_ptr clone() const;
    std::string get_name() const;
    
    /* Utils */
    structure_t encode_genes(const decision_vector& genes) const;

protected:
    /* Pagmo stuff */
    void objfun_impl(fitness_vector&, const decision_vector&) const;
    void compute_constraints_impl(fitness_vector&, const decision_vector&) const;
    std::string human_readable_extra() const;
        
private:
    /* Serialization related stuff (still Pagmo stuff) */
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int)
    {
        ar & boost::serialization::base_object<base>(*this);
        ar & joints_;
        ar & gene_mask_;
        ar & m_max_;
        ar & mandatory_vertices_;
    }
    
    /* Member Data */
    /// List of available structure joints
    const std::vector<joint_t>& joints_;             

    /// Genes mask
    const gene_mask_t& gene_mask_;

    /// Maximum mass
    const int m_max_;         

    /// Indices of mandatory vertices (the ones with either an applied
    /// load or a specified BC)
    std::vector<eva::index_t> mandatory_vertices_;
};


}} // namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::continuous_structural_problem)
# endif // __PAGMO_CONTINUOUS_STRUCTURAL_PROBLEM__
