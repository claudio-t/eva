# ifndef __PAGMO_STRUCTURAL_PROBLEM__
# define __PAGMO_STRUCTURAL_PROBLEM__
/**
 * @file structural_problem.hpp
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

namespace eva {

/// Base class for a structural optimization problem.
/// Functions that have to be reimplemented are:
/// -] encode_genes
/// -] clone
/// -] objfun_impl
/// -] set_bounds
/// -] compute_constraints_impl
/// -] compute_mass
/// -]
class structural_problem : public pagmo::problem::base
{
public:
    using structure_t = eva::frame2d;
    using joint_t     = typename eva::joint_of  <structure_t>::type;
    using element_t   = typename eva::element_of<structure_t>::type;
    // using gene_mask_t = std::vector< std::vector<size_t> >;
    
    /// Constructor that initializes class members
    structural_problem(
        unsigned int                dim,
        const std::vector<joint_t>& joints,
        const double                max_mass,
        );
    
    /// Encode genes into the proper structure
    virtual structure_t encode_genes(const decision_vector& genes) const = 0;

    /// Clone method
    virtual base_ptr clone() const = 0;

    /// Returns the problem name
    virtual std::string get_name() const;

    /// Returns extra information
    virtual std::string human_readable_extra() const;
    
    
protected:
    /// Objective function implementation; must be reimplemented.
    virtual void objfun_impl(fitness_vector&, const decision_vector&) const = 0;

    
    /// Defines problem variables bounds; must be reimplemented.
    virtual void set_bounds() = 0;

    /// Evaluates constraints: the default implementation requires  
    /// the structure to be a single connected component containing
    /// all mandatory vertices; moreover the structure mass is
    /// required to be lesser than max_mass_ (thus a suitable
    /// implementation of compute_mass function must be provided).    
    virtual void compute_constraints_impl(fitness_vector&, const decision_vector&) const;
    
    /// Computes the structure mass; it is called by
    /// compute_constraints_impl in order to enforce the mass constraint.
    virtual double compute_mass(const structure_t& structure) const = 0;
    
    /// Save indices of joints with either an applied load or a BC
    virtual void store_mandatory_joints();

    
private:
    friend class boost::serialization::access;

    /// Object serialization
    template <class Archive>
    void serialize(Archive &ar, const unsigned int)
    {
        ar & joints_;
        ar & max_mass_;
        ar & mandatory_joints_;
    }
    
    /// List of available structure joints
    const std::vector<joint_t>& joints_;
    
    /// Maximum mass
    const double max_mass_;         

    /// Indices of mandatory joints (the ones with either an applied
    /// load or a specified BC)
    std::vector<eva::index_t> mandatory_joints_;
};

} // namespaces

# endif // __PAGMO_STRUCTURAL_PROBLEM__
