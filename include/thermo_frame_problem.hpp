# ifndef __THERMO_FRAME_PROBLEM__
# define __THERMO_FRAME_PROBLEM__


// eva
# include "truss.hpp"
# include "frame.hpp"
# include "thermo.hpp"
# include "vtk.hpp"
// boost
# include <boost/container/flat_map.hpp>
# include "boost_flat_container_serialization.hpp"
// pagmo
# include <pagmo/src/config.h>
# include <pagmo/src/serialization.h>
# include <pagmo/src/types.h>
# include <pagmo/src/problem/base.h>


namespace pagmo { namespace problem {

class __PAGMO_VISIBLE thermo_frame : public base
{
public:

    using struct_kind_t = eva::frame_kind<2>;
    using thermo_kind_t = eva::thermo_kind<struct_kind_t>;
    // using thermo_structure_t = eva::thermo_structure<struct_kind_t>;

    using structure_t = eva::thermo_structure<struct_kind_t>;
    using joint_t     = typename eva::joint_of  <structure_t>::type;
    using element_t   = typename eva::element_of<structure_t>::type;

    using topology_t     = Eigen::SparseMatrix<bool>;
    using boundary_map_t = boost::container::flat_map<size_t, int>;
    

    struct bounds;
    // struct physical_properties;

    /// Defualt constructor
    thermo_frame() : base(1) {}
    
    /// Constructor that initializes class members
    thermo_frame(
        const structure_t & base_frame,
        // const topology_t  & frame_topo,
        const boundary_map_t & boundary_map,
        const int problem_dim,
        const bounds bnds
        // const physical_properties phys_props
        );
    
    /// Encode genes into the proper structure
    structure_t encode_genes(const decision_vector & genes) const;

    
    /// Return flags indicating wheter the node lies or not on
    /// a boundary row or column
    std::pair<bool, bool> is_boundary(const size_t vertex_idx) const;

    /// Clone problem
    base_ptr clone() const;

    /// Get problem name
    std::string get_name() const;

    /// Extra human readable info for the problem.
    std::string human_readable_extra() const;
    
    /// Simple POD containing the problem bounds
    struct bounds
    {
        eva::real x_jitter;
        eva::real y_jitter;

        eva::real min_radius;
        eva::real max_radius;

        // eva::real min_volume;
        // eva::real max_volume;
        // eva::real tol;

        template <typename Archive>
        void serialize(Archive & ar, const unsigned int /* file_version */)
        {
            ar & x_jitter & y_jitter & min_radius & max_radius;// & min_volume & max_volume & tol;
        }
    };

    /// Simple POD containing the problem physical properties
    // struct physical_properties
    // {
    //     eva::real E;
    //     eva::real k;

    //     template <typename Archive>
    //     void serialize(Archive & ar, const unsigned int /* file_version */)
    //     {
    //         ar & E & k;
    //     }        
    // };
    
    
protected:
    /// Evaluates the fitness function f:
    /// x -> f(x) = (compliance, )
    void objfun_impl(
        fitness_vector & f,
        const decision_vector & x) const;

    /// Computes constraints:
    /// x -> h(x) = (volume - max_volume, ) = 0
    // void compute_constraints_impl(
    //     constraint_vector & constraints,
    //     const decision_vector & x) const;
    
private:
    /// Serialization friend function
    friend class boost::serialization::access;
    template <class Archive>
        void serialize(Archive & ar, const unsigned int);

    
    /// Base frame (i.e. joints with loads applied
    /// having the external boundary completely connected)
    structure_t base_frame_;

    /// Set of allowed connections between frame joints
    // topology_t  frame_topo_;

    /// Boundary map
    boundary_map_t boundary_map_;

    /// Problem bounds
    bounds bounds_; ///< Frame maximum allowed volume

    /// Problem physical properties
    // physical_properties phys_props_;
};


template <class Archive>
void thermo_frame::serialize(Archive & ar, const unsigned int)
{
    ar & boost::serialization::base_object<base>(*this);
    ar & base_frame_;
    // ar & frame_topo_;
    ar & boundary_map_;
    ar & bounds_;
    // ar & phys_props_;
}




}}// end namespace pagmo & problem
BOOST_CLASS_EXPORT_KEY(pagmo::problem::thermo_frame)
#endif
