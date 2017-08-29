# ifndef __THERMO_FRAME_PROBLEM_UTILS__
# define __THERMO_FRAME_PROBLEM_UTILS__

# include "thermo_frame_problem.hpp"
# include <pagmo/src/archipelago.h>
# include <boost/program_options.hpp>
# include <cassert>
# include <random>

namespace utils { namespace thermo_frame_problem {

//---------------------------------  Aliases ------------------------------------//
using eva::real;
namespace po = boost::program_options;
namespace pb = pagmo::problem;

using problem_t =  pb::thermo_frame;

using struct_kind_t = problem_t::struct_kind_t;
using thermo_kind_t = problem_t::thermo_kind_t;
using structure_t   = problem_t::structure_t;
using joint_t       = problem_t::joint_t;
using element_t     = problem_t::element_t;

using structure_t = pb::thermo_frame::structure_t;
using joint_t     = pb::thermo_frame::joint_t;
using element_t   = pb::thermo_frame::element_t;

using joint_grid_t   = std::vector< std::vector<joint_t> >;
using boundary_map_t = pb::thermo_frame::boundary_map_t;
using topology_t     = pb::thermo_frame::topology_t;
using bounds_t       = pb::thermo_frame::bounds;


// ---------------------------- Declarations -----------------------------------//
/// Builds the map containing containing the problem configuration
/// by reaading the specified file, throwing the proper exceptions 
/// when  a field is missing or has an invalid value.
/// @param[in] Filename name of the file containing the configuration
/// @returns A map containing the options red
po::variables_map
read_config_file(const std::string & filename);


/// Creates a grid of node properties setting coordinates.
/// @param[in] length The grid length
/// @param[in] height The grid height
/// @param[in] m Number of grid rows
/// @param[in] n Number of grid columns
/// @param[in] jitter Whether to apply or not a random jitter to the
///                   node coordinates
/// @returns A 2D grid of node properties.
joint_grid_t make_grid(
    real length, real height, size_t m,
    size_t n, bool jitter = false);

/// Adds thermal and mechanical boundary conditions, filling the 
/// associated fields of the proper nodes
void
add_bcs(joint_grid_t & grid, const po::variables_map & config);


/// Adds thermal and mechanical loads, filling the 
/// associated fields of the proper nodes
void
add_loads(joint_grid_t & grid, const po::variables_map & config);


/// Computes the number of DOFs associated with the structure
/// given its topology and boundary map
/// @param[in] topology The structure topology
/// @param[in] boundary_map Map indicating boundary nodes
/// @returns The total number of DOFs
size_t
compute_dof_nr(
    const topology_t & topology,
    const boundary_map_t & boundary_map);
// size_t
// compute_dof_nr(const size_t m , const size_t n,
//                const std::string & rule);


/// Creates a map which indicates if a not lies on
/// the grid boundary or not. Possible values
/// are:
/// -] 0 -> internal node
/// -] 1 -> boundary row node
/// -] 2 -> boundary column node
/// -] 3 -> corner node
/// @param[in] grid The grid
boundary_map_t
make_boundary_map(const joint_grid_t & grid);


/// Creates a topology connecting the provided grid nodes
/// using the given rule.
/// @param[in] grid The grid to connect
/// @tparam[in] rule A functor that overloads operator
///                  (int i1, int j1, inti2, int j2) and 
///                  and evaluates to true node (i1, j1)
///                  has to be connected to (12, j2)
/// @returns A sparse matrix representing the generated
///          topology
template <typename Functor>
topology_t make_topology(
    const joint_grid_t & grid,
    const Functor & rule,
    long reserve_hint = 4);


/// Generates a thermo frame given a grid containing
/// the nodes positions and properties and its topology.
/// @param[in] grid A grid of node properties
/// @param[in] topology The structure topology
structure_t
make_frame(
    const joint_grid_t & grid,
    const topology_t & topology,
    const po::variables_map & config);


/// Builds the problem object.
/// @param[in] structure The Base frame containing the structure nodes
/// @param[in] boundary_map The map which indicates boundary nodes
/// @param[in] topology The network representing the allowed set of
///                     elements
/// @param[in] config The map representing the config file content
/// @returns An instance of the problem
problem_t
make_problem(
    const structure_t & structure,
    const boundary_map_t & boundary_map,
    const topology_t & topology,
    const po::variables_map & config);



struct neumann_rule
{
    template <typename Joints>
    bool operator()(Joints & joints, long i1, long j1, long i2, long j2) const
    {
        if (std::abs(i1 - i2) + std::abs(j1 - j2) <= 1) return true;

        return false;
    }
};


struct moor_rule
{
    template <typename Joints>
    bool operator()(Joints & joints, long i1, long j1, long i2, long j2) const
    {
        if (std::abs(i1 - i2) <= 1 && std::abs(j1 - j2) <= 1) return true;

        return false;
    }
};


int
get_champion_island_idx(const pagmo::archipelago & archi);

pagmo::population::champion_type
get_champion(const pagmo::archipelago & archi);

void
display_displaced(const structure_t & frame);

void
export_to_vtu(const structure_t & frame,
              const std::string & filename);

struct problem_components
{
    structure_t    base_frame;
    boundary_map_t boundary_map;
    size_t         dof_nr;
    bounds_t       bounds;
    /// Serialization friend function
    friend class boost::serialization::access;
    template <class Archive>
        void serialize(Archive & ar, const unsigned int)
    {
        ar & base_frame & boundary_map & dof_nr & bounds;
    }
    
};


// void
// display_and_save(
//     const pagmo::population::champion_type & champion,
//     const problem_t & problem, const std::string & filename);


// ----------------------------- Definitions --------------------------------//


template <typename Functor>
topology_t make_topology(
    const joint_grid_t & joints,
    const Functor & rule,
    long reserve_hint)
{
    // Get number of joints
    auto m = joints.size();
    auto n = joints.front().size();
    auto nr_joints = m * n;
    
    // Prealloc connection matrix
    auto ret = topology_t(nr_joints, nr_joints);
    ret.reserve(Eigen::VectorXi::Constant(nr_joints, 4));
    
    
    // Insert edges    
    for (auto it1 = 0u; it1 < nr_joints; ++it1)
        for (auto it2 = it1 + 1; it2 < nr_joints; ++it2)
            // If on upper diagonal & condition is true
            if (rule(joints, it1 / n, it1 % n, it2 / n, it2 % n))
                ret.insert(it1, it2) = true;

    // Compress matrix & return
    ret.makeCompressed();
    return ret;
}
    


}}// end namespace
#endif
