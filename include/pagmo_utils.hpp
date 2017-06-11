# ifndef __EVA_PAGMO_UTILS__
# define __EVA_PAGMO_UTILS__

#include "core.hpp"

namespace utils {

//---------------------------------  Aliases ------------------------------------//
using eva::real;

template <typename Joint>
using joint_grid = std::vector< std::vector<Joint> >;
using topology_t = Eigen::SparseMatrix<bool>;


// ---------------------------- Declarations -----------------------------------//
template <typename Joint>
joint_grid<Joint> make_grid(real length, real height, long m, long n, bool jitter = false);

template <typename Joint, typename Functor>
topology_t make_topology(
    const joint_grid<Joint>& grid,
    const Functor& rule,
    long reserve_hint = 4);
    


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


// ----------------------------- Definitions --------------------------------//


template <typename Joint>
joint_grid<Joint> make_grid(real length, real height, long m, long n, bool jitter)
{
    if (m < 2) throw std::invalid_argument("Parameter m must be >= 2");
    if (n < 2) throw std::invalid_argument("Parameter n must be >= 2");
    
    // Prealloc joint rows
    auto joints = joint_grid<Joint>(m);
    joints.reserve(m);

    // Compute grid x and y spacings
    auto delta_x = length / (n - 1);
    auto delta_y = height / (m - 1);

    // // Add some randomess
    std::mt19937 rng;
    rng.seed(std::random_device()());
    auto dx = std::uniform_real_distribution<float>(-delta_x / 2, delta_x / 2);
    auto dy = std::uniform_real_distribution<float>(-delta_y / 2, delta_y / 2);

    
    for (auto i = 0u; i < m; ++i)
    {        
        // Reserve space
        auto n_nodes = n;//(interleaved_row ? n-1 : n);
        joints[i].resize(n_nodes);
        
        for (auto j = 0u; j < n_nodes; ++j)
        {
            // Init coordinates
            auto x = j*delta_x;
            auto y = i*delta_y;

            if (jitter)
            {                
            auto is_boundary_row = i == 0 || i == m-1;
            auto is_boundary_col = j == 0 || j == n-1;

            // Add random jitter
            if (!is_boundary_row) y += dy(rng);
            if (!is_boundary_col) x += dx(rng);
            }
            
            
            // Set joint coords
            joints[i][j].coords << x, y;
        }
    }

    return joints;
}


template <typename Joint, typename Functor>
topology_t make_topology(
    const joint_grid<Joint> & joints,
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



}



#endif
