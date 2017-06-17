# include "thermo_frame_problem.hpp"
# include <boost/program_options.hpp>
# include "graphviz.hpp"
# include <limits>


namespace pagmo { namespace problem {

thermo_frame::thermo_frame(
    const structure_t & base_frame,
    // const thermo_frame::topology_t & frame_topo,
    const thermo_frame::boundary_map_t & boundary_map,
    const int problem_dim,
    const thermo_frame::bounds bnds
    // const physical_properties phys_props
    )
    : base(problem_dim, // Global dimension
           int(0),   // Integer dimension
           int(3),   // Fitness dimension
           int(0),   // Global constraints dimension
           int(0),   // Inequality constraints dimension
           0.)       // Constraints tolerance
    , base_frame_(base_frame)
    // , frame_topo_(frame_topo)
    , boundary_map_(boundary_map)
    , bounds_(bnds)
    // , phys_props_(phys_props)
{
    // Set problem bounds
    auto lbs = std::vector<double>(problem_dim);
    auto ubs = std::vector<double>(problem_dim);

    // Iterator that spans each DOF/gene
    auto bound_it = 0u;

    // First bounds are related to node jitters
    for (const auto v : boost::make_iterator_range(vertices(base_frame)))
    {
        // Get boundary info
        bool is_boundary_row, is_boundary_col;
        std::tie(is_boundary_row, is_boundary_col) = is_boundary(v);

        // If node DOES NOT lie on a boundary column
        // then it is free to move on the x-axis
        if (!is_boundary_col)
        {
            lbs[bound_it] = -bounds_.x_jitter;
            ubs[bound_it] =  bounds_.x_jitter;
            ++bound_it;
        }

        // If node DOES NOT lie on a boundary row
        // then it is free to move on the y-axis        
        if (!is_boundary_row)
        {
            lbs[bound_it] = -bounds_.y_jitter;
            ubs[bound_it] =  bounds_.y_jitter;
            ++bound_it;
        }
    }

    // Remaining bounds are associated to element thicknesses
    // Add elements to the structure
    for (auto e : make_iterator_range(edges(base_frame_)))
    {
        auto src = source(e, base_frame_);
        auto trg = target(e, base_frame_);
        
        // If the edge connects two boundary nodes
        // the edge thickness lower bound is the min radius
        if (boundary_map_[src] > 0 && boundary_map_[trg] > 0)
        {
            lbs[bound_it] = bounds_.min_radius;
            ubs[bound_it] = bounds_.max_radius;
        }
        // Otherwise the thickness can vanish
        else
        {
            lbs[bound_it] = bounds_.min_radius;
            ubs[bound_it] = bounds_.max_radius;
        }
        ++bound_it;        
    }
    
    // for (int k=0; k < frame_topo_.outerSize(); ++k)
    //     for (topology_t::InnerIterator it(frame_topo_, k); it; ++it)
    //         if (it.value())
    //         {
    //             // If the edge connects two boundary nodes
    //             // the edge thickness lower bound is the min radius
    //             if (boundary_map_[it.row()] > 0 && boundary_map_[it.col()] > 0)
    //             {
    //                 lbs[bound_it] = bounds_.min_radius;
    //                 ubs[bound_it] = bounds_.max_radius;
    //             }
    //             // Otherwise the thickness can vanish
    //             else
    //             {
    //                 lbs[bound_it] = 0.;//bounds_.min_radius;
    //                 ubs[bound_it] = bounds_.max_radius;
    //             }
    //             ++bound_it;
    //         }
    if (bound_it != problem_dim) throw std::runtime_error("Problem dim is wrong!");
    
    // SET BOUNDS
    this->set_lb(lbs);
    this->set_ub(ubs);
}


thermo_frame::structure_t
thermo_frame::encode_genes(const decision_vector & genes) const
{
    if (genes.size() != this->get_dimension())
        throw std::runtime_error("Decision vector size must match problem size!");
    
    
    // Copy the base structure
    auto frame = base_frame_;

    // Init gene iterator
    auto gene_it = 0u;
    
    for (const auto v : boost::make_iterator_range(vertices(frame)))
    {
        // Get boundary info
        bool is_boundary_row, is_boundary_col;
        std::tie(is_boundary_row, is_boundary_col) = is_boundary(v);
        
        // Apply jitter
        eva::fixed_vector<2> jitter = eva::fixed_vector<2>::Zero();
        
        if (!is_boundary_col) jitter[0] = genes[gene_it++]; // dx
        if (!is_boundary_row) jitter[1] = genes[gene_it++]; // dy
        
        frame[v].coords += jitter;
    }


    for (auto e : make_iterator_range(edges(frame)))
    {
         constexpr double pi = boost::math::constants::pi<double>();
         const auto r = genes[gene_it++];
        
         auto & element = frame[e];
         element.A = pi * r*r;
         element.I = pi/2 * r*r*r*r;
    }
    
    // The remaining card(topology) genes are the element thicknesses
    // for (int k=0; k < frame_topo_.outerSize(); ++k)
    //     for (topology_t::InnerIterator it(frame_topo_, k); it; ++it)
    //         if (it.value())
    //         {
    //             // Init element & insert it
    //             constexpr double pi = boost::math::constants::pi<double>();
    //             const auto r = genes[gene_it];//2.e-3;
                
    //             auto element = element_t();
    //             element.E = phys_props_.E;
    //             element.A = pi * r*r;
    //             element.I = pi/2 * r*r*r*r;
    //             element.k = phys_props_.k;
                    
    //             add_edge(it.row(), it.col(), element, frame);
    
                
    //             // Increment gene iterator
    //             ++gene_it;
    //         }

    // Return frame
    return frame;
}


std::pair<bool, bool>
thermo_frame::is_boundary(const size_t vertex_idx) const
{
    // bool is_boundary_row = false;
    // bool is_boundary_col = false;

    switch (boundary_map_.at(vertex_idx))
    {
    default:
        throw std::invalid_argument("Boundary code not supported!");
        return {false, false};
        
    case 0: return {false, false};   
        // is_boundary_row = false;
        // is_boundary_col = false;
        // break;
            
    case 1: return {true, false};
        // is_boundary_row = true;
        // is_boundary_col = false;
        // break;

    case 2: return {false, true};
        // is_boundary_row = false;
        // is_boundary_col = true;
        // break;

    case 3: return {true, true};
        // is_boundary_row = true;
        // is_boundary_col = true;
        // break;
    }

    // return {is_boundary_row, is_boundary_col};
}


void thermo_frame::objfun_impl(
    fitness_vector & f,
    const decision_vector & x) const
{
    // f = (compliance, )

    // Build structure
    // FIXME!!!! build structure either here or in compute_constraints_impl
    auto frame = encode_genes(x);

    // Solve structural problem & compute compliance
    using struct_solver_t = eva::sparse_solver_params<>;
    auto struct_res  = solve(frame, struct_kind_t(), struct_solver_t());
    auto compliance  = eva::compute_compliance(frame, struct_res);

    // Compute volume
    auto volume = compute_mass(frame, 1.);
    
    // Solve thermal problem & compute max temperature
    using thermo_solver_t = eva::dense_solver_params<thermo_kind_t::default_dense_solver_t>;
    auto thermo_results = solve(frame, thermo_kind_t(), thermo_solver_t());
    auto max_t = eva::get_max_temperature(frame, thermo_results);


    // Set fitness tuple
    f[0] = compliance;
    f[1] = volume;
    f[2] = max_t;
    
    if (compliance < 0. || volume < 0.)
    {
        // Something went wrong
        std::cerr << "Something's wrong: negative compliance." 
                  << "Saving structure as 'wrong.dot'...\n";

        eva::write_graphviz(frame, "wrong.dot");
        
        // Set fitness tuple to unfeasible-like value
        auto err_val = std::numeric_limits<double>::max();
        f[0] = err_val;
        f[1] = err_val;
        f[2] = err_val;
    }
}


// void thermo_frame::compute_constraints_impl(
//         constraint_vector & constraints,
//         const decision_vector & x) const
// {
//     // Build structure
//     // FIXME!!!! build structure either here or in objfun_impl
//     auto frame = encode_genes(x);

//     // Compute volume (use rho=1)
//     auto volume = compute_mass(frame, 1.);

//     // Write constraints
//     // NB: h(x) <= 0 is OK
//     constraints[0] = volume - bounds_.max_volume;
//     constraints[1] = bounds_.min_volume - volume;
// }



// thermo_frame::structure_t
// thermo_frame::make_base_frame(
//     const thermo_frame::joint_grid_t & grid,
//     const thermo_frame::topology_t & topology,
//     const thermo_frame::config_map_t & config)
// {
//     auto E = config["E"].as<eva::real>(); assert(E >= 0.);
//     auto k = config["k"].as<eva::real>(); assert(k >= 0.);
    
//     auto structure = thermo_frame::structure_t();

//     const auto m = grid.size();
//     const auto n = grid.front().size();
    
//     for (auto i = 0u; i < m; ++i)
//         for (auto j = 0u; j < n; ++j)
//             add_vertex(grid[i][j], structure);

//     // Add elements to the structure
//     for (int kk = 0; kk < topology.outerSize(); ++kk)
//         for (thermo_frame::topology_t::InnerIterator it(topology, kk); it; ++it)
//             if (it.value())
//             {
//                 auto r = 2.e-3;  
//                 constexpr double pi = boost::math::constants::pi<double>();

//                 auto element = thermo_frame::element_t();
//                 element.E = E;
//                 element.A = pi * r*r;
//                 element.I = pi/2 * r*r*r*r;
//                 element.k = k;
                
//                 add_edge(it.row(), it.col(), element, structure);
//             }
    
//     return structure;
// }

// thermo_frame::boundary_map_t
// thermo_frame::make_boundary_map(const thermo_frame::joint_grid_t & grid)
// {
//    // Build structure and tag nodes indicating on 
//     // which boundary they lie:
//     // -] 0 -> internal node
//     // -] 1 -> lies on a boundary row
//     // -] 2 -> lies on a boundary column
//     // -] 3 -> 1+2 hold
//     auto vertex_idx   = 0u;
//     auto boundary_map = boost::container::flat_map<size_t, int>();
    
//     const auto m = grid.size();
//     const auto n = grid.front().size();  
    
//     for (auto i = 0u; i < m; ++i)
//         for (auto j = 0u; j < n; ++j)
//         {
//             // Init value as 0 by default
//             auto & value = boundary_map[vertex_idx];
//             value = 0;
            
//             if (i == 0 || i == m-1) value += 1;
//             if (j == 0 || j == n-1) value += 2;

//             // Also keep track of nodes with loads applied ?
//             ++vertex_idx;
//         }
//     return boundary_map;
// }



// void thermo_frame::add_bcs(
//     thermo_frame::joint_grid_t & grid,
//     const thermo_frame::config_map_t & config)
// {    
//     const auto m = grid.size();
//     const auto n = grid.front().size();
    
//     // Fix first and last columns joints (set bcs)
//     for (auto i = 0u; i < m; ++i)
//     {
//         //                      x , y , theta
//         grid[i][ 0 ].bcs << 0., 0., 0.;
//         grid[i][n-1].bcs << 0., 0., 0.;
//     }
    
//     // Apply homogenous Dirchlet BC on the outer elements
//     for (auto i = 0u; i < m; ++i)
//     {
//         grid[i][ 0 ].T_bc = 0.;
//         grid[i][n-1].T_bc = 0.;
//     }
// }



// void thermo_frame::add_loads(
//     thermo_frame::joint_grid_t & grid,
//     const thermo_frame::config_map_t & config)
// { 
//     const auto m = grid.size();
//     const auto n = grid.front().size();
    
//     // Apply mechanical load on joint [lowest_row][3] (lowest row)
//     grid[0][3].load << 0, 16 * 1.e4;

//     // Apply thermal load on joint [highest_row][2]
//     grid[m-1][2].flux_bc = 0.5 * 1.3;
// }


// thermo_frame::joint_grid_t make_grid(eva::real length, eva::real height, long m, long n)
// {
//     assert(Joint().coords.size() == 2);
        
//     if (m < 2) throw std::invalid_argument("Parameter m must be >= 2");
//     if (n < 2) throw std::invalid_argument("Parameter n must be >= 2");
    
//     // Prealloc joint rows
//     auto joints = thermo_frame::joint_grid_t(m);
//     joints.reserve(m);

//     // Compute grid x and y spacings
//     auto delta_x = length / (n - 1);
//     auto delta_y = height / (m - 1);

//     // // Add some randomess
//     // std::mt19937 rng;
//     // rng.seed(std::random_device()());
//     // auto dx = std::uniform_real_distribution<float>(-delta_x / 2, delta_x / 2);
//     // auto dy = std::uniform_real_distribution<float>(-delta_y / 2, delta_y / 2);

//     for (auto i = 0u; i < m; ++i)
//     {        
//         // Reserve space
//         auto n_nodes = n;
//         joints[i].resize(n_nodes);
        
//         for (auto j = 0u; j < n_nodes; ++j)
//         {
//             // Init coordinates
//             auto x = j*delta_x;
//             auto y = i*delta_y;

//             // if (jitter)
//             // {                
//             // auto is_boundary_row = i == 0 || i == m-1;
//             // auto is_boundary_col = j == 0 || j == n-1;

//             // // Add random jitter
//             // if (!is_boundary_row) y += 0.8*dy(rng);
//             // if (!is_boundary_col) x += 0.8*dx(rng);
//             // }
            
            
//             // Set joint coords
//             joints[i][j].coords << x, y;
//         }
//     }

//     return joints;
// }



/// Clone method.
base_ptr thermo_frame::clone() const
{
    return base_ptr(new thermo_frame(*this));
}


std::string thermo_frame::get_name() const
{
    return "Clamped beam problem";
}

std::string thermo_frame::human_readable_extra() const
{
    std::ostringstream oss;
    oss << "\n\tPrint something...";
    return oss.str();
}


}} //end namespaces
