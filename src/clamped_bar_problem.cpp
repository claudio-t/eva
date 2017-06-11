# include "clamped_bar_problem.hpp"


namespace pagmo { namespace problem {

clamped_bar::clamped_bar(
    const structure_t & base_frame,
    const clamped_bar::topology_t & frame_topo,
    const clamped_bar::boundary_map_t & boundary_map,
    const int problem_dim,
    // const double length,
    // const double height,
    // int num_rows,
    // int num_cols,
    const clamped_bar::bounds bnds,
    const physical_properties phys_props
    )
    : base(problem_dim, // Global dimension
           int(0),   // Integer dimension
           int(2),   // Fitness dimension
           int(1),   // Global constraints dimension
           int(1),   // Inequality constraints dimension
           0.)       // Constraints tolerance
    , base_frame_(base_frame)
    , frame_topo_(frame_topo)
    , boundary_map_(boundary_map)
    , bounds_(bnds)
    , phys_props_(phys_props)
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
    for (int k=0; k < frame_topo_.outerSize(); ++k)
        for (topology_t::InnerIterator it(frame_topo_, k); it; ++it)
            if (it.value())
            {
                if (boundary_map_[it.row()] > 0 && boundary_map_[it.col()])
                {
                    lbs[bound_it] = bounds_.min_radius;
                    ubs[bound_it] = bounds_.max_radius;
                }
                else
                {
                    lbs[bound_it] = 0.;//bounds_.min_radius;
                    ubs[bound_it] = bounds_.max_radius;
                }
                ++bound_it;
            }
    if (bound_it != problem_dim) throw std::runtime_error("Problem dim is wrong!");
    // for (; bound_it < problem_dim; ++bound_it)
    // {
    //     lbs[bound_it] = 0.;//bounds_.min_radius;
    //     ubs[bound_it] = bounds_.max_radius;       
    // }
    
}


clamped_bar::structure_t
clamped_bar::encode_genes(const decision_vector & genes) const
{
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
        auto jitter = eva::fixed_vector<2>();
        
        if (!is_boundary_col) jitter[0] = genes[gene_it++]; // dx
        if (!is_boundary_row) jitter[1] = genes[gene_it++]; // dy
        
        frame[v].coords += jitter;
    }

    
    // The remaining card(topology) genes are the element thicknesses
    for (int k=0; k < frame_topo_.outerSize(); ++k)
        for (topology_t::InnerIterator it(frame_topo_, k); it; ++it)
            if (it.value())
            {
                // Init element & insert it
                constexpr double pi = boost::math::constants::pi<double>();
                const auto r = genes[gene_it];//2.e-3;

                // If radius is big enough add it
                if (r >= bounds_.min_radius)
                {
                    auto element = element_t();
                    element.E = 69.e9; // Steel ASTM-A36
                    element.A = pi * r*r;
                    element.I = pi/2 * r*r*r*r;
                    element.k = 1.0;
                    
                    add_edge(it.row(), it.col(), element, frame);
                }
                
                // Increment gene iterator
                ++gene_it;
            }

    // Return frame
    return frame;
}


std::pair<bool, bool>
clamped_bar::is_boundary(const size_t vertex_idx) const
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


void clamped_bar::objfun_impl(
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
    // using thermo_solver_t = eva::dense_solver_params<thermo_kind_t::default_dense_solver_t>;
    // auto thermo_results = solve(frame, thermo_kind_t(), thermo_solver_t());
    // auto max_t = get_max_temperature(thermo_results);

    // Assemble result tuple
    f[0] = compliance;
    f[1] = volume;
    // f[1] = ...
}


void clamped_bar::compute_constraints_impl(
        constraint_vector & constraints,
        const decision_vector & x) const
{
    // Build structure
    // FIXME!!!! build structure either here or in objfun_impl
    auto frame = encode_genes(x);

    // Compute volume (use rho=1)
    auto volume = compute_mass(frame, 1.);

    // Write constraint
    // NB: h(x) <= 0 is OK
    constraints[0] = volume - bounds_.min_volume;
}




/// Clone method.
base_ptr clamped_bar::clone() const
{
    return base_ptr(new clamped_bar(*this));
}


std::string clamped_bar::get_name() const
{
    return "Clamped beam problem";
}

std::string clamped_bar::human_readable_extra() const
{
    std::ostringstream oss;
    oss << "\n\tPrint something...";
    return oss.str();
}


}} //end namespaces
