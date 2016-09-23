// eva
# include "continuous_structural_problem.hpp"
# include "vtk.hpp"
// pagmo
# include <pagmo/src/population.h>
# include <pagmo/src/algorithm/sga.h>
# include <pagmo/src/algorithm/cstrs_co_evolution.h>
# include <pagmo/src/algorithm/gsl_fr.h>

int main(int argc, char * argv [])
{
    // Config types
    using structure_t = eva::frame2d;
    using joint_t     = typename eva::joint_of  <structure_t>::type;
    using element_t   = typename eva::element_of<structure_t>::type;
    
    // Setup structure points
    auto joints = std::vector<joint_t>();
    auto imax = 14u;
    auto jmax = 14u;
    joints.reserve(imax*jmax);
    
    for (auto i = 0u; i < imax; ++i)
        for (auto j = 0u; j < jmax; ++j)
        {
            // Setup properties
            auto joint = joint_t ();
            joint.coords << i, j;
            if (j == 0)
            {
                if (i == 0 || i == imax-1)
                    joint.bcs << 0., 0., 0.;
                else
                    joint.load << 0., 1.;
            }
            // Add joint
            joints.emplace_back(std::move(joint));
        }

    // Setup feasible connections //
    auto max_dist = 1.;
    auto rule = [&](auto& j1, auto& j2){

        auto diff = eva::fixed_vector<2>(j1.coords-j2.coords);
        auto dist = diff.lpNorm<Eigen::Infinity>();

        if (dist > max_dist) return false;

        // No co-linear links with length > 1.
        if (dist > 1. && (diff[0] == 0. || diff[1] == 0.)) return false;

        // No 45deg links with length > 1.
        if (dist > 1. && diff[1] / diff[0] == 1.) return false;

        return true;
    };
    auto gene_mask = eva::utils::make_topology(joints, rule);
    
    // Build problem
    auto prb_dim = 0u;
    for (const auto& el: gene_mask) prb_dim += el.size();

    // Setup problem
    pagmo::problem::structural_problem prb_base(prb_dim, joints, gene_mask, 20);
    pagmo::problem::death_penalty prb(
        prb_base,
        // pagmo::problem::structural_problem (prb_dim, joints, gene_mask, 20),
        pagmo::problem::death_penalty::SIMPLE
        );

    // Setup population
    pagmo::population pop(prb, 1000);
    
    // Setup algorithm
    // pagmo::algorithm::sga alg(500, 0.9, 0.1);
    auto alg = pagmo::algorithm::gsl_fr();
    
    
    // pagmo::algorithm::cstrs_co_evolution alg(pagmo::algorithm::sga(500),
    //                                          pagmo::algorithm::sga(500),
    //                                          50, 500);
    alg.set_screen_output(true);

    // Solve
    auto best_structure = structure_t();
    
    // for (auto i = 0u; i < 10u; ++i)
    // {
        alg.evolve(pop);

        // Print champion
        std::cout << "Done .. the champion is ... \n" << std::endl;
        auto champion = pop.champion();
    
        auto nr_links = 0u;
        for (auto el : champion.x)
        {
            if (el) nr_links++; 
            std::cout << el << " ";
        }
        std::cout << "\n";
        // std::cout << "Nr of links = " << nr_links << std::endl;
        // std::cout << "Feasibility : " << prb.feasibility_x(pop.champion().x) << std::endl;
        best_structure = prb_base.encode_genes(champion.x);
        display(best_structure);
    // }

    // Save to vtu for paraview
    write_vtu(best_structure, solve(best_structure), "best.vtu");
    
    return 0;
}



    
 

    
    
    
