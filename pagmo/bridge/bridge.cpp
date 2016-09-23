// eva
# include "continuous_structural_problem.hpp"
# include "vtk.hpp"
// pagmo
# include <pagmo/src/population.h>
# include <pagmo/src/algorithms.h>
# include <pagmo/src/archipelago.h>
# include <pagmo/src/topologies.h>

int get_champion_island_idx(const pagmo::archipelago& archi)
{
    int idx = 0;
    auto min_glob = archi.get_island(0)->get_population().champion().f[0];
    auto problem  = archi.get_island(0)->get_population().problem().clone();
    
	for (size_t i = 1;i<archi.get_size(); ++i)
    {
		auto champ_loc = archi.get_island(i)->get_population().champion();
		if (champ_loc.f[0] < min_glob && problem->feasibility_x(champ_loc.x))
        {
			min_glob = champ_loc.f[0];
			idx = i;
		}
	}
	return idx;
}

int main(int argc, char * argv [])
{
    // Config types //
    using structure_t = eva::frame2d;
    using joint_t     = typename eva::joint_of  <structure_t>::type;
    using element_t   = typename eva::element_of<structure_t>::type;
    
    // Setup structure points //
    auto joints = std::vector<joint_t>();
    auto imax = 20u;
    auto jmax = 8u;
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
                    joint.load << 0., -1.e3;
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
    
    // Build problem //
    auto prb_dim = 0u;
    for (const auto& el: gene_mask) prb_dim += el.size();
    auto prob = pagmo::problem::continuous_structural_problem(prb_dim, joints, gene_mask, 50.e-2);//11.e-2);
    
    // Setup algorithm //
    auto algo_original = pagmo::algorithm::jde(/*#generations=*/100);
    auto algo_coevo    = pagmo::algorithm::jde(/*#generations=*/100);
    // auto algo = pagmo::algorithm::cstrs_co_evolution(algo_original, algo_coevo,
    //                                                  /*#coevo individuals=*/100,
    //                                                  /*#tot generations=*/3);
    auto algo = pagmo::algorithm::cstrs_self_adaptive(algo_original, /*#generations=*/600);
    
    algo.set_screen_output(true);


    // Setup archipelago //
    auto topo  = pagmo::topology::ring();
    auto archi = pagmo::archipelago(/*algorithm=*/algo,
                                    /*problem=*/prob,
                                    /*#islands=*/6,
                                    /*#individuals*/100);
    
    // Solve //
    for (auto i = 0u; i < 5; ++i)
    {
        std::cout << "Starting generation " << i << "\n";
        archi.evolve(/*#migrations=*/1);
        archi.join();
    }

    // Retrieve champion //
    auto champ_isl_idx = get_champion_island_idx(archi);
    auto champion = archi.get_island(champ_isl_idx)->get_population().champion();

    // Print infos //
    std::cout << "Done .. the champion is ... \n" << std::endl;
    std::cout << champion.x << std::endl;
    
    std::cout << "Fitness = " << champion.f << std::endl;
    std::cout << "Feasibility = " << champion.c << std::endl;

    // Build best structure //
    auto best_structure = prob.encode_genes(champion.x);
    auto results = solve(best_structure);
    
    // Display&Save //
    display(best_structure);
    write_vtu(best_structure, results, "best.vtu");
    
    return 0;
}


// Compliance = 0.207616
//                  auto best_genes = std::vector<double> {
//     0.992716, 0.800316, 0.881268, 0.013468, 0.0107477, 0.0273263, 0.03187, 0.010434,
//     0.0028074, 0.326523, 0.997589, 0.976833, 0.033075, 0.924048, 0.00154099, 0.105252,
//     0.0833592, 0.0623852, 0.454175, 0.829388, 0.0032706, 0.00736974, 0.990774, 0.0115424,
//     0.0275509, 0.0315993, 0.00276079, 0.840914, 0.00196601
// };

    
 

    
    
    
