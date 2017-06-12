# include "pagmo_utils.hpp"

namespace utils {


int get_champion_island_idx(const pagmo::archipelago & archi)
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

}//end namespace
