    
    // Make random number generator(s) //
    auto generator = std::default_random_engine();
    auto realdistr = std::uniform_real_distribution<eva::real>();
    auto intdistr  = std::uniform_int_distribution<int>(0,1);
    
    auto make_die = [&](auto&& d) { return [&](){return d(generator); }; };
    auto real_die = make_die(realdistr);
    auto int_die  = make_die(intdistr);
    
    // Add nodes
    auto coords = std::vector<eva::real> {0., 4.}
    auto n = 5u; //nuber of nodes
    for (auto i = 0u; i < n; ++i)
        //           |         coords        |   load   |  bcs  |
        add_vertex({ {real_die(), real_die()}, {1., 1.}, {0., 0.} }, truss);
    
    // Add egdes
    for(auto&& i : make_iterator_range(vertices(truss)))
        for(auto&& j : make_iterator_range(vertices(truss)))
            if (int_die()) add_edge(i, j, { 1., 2. }, truss);
