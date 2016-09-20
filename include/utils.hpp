# ifndef __EVA_UTILS__
# define __EVA_UTILS__

# include <boost/graph/depth_first_search.hpp>
# include <boost/graph/breadth_first_search.hpp>

namespace eva { namespace utils {


/// Creates a connectivity map adding a link between two points
/// whenever the rule is satisfied.
template <typename J, typename F>
std::vector< std::vector<index_t> >
make_topology(const std::vector<J>& joints, const F& functor_rule);


/// Checks for an existing path between a source and a target node.
template <typename S>
bool
check_path(typename S::vertex_descriptor src,
           typename S::vertex_descriptor trg,
           const S& structure);


/// Checks if the graph representing the structure is formed by single
/// connected component
template <typename S>
bool is_connected(const S& structure);


/// Provides a view of a subrange of contiguous elements. The object
/// from which the subrange gets extracted must provide a suitable
/// overload for the (random) access operator[].
template <typename Iterator>
class subrange_view;


// ---------------------------------------------------------------------------------------------- //

template <typename J, typename F>
std::vector< std::vector<index_t> >
make_topology(const std::vector<J>& joints, const F& rule)
{
    // Init empty ret var
    auto nr_joints = joints.size();
    auto ret = std::vector< std::vector<index_t> >(nr_joints);
    
    for (auto i = 0u; i < nr_joints; ++i)
        for (auto j = 0u; j < nr_joints; ++j)
            if (i != j && rule(joints[i], joints[j]))
                ret[i].emplace_back(j);
    
    return ret;
}

template <typename S>
bool is_connected(const S& structure)
{
    // Check for non empty structure
    if (num_vertices(structure) == 0)
        throw std::runtime_error("Structure must be non-empty!");
    
    // Perform a BFS search
    using vertex_t    = typename S::vertex_descriptor;
    using color_map_t = std::vector<boost::default_color_type>;

    // Allocate colors storage
    auto colors = color_map_t(num_vertices(structure));

    // Build color map using colors as storage
    auto colors_map = boost::make_iterator_property_map(
        colors.begin(),
        boost::get(boost::vertex_index, structure)
        );
    
    // Do BFS
    auto root_v = *vertices(structure).first; // Search starts here 
    auto queue  = boost::queue<vertex_t> (); // Queue buffer
    
    boost::breadth_first_search(structure, root_v, queue,
                                boost::default_bfs_visitor(),
                                colors_map);

    // If there exists an unreachable node (color != black)
    // then the graph is not connected
    auto black = boost::color_traits<boost::default_color_type>::black();
    for (auto color : colors) if (color != black) return false;
    
    // Otherwise the graph is connected
    return true;
}


namespace details {
/// Custom exception that is thrown when bfs finds a given vertex
class vertex_reached_exception : public std::exception {
private:
public:
		/// Throws everytime
		~vertex_reached_exception() throw() {};
};

/**
 * @brief Visitor that throws a vertex_reached_exception when the
 * target vertex gets found
 * @tparam S The structure type
 */
template <class S>
class ex_visitor : public boost::default_bfs_visitor {
public:
    /// Alias for the structure vertex descriptor
    using vd_t = typename S::vertex_descriptor;
    
	/// Function that gets called by the bfs algorithm at each vertex
	void discover_vertex(vd_t v, const S& g) const {
        // std::cout << "Hi, I'm vertex " << v << "\n";
		if (v == target_)
        {
            // std::cout << "Found!\n";
            throw vertex_reached_exception();
        }
	}
    
	/// Constructs the object and sets the target vertex
	ex_visitor(vd_t target) : target_(target) {}

private:
	vd_t  target_; ///< Target vertex
	ex_visitor();  ///< Private default destructor
};



} // end namespace details


template <typename S>
bool
check_path(typename S::vertex_descriptor src,
           typename S::vertex_descriptor trg,
           const S& structure)
{
    using vertex_t    = typename S::vertex_descriptor;
    using color_map_t = std::vector<boost::default_color_type>;

    // Allocate colors storage
    auto colors = color_map_t(num_vertices(structure));

    // Build color map using colors as storage
    auto colors_map = boost::make_iterator_property_map(
        colors.begin(),
        boost::get(boost::vertex_index, structure)
        );
    
    // Do BFS
    auto path_exists = false;
    try
    {
        boost::queue<vertex_t> queue;
        boost::breadth_first_search(structure, src, queue,
                                    // boost::default_bfs_visitor(),
                                    details::ex_visitor<S>(trg),
                                    colors_map);
    }
    catch(details::vertex_reached_exception& ex)
    {
        path_exists = true;
    }
    return path_exists;
}


template <typename Iterator>
class subrange_view
{
public:
    /// Constructs the subrange in between begin and end iterators
    subrange_view(Iterator begin, Iterator end)
        : begin_(begin), end_(end) {}

    /// Returns iterator to the subrange beginning
    inline constexpr Iterator begin() { return begin_; }

    /// Returns iterator to the subrange end
    inline constexpr Iterator end() { return end_; }

    /// Returns the subrange width
    inline constexpr size_t size() { return end_ - begin_; }
    
    /// Provide access to elements belonging to the subrange (no check
    /// is performed)
    decltype(auto) operator[](size_t idx) { return begin_[idx]; }
    

private:
    Iterator begin_; ///< Iterator begin
    Iterator end_;   ///< Iterator end 
};


} } // end namespace eva & utils

# endif //__EVA_UTILS__
