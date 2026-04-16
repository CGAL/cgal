namespace CGAL {

/*!
\ingroup PkgFeatureGraphRef

\brief The feature graph structure represent polylines using points in space.

\details The feature graph is index-based,
meaning that the lines are described
using the indices of the points in the point range.

\tparam Point_d is the type of points.

*/

template <typename Point_d>
struct Feature_graph
{
public:

/// \name Types
/// @{

/// The point type
typedef Point_d Point;
/// The type for multiple points
typedef unspecified_type Point_range;
/// The point index type
typedef std::size_t Point_index;
/// The line of point index type
typedef unspecified_type Line;
/// The type for multiple lines
typedef unspecified_type Line_range;

/// @}

  Point_range points;
  Line_range lines;
};

// TODO : add infos in run time (not compile time) with property maps

// Feature_graph feature_graph;
// const Line& single_line_access = feature_graph.get_line(line_index);
// const auto& point_map = feature_graph.get_point_map();
// for (const Line &line : feature_graph.lines) // for each on lines
//   for (const Point_index& point_index : line)
//    const Point &point = feature_graph.get_point(point_index); // or
//    const Point &point = get(point_map, point_index);
//    const InfoA& info_a = get(info_a_map, point_index); // outside of feature_graph (via property maps ?)
//    const InfoB& info_b = get(info_b_map, point_index); // outside of feature_graph (via property maps ?)

// corner iterator ? Iterating on corner is not trivial, need a marker to avoid duplicates

} /* namespace CGAL */