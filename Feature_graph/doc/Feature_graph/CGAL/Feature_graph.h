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

/*!
The point index type.
*/
typedef Point_d Point;
/*!
The point index type.
\cgalModels{Index,LessThanComparable,Hashable}
*/
typedef unspecified_type Point_index;

/*!
\brief The type of line containing point incides.
A model of <a href="https://www.boost.org/libs/range/doc/html/range/concepts/bidirectional_range.html">BidirectionalRange</a> with value type `Point_index`.
*/
typedef unspecified_type Line;

/*!
The line index type.
\cgalModels{Index,LessThanComparable,Hashable}
*/
typedef unspecified_type Line_index;

/*!
\brief The range over all points.
A model of <a href="https://www.boost.org/libs/range/doc/html/range/concepts/bidirectional_range.html">BidirectionalRange</a> with value type `Point`.
*/
typedef unspecified_type Point_range;
/*!
\brief The range over all lines.
A model of <a href="https://www.boost.org/libs/range/doc/html/range/concepts/bidirectional_range.html">BidirectionalRange</a> with value type `Line`.
*/
typedef unspecified_type Line_range;

/// @}

/// \name Accessors
/// @{

/*!
Return the iterator range of the points of the feature graph.
*/
Point_range points() const;
/*!
Return the iterator range of the lines of the feature graph.
*/
Line_range lines() const;

/*!
Return the point from its index.
\param point_index the index of the point to return.
*/
Point point(Point_index point_index) const;
/*!
Return the line from its index.
\param line_index the index of the line to return.
*/
Line line(Line_index line_index) const;

/// @}

private:
  Point_range points;
  Line_range lines;
};

} /* namespace CGAL */