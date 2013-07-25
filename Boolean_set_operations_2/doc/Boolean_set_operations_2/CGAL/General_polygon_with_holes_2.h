
namespace CGAL {

/*!
\ingroup PkgBooleanSetOperations2

\cgalModels `GeneralPolygonWithHoles_2`

*/
template< typename Polygon >
class General_polygon_with_holes_2 {
public:

/// \name Definition 
/// The class `General_polygon_with_holes_2` models the concept
/// `GeneralPolygonWithHoles_2`. It represents a general polygon with
/// holes. It is parameterized with a type `Polygon` used to define
/// the exposed type `General_polygon_2`. This type represents the
/// outer boundary of the general polygon and the outer boundaries of
/// each hole.
/// @{

/*!

*/ 
typedef Polygon General_polygon_2; 

/// @}

}; /* end General_polygon_with_holes_2 */

/*!
This operator imports a General_polygon_with_holes_2 from the input stream `in`.

An ASCII and a binary format exist. The stream detects the format
automatically and can read both.

The format consists of the number of curves of the outer boundary
followed by the curves themselves in counterclockwise order, followed
by the number of holes, and for each hole, the number of curves on its
outer boundary is followed by the curves themselves in clockwise
order.

\relates General_polygon_with_holes_2
*/
template <class Polygon>
std::istream& operator>>(std::istream& in, CGAL::General_polygon_with_holes_2<Polygon>& P);


/*!
This operator exports a General_polygon_with_holes_2 to the output stream `out`.

An ASCII and a binary format exist. The format can be selected with
the \cgal modifiers for streams, `set_ascii_mode(0` and `set_binary_mode()`
respectively. The modifier `set_pretty_mode()` can be used to allow for (a
few) structuring comments in the output. Otherwise, the output would
be free of comments. The default for writing is ASCII without
comments.

The number of curves of the outer boundary is exported followed by the
curves themselves in counterclockwise order. Then, the number of holes
is exported, and for each hole, the number of curves on its outer
boundary is exported followed by the curves themselves in clockwise
order.

\relates General_polygon_with_holes_2
*/
template <class Polygon>
std::ostream& operator<<(std::ostream& out, CGAL::General_polygon_with_holes_2<Polygon>& P);

} /* end namespace CGAL */
