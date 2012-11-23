
namespace CGAL {

/*!
\ingroup PkgBooleanSetOperations2

The class `Polygon_with_holes_2` models the concept `GeneralPolygonWithHoles_2`. 
It represents a linear polygon with holes. It is parameterized with two 
types (`Kernel` and `Container`) that are used to instantiate 
the type `Polygon_2<Kernel,Container>`. The latter is used to 
represents the outer boundary and the boundary of the holes (if any exist). 

\cgalModels `GeneralPolygonWithHoles_2`

*/
template< typename Kernel, typename Container >
class Polygon_with_holes_2 {
public:

/// @}

}; /* end Polygon_with_holes_2 */

/*!
This operator imports a polygon with holes from the input stream `in`.

An ASCII and a binary format exist. The stream detects the format
automatically and can read both.

The format consists of the number of points of the outer boundary followed 
by the points themselves in counterclockwise order, followed by the number of holes,
and for each hole, the number of points of the outer boundary is followed 
by the points themselves in clockwise order.

\relates Polygon_with_holes_2
*/
template <class Kernel, Class Container>
std::istream& operator>>(std::istream& in, CGAL::Polygon_with_holes_2<Kernel, Container>& P);


/*!
This operator exports a polygon with holes to the output stream `out`.

An ASCII and a binary format exist. The format can be selected with
the \cgal modifiers for streams, `set_ascii_mode()` and `set_binary_mode()`
respectively. The modifier `set_pretty_mode()` can be used to allow for (a
few) structuring comments in the output. Otherwise, the output would
be free of comments. The default for writing is ASCII without
comments.

The number of points of the outer boundary is exported followed by the
points themselves in counterclockwise order. Then, the number of holes
is exported, and for each hole, the number of points on its outer
boundary is exported followed by the points themselves in clockwise
order.

\relates Polygon_with_holes_2
*/
template <class Polygon>
std::ostream& operator<<(std::ostream& out, CGAL::Polygon_with_holes_2<Kernel, Polygon>& P);

} /* end namespace CGAL */
