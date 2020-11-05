namespace CGAL {
namespace Qt {
/*!
\ingroup PkgGraphicsViewInputClasses

An object of type `GraphicsViewCircleInput` creates a circle, defined by either
center and radius, or two or three points on the circle. A new
vertex is inserted every time the left mouse button is released.
The `Del` key removes the last entered point. The `Esc`
key removes all entered points.



\tparam K must be a model of `Kernel`.

*/
template< typename K >
class GraphicsViewCircleInput : Qt::GraphicsViewInput {
public:

/// \name Creation
/// @{

/*!
\param p is a parent object.
\param s is the scene where the circle is generated.
\param pointsOnCircle is the
number of points on the circle to be generated, that is the default value `1`
corresponds to the case center/radius.
*/
GraphicsViewCircleInput<T>(QObject *p, QGraphicsScene* s, int pointsOnCircle = 1);

/// @}

/// \name Signals
/// @{

/*!
The object `o` contains a `std::pair<K::Point_2, K::FT>`
for center and radius, or a `std::pair<K::Point_2, K::Point_2>` for two points defining
the circle, or `CGAL::array<K::Point_2, 3>` for three points defining
the circle.
*/
void generate(Object o);

/// @}

}; /* end GraphicsViewCircleInput */
} /* end namespace Qt */
} /* end namespace CGAL */
