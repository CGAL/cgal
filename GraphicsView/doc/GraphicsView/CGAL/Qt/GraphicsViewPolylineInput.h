namespace CGAL {
namespace Qt {
/*!
\ingroup PkgGraphicsViewInputClasses

An object of type `GraphicsViewPolylineInput` creates a list of points. A new
vertex is inserted every time the left mouse button is pressed.
The list of points is emitted on a right click or when the number of
points specified in the constructor is reached. You can use the 'Del'
or 'Backspace' key if you want to remove your last entered point in the polygon,
and the 'Esc' key if you want to remove all points.

The tool can serve at the same time for entering a single point,
a polyline with a given number of points, and for open as well as closed
polylines.

For polylines the segment between the last entered point and the current
mouse position is only drawn correctly when mouse tracking is enabled
in the graphics view. The same holds for closed polygons.

\tparam K must be a model of `Kernel`.

*/
template< typename K >
class GraphicsViewPolylineInput : public Qt::GraphicsViewInput {
public:

/// \name Creation
/// @{

/*!
\param p is a parent object.
\param s is the scene where the polyline is generated.
\param n is the number of points of the polyline to be generated. If `n = 0`,
that is the default value, the number of points of the polyline is not
limited.
\param closed determines if a closed polygon is drawn during the input.
*/
GraphicsViewPolylineInput<T>(QObject *p, QGraphicsScene* s, int n = 0,
bool closed = true);

/// @}

/// \name Signals
/// @{

/*!
The object `o` contains a `std:list<K::Point_2>`.
*/
void generate(CGAL::Object o);

/// @}

}; /* end GraphicsViewPolylineInput */
} /* end namespace Qt */
} /* end namespace CGAL */
