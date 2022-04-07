namespace CGAL {
namespace Qt {
/*!
\ingroup PkgGraphicsViewGraphicsItemClasses

An object of type `CircularArcGraphicsItem` is a graphics item that encapsulates a circular arc.

\tparam CK  must be a \cgal `CircularKernel`.

*/
template< typename CK >
class CircularArcGraphicsItem : public Qt::GraphicsItem {
public:

/// \name Creation
/// @{

/*!
Constructs a graphics
item for a circular arc.
*/
CircularArcGraphicsItem<CK>();

/// @}

/// \name Operations
/// @{

/*!
Returns the pen used to draw edges.
*/
QPen edgesPen()() const;

/*!
Sets the pen used to draw edges.
*/
void setEdgesPen()(const QPen& p);

/*!
Sets the circular arc.
*/
void setArc(const Circular_arc_2<CK>& ca2);

/// @}

}; /* end CircularArcGraphicsItem */
} /* end namespace Qt */
} /* end namespace CGAL */
