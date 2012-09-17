namespace CGAL {
namespace Qt {
/*!
\ingroup PkgGraphicsView

An object of type `GraphicsViewCircularArcInput` creates a circular arc, defined by three points on a circle. A new 
vertex is inserted every time the left mouse button is released. 
The `Del` key removes the last entered point. The `Esc` 
key removes all entered points. 

\tparam K must be a \cgal `CircularKernel`. 

*/
template< typename K >
class GraphicsViewCircularArcInput : public Qt::GraphicsViewInput {
public:

/// \name Creation 
/// @{

/*! 
`p` is a parent object. `s` is the scene where the circular arc is generated. 
*/ 
GraphicsViewCircularArcInput<T>(QObject *p, QGraphicsScene* s); 

/// @} 

/// \name Signals 
/// @{

/*! 
The object `o` contains a `Circular_arc_2<K>`. 
*/ 
void generate(Object o); 

/// @}

}; /* end GraphicsViewCircularArcInput */
} /* end namespace Qt */
} /* end namespace CGAL */
