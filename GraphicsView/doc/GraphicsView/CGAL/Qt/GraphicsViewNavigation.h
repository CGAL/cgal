namespace CGAL {
namespace Qt {
/*!
\ingroup PkgGraphicsViewMiscClasses

An object of type `GraphicsViewNavigation` can be added as event filter to a `Qt::QGraphicsView` and its viewport. 

Dragging the left mouse button while holding the 'Ctrl' key defines a zoom rectangle. 
Dragging the right mouse button while holding the 'Ctrl' key translates the scene. 
'Ctrl-Shift' and a click of the right mouse button translates what is under the mouse to the center. 

*/

class GraphicsViewNavigation {
public:

/// \name Operations 
/// @{

/*!
The event filter. 
*/ 
bool eventFilter(QObject *obj, QEvent *event); 

/// @} 

/// \name Signals 
/// @{

/*!
Emits the real world mouse coordinates. 
*/ 
void mouseCoordinates(QPointF p); 

/// @}

}; /* end GraphicsViewNavigation */
} /* end namespace Qt */
} /* end namespace CGAL */
