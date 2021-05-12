namespace CGAL {
namespace Qt {
/*!
\ingroup PkgGraphicsViewGraphicsItemClasses

The \cgal graphics items hold pointers to \cgal data structures. When a data structure changes,
a signal is emitted. An object derived from type `GraphicsItem` must provide an implementation of the
virtual slot `modelChanged()`. This typically triggers redrawing, as well as recomputation
of the bounding box which in turn may trigger redrawing of overlapping graphics items.

*/
class GraphicsItem : public QGraphicsItem {
public:

/// \name Slots
/// @{

/*!
This slot must be provided by derived classes.
*/
virtual void modelChanged();

/// @}

}; /* end GraphicsItem */
} /* end namespace Qt */
} /* end namespace CGAL */
