namespace CGAL {
namespace Qt {
/*!
\ingroup PkgGraphicsViewGraphicsItemClasses

An object of type `VoronoiGraphicsItem` is a graphics item that 
encapsulates a Delaunay triangulation in order to draw its dual, the 
Voronoi diagram. 

\tparam DT must be a 2D Delaunay triangulation class. 

*/
template< typename DT >
class VoronoiGraphicsItem : public Qt::GraphicsItem {
public:

/// \name Creation 
/// @{

/*!
Constructs a graphics item for the dual of the 
Delaunay triangulation `dt`. 
*/ 
VoronoiGraphicsItem<DT>(DT* dt); 

/// @} 

/// \name Operations 
/// @{

/*!
Returns the pen used to draw edges. 
*/ 
QPen edgesPen()() const; 

/*!
Set the pen used to draw edges. 
*/ 
void setEdgesPen()(const QPen& p); 

/*!
Returns `true`, iff edges are drawn. 
*/ 
bool visibleEdges(); 

/*!
Set the property. 
*/ 
bool setVisibleEdges(bool b); 

/// @}

}; /* end VoronoiGraphicsItem */
} /* end namespace Qt */
} /* end namespace CGAL */
