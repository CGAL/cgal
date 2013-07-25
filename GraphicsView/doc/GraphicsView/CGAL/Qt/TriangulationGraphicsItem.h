namespace CGAL {
namespace Qt {
/*!
\ingroup PkgGraphicsViewGraphicsItemClasses

An object of type `TriangulationGraphicsItem` is a graphics item that encapsulates a 2D triangulation. 

\tparam T must be a 2D triangulation class. 

*/
template< typename T >
class TriangulationGraphicsItem : public Qt::GraphicsItem {
public:

/// \name Creation 
/// @{

/*!
Constructs a graphics 
item for a triangulation pointed by `t`. 
*/ 
TriangulationGraphicsItem<T>(T* t); 

/// @} 

/// \name Operations 
/// @{

/*!
Returns the pen used to draw vertices. 
*/ 
QPen verticesPen()() const; 

/*!
Sets the pen used to draw vertices. 
*/ 
void setVerticesPen()(const QPen& p); 

/*!
Returns the pen used to draw edges. 
*/ 
QPen edgesPen()() const; 

/*!
Sets the pen used to draw edges. 
*/ 
void setEdgesPen()(const QPen& p); 

/*!
Returns `true`, iff vertices are drawn. 
*/ 
bool visibleVertices(); 

/*!
Sets the property. 
*/ 
void setVisibleVertices(bool b); 

/*!
Returns `true`, iff edges are drawn. 
*/ 
bool visibleEdges(); 

/*!
Sets the property. 
*/ 
void setVisibleEdges(bool b); 

/// @}

}; /* end TriangulationGraphicsItem */
} /* end namespace Qt */
} /* end namespace CGAL */
