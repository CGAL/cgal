namespace CGAL {
namespace Qt {
/*!
\ingroup PkgGraphicsViewGraphicsItemClasses

An object of type `ConstrainedTriangulationGraphicsItem` is a graphics item that encapsulates a constrained triangulation. 


\tparam CT must be a \cgal 2D constrained triangulation class. 

*/
template< typename CT >
class ConstrainedTriangulationGraphicsItem : public Qt::TriangulationGraphicsItem {
public:

/// \name Creation 
/// @{

/*!
Constructs 
a graphics item for triangulation `ct`. 
*/ 
ConstrainedTriangulationGraphicsItem<CT>(CT* ct); 

/// @} 

/// \name Operations 
/// @{

/*!
Returns the pen used to draw constraints. 
*/ 
QPen constraintsPen()() const; 

/*!
Sets the pen used to draw constraints. 
*/ 
void setConstraintsPen()(const QPen& p); 

/// @}

}; /* end ConstrainedTriangulationGraphicsItem */
} /* end namespace Qt */
} /* end namespace CGAL */
