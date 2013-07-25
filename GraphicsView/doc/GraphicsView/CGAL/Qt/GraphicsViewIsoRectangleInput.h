namespace CGAL {
namespace Qt {
/*!
\ingroup PkgGraphicsViewInputClasses

An object of type `GraphicsViewIsoRectangleInput` creates an axis parallel rectangle. 

\tparam K must be a model of `Kernel`. 

*/
template< typename K >
class GraphicsViewIsoRectangleInput : public GraphicsViewInput {
public:

/// \name Creation 
/// @{

/*!
\param p is a parent object. 
\param s is the scene where the iso rectangle is generated. 
*/ 
GraphicsViewIsoRectangleInput<T>(QObject *p, QGraphicsScene* s); 

/// @} 

/// \name Signals 
/// @{

/*!
The object `o` contains a `K::Iso_rectangle_2`. 
*/ 
void generate(Object o); 

/// @}

}; /* end GraphicsViewIsoRectangleInput */
} /* end namespace Qt */
} /* end namespace CGAL */
