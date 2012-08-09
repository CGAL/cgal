namespace CGAL {
namespace Qt {
/*!
\ingroup PkgGraphicsView

An object of type `GraphicsViewIsoRectangleInput` creates an axis parallel rectangle. 

Parameters 
-------------- 

The template parameter of `GraphicsViewIsoRectangleInput` must be a \cgal `Kernel`. 

*/
template< typename K >
class GraphicsViewIsoRectangleInput : public GraphicsViewInput {
public:

/// \name Creation 
/// @{

/*! 
`p` is a parent object. `s` is the scene where the iso rectangle is generated. 
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
