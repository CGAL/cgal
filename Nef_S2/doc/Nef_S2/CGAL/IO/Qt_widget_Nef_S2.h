
namespace CGAL {

/*!
\ingroup PkgNefS2

The class `Qt_widget_Nef_S2` uses the OpenGL interface of Qt in order to 
display a 
`Nef_polyhedron_S2`. Its purpose is to provide an easy to use viewer for 
`Nef_polyhedron_S2`. There are no means provided to enhance the 
functionality of the viewer. 

In addition to the functions inherited from the Qt class `QGLWidget`, 
`Qt_widget_Nef_S2` only has a single public 
constructor. For the usage of `Qt_widget_Nef_S2` see the example 
below. 

### Parameters ###

The template parameter expects an instantiation of `Nef_polyhedron_S2<Traits>`. 

\sa `CGAL::Nef_polyhedron_S2<Traits>` 

### Example ###

This example creates some random `Sphere_segments` and distributes them on 
two `Nef_polyhedron_2`. The two Nef polyhedra are combined by a symmetric 
diffrence and the result is displayed in a Qt widget. 

\cgalExample{Nef_S2/nef_S2.cpp} 

*/
template< typename Nef_polyhedron_S2 >
class Qt_widget_Nef_S2 {
public:

/// \name Creation 
/// @{

/*! 
Creates a widget `W` for displaying the `Nef_polyhedron_S2` N. 
*/ 
Qt_widget_Nef_3(const Nef_polyhedron_S2& N); 

/// @}

}; /* end Qt_widget_Nef_S2 */
} /* end namespace CGAL */
