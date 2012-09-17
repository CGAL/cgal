
namespace CGAL {

/*!
\ingroup PkgNef3

The class `Qt_widget_Nef_3` uses the OpenGL interface of Qt to display a 
`Nef_polyhedron_3`. Its purpose is to provide an easy to use viewer for 
`Nef_polyhedron_3`. There are no means provided to enhance the 
functionality of the viewer. 

In addition to the functions inherited from the Qt class `OGLWidget` via, 
`Qt_widget_Nef_3` only has a single public 
constructor. For the usage of `Qt_widget_Nef_3` see the example 
below. 

### Parameters ###

The template parameter expects an instantiation of `Nef_polyhedron_3<Traits>`. 

\sa `CGAL::Nef_polyhedron_3<Traits>` 

### Example ###

This example reads a 3D Nef polyhedron from standard input and displays it 
in a Qt widget. 

\cgalexample{visualization_SNC.cpp} 

*/
template< typename Nef_polyhedron_3 >
class Qt_widget_Nef_3 {
public:

/// \name Creation 
/// @{

/*! 
Creates a widget `W` for displaying the 3D Nef polyhedron N. 
*/ 
Qt_widget_Nef_3(const Nef_polyhedron_3& N); 

/// @}

}; /* end Qt_widget_Nef_3 */
} /* end namespace CGAL */
