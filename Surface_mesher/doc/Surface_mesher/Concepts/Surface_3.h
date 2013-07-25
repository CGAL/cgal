
/*!
\ingroup PkgSurfaceMesher3Concepts
\cgalConcept

The concept `Surface_3` describes the types of surfaces to be meshed. 
The surface types 
are required to be copy constructible 
and assignable. 

\cgalHasModel `CGAL::Implicit_surface_3<Traits, Function>` 

\sa `CGAL::make_surface_mesh()`
\sa `SurfaceMeshTraits_3` 
\sa `CGAL::Surface_mesh_traits_generator_3<Surface>` 

*/

class Surface_3 {
public:

/// \name Types 
/// In addition, `Surface_3` is required 
/// - either to provide a nested type: `Surface_mesher_traits_3`
/// - or to be a surface type for which a specialization of the traits generator
///   `CGAL::Surface_mesh_traits_generator_3<Surface_3>` exists. 
/// @{

/*!
a model of `SurfaceMesherTraits_3` 
*/ 
typedef unspecified_type Surface_mesher_traits_3; 

/// @}

}; /* end Surface_3 */

