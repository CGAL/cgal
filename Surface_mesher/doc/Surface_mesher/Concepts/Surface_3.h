
/*!
\ingroup PkgSurfaceMesher3Concepts
\cgalconcept

The concept `Surface_3` describes the types of surfaces to be meshed. 
The surface types 
are required to be copy constructible 
and assignable. 

\hasModel Implicit_surface_3<Traits, Function> 

\sa `make_surface_mesh`
\sa `SurfaceMeshTraits_3` 
\sa `Surface_mesh_traits_generator_3<Surface>` 

*/

class Surface_3 {
public:

/// \name Types 
/// In addition, surface types are required 
/// - either to provide a nested type: `Surface_mesher_traits_3`
/// - or to be a surface type for which a specialization of the traits generator
///   `Surface_mesh_traits_generator_3<Surface>` exists. 
/// @{

/*! 
a model of `SurfaceMesherTraits_3` 
*/ 
typedef Hidden_type Surface_mesher_traits_3; 

/// @}

}; /* end Surface_3 */

