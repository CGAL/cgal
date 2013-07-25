
/*!
\ingroup PkgSkinSurface3Concepts
\cgalConcept

The concept `SkinSurface_3` defines a skin surface and provides an 
interface for the dedicated skin surface mesher. 
The concept requires a constructor from an iterator range of 
weighted points and a shrink factor. By default the input balls are 
grown in such that the skin surface wraps around the input balls. 

\cgalHasModel `CGAL::Skin_surface_3<SkinSurfaceTraits_3>` 
\cgalHasModel `CGAL::Union_of_balls_3<SkinSurfaceTraits_3>`

*/

class SkinSurface_3 {
public:

/// \name Types 
/// @{

/*!
The geometric traits used for the 
construction of the regular triangulation. 
*/ 
typedef unspecified_type Geometric_traits; 

/*!
The `Weighted_point` type defined 
by the `Geometric_traits`. 
*/ 
typedef unspecified_type Weighted_point; 

/*!
The `Bare_point` type defined by the 
`Geometric_traits`. 
*/ 
typedef unspecified_type Bare_point; 

/*!
The `FT` type defined by the 
`Geometric_traits`. This is the number type used by the 
weighted points. 
*/ 
typedef unspecified_type FT; 

/// @} 

/// \name Creation 
/// @{

/*!
The mandatory arguments to the constructor are 
an iterator range `[begin,end)` of weighted points and a 
shrink factor between 0 and 1. 
*/ 
 template < class WP_iterator >
Skin_surface_3 ( WP_iterator begin, WP_iterator end, 
RT shrink_factor); 

/// @} 

/// \name Operations 
/// @{

/*!
Constructs a coarse mesh in `p`.
\tparam Polyhedron must be an instance of `Polyhedron_3`.
\pre `Polyhedron::HDS` can be used as the template argument of the `CGAL::Polyhedron_incremental_builder_3<HDS>`. 
*/ 
template <class Polyhedron> void 
mesh_skin_surface_3 (Polyhedron &p); 

/*!
Subdivides the skin surface using `nSubdiv` 1-4 split 
operations (each triangle is split into four sub-triangles) and 
the new vertices are moved towards the skin surface. 

\tparam Polyhedron must be an instance of `Polyhedron_3`.
*/ 
template <class Polyhedron> void 
subdivide_skin_surface_mesh_3 (Polyhedron &p, int nSubdiv = 
1); 

/// @}

}; /* end SkinSurface_3 */

