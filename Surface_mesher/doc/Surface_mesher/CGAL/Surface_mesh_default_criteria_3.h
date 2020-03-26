
namespace CGAL {

/*!
\ingroup PkgSurfaceMesher3Classes

The class `Surface_mesh_default_criteria_3` implements the most commonly used combination
of meshing criteria. It involves mainly three criteria which are
in order:
<UL>
<LI>a lower bound on the minimum angle in degrees of the surface mesh facets.

<LI>an upper bound on the radius of surface Delaunay balls.
A surface Delaunay ball is a ball circumscribing a facet,
centered on the surface and empty of vertices.
Such a ball exists for each facet
of the current surface mesh.
Indeed the current surface mesh
is the Delaunay triangulation of the current sampling restricted to
the surface
which is just the set of facets in the three dimensional Delaunay triangulation of
the sampling that have a Delaunay surface ball.

<LI>an upper bound on the center-center distances of the surface mesh facets.
The center-center distance of a surface mesh facet
is the distance between the facet circumcenter and the
center of its surface Delaunay ball.
</UL>

\cgalModels `SurfaceMeshFacetsCriteria_3`

\sa `make_surface_mesh`

*/
template< typename Tr >
class Surface_mesh_default_criteria_3 {
public:

/// \name Types
/// @{

/*!
The numerical type.
*/
typedef Tr::FT FT;

/// @}

/// \name Creation
/// @{

/*!
Returns a `Surface_mesh_default_criteria_3` with `angle_bound`, `radius_bound`,
`distance_bound` as bounds for the minimum facet angle in degrees,
the radius of the surface Delaunay balls
and the center-center distances respectively.
*/
Surface_mesh_default_criteria_3(FT angle_bound,
FT radius_bound,
FT distance_bound);

/// @}

}; /* end Surface_mesh_default_criteria_3 */
} /* end namespace CGAL */
