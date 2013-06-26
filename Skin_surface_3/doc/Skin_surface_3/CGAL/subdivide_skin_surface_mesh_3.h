namespace CGAL {

/*!
\ingroup PkgSkinSurface3


subdivides a skin surface mesh constructed 
by the function `mesh_skin_surface_3()`
using `nSubdiv` 1-4 split
operations (each triangle is split into four sub-triangles) and
the new vertices are moved towards the skin surface.


\tparam SkinSurface_3  must be a model of the concept `SkinSurface_3`.
\tparam Polyhedron must be an instance of `Polyhedron_3`.
 
*/
template <class SkinSurface_3, class Polyhedron>
void subdivide_skin_surface_mesh_3 (const
SkinSurface_3 &skin_surface, Polyhedron &p, int nSubdiv = 1);

} /* namespace CGAL */

