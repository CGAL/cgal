namespace CGAL {

/*!
\ingroup PkgSkinSurface3

Subdivides the skin surface using `nSubdiv` 1-4 split
operations (each triangle is split into four sub-triangles) and
the new vertices are moved towards the skin surface.

The function `subdivide_skin_surface_mesh_3` subdivides a skin surface mesh constructed 
by the function `mesh_skin_surface_3<SkinSurface_3, Polyhedron_3>`. 
*/
void
subdivide_skin_surface_mesh_3<SkinSurface_3,Polyhedron_3> (const
SkinSurface_3 &skin_surface, Polyhedron_3 &p, int nSubdiv =
1);

} /* namespace CGAL */

