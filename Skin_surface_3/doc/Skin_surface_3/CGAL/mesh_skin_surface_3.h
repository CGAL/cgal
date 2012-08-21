namespace CGAL {

/*!
\ingroup PkgSkinSurface3

Constructs a mesh of the `skin_surface` in `p`.

The function `mesh_skin_surface_3` constructs a mesh isotopic to the skin 
surface based on the algorithm in \cite cgal:kv-mssct-05. It takes 
as input a `SkinSurface_3` object, which is a model of the 
`SkinSurface_3` concept and outputs the mesh in a 
`Polyhedron_3` object. 

\pre `SkinSurface_3` is
a model of the concept `SkinSurface_3` and `Polyhedron_3::HDS` can be
used as the template argument of the
`Polyhedron_incremental_builder_3<HDS>`.

*/
template <typename SkinSurface_3, typename Polyhedron_3>
void mesh_skin_surface_3
(const SkinSurface_3 &skin_surface, Polyhedron_3 &p);

} /* namespace CGAL */

