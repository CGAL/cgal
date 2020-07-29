
namespace CGAL {

/*!
\ingroup PkgSkinSurface3Ref

constructs a mesh of the `skin_surface` in `p`.

The function `mesh_skin_surface_3()` constructs a mesh isotopic to the skin
surface based on the algorithm in \cgalCite{cgal:kv-mssct-05}. It takes
as input a `SkinSurface_3` object, which is a model of the
`SkinSurface_3` concept and outputs the mesh in a
`Polyhedron_3` object.

\tparam SkinSurface_3  must be a model of the concept `SkinSurface_3`.
\tparam Polyhedron must be an instance of `Polyhedron_3`.

\pre `Polyhedron::HDS` can be used as the template argument
of the `Polyhedron_incremental_builder_3<HDS>`.

*/
template <typename SkinSurface_3, typename Polyhedron>
void mesh_skin_surface_3
(const SkinSurface_3 &skin_surface, Polyhedron &p);

} /* namespace CGAL */

