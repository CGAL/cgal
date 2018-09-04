
namespace CGAL {

/*!
\ingroup PkgSkinSurface3Ref

An items class for the `Polyhedron_3` that stores information in faces instead of an auxiliary data structure. Using it accelerates the subdivision function `subdivide_skin_surface_mesh_3()`.
*/

template<class SkinSurface3>
struct Skin_surface_polyhedral_items_3: public Polyhedron_items_3
{};

} // namespace CGAL



