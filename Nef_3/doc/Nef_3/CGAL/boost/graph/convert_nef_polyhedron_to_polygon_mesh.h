namespace CGAL {

/// \ingroup PkgNef3IOFunctions
/// Converts an objet of type `Nef_polyhedron_3` into a polygon mesh model of `MutableFaceGraph`.
/// Note that contrary to `Nef_polyhedron_3::convert_to_polyhedron()`, the output is not triangulated
/// (but faces with more than one connected component of the boundary).
/// The polygon mesh can be triangulated by setting `triangulate_all_faces` to `true` or by calling the function `triangulate_faces()`.
/// \pre `Polygon_mesh` must have an internal point property map with value type being `Nef_polyhedron_3::Point_3`.
/// \pre `nef.simple()`
  
  template <class Nef_polyhedron, class Polygon_mesh>
  void convert_nef_polyhedron_to_polygon_mesh(const Nef_polyhedron& nef, Polygon_mesh& pm, bool triangulate_all_faces = false);
}
