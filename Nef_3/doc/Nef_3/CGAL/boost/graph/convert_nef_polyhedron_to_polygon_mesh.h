namespace CGAL {

/// \ingroup PkgNef3IOFunctions
/// Converts an objet of type `Nef_polyhedron_3` into a polygon mesh model of `MutableFaceGraph`.
/// Note that contrary to `Nef_polyhedron_3::convert_to_polyhedron()`, the output is not triangulated
/// (but faces with more than one connected component of the boundary).
/// The polygon mesh can be triangulated by setting `triangulate_all_faces` to `true` or by calling the function `triangulate_faces()`.
/// @tparam Nef_polyhedron an object of type `Nef_polyhedron_3`.
/// @tparam Polygon_mesh a model of `MutableFaceGraph`.
/// 
/// @param nef the input.
/// @param pm the output.
/// @param triangulate_all_faces the bool for triangulating the faces.
/// 
/// \pre `Polygon_mesh` must have an internal point property map with value type being `Nef_polyhedron_3::Point_3`.
/// \pre `nef.simple()`
  
  template <class Nef_polyhedron, class Polygon_mesh>
  void convert_nef_polyhedron_to_polygon_mesh(const Nef_polyhedron& nef, Polygon_mesh& pm, bool triangulate_all_faces = false);
  
  /// \ingroup PkgNef3IOFunctions
  /// Converts an objet of type `Nef_polyhedron_3` into a polygon soup.
  /// The polygons can be triangulated by setting `triangulate_all_faces` to `true`.
  /// @tparam Output_kernel the Kernel from which the point type of points come from.
  /// @tparam Nef_polyhedron an object of type `Nef_polyhedron_3`.
  /// @tparam Polygon_mesh a model of `MutableFaceGraph`.
  /// @tparam PolygonRange a container of collections of indices of the points in their range.
  /// @tparam PointRange a container of `Output_kernel::Point_3`.
  /// 
  /// @param nef the input.
  /// @param pm the output.
  /// @param triangulate_all_faces the bool for triangulating the faces.
  /// \pre `nef.simple()`
  template <class Output_kernel, class Nef_polyhedron, typename PolygonRange, typename PointRange>
  void convert_nef_polyhedron_to_polygon_soup(const Nef_polyhedron& nef,
                                                    PointRange& points,
                                                    PolygonRange& polygons,
                                                    bool triangulate_all_faces = false);
}
