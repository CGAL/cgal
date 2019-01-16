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
/// @param triangulate_all_faces indicates whether all the faces must be triangulated.
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
  /// @tparam PointRange a model of the concepts `RandomAccessContainer` and 
  /// `BackInsertionSequence` whose `value_type` is the point type
  /// @tparam PolygonRange a model of the concept `RandomAccessContainer` whose 
  /// `value_type` is a model of the concept `RandomAccessContainer` whose `value_type` is `std::size_t`.
  /// 
  /// @param nef the input.
  /// @param pm the output.
  /// @param triangulate_all_faces indicates whether all polygons must be triangulated.
  template <class Output_kernel, class Nef_polyhedron, typename PolygonRange, typename PointRange>
  void convert_nef_polyhedron_to_polygon_soup(const Nef_polyhedron& nef,
                                                    PointRange& points,
                                                    PolygonRange& polygons,
                                                    bool triangulate_all_faces = false);
}
