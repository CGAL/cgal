/// \ingroup PkgPolygonMeshProcessingConcepts
/// \cgalConcept
///
/// The concept `PMPOrientationVisitor` defines the requirements for the visitor
/// used in `CGAL::Polygon_mesh_processing::orient_polygon_soup()` to track
/// the non-manifold simplices and the modifications done to polygons during the
/// orientation process.
///
/// \cgalRefines `CopyConstructible`
/// \cgalHasModel `CGAL::Polygon_mesh_processing::Default_orientation_visitor`.

class PMPOrientationVisitor{
public:

/// @name Functions used to report non-manifold simplices
/// @{
  /// called  after the filling of the edge map. If an edge appears in more than 2 polygons,
  /// it is marked as non-manifold.
  /// `id1` and `id2` are the two ids of the edge endpoints, in the input points range.
  void non_manifold_edge(const std::size_t & id1, const std::size_t& id2);
/// @}

/// @name Functions used to report modifications done to the polygons
/// @{
  ///called after the filling of the edge map. When a polygon's orientation is not
  ///compatible with its neighbors, and there is no non-manifold edge involved, this
  /// polygon's orientation is reversed.
  /// `id` is the index of the reversed polygon in the input polygons range.
  void polygon_orientation_reversed(const std::size_t& id) ;
  /// called after the orientation to fix the non-manifold vertices.
  /// If a vertex is an endpoint of a non-manifold edge, it is considered non-manifold.
  /// If the link of a vertex is neither a cycle nor a chain, but several cycles and chains,
  /// then this vertex is considered non-manifold.
  /// `v1` is the index of the input point, and `v2` is the index of the new vertex.
  void duplicated_vertex(const std::size_t& v1, const std::size_t& v2);
  ///called after a vertex has been duplicated. The new vertex id will replace the old one
  /// in all the polygons that need it.
  /// `polygon_id` is the id of the impacted polygon, `old_vertex` the index of the vertex
  ///  that will be replaced in the polygon, and `new_vertex` the index of the new vertex.
  void point_id_in_polygon_updated(const std::size_t& polygon_id, const std::size_t& old_vertex, const std::size_t& new_vertex);
/// @}

};



