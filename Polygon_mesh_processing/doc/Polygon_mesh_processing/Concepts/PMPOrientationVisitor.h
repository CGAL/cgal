/// \ingroup PkgPolygonMeshProcessingConcepts
/// \cgalConcept
///
/// The concept `PMPOrientationVisitor` defines the requirements for the visitor
/// used in `CGAL::Polygon_mesh_processing::orient_polygon_soup()` to be notified
/// of the presence of non-manifold simplices and of the modifications done to polygons
/// during the orientation process.
///
/// \cgalRefines `CopyConstructible`
/// \cgalHasModel `CGAL::Polygon_mesh_processing::Default_orientation_visitor`.

class PMPOrientationVisitor{
public:

/// @name Functions used to report non-manifold simplices
/// @{
  /// is called each time an edge appears in more than 2 polygons.
  /// `id1` and `id2` are the ids, from the input point range, of the endpoints of the edge.
  void non_manifold_edge(const std::size_t & id1, const std::size_t& id2);
  ///is called each time a vertex is found non-manifold.
  /// `id`is the id of the non manifold vertex.
  void non_manifold_vertex(const std::size_t & id);
/// @}

/// @name Functions used to report modifications done to the polygons
/// @{
  /// is called when the orientation of a polygon is not compatible
  /// with its neighbors and is reversed. `id` is the index of the polygon in the input polygon range.
  void polygon_orientation_reversed(const std::size_t& id) ;
  /// is called for each non-manifold vertex (part of a non-manifold edge or not).
  /// Non-manifoldness is resolved in the algorithm by duplicating the point.
  /// `input_id` is the index of the input point, and `new_id` is the index of the new point.
  /// Note that a point might be duplicated several times.
  void duplicated_vertex(const std::size_t& input_id, const std::size_t& new_id);
  /// is called when a the point with id `input_id` in polygon with id `pid` is replaced
  /// by the point with id `new_id`. This functions is called after all calls to `duplicated_vertex()` are done.
  void point_id_in_polygon_updated(const std::size_t& pid, const std::size_t& input_id, const std::size_t& new_id);
/// @}

};



