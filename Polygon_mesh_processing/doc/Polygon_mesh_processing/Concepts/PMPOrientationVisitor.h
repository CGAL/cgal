/// \ingroup PkgPolygonMeshProcessingConcepts
/// \cgalConcept
///
/// The concept `PMPOrientationVisitor` defines the requirements for the visitor
/// used in `CGAL::Polygon_mesh_processing::orient_polygon_soup()` to track
/// the detection of non-manifold simplices and the modifications performed to polygons
/// during the orientation process.
///
/// \cgalRefines `CopyConstructible`
/// \cgalHasModel `CGAL::Polygon_mesh_processing::Default_orientation_visitor`.

class PMPOrientationVisitor{
public:

/// @name Functions used to report non-manifold simplices.
/// @{

/// called each time an edge appears in more than two polygons.
/// `id1` and `id2` are the ids of the endpoints of the edge in the input point range.
void non_manifold_edge(const std::size_t & id1, const std::size_t& id2);

/// called each time a non-manifold vertex is detected.
/// `id` is the id of the non-manifold vertex.
void non_manifold_vertex(const std::size_t & id);

/// @}

/// @name Functions used to report modifications done to the polygons.
/// @{

/// called when the orientation of a polygon is not compatible with its neighbors
/// and the polygon is reversed. `id` is the index of the polygon in the input polygon range.
void polygon_orientation_reversed(const std::size_t& id) ;

/// called for each non-manifold vertex (whether part of a non-manifold edge or not).
/// Non-manifoldness is resolved in the algorithm by duplicating the point.
/// `input_id` is the index of the input point, and `new_id` is the index of the new point.
/// Note that a point might be duplicated several times.
void duplicated_vertex(const std::size_t& input_id, const std::size_t& new_id);

/// called when the point with id `input_id` in polygon with id `pid` is replaced
/// by the point with id `new_id`. This function is called after all calls to `duplicated_vertex()`.
void point_id_in_polygon_updated(const std::size_t& pid, const std::size_t& input_id, const std::size_t& new_id);

/// @}
};


