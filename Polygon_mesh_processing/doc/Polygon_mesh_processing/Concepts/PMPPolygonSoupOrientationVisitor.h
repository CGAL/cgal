/// \ingroup PkgPolygonMeshProcessingConcepts
/// \cgalConcept
///
/// The concept `PMPPolygonSoupOrientationVisitor` defines the requirements for the visitor
/// used in `CGAL::Polygon_mesh_processing::orient_polygon_soup()` to track
/// the detection of non-manifold simplices and the modifications performed to polygons
/// during the orientation process.
///
/// \cgalRefines{CopyConstructible}
/// \cgalHasModelsBegin
/// \cgalHasModels{CGAL::Polygon_mesh_processing::Default_orientation_visitor}
/// \cgalHasModelsEnd

class PMPPolygonSoupOrientationVisitor{
public:

/// @name Functions used to report non-manifold simplices.
/// @{

/// called each time an edge appears in more than two polygons.
/// `id1` and `id2` are the vertex ids of the endpoints of the edge.
/// `nb_polygons` indicates the number of polygons containing that edge.
void non_manifold_edge(std::size_t id1, std::size_t id2, std::size_t nb_polygons);

/// called each time a non-manifold vertex is detected.
/// `vid` is the id of the vertex that is non-manifold.
/// `nb_link_ccs` is the number of edge connected components of vertices
/// in the link of the corresponding vertex.
void non_manifold_vertex(std::size_t vid, std::size_t nb_link_ccs);

/// called during the detection of a non-manifold vertex, one time
/// per incident edge connected component of vertices in the link of the vertex with id `vid`.
/// `id` is the id of the vertex that is non-manifold.
/// `polygon_ids` contains the ids of the polygon in such a connected component.
/// This function is called a number of times exactly equal to the parameter `nb_link_ccs`
/// for the same vertex in `non_manifold_vertex()`. Note that the aforementioned function is
/// called after all the calls to this function are done.
void link_connected_polygons(std::size_t vid, const std::vector<std::size_t>& polygon_ids){}

/// @}

/// @name Functions used to report modifications done to the polygons.
/// @{

/// called when the orientation of a polygon is not compatible with its neighbors
/// and the polygon is reversed. `id` is the index of the polygon in the input polygon range.
void polygon_orientation_reversed(std::size_t id);

/// called for each non-manifold vertex (whether part of a non-manifold edge or not).
/// Non-manifoldness is resolved in the algorithm by duplicating the vertex.
/// `input_id` is the index of the input vertex, and `new_id` is the index of the new vertex.
/// Note that a vertex might be duplicated several times.
void duplicated_vertex(std::size_t input_id, std::size_t new_id);

/// called when the vertex with id `input_id` in polygon with id `pid` is replaced
/// by the vertex with id `new_id`. This function is called after all calls to `duplicated_vertex()`.
void vertex_id_in_polygon_replaced(std::size_t pid, std::size_t input_id, std::size_t new_id);

/// @}
};


