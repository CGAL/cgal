/// \ingroup PkgPolygonMeshProcessingConcepts
/// \cgalConcept
///
/// The concept `PMPSizingField` defines the requirements for the sizing field
/// used in `CGAL::Polygon_mesh_processing::isotropic_remeshing()` to define
/// the target length for every individual edge during the remeshing process.
///
/// \cgalHasModelsBegin
/// \cgalHasModels{CGAL::Polygon_mesh_processing::Uniform_sizing_field}
/// \cgalHasModels{CGAL::Polygon_mesh_processing::Adaptive_sizing_field}
/// \cgalHasModelsEnd
///


class PMPSizingField{
public:

/// @name Types
///  These types are used for the documentation of the functions of the concept and not needed implementation wise.
/// @{

/// Vertex descriptor type
typedef boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

/// Halfedge descriptor type
typedef boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

/// 3D point type matching the value type of the vertex property map passed to `CGAL::Polygon_mesh_processing::isotropic_remeshing()`
typedef unspecified_type Point_3;

/// Polygon mesh type matching the type passed to `CGAL::Polygon_mesh_processing::isotropic_remeshing()`
typedef unspecified_type PolygonMesh;

/// Number type matching the `FT` type of the geometric traits passed to `CGAL::Polygon_mesh_processing::isotropic_remeshing()`
typedef unspecified_type FT;

/// @}

/// @name Functions
/// @{

/// returns the sizing value at `v` (used during tangential relaxation).
FT at(const vertex_descriptor v, const PolygonMesh& pmesh) const;

/// returns the ratio of the current edge squared length and the local target edge squared length between
/// the points of `va` and `vb` in case the current edge is too long, and `std::nullopt` otherwise
/// (used for triggering edge splits and preventing some edge collapses).
std::optional<FT> is_too_long(const vertex_descriptor va,
                              const vertex_descriptor vb,
                              const PolygonMesh& pmesh) const;

/// returns the ratio of the squared length of `h` and the
/// local target edge squared length if it is too short, and `std::nullopt` otherwise
/// (used for triggering edge collapses).
std::optional<FT> is_too_short(const halfedge_descriptor h,
                               const PolygonMesh& pmesh) const;

/// returns the position of the new vertex created when splitting the edge of `h`.
Point_3 split_placement(const halfedge_descriptor h,
                        const PolygonMesh& pmesh) const;

/// function called after the addition of the split vertex `v` in `pmesh`.
/// This function can be used for example to update a pre-computed sizing field.
void register_split_vertex(const vertex_descriptor v,
                           const PolygonMesh& pmesh);

/// @}
};


