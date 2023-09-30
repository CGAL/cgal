/// \ingroup PkgPolygonMeshProcessingConcepts
/// \cgalConcept
///
/// The concept `PMPSizingField` defines the requirements for the sizing field
/// used in `CGAL::Polygon_mesh_processing::isotropic_remeshing()` to define
/// the target length for every individual edge during the remeshing process.

class PMPSizingField{
public:

/// @name Types
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

/// a function controlling edge split and edge collapse,
/// returning the ratio of the current edge length and the local target edge length between
/// the points of `va` and `vb` in case the current edge is too long, and `std::nullopt` otherwise.
std::optional<FT> is_too_long(const vertex_descriptor va,
                                const vertex_descriptor vb) const;

/// a function controlling edge collapse by returning the ratio of the squared length of `h` and the
/// local target edge length if it is too short, and `std::nullopt` otherwise.
std::optional<FT> is_too_short(const halfedge_descriptor h,
                                 const PolygonMesh& pmesh) const;

/// a function returning the location of the split point of the edge of `h`.
Point_3 split_placement(const halfedge_descriptor h,
                        const PolygonMesh& pmesh) const;
/// @}
};


