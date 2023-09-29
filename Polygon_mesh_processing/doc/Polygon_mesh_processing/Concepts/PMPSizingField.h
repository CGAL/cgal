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
typedef unspecified_type halfedge_descriptor;

/// 3D point type
typedef unspecified_type Point_3;

/// Polygon mesh type
typedef unspecified_type PolygonMesh;

/// Numerical type
typedef unspecified_type FT;

/// @}

/// @name Functions
/// @{

/// called to check whether the halfedge `h` is longer than the target edge size
/// and as such should be split. If the halfedge is longer, it returns the squared
/// length of the edge.
std::optional<FT> is_too_long(const halfedge_descriptor h,
                                const PolygonMesh& pmesh) const;

/// called to check whether the halfedge with end vertices `va` and `vb` is longer
/// than the target edge size and as such should be split. If the halfedge is longer,
/// it returns the squared length of the edge.
std::optional<FT> is_too_long(const vertex_descriptor va,
                                const vertex_descriptor vb) const;

/// called to check whether the halfedge `h` should be collapsed in case it is
/// shorter than the target edge size.
std::optional<FT> is_too_short(const halfedge_descriptor h,
                                 const PolygonMesh& pmesh) const;

/// called to define the location of the halfedge `h` split in case `is_too_long()`
/// returns a value.
Point_3 split_placement(const halfedge_descriptor h,
                        const PolygonMesh& pmesh) const;
/// @}
};


