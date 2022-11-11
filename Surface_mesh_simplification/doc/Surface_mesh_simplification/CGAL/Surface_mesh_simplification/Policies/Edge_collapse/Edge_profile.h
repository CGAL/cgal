
namespace CGAL {
namespace Surface_mesh_simplification {

/*!
\ingroup PkgSurfaceMeshSimplificationRef

The class `Edge_profile` regroups useful information about an edge, such as its incident vertices
and faces. Its purpose is to pass this information onto the chosen models of `GetCost` and `GetPlacement`,
and as such an object of this type appears as a parameter in the `operator()` of both concepts.

The template parameters of this class must be consistent with the types used in the main call
to `edge_collapse()`: if you have specified a vertex point map or a geometric traits
(via `CGAL::parameters::vertex_point_map()` and `CGAL::parameters::geom_traits()` respectively),
then the template parameters must be identical.

Note however that if you wish to define your own models of the `GetCost` or `GetPlacement` concepts, you can
simply template the call to `operator()` by an "EdgeProfile" template parameter, which will
avoid having to explicit the template parameters of the `Edge_profile` class: indeed, the profile appears
as an argument of the `operator()` functions of these two concepts and automatic template deduction
can therefore be used.

\tparam TriangleMesh is the type of surface mesh being simplified
\tparam VertexPointMap is the type of the map that associates geometric positions to vertices of the mesh.
                       If you have specified a vertex point map in the named parameters in the call to `edge_collapse()`,
                       it must the same map type. Otherwise, you need not specify this template
                       parameter and it will be automatically deduced.
\tparam GeomTraits is the type of the kernel that is used to perform predicates and constructions
                   on geometric objects. If you have specified a traits class in the named parameters
                   in the call to `edge_collapse()`, it must be the same traits class type. Otherwise, you need not
                   specify this template parameter and it will be automatically deduced.

\sa `GetCost`
\sa `GetPlacement`

*/
template< typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
class Edge_profile
{
public:
  /// \name Types
  /// @{

  /*!
  The type of the surface mesh to simplify.
  */
  typedef TriangleMesh Triangle_mesh;

  /*!
  The type of a property map that maps vertices on points.
  */
  typedef VertexPointMap Vertex_point_map;

  /*!
  The type of a kernel-like objects used for predicates and constructions.
  */
  typedef GeomTraits Geom_traits;

  /*!
  A \bgl vertex descriptor representing a vertex of the surface mesh.
  */
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;

  /*!
  A \bgl halfedge descriptor representing a halfedge of the surface mesh.
  */
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  /*!
  A \bgl edge descriptor representing an edge of the surface mesh.
  */
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor edge_descriptor;

  /*!
  The unsigned integer type used for representing the number of edges in the graph.
  */
  typedef typename boost::graph_traits<TriangleMesh>::edges_size_type edges_size_type;

  /*!
  The point type for the surface mesh vertex.
  */
  typedef typename boost::property_traits<VertexPointMap>::value_type Point;

  /*!
  The coordinate type of points.
  */
  typedef typename GeomTraits::FT FT;
  /// @}

  /// \name Access Functions
  /// @{

  /*!
  One of vertices of the edge to be collapsed.
  */
  vertex_descriptor v0() const;

  /*!
  The other vertex of the edge to be collapsed.
  */
  vertex_descriptor v1() const;

  /*!
  One of the directed edges corresponding to the halfedge being collapsed.
  */
  halfedge_descriptor v0_v1() const;

  /*!
  The other directed edge corresponding to the halfedge being collapsed.
  */
  halfedge_descriptor v1_v0() const;

  /*!
  The point of vertex `v0`.
  */
  const Point& p0() const;

  /*!
  The point of vertex `v1`.
  */
  const Point& p1() const;

  /*!
  If ` v0v1` belongs to a finite face (is not a border edge)
  the third vertex of that triangular face, a <I>null descriptor</I> otherwise.
  */
  vertex_descriptor vL() const;

  /*!
  If ` v0v1` belongs to a finite face (is not a border edge)
  the directed edge from ` v1` to ` vL`, a <I>null descriptor</I> otherwise.
  */
  halfedge_descriptor v1_vL() const;

  /*!
  If ` v0v1` belongs to a finite face (is not a border edge)
  the directed edge from ` vL` to ` v0`, a <I>null descriptor</I> otherwise.
  */
  halfedge_descriptor vL_v0() const;

  /*!
  If ` v1v0` belongs to a finite face (is not a border edge)
  the third vertex of that triangular face, a <I>null descriptor</I> otherwise.
  */
  vertex_descriptor vR() const;

  /*!
  If ` v1v0` belongs to a finite face (is not a border edge)
  the directed edge from ` v0` to ` vR`, a <I>null descriptor</I> otherwise.
  */
  halfedge_descriptor v0_vR() const;

  /*!
  If ` v1v0` belongs to a finite face (is not a border edge)
  the directed edge from ` vR` to ` v1`, a <I>null descriptor</I> otherwise.
  */
  halfedge_descriptor vR_v1() const;

  /*!
  The unique sequence of the vertices around ` v0v1` in topological order (<I>ccw</I> or <I>cw</I> depending
  on the relative ordering of `v0` and `v1` in the profile).
  The sequence may have duplicates, but when this happens the edge is not collapsible.
  */
  std::vector<vertex_descriptor> link() const;

  /*!
  The unique collection of the border directed edges incident upon ` v0` and ` v1`.
  */
  std::vector<halfedge_descriptor> border_edges() const;

  /*!
  Indicates if `v0v1` belongs to a finite face of the surface mesh (i.e, `v0v1` is not a border edge).
  */
  bool left_face_exists() const;

  /*!
  Indicates if `v0v1` belongs to a finite face of the surface mesh (i.e, `v1v0` is not a border edge).
  */
  bool right_face_exists() const;

  /*!
  Returns the surface mesh the edge belongs to.
  */
  const Triangle_mesh& surface_mesh() const;

  /*!
  Returns the vertex point property map.
  */
  const Vertex_point_map& vertex_point_map() const;

  /*!
  Returns the geometric traits class.
  */
  const Geom_traits& geom_traits() const;

  /// @}
};

} // namespace Surface_mesh_simplification
} // namespace CGAL
