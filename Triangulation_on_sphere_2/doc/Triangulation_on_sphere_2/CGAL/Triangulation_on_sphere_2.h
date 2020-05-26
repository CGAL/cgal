namespace CGAL {

/*!
\ingroup PkgTriangulationOnSphere2TriangulationClasses

The class `Triangulation_on_sphere_2` is the basic class designed to handle triangulations
of set of points \f$ \mathcal{S}\f$ on the sphere : it has vertices at the points of \f$ \mathcal{S}\f$
and its domain covers the convex hull of \f$ \mathcal{S}\f$.

\warning This triangulation supports neither the insertion nor the removal of vertices,
see `CGAL::Delaunay_triangulation_on_sphere_2` for such purposes.

This triangulation class is very similar to `CGAL::Triangulation_2` as both classes represent
triangulations of 2-manifold domain without boundary. A significant difference is that
in the case of Euclidean 2D triangulation, it is necessary to introduce so-called <i>infinite
faces</i> to complete the convex hull into an actual 2-manifold without boundary that the triangulation
data structure can represent. This is not necessary for triangulations of the sphere,
that is already perfectly adapted to the triangulation data structure.

There is an exception to the previous statement: in the degenerate configuration
where all points of \f$ \mathcal{S}\f$ lie on the same hemisphere, the triangulation has a border.
Internally, the triangulation data structure must however remain a 2-manifold at all time,
and to ensure this fictitious faces called <i>ghost faces</i> are added. In contrast, faces that
not ghost-faces are called <i>solid</i> faces.

\tparam Traits is the geometric traits, which must be a model of the concept `TriangulationOnSphereTraits_2`.

\tparam TDS is the triangulation data structure, which must be a model of the concept `TriangulationDataStructure_2`.
        By default, the triangulation data structure is instantiated by
        `Triangulation_data_structure_2 < Triangulation_on_sphere_vertex_base_2<Gt>, Triangulation_on_sphere_face_base_2<Gt> >`.

\sa `CGAL::Delaunay_triangulation_on_sphere_2<Traits, TDS>`
*/
template< typename Traits, typename TDS >
class Triangulation_on_sphere_2
  : public Triangulation_cw_ccw_2
{
public:

  /// \name Types
  /// @{

  /*!
  The traits class.
  */
  typedef Traits Geom_traits;

  /*!
  The triangulation data structure type.
  */
  typedef TDS Triangulation_data_structure;

  /*!
  The number type.
  */
  typedef Traits::FT FT;

  /*!
  The point type representing a point on the sphere.
  */
  typedef Traits::Point_on_sphere_2 Point;

  /*!
  The 3D point type.
  */
  typedef Traits::Point_3 Point_3;

public:
  /*!
  The vertex type.
  */
  typedef TDS::Vertex Vertex;

  /*!
  The edge type.
  */
  typedef TDS::Edge Edge;

  /*!
  The face type.
  */
  typedef TDS::Face Face;

public:
  /// \name Handles, Iterators, and Circulators
  ///
  /// The vertices and faces of the triangulation are accessed through
  /// handles, iterators and circulators. The handles are models
  /// of the concept `Handle` which basically offers the two dereference
  /// operators `.` and `->`. The handles are also model of the concepts
  ///  `LessThanComparable` and `Hashable`, that is they can be used as keys
  /// in containers such as `std::map` and `boost::unordered_map`.
  /// The iterators and circulators are all bidirectional and non-mutable.
  /// The circulators and iterators are convertible to handles with the same value type,
  /// so that whenever a handle appear in the parameter list of a function,
  /// an appropriate iterator or circulator can be passed as well.
  ///
  /// The edges of the triangulation can also be visited through
  /// iterators and circulators, the edge circulators and iterators are
  /// also bidirectional and non mutable.
  ///

  /*!
  Handle to a vertex.
  */
  typedef TDS::Vertex_handle Vertex_handle;

  /*!
  Handle to a face.
  */
  typedef TDS::Face_handle Face_handle;

  /*!
  Iterator over all vertices.
  */
  typedef TDS::Vertex_iterator All_vertices_iterator;

  /*!
  Iterator over all edges.
  */
  typedef TDS::Edge_iterator All_edges_iterator;

  /*!
  Iterator over all faces.
  */
  typedef TDS::Face_iterator All_faces_iterator;

  /*!
  Iterator over all solid edges.
  */
  typedef unspecified_type Solid_edges_iterator;

  /*!
  Iterator over all solid faces.
  */
  typedef unspecified_type Solid_faces_iterator;

  /*!
  Circulator over all vertices incident to a given vertex.
  */
  typedef unspecified_type Vertex_circulator;

  /*!
  Circulator over all edges incident to a given vertex.
  */
  typedef unspecified_type Edge_circulator;

  /*!
  Circulator over all faces incident to a given vertex.
  */
  typedef unspecified_type Face_circulator;

  /// @}

public:
  /// \name Creation
  /// @{

  /*!
  Introduces an empty triangulation.
  */
  Triangulation_on_sphere_2(const Traits& gt = Traits());

  /*!
  Introduces an empty triangulation and sets the center and radius to `c` and `r` respectively.
  */
  Triangulation_on_sphere_2(const Point_3& c, const FT r);

  /*!
  Copy constructor. All the vertices and faces are duplicated.
  After the copy, `*this` and `tr` refer to different triangulations: if `tr` is modified, `*this` is not.
  */
  Triangulation_on_sphere_2(const Triangulation_on_sphere_2& tr);

  /*!
  Assignment. All the vertices and faces are duplicated.
  After the assignment, `*this` and `tr` refer to different triangulations: if `tr` is modified, `*this` is not.
  */
  Triangulation_on_sphere_2 operator=(const Triangulation_on_sphere_2<Traits,TDS>& tr);

  /*!
  The triangulations `tr` and `*this` are swapped.
  */
  void swap(Triangulation_on_sphere_2& tr);

  /*!
  Deletes all faces and vertices, resulting in an empty triangulation.
  */
  void clear();

  /// @}

public:
  /// \name Access Functions
  /// @{

  /*!
  Returns the dimension of the convex hull.
  */
  int dimension() const;

  /*!
  Returns the number of vertices.
  */
  size_type number_of_vertices() const;

  /*!
  Returns the number of faces.
  */
  size_type number_of_faces() const;

  /*!
  Returns the number of ghost_faces.
  */
  size_type number_of_ghost_faces() const;

  /*!
  Returns a const reference to the triangulation traits object.
  */
  const Geom_traits& geom_traits() const;

  /*!
  Returns a const reference to the triangulation data structure.
  */
  const Triangulation_data_structure& tds() const;

  /*!
  Returns the geometric position of the vertex `*v`.
  */
  const Point& point(const Vertex_handle v);

  /*!
  Returns the geometric position of the `i`-th vertex of the face `*f`.
  */
  const Point& point(const Face_handle f, const int i);

  /// @}

  /// \name All Face, Edge and Vertex Iterators
  ///
  /// The following iterators allow respectively to visit all faces, edges and
  /// vertices of the triangulation. These iterators are non mutable,
  /// bidirectional and their value types are respectively `Face`,
  /// `Edge` and `Vertex`. They are all invalidated by any change in the triangulation.
  ///
  /// @{

  /*!
  Starts at an arbitrary vertex
  */
  All_vertices_iterator all_vertices_begin() const;

  /*!
  Past-the-end iterator
  */
  All_vertices_iterator all_vertices_end() const;

  /*!
  Starts at an arbitrary edge
  */
  All_edges_iterator all_edges_begin() const;

  /*!
  Past-the-end iterator
  */
  All_edges_iterator all_edges_end() const;

  /*!
  Starts at an arbitrary face
  */
  All_faces_iterator all_faces_begin() const;

  /*!
  Past-the-end iterator
  */
  All_faces_iterator all_faces_end() const;

  /*!
  Starts at an arbitrary face
  */
  Solid_faces_iterator solid_faces_begin() const;

  /*!
  Past-the-end iterator
  */
  Solid_faces_iterator solid_faces_end() const;

  /*!
  Starts at an arbitrary face
  */
  Solid_edges_iterator solid_edges_begin() const;

  /*!
  Past-the-end iterator
  */
  Solid_edges_iterator solid_edges_end() const;

  /// @}

  /// \name Face, Edge and Vertex Circulators
  ///
  /// The triangulation also provides circulators that allows to visit respectively
  /// all faces or edges incident to a given vertex or all vertices adjacent to a given vertex.
  /// These circulators are non-mutable and bidirectional. The `operator++` moves the
  /// circulator counterclockwise around the vertex while the `operator-` moves clockwise.
  /// A face circulator is invalidated by any modification of the face pointed to.
  /// An edge or a vertex circulator are invalidated by any modification of one of the two
  /// faces incident to the edge pointed to.
  ///
  /// @{

  /*!
  Starts at an arbitrary vertex incident to `v`.
  */
  Vertex_circulator incident_vertices(Vertex_handle v) const;

  /*!
  Starts at the first vertex of `f` adjacent to `v` in counterclockwise order around `v`.
  \pre Face `f` is incident to vertex `v`.
  */
  Vertex_circulator incident_vertices(Vertex_handle v, Face_handle f);

  /*!
  Starts at an arbitrary edge incident to `v`.
  */
  Edge_circulator incident_edges(Vertex_handle v) const;

  /*!
  Starts at the first edge of `f` incident to `v`, in counterclockwise order around `v`.
  \pre Face `f` is incident to vertex `v`.
  */
  Edge_circulator incident_edges(Vertex_handle v, Face_handle f) const;

  /*!
  Starts at an arbitrary face incident to `v`.
  */
  Face_circulator incident_faces(Vertex_handle v) const;

  /*!
  Starts at face `f`.
  \pre Face `f` is incident to vertex `v`.
  */
  Face_circulator incident_faces(Vertex_handle v, Face_handle f) const;

  /// @}

  /// \name Combinatorial Predicates
  ///
  /// The class `Triangulation_on_sphere_2` provides methods to test the
  /// presence in the triangulation of a particular feature (edge or face).
  ///
  /// @{

  /*!
  Returns `true` if there exists an edge having `va` and `vb` as vertices.
  */
  bool is_edge(Vertex_handle va, Vertex_handle vb);

  /*!
  Returns `true` if there exists an edge having `va` and `vb` as vertices.
  If `true` is returned, the edge with vertices `va` and `vb` is the edge `e=(fr,i)`
  where `fr` is a handle to the face incident to `e` and on the right side of `e` oriented from `va` to `vb`.
  */
  bool is_edge(Vertex_handle va, Vertex_handle vb, Face_handle& fr, int & i);

  /*!
  Returns `true` if there exists a face having `v1`, `v2` and `v3` as vertices.
  */
  bool is_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3);

  /*!
  Returns `true` if there exists a face having `v1`, `v2` and `v3` as vertices.
  If `true` is returned, fr is a handle to the face with `v1`, `v2` and `v3` as vertices.
  */
  bool is_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3, Face_handle &fr);

  /// @}

  /// \name Queries
  ///
  /// The class `Triangulation_on_sphere_2` provides methods to locate a given
  /// point with respect to a triangulation. It also provides methods to locate
  /// a point with respect to a given face of the triangulation.
  ///
  /// @{

  /*!
  Specifies which case occurs when locating a point in the triangulation.
  @todo contour type is awkward
  */
  enum Locate_type { VERTEX=0, /*!< when the located point coincides with a vertex of the triangulation */
                     EDGE, /*!< when the point is in the relative interior of an edge */
                     FACE, /*!< when the point is in the interior of a facet */
                     OUTSIDE_CONVEX_HULL, /*!< when the point is outside the convex hull but in the affine hull of the current triangulation */
                     OUTSIDE_AFFINE_HULL, /*!< when the point is outside the affine hull of the current triangulation. */
                     CONTOUR, /*!< when the face that contains the point is a ghost face, but it is in conflict with a neighboring face*/
                     NOT_ON_SPHERE, /*!< when the point is not on the sphere */
                     TOO_CLOSE /*!< when the point is too close to a vertex of the triangulation */
                   };

  /*!
  Locates the point in the triangulation, and returns information on this location.

  If the point is outside the affine hull, or if the point is (according to the traits) not
  on the sphere or too close to an existing vertex, then the returned `Face_handle` is `nullptr`.

  Otherwise, a face such that the orientation test of the three vertices of the face and `query`
  is positive is returned.

  The optional `Face_handle` argument, if provided, is used as a hint
  of where the locate process should start its search.
  */
  Face_handle locate(const Point& query, Face_handle f = Face_handle()) const;

  /*!
  Same as above. Additionally, the parameters `lt` and `li` describe where the query point is located.
  The variable `lt` is set to the locate type of the query.
  If `lt==VERTEX`, the variable `li` is set to the index of the vertex,
  and if `lt==EDGE` `li` is set to the index of the vertex opposite to the edge.
  Be careful that `li` has no meaning when the query type is `FACE`, `OUTSIDE_CONVEX_HULL`,
  or `OUTSIDE_AFFINE_HULL` or when the triangulation is \f$ 0\f$-dimensional.
  */
  Face_handle locate(const Point& query, Locate_type& lt, int& li, Face_handle h = Face_handle()) const;

  /// @}

  /// \name Miscellaneous
  ///
  /// @{

  /*!
  tests the validity of the triangulation as a `Triangulation_on_sphere_2`
  This method is mainly useful for debugging Delaunay triangulation algorithms designed by the user.
  */
  bool is_valid(bool verbose = false, int level = 0) const;

  /// @}
};

} // namespace CGAL
