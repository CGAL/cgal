namespace CGAL {

/*!
\ingroup PkgTriangulationOnSphere2TriangulationClasses

The class `Triangulation_on_sphere_2` is the basic class designed to represent a triangulation
of a point set on a sphere: its vertices coincide with the points of the set.

\warning This triangulation supports neither the insertion nor the removal of vertices,
see `CGAL::Delaunay_triangulation_on_sphere_2` for such purposes.

This triangulation class is very similar to `CGAL::Triangulation_2` as both classes represent
triangulations of a 2-manifold domain without boundary. A significant difference is that
in the case of Euclidean 2D triangulation, it is necessary to introduce so-called <i>infinite
faces</i> to complete the convex hull into an actual 2-manifold without boundary that the triangulation
data structure can represent. This is not necessary for triangulations on the sphere,
which are already perfectly adapted to the triangulation data structure.

There is an exception to the previous statement: in the degenerate configuration
where all points lie on the same hemisphere, the triangulation has a border.
Internally, the triangulation data structure must however remain a 2-manifold at all time,
and to ensure this property, fictitious faces called <i>ghost faces</i> are added. We call faces that
are not ghost faces <em>solid faces</em>, and edges of such faces <em>solid edges</em>.

\tparam Traits is the geometric traits; it must be a model of the concept `TriangulationOnSphereTraits_2`.

\tparam TDS is the triangulation data structure; it must be a model of the concept `TriangulationDataStructure_2`,
        whose vertex base must be a model of `TriangulationOnSphereVertexBase_2` and whose face base
        must be a model of `TriangulationOnSphereFaceBase_2`. By default, the triangulation data structure
        is instantiated by `Triangulation_data_structure_2<Triangulation_on_sphere_vertex_base_2<Gt>, Triangulation_on_sphere_face_base_2<Gt> >`.

\sa `CGAL::Delaunay_triangulation_on_sphere_2<Traits, TDS>`
*/
template <typename Traits, typename TDS>
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
  Size type (an unsigned integral type).
  */
  typedef Triangulation_data_structure::size_type size_type;

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

  /*!
  The 3D segment type.
  */
  typedef Traits::Segment_3 Segment_3;

  /*!
  The 3D triangle type.
  */
  typedef Traits::Triangle_3 Triangle_3;

  /*!
  An arc of a great circle, used to represent a curved segment (Voronoi or Delaunay edge).
  */
  typedef Traits::Arc_on_sphere_2 Arc_on_sphere_2;

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
  /// operators `*` and `->`. The handles are also model of the concepts
  ///  `LessThanComparable` and `Hashable`, that is they can be used as keys
  /// in containers such as `std::map` and `std::unordered_map`.
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
  typedef TDS::Vertex_iterator Vertices_iterator;

  /*!
  Range type for iterating over all vertices, with a nested
  type `iterator` that has as value type `Vertex_handle`.
  */
  typedef Iterator_range<unspecified_type> Vertex_handles;

  /*!
  Iterator over all edges.
  */
  typedef TDS::Edge_iterator All_edges_iterator;

  /*!
  Range type for iterating over all edges (including non-solid ones).
  */
  typedef Iterator_range<All_edges_iterator> All_edges;

  /*!
  Iterator over all faces.
  */
  typedef TDS::Face_iterator All_faces_iterator;

  /*!
  Range type for iterating over all faces (including ghost faces),  with a nested
  type `iterator` that has as value type `Face_handle`.
  */
  typedef Iterator_range<All_faces_iterator> All_face_handles;

  /*!
  Iterator over all solid edges.
  */
  typedef unspecified_type Solid_edges_iterator;

  /*!
  Range type for iterating over all solid edges.
  */
  typedef Iterator_range<Solid_edges_iterator> Solid_edges;

  /*!
  Iterator over all solid faces.
  */
  typedef unspecified_type Solid_faces_iterator;

  /*!
  Range type for iterating over solid faces, with a nested
  type `iterator` that has as value type `Face_handle`.
  */
  typedef Iterator_range<unspecified_type> Solid_face_handles;

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

  /*!
  Iterator over the points corresponding the vertices of the triangulation.
  */
  typedef unspecified_type Point_iterator;

  /*!
  Range type for iterating over the points of the finite vertices.
  */
  typedef Iterator_range<Point_iterator> Points;


  /// @}

public:
  /// \name Creation
  /// @{

  /*!
  constructs an empty triangulation.

  \note The values for the center and radius must be either already set in the traits
  (if `gt` is passed) or must be set after the construction, using the function `set_center_and_radius()`.
  */
  Triangulation_on_sphere_2(const Traits& gt = Traits());

  /*!
  constructs an empty triangulation and sets the center and radius to `c` and `r` respectively.
  */
  Triangulation_on_sphere_2(const Point_3& c, const FT r);

  /*!
  Copy constructor. All the vertices and faces are duplicated.
  After the copy, `*this` and `tr` refer to different triangulations:
  if `tr` is modified, `*this` is not.
  */
  Triangulation_on_sphere_2(const Triangulation_on_sphere_2& tr);

  /*!
  Assignment operator. This performs a deep copy of the triangulation, duplicating both vertices and faces.
  */
  Triangulation_on_sphere_2& operator=(Triangulation_on_sphere_2<Traits,TDS> tr);

  /*!
  The triangulations `tr` and `*this` are swapped.
  */
  void swap(Triangulation_on_sphere_2& tr);

  /*!
  deletes all faces and vertices, resulting in an empty triangulation.
  */
  void clear();

  /// @}

public:
  /// \name Ghost Predicates
  /// @{

  /*!
  returns `true` if `f` is a ghost face, and `false` otherwise.

  \pre `dimension() == 2`
  */
  bool is_ghost(const Face_handle f) const;

  /*!
  returns `true` if `e` is a ghost edge, that is if both its incident faces are ghost faces,
  and `false` otherwise.

  \pre `dimension() == 2`
  */
  bool is_ghost(const Edge& e) const;

 /// @}

public:
  /// \name Modifying the domain
  ///
  /// The following functions can be used to modify the geometry of the spherical domain.
  ///
  /// \warning The triangulation is cleared in the process.
  ///
  /// @{

  /*! */
  void set_center(const Point_3& c);

  /*! */
  void set_radius(const FT radius);

  /*! */
  void set_center_and_radius(const Point_3& c, const FT radius);

  /// @}

public:
  /// \name Access Functions
  /// @{

  /*!
  returns:
  - `-2` if the triangulation is empty
  - `-1` if the triangulation contains a single vertex
  - `0` if the triangulation contains exactly two vertices
  - `1` if the triangulation contains three (or more) coplanar vertices
  - `2` if the triangulation contains at least four non-coplanar vertices

  Note that a triangulation of dimension `1` is just a polygon drawn on a circle. The polygon is
  not triangulated itself. Thus the triangulation of dimension one consists of one polygon and
  has no faces.
  */
  int dimension() const;

  /*!
  returns a const reference to the triangulation traits object.
  */
  const Geom_traits& geom_traits() const;

  /*!
  returns a const reference to the triangulation data structure.
  */
  const Triangulation_data_structure& tds() const;

  /*!
  returns the number of vertices.
  */
  size_type number_of_vertices() const;

  /*!
  returns the number of faces. Note that this includes ghost faces.
  */
  size_type number_of_faces() const;

  /*!
  returns the number of ghost faces.
  */
  size_type number_of_ghost_faces() const;

  /*!
  returns the number of solid faces.
  */
  size_type number_of_solid_faces() const;

  /// @}

  /// \name Geometric Access Functions
  /// @{

  /*!
  returns the geometric position of the vertex `v`.
  */
  const Point& point(const Vertex_handle v);

  /*!
  returns the geometric position of the `i`-th vertex of the face `f`.
  */
  const Point& point(const Face_handle f, const int i);

  /*!
  returns the 3D line segment formed by the vertices of the edge `e`.
  \pre `t.dimension()` \f$ \geq1\f$.
  */
  Segment_3 segment(const Edge& e) const;

  /*!
  returns the 3D line segment formed by the vertices of the edge `(f, i)`.
  \pre `t.dimension()` \f$ \geq1\f$.
  */
  Segment_3 segment(const Face_handle f, int i) const;

  /*!
  returns the great circle arc formed by the vertices of the edge `e`.
  \pre `t.dimension()` \f$ \geq1\f$.
  */
  Arc_on_sphere_2 segment_on_sphere(const Edge& e) const;

  /*!
  returns the great circle arc formed by the vertices of the edge `(f, i)`.
  \pre `t.dimension()` \f$ \geq1\f$.
  */
  Arc_on_sphere_2 segment_on_sphere(const Face_handle f, int i) const;

  /*!
  returns the 3D triangle formed by the three vertices of the face `f`.
  \pre `t.dimension()` \f$ \geq2\f$.
  */
  Triangle_3 triangle(const Face_handle f) const;

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
  Starts at an arbitrary vertex.
  */
  Vertices_iterator vertices_begin() const;

  /*!
  Past-the-end iterator
  */
  Vertices_iterator vertices_end() const;

  /*!
  returns a range of iterators over all vertices.
  \note While the value type of `Vertices_iterator` is `Vertex`, the value type of
        `Vertex_handles::iterator` is `Vertex_handle`.
  */
  Vertex_handles vertex_handles() const;

  /*!
  Starts at an arbitrary edge.
  */
  All_edges_iterator all_edges_begin() const;

  /*!
  Past-the-end iterator
  */
  All_edges_iterator all_edges_end() const;

  /*!
  returns a range of iterators over all edges.
  */
  All_edges all_edges() const;

  /*!
  Starts at an arbitrary face.
  */
  All_faces_iterator all_faces_begin() const;

  /*!
  Past-the-end iterator
  */
  All_faces_iterator all_faces_end() const;

  /*!
  returns a range of iterators over all faces.
  \note While the value type of `All_faces_iterator` is `Face`, the value type of
        `All_face_handles::iterator` is `Face_handle`.
  */
  All_face_handles all_face_handles() const;

  /*!
  Starts at an arbitrary face.
  */
  Solid_faces_iterator solid_faces_begin() const;

  /*!
  Past-the-end iterator
  */
  Solid_faces_iterator solid_faces_end() const;

  /*!
  returns a range of iterators over all solid faces.
  \note While the value type of `Solid_faces_iterator` is `Face`, the value type of
        `Solid_face_handles::iterator` is `Face_handle`.
  */
  Solid_face_handles solid_face_handles() const;

  /*!
  Starts at an arbitrary face.
  */
  Solid_edges_iterator solid_edges_begin() const;

  /*!
  Past-the-end iterator
  */
  Solid_edges_iterator solid_edges_end() const;

  /*!
  returns a range of iterators over all solid edges.
  */
  Solid_edges solid_edges() const;

  /*!
  returns a range of iterators over all the points of the triangulations.
  */
  Points points() const;

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
  Starts at an arbitrary vertex adjacent to the vertex `v`.
  */
  Vertex_circulator adjacent_vertices(Vertex_handle v) const;

  /*!
  Starts at the first vertex of `f` adjacent to `v` in counterclockwise order around `v`.
  \pre The face `f` is incident to the vertex `v`.
  */
  Vertex_circulator adjacent_vertices(Vertex_handle v, Face_handle f);

  /*!
  Starts at an arbitrary edge incident to the vertex `v`.
  */
  Edge_circulator incident_edges(Vertex_handle v) const;

  /*!
  Starts at the first edge of `f` incident to `v`, in counterclockwise order around `v`.
  \pre The face `f` is incident to the vertex `v`.
  */
  Edge_circulator incident_edges(Vertex_handle v, Face_handle f) const;

  /*!
  Starts at an arbitrary face incident to the vertex `v`.

  Note that this may contain ghost faces.
  */
  Face_circulator incident_faces(Vertex_handle v) const;

  /*!
  Starts at face `f`.
  \pre The face `f` is incident to the vertex `v`.
  */
  Face_circulator incident_faces(Vertex_handle v, Face_handle f) const;

  /// @}

  /// \name Combinatorial Predicates
  ///
  /// The class `Triangulation_on_sphere_2` provides methods to test the
  /// presence of a particular feature (edge or face) in the triangulation.
  ///
  /// @{

  /*!
  returns `true` if there exists an edge (ghost or solid) having `va` and `vb` as vertices.
  */
  bool is_edge(Vertex_handle va, Vertex_handle vb);

  /*!
  returns `true` if there exists an edge (ghost or solid) having `va` and `vb` as vertices.
  If `true` is returned, the edge with vertices `va` and `vb` is the edge `e=(fr,i)`
  where `fr` is a handle to the face incident to `e` and on the right side of `e` oriented from `va` to `vb`.
  */
  bool is_edge(Vertex_handle va, Vertex_handle vb, Face_handle& fr, int& i);

  /*!
  returns `true` if there exists a face (ghost or solid) having `v1`, `v2` and `v3` as vertices.
  */
  bool is_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3);

  /*!
  returns `true` if there exists a face (ghost or solid) having `v1`, `v2` and `v3` as vertices.
  If `true` is returned, `fr` is a handle to the face with `v1`, `v2` and `v3` as vertices.
  */
  bool is_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3, Face_handle& fr);

  /// @}

  /// \name Queries
  ///
  /// The class `Triangulation_on_sphere_2` provides methods to locate a given
  /// point with respect to a triangulation. It also provides methods to locate
  /// a point with respect to a given face of the triangulation.
  ///
  /// @{

  /*!
  specifies which case occurs when locating a query point in the triangulation.
  */
  enum Locate_type { VERTEX=0, /*!< when the point coincides with a vertex of the triangulation */
                     EDGE, /*!< when the point is in the relative interior of an edge */
                     FACE, /*!< when the point is in the interior of a face */
                     OUTSIDE_CONVEX_HULL, /*!< when the point is on the same 3D plane as the existing vertices, but not on an existing edge */
                     OUTSIDE_AFFINE_HULL, /*!< when the insertion of the point would increase the dimension of the triangulation. */
                     NOT_ON_SPHERE, /*!< when the point is not on the sphere */
                     TOO_CLOSE /*!< when the point is too close to a vertex of the triangulation (see `TriangulationOnSphereTraits_2` for more details) */
                   };

  /*!
  locates the point `query` in the triangulation, and returns information on this location.

  If the point is (according to the traits) not on the sphere or is too close to an existing vertex,
  or if the dimension of the triangulation is not 2, or if the point is outside the affine hull,
  then the returned `Face_handle` is `nullptr`.

  Otherwise, a face intersected by the ray spawned from the center of the sphere
  and passing through the point is returned.

  The optional `Face_handle` argument, if provided, is used as a hint
  of where the locate process should start its search.
  */
  Face_handle locate(const Point& query, Face_handle f = Face_handle()) const;

  /*!
  Same as above. Additionally, the parameters `lt` and `li` describe where the query point is located.
  The variable `lt` is set to the locate type of the query.
  If `lt==VERTEX`, the variable `li` is set to the index of the vertex,
  and if `lt==EDGE` `li` is set to the index of the vertex opposite to the edge.
  Note that `li` has no meaning when the query type is `FACE`, `OUTSIDE_CONVEX_HULL`,
  or `OUTSIDE_AFFINE_HULL`.
  */
  Face_handle locate(const Point& query, Locate_type& lt, int& li, Face_handle h = Face_handle()) const;

  /// @}

  /// \name Miscellaneous
  ///
  /// @{

  /*!
  tests the validity of the triangulation as a `Triangulation_on_sphere_2`.

  This function tests the validity of the underlying data structure (using the function
  `TriangulationDataStructure_2::is_valid()`), and the validity of the triangulation itself
  (consistency between the dimension and the number of simplices, face orientation checks, etc.).

  This method is mainly useful for debugging Delaunay triangulation algorithms.
  */
  bool is_valid(bool verbose = false, int level = 0) const;

  /// @}
};

} // namespace CGAL
