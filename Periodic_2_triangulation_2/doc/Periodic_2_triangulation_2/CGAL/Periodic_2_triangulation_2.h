// Copyright (c) 1997-2013 INRIA Sophia-Antipolis (France).
// All rights reserved.

namespace CGAL
{

/*!
\ingroup PkgPeriodic2Triangulation2MainClasses

The class `Periodic_2_triangulation_2` represents a 2-dimensional
triangulation of a point set in \f$ \mathbb T_c^2\f$.

\tparam Traits is the geometric traits, it
has to be instantiated by a model of the concept
`Periodic_2TriangulationTraits_2`.

\tparam Tds is the triangulation data data structure and must be a model of `TriangulationDataStructure_2`
whose vertex and face are models of `Periodic_2TriangulationVertexBase_2` and `Periodic_2TriangulationFaceBase_2`.
It defaults to:
\code
CGAL::Triangulation_data_structure_2<
  CGAL::Periodic_2_triangulation_vertex_base_2<Gt>,
  CGAL::Periodic_2_triangulation_face_base_2<Gt> > >
\endcode

\cgalHeading{Traversal of the Triangulation}

The periodic triangulation class provides several iterators and circulators that allow one to traverse it.

\cgalHeading{I/O}

The I/O operators are defined for `iostream`. The format for the
iostream is an internal format.

The information in the `iostream` is:

- the original domain
- the number of sheets of the covering space as in `number_of_sheets()`
- the number of vertices
- the non-combinatorial information of vertices (point resp. point-offset pairs, etc.)
- the number of faces
- the indices of the vertices of each face
- the indices of the neighbors of each face, where the index
  corresponds to the preceding list of faces
- the offsets corresponding to the vertices of the faces
- the non-combinatorial information of each face

<p></p><!-- work around for a doxygen bug -->

\cgalHeading{Implementation}

Locate is implemented by a randomized walk from a vertex of the
face given as optional parameter (or from an arbitrary vertex of if no
optional parameter is given).

Insertion of a point is done by locating a face that contains the
point, and then splitting this face. Apart from the location,
insertion takes a time \f$ O(1)\f$.

Removal of a vertex is more difficult than in the Euclidean space,
since the star of a vertex may not be disjoint from the star of a
virtual copy of that vertex. Therefore generic removal of vertices is
not implemented. Several more constrained cases are implemented: the
removal of the last vertex in the triangulation and the removal of a
vertex of degree 3.

The face, edge, and vertex iterators on features are derived from
their counterparts visiting all (non-virtual and virtual) features
which are themselves derived from the corresponding iterators of the
triangulation data structure.

\sa `CGAL::Periodic_2_triangulation_2<Traits,Tds>`
\sa `CGAL::Periodic_2_triangulation_hierarchy_2<Tr>`
\sa `CGAL::Triangulation_2<Traits, Tds>`
*/
template< typename Traits, typename Tds >
class Periodic_2_triangulation_2 : public Triangulation_cw_ccw_2
{
public:

/// The enum Iterator_type is defined by `Periodic_2_triangulation_2` to
/// specify the behavior of geometric iterators.
  enum Iterator_type
  {
    /// Return all geometric primitives as they are
    /// stored internally in `Triangulation_data_structure_2`.
    STORED = 0,
    /// Return only one representative of each geometric
    /// primitive even if the triangulation is computed in a multiply
    /// sheeted covering space. Choose the representative whose maximum
    /// offset is minimal but non-negative in each direction of space.
    UNIQUE,
    /// Same as `STORED` but return
    /// additionally all primitives whose intersection with the original
    /// domain of the current covering space is non-empty.
    STORED_COVER_DOMAIN,
    /// Same as `UNIQUE` but return
    /// additionally all primitives whose intersection with the original
    /// domain is non-empty.
    UNIQUE_COVER_DOMAIN
  };

/// The enum `@` is defined by `Periodic_2_triangulation_2` to
/// specify which case occurs when locating a point in the
/// triangulation. If the triangulation does not contain any points
/// `EMPTY` is returned.
  enum Locate_type
  {
    /// when the located point coincides with a vertex of the triangulation
    VERTEX = 0,
    /// when the point is in the relative interior of an edge
    EDGE,
    /// when the point is in the interior of a facet
    FACE,
    /// when the triangulation is empty
    EMPTY
  };

/// \name Types
/// @{

  /*!
  the traits class.
  */
  typedef Traits Geom_traits;

  /*!
  the triangulation data structure type.
  */
  typedef Tds Triangulation_data_structure;

  /*!
  the offset type.
  */
  typedef Geom_traits::Periodic_2_offset_2 Offset;

  /*!
  the iso rectangle type.
  */
  typedef Geom_traits::Iso_rectangle_2 Iso_rectangle;

  /*!
  integer tuple to store the number of sheets in each direction of space.
  */
  typedef array<int, 2> Covering_sheets;

  /*!
  the point type.
  */
  typedef Geom_traits::Point_2 Point;

  /*!
  the segment type.
  */
  typedef Geom_traits::Segment_2 Segment;

  /*!
  the triangle type.
  */
  typedef Geom_traits::Triangle_2 Triangle;

  /*!
  represents a point-offset pair. The point in the pair lies in the
  original domain.
  */
  typedef std::pair< Point, Offset > Periodic_point;

  /*!
  a pair of periodic points representing a segment in the periodic
  domain.
  */
  typedef array< Periodic_point, 2> Periodic_segment;

  /*!
  a triple of periodic points representing a triangle in the periodic
  domain.
  */
  typedef array< Periodic_point, 3>
  Periodic_triangle;

  /*!
  the vertex type.
  */
  typedef Tds::Vertex Vertex;

  /*!
  the face type.
  */
  typedef Tds::Face Face;

  /*!
  the edge type.
  */
  typedef Tds::Edge Edge;

  /*!
  size type (an unsigned integral type).
  */
  typedef Tds::size_type size_type;

  /*!
  difference type (a signed integral type).
  */
  typedef Tds::difference_type difference_type;

/// @}

  /*!
  \name Handles, Iterators and Circulators

  The vertices and faces of the triangulations are accessed through
  `handles`, `iterators` and `circulators`. The handles are \cgalModels of
  the concept `Handle` which basically offers the two dereference
  operators and `->`. The iterators and circulators are all
  bidirectional and non-mutable. The circulators and iterators are
  convertible to handles with the same value type, so that whenever a
  handle appear in the parameter list of a function, an appropriate
  iterator or circulator can be passed as well.

  The edges of the triangulation can also be visited through iterators
  and circulators, the edge circulators and iterators are also
  bidirectional and non-mutable.
  */
/// @{

  /*!
  handle to a vertex.
  */
  typedef Tds::Vertex_handle Vertex_handle;

  /*!
  handle to a face.
  */
  typedef Tds::Face_handle Face_handle;

  /*!
  iterator over all faces.
  */
  typedef Tds::Face_iterator Face_iterator;

  /*!
  iterator over all edges.
  */
  typedef Tds::Edge_iterator Edge_iterator;

  /*!
  iterator over all vertices.
  */
  typedef Tds::Vertex_iterator Vertex_iterator;

  /*!
  iterator over the vertices whose
  corresponding points lie in the original domain, i.e. for each set
  of periodic copies the `Unique_vertex_iterator` iterates over
  exactly one representative.
  */
  typedef unspecified_type Unique_vertex_iterator;

  /*
  \cgalAdvancedBegin
  For compatibility with `Triangulation_2`.
  \cgalAdvancedEnd
  */
  typedef Face_iterator Finite_faces_iterator;

  /*
  \cgalAdvancedBegin
  For compatibility with `Triangulation_2`.
  \cgalAdvancedEnd
  */
  typedef Edge_iterator Finite_edges_iterator;

  /*
  \cgalAdvancedBegin
  For compatibility with `Triangulation_2`.
  \cgalAdvancedEnd
  */
  typedef Vertex_iterator Finite_vertices_iterator;

  /*
  \cgalAdvancedBegin
  For compatibility with `Triangulation_2`.
  \cgalAdvancedEnd
  */
  typedef Face_iterator All_faces_iterator;

  /*!
  circulator over all faces incident to a given vertex.
  */
  typedef unspecified_type Face_circulator;

  /*!
  circulator over all edges incident to a given vertex.
  */
  typedef unspecified_type Edge_circulator;

  /*!
  circulator over all vertices adjacent to a given vertex.
  */
  typedef unspecified_type Vertex_circulator;

/// @}

/// \name Geometric iterators:
/// @{

  /*!
  iterator over the triangles
  corresponding to faces of the triangulation.
  */
  typedef unspecified_type Periodic_triangle_iterator;

  /*!
  iterator over the segments
  corresponding to edges of the triangulation.
  */
  typedef unspecified_type Periodic_segment_iterator;

  /*!
  iterator over the points
  corresponding to vertices of the triangulation.
  */
  typedef unspecified_type Periodic_point_iterator;

/// @}

/// \name Creation
/// @{

  /*!
  Introduces an empty triangulation `t` with
  `domain` as original domain. \pre `domain` is a square.
  */
  Triangulation_2(const Iso_rectangle & domain =
                    Iso_rectangle(0, 0, 1, 1), const Geom_traits & traits =
                    Geom_traits());

  /*!
  Copy constructor. All the vertices and faces are duplicated.
  After the copy, `this` and `tr`
  refer to different triangulations:
  if `tr` is modified, `this` is not.
  */
  Triangulation_2(const Triangulation_2& tr);

  /*!
  Assignment. All the vertices and faces are duplicated.
  After the assignment, `this` and `tr`
  refer to different triangulations:
  if `tr` is modified, `this` is not.
  */
  Triangulation_2 operator=(const Triangulation_2<Traits, Tds>& tr);

  /*!
  The triangulations `tr` and `this` are swapped.
  `t.swap(tr)` should be preferred to `this` = `tr` or to
  `t(tr)` if `tr` is deleted after that.
  */
  void swap(Triangulation_2& tr);

  /*!
  Deletes all faces and vertices
  resulting in an empty triangulation.
  */
  void clear();

/// @}

/// \name Access Functions
/// The responsibility of keeping a valid triangulation belongs to the
/// user when using advanced operations allowing a direct manipulation
/// of the `tds`.
/// @{

  /*!
  Returns a const reference to the triangulation traits object.
  */
  const Geom_traits& geom_traits() const;

  /*!
  Returns a const reference to the triangulation data structure.
  */
  const Triangulation_data_structure_2 & tds() const;

  /*!
  Returns the original domain.
  */
  Iso_rectangle domain() const;

  /*!
  Returns the number of sheets of the covering space the triangulation is
  currently computed in.
  */
  Covering_sheets number_of_sheets() const;

  /*!
  Returns the dimension of the convex hull. The dimension is zero if
  the triangulation is empty and two otherwise.
  */
  int dimension() const;

  /*!
  Returns the number of vertices. Counts all vertices that are
  representatives of the same point in \f$ \mathbb T_c^2\f$ as one vertex.
  */
  size_type number_of_vertices() const;

  /*!
  Returns the number of faces. Counts all faces that are
  representatives of the same triangle in \f$ \mathbb T_c^2\f$ as one face.
  */
  size_type number_of_faces() const;

  /*!
  Returns the number of vertices in the data structure. This is the
  same as the number of sheets times `number_of_vertices()`.
  */
  size_type number_of_stored_vertices() const;

  /*!
  Returns the number of faces in the data structure. This is the
  same as the number of sheets times `number_of_faces()`.
  */
  size_type number_of_stored_faces() const;

/// @}

/// \name Non const access
/// \cgalAdvancedBegin
/// This method is mainly a help for users implementing
/// their own triangulation algorithms.
/// \cgalAdvancedEnd
/// @{

  /*!
  \cgalAdvancedFunction
  \cgalAdvancedBegin
  Returns a reference to the triangulation data structure.
  \cgalAdvancedEnd
  */
  Triangulation_data_structure_2 & tds();

/// @}

/// \name Non-constant-time access functions
/// @{

  /*!
  Returns the number of edges. Counts all edges that are
  representatives of the same segment in \f$ \mathbb T_c^2\f$ as one edge.
  */
  size_type number_of_edges() const;

  /*!
  Returns the number of edges in the data structure. This is the same
  as the number of sheets times `number_of_edges()`.
  */
  size_type number_of_stored_edges() const;

/// @}

/// \name Non-constant-time queries and conversions
/// \cgalAdvancedBegin
/// It is not recommended to interfere with the built-in
/// covering management. Especially a premature conversion to the
/// 1-sheeted covering space might lead to problems when modifying the
/// triangulation later.
/// \cgalAdvancedEnd
/// @{

  /*!
  \cgalAdvancedFunction
  \cgalAdvancedBegin
  The current triangulation remains a triangulation in the 1-sheeted
  covering space even after adding points if this method returns
  `true`. This test relies on a heuristic, i.e. if it answers
  `false` nothing is known. This test runs in constant-time when
  not computing in the 1-sheeted covering space. (This test uses the length
  of the longest edge in the triangulation as a
  criterion \cgalCite{cgal:ct-c3pt-09}.)
  \cgalAdvancedEnd
  */
  bool is_extensible_triangulation_in_1_sheet_h1() const;

  /*!
  \cgalAdvancedFunction
  \cgalAdvancedBegin
  The same as `is_extensible_triangulation_in_1_sheet_h1()` but with
  a more precise heuristic, i.e. it might answer `true` in cases in which
  `is_extensible_triangulation_in_1_sheet_h1()` would not. However, it is
  much less time efficient when not computing in the 1-sheeted covering
  space. (This test uses the diameter of the largest empty circle in the
  input point set as a criterion \cgalCite{cgal:ct-c3pt-09}.)
  \cgalAdvancedEnd
  */
  bool is_extensible_triangulation_in_1_sheet_h2() const;

  /*!
  \cgalAdvancedFunction
  \cgalAdvancedBegin
  Returns `true` if the current triangulation would still be a
  triangulation in the 1-sheeted covering space, returns `false` otherwise.
  \cgalAdvancedEnd
  */
  bool is_triangulation_in_1_sheet() const;

  /*!
  \cgalAdvancedFunction
  \cgalAdvancedBegin
  Converts the current triangulation into the same periodic
  triangulation in the 1-sheeted covering space.
  \pre `is_triangulation_in_1_sheet()`
  \cgalAdvancedEnd
  */
  void convert_to_1_sheeted_covering();

  /*!
  \cgalAdvancedFunction
  \cgalAdvancedBegin
  Converts the current triangulation into the same periodic
  triangulation in the 9-sheeted covering space.
  \cgalAdvancedEnd
  */
  void convert_to_9_sheeted_covering();

/// @}

/// \name Geometric access functions
/// @{

  /*!
  Returns the periodic point given by vertex `v`. If `this` is
  represented in the 1-sheeted covering space, the offset is always
  zero. Otherwise `v` can correspond to a periodic copy outside the
  `domain` of an input point.
  */
  Periodic_point periodic_point(const Vertex_handle v) const;

  /*!
  If `this` is represented in the 1-sheeted covering space, this
  function returns the periodic point given by the \f$ i\f$-th vertex of
  face `f`, that is the point in the original domain and the
  offset of the vertex in `f`. If `this` is represented in the
  9-sheeted covering space, this offset is possibly added to another
  offset determining the periodic copy. \pre \f$ i \in\{0,1,2\}\f$
  */
  Periodic_point periodic_point(const Face_handle f, int i)
  const;

  /*!
  Returns the periodic segment formed by the two point-offset pairs
  corresponding to the two vertices of edge `(f,i)`.
  \pre \f$ i \in\{0,1,2\}\f$
  */
  Periodic_segment periodic_segment(const Face_handle f, int
                                    i) const;

  /*!
  Same as the previous method for edge `e`.
  */
  Periodic_segment periodic_segment(const Edge & e) const;

  /*!
  Returns the periodic triangle formed by the three point-offset pairs
  corresponding to the three vertices of facet `f`.
  */
  Periodic_triangle periodic_triangle(const Face_handle f) const;

/// @}

/// \name
/// Note that a traits class providing exact constructions should be
/// used in order to guarantee the following operations to be exact
/// (as opposed to computing the triangulation only, which requires
/// only exact predicates).
/// @{

  /*!
  Converts the `Periodic_point` `pp` (point-offset pair) to the
  corresponding `Point` in \f$ \mathbb R^3\f$.
  */
  Point point(const Periodic_point & pp ) const;

  /*!
  Converts the `Periodic_segment` `s` to a `Segment`.
  */
  Segment segment(const Periodic_segment & s) const;

  /*!
  Converts the `Periodic_triangle` `this` to a `Triangle`.
  */
  Triangle triangle(const Periodic_triangle & t) const;

  /*!
  Equivalent to
  the call `t.segment(t.periodic_segment(f,i));`
  */
  Segment segment(Face_handle f, int i) const;

  /*!
  Equivalent to the
  call `t.segment(t.periodic_segment(e));`
  */
  Segment segment(const Edge& e) const;

  /*!
  Equivalent to the call `t.segment(t.periodic_segment(ec->first, ec->second));`
  */
  Segment segment(const Edge_circulator& ec) const;

  /*!
  Equivalent to the call `t.segment(t.periodic_segment(ei->first, ei->second));`
  */
  Segment
  segment(const Edge_iterator& ei) const;

  /*!
  Equivalent to the call `t.triangle(t.periodic_triangle(f));`
  */
  Triangle triangle(Face_handle f) const;

/// @}

/// \name Predicates
/// The class `Periodic_2_triangulation_2` provides methods to test
/// the presence in the triangulation of a particular feature (edge or
/// face).
/// @{

  /*!
  `true` if there is an edge having `va` and `vb` as
  vertices.
  */
  bool is_edge(Vertex_handle va, Vertex_handle vb);

  /*!
  as above. In addition, if `true` is returned, the edge with
  vertices `va` and `vb` is the edge `e=(fr,i)` where
  `fr` is a handle to the face incident to `e` and
  on the right side of `e` oriented from `va` to `vb`.
  */
  bool is_edge(Vertex_handle va, Vertex_handle vb, Face_handle& fr,
               int & i);

  /*!
  `true` if there is a face having `v1`, `v2` and `v3`
  as vertices.
  */
  bool is_face(Vertex_handle v1, Vertex_handle v2,
               Vertex_handle v3);

  /*!
  as above. In addition, if `true` is returned, `fr` is a
  handle to the face with `v1`, `v2` and `v3` as
  vertices.
  */
  bool is_face(Vertex_handle v1, Vertex_handle v2,
               Vertex_handle v3, Face_handle &fr);

/// @}

/// \name Queries
/// The class `Periodic_2_triangulation_2` provides methods to locate
/// a given point with respect to a triangulation. It also provides
/// methods to locate a point with respect to a given face of the
/// triangulation.
/// @{

  /*!
  If the triangulation is not empty, a face
  that contains the query in its interior or on its
  boundary is returned.

  If the triangulation is empty, the default constructed `Face_handle` is returned.

  */
  Face_handle
  locate(const Point& query,
         Face_handle f = Face_handle()) const;

  /*!
  Same as above. Additionally, the parameters `lt`
  and `li`
  describe where the query point is located.
  The variable `lt` is set to the locate type of the query.
  If `lt==VERTEX`
  the variable `li`
  is set to the index of the vertex, and if `lt==EDGE`
  `li`
  is set to the index
  of the vertex opposite to the
  edge.
  Be careful that `li`
  has no meaning when the query type is `FACE` or when the
  triangulation is \f$ 0\f$-dimensional.
  */
  Face_handle
  locate(const Point& query,
         Locate_type& lt,
         int& li,
         Face_handle h = Face_handle() ) const;

  /*!
  Returns on which side of the oriented boundary of `f`
  the point `p` lies.
  */
  Oriented_side
  oriented_side(Face_handle f,
                const Point& p) const;

/// @}

/// \name Face, Edge and Vertex Iterators
/// The following iterators allow the user to visit faces, edges and
/// vertices of the stored triangulation, i.e. in case of computing in
/// a multiply sheeted covering space all stored periodic copies of
/// each item are returned. These iterators are non-mutable,
/// bidirectional and their value types are respectively `Face`,
/// `Edge` and `Vertex`. They are all invalidated by any change in the
/// triangulation.
/// @{

  /*!
  Starts at an arbitrary vertex
  */
  Vertex_iterator vertices_begin() const;

  /*!
  Past-the-end iterator
  */
  Vertex_iterator vertices_end() const;

  /*!
  Starts at an arbitrary edge
  */
  Edge_iterator edges_begin() const;

  /*!
  Past-the-end iterator
  */
  Edge_iterator edges_end() const;

  /*!
  Starts at an arbitrary face
  */
  Face_iterator faces_begin() const;

  /*!
  Past-the-end iterator
  */
  Face_iterator faces_end() const;

/// @}

/// \name Geometric iterators
/// The following iterators allow the user to obtain geometric
/// primitives corresponding to faces, edges, and vertices of the
/// triangulation. These iterators are non-mutable, bidirectional and
/// their value types are respectively `Periodic_triangle`,
/// `Periodic_segment` and `Periodic_point`. They are all invalidated
/// by any change in the triangulation. If the periodic triangulation
/// is not computed in the 1-sheeted covering space, these iterators
/// can be used to retain only the geometric primitives in the
/// original domain. This can be controlled using the enum
/// `Iterator_type`, see \ref
/// ::CGAL::Periodic_2_triangulation_2::Iterator_type.
///
/// \anchor P2Triangulation2figgeom_iterators
/// \image html 3pts_stored.png
/// \image html 3pts_stored_cover_domain.png
/// \image html 3pts_unique.png
/// \image html 3pts_unique_cover_domain.png
/// <center><b>The four different modes of the geometric iterators:
/// `STORED`, `STORED_COVER_DOMAIN`, `UNIQUE`,
/// `UNIQUE_COVER_DOMAIN`. Note that in case of computing in the
/// 1-sheeted covering space, stored and unique give the same
/// result.</b></center>
/// @{

  /*!
  Iterates over the points of the triangulation. Its behavior is
  defined by the `Iterator_type` `it` as described on
  \ref ::CGAL::Periodic_2_triangulation_2::Iterator_type.
  */
  Periodic_point_iterator periodic_points_begin(Iterator_type it =
        STORED) const;

  /*!
  Past-the-end iterator. Note that to match another
  `Periodic_point_iterator` both must have the same
  `Iterator_type` `it`.
  */
  Periodic_point_iterator periodic_points_end(Iterator_type it =
        STORED) const;

  /*!
  Iterates over the segments of the triangulation. Its behavior is
  defined by the `Iterator_type` `it` as described on
  \ref ::CGAL::Periodic_2_triangulation_2::Iterator_type.
  */
  Periodic_segment_iterator periodic_segments_begin(Iterator_type it =
        STORED) const;

  /*!
  Past-the-end iterator. Note that to match another
  `Periodic_segment_iterator` both must have the same
  `Iterator_type` `it`.
  */
  Periodic_segment_iterator periodic_segments_end(Iterator_type it =
        STORED) const;

  /*!
  Iterates over the triangles of the triangulation. Its behavior is
  defined by the `Iterator_type` `it` as described on
  \ref ::CGAL::Periodic_2_triangulation_2::Iterator_type.
  */
  Periodic_triangle_iterator periodic_triangles_begin(Iterator_type it =
        STORED) const;

  /*!
  Past-the-end iterator. Note that to match another
  `Periodic_triangle_iterator` both must have the same
  `Iterator_type` `it`.
  */
  Periodic_triangle_iterator periodic_triangles_end(Iterator_type it =
        STORED) const;

/// @}

/// \name Face, Edge and Vertex Circulators
/// The triangulation also provides circulators that allows to visit
/// respectively all faces or edges incident to a given vertex or all
/// vertices adjacent to a given vertex. These circulators are
/// non-mutable and bidirectional. The `operator++` moves the
/// circulator counterclockwise around the vertex while the
/// `operator-` moves clockwise. A face circulator is invalidated by
/// any modification of the face pointed to. An edge or a vertex
/// circulator are invalidated by any modification of one of the two
/// faces incident to the edge pointed to.
/// @{

  /*!
  Starts at an arbitrary face incident
  to `v`.
  */
  Face_circulator incident_faces(Vertex_handle v) const;

  /*!
  Starts at face `f`.
  \pre Face `f` is incident to vertex `v`.
  */
  Face_circulator incident_faces(Vertex_handle v, Face_handle f) const;

  /*!
  Starts at an arbitrary edge incident
  to `v`.
  */
  Edge_circulator incident_edges(Vertex_handle v) const;

  /*!
  Starts at the first edge of `f` incident to
  `v`, in counterclockwise order around `v`.
  \pre Face `f` is incident to vertex `v`.
  */
  Edge_circulator incident_edges(Vertex_handle v, Face_handle f) const;

  /*!
  Starts at an arbitrary vertex adjacent to `v`.
  */
  Vertex_circulator adjacent_vertices(Vertex_handle v) const;

  /*!
  Starts at the first vertex of `f` adjacent to `v`
  in counterclockwise order around `v`.
  \pre Face `f` is incident to vertex `v`.
  */
  Vertex_circulator adjacent_vertices(Vertex_handle v, Face_handle f) ;

/// @}

/// \name Traversal between adjacent faces
/// @{

  /*!
  returns the vertex of the \f$ i^{th}\f$ neighbor of `f` that is
  opposite to `f`.
  \pre \f$ 0 \leq i \leq 2\f$.
  */
  Vertex_handle mirror_vertex(Face_handle f, int i) const;

  /*!
  returns the index of `f` in its \f$ i^{th}\f$ neighbor.
  \pre \f$0 \leq i \leq 2\f$.
  */
  int mirror_index(Face_handle f, int i) const;

/// @}

// \name Modifiers
// The following operations are guaranteed to lead to a valid
// triangulation when they are applied on a valid triangulation.
// @{

  /*
  Inserts point `p` in the triangulation and returns the
  corresponding vertex.
  If point `p` coincides with an already
  existing vertex, this vertex is returned and the triangulation
  remains unchanged.
  If point `p` is on an edge, the two
  incident faces are split in two, see
  Figure \ref Triangulation_ref_Fig_insert1.
  If point `p` is
  strictly inside a face of the triangulation, the face is split in
  three, see Figure \ref Triangulation_ref_Fig_insert2.
  If the
  triangulation is empty, the triangulation with a single vertex at
  point `p` is created.
  The last argument `f` is an
  indication to the underlying locate algorithm of where to start.
  \pre `p` lies in the original domain.
  */
  Vertex_handle insert(const Point& p, Face_handle f =
                         Face_handle());

  /*
  Same as above except that the location
  of the point `p` to be inserted is assumed to be given by
  `(lt,loc,i)` (see the description of the `locate` method
  above.)
  */
  Vertex_handle insert(const Point& p, Locate_type lt,
                       Face_handle loc, int li );

  /*
  Equivalent to
  `insert(p)`.
  */
  Vertex_handle push_back(const Point& p);

  /*
  Inserts the points in the range
  \f$ \left[\right.\f$`first`, `last`\f$ \left.\right)\f$. Returns the
  number of inserted points. \pre The `value_type` of
  `InputIterator` is `Point` and all points lie in the
  original domain.
  */
  template < class InputIterator > int insert(InputIterator
      first, InputIterator last);

// @}

/// \anchor Triangulation_ref_Fig_insert1
/// \image html insert1.png "Insertion of a point on an edge."
/// \anchor Triangulation_ref_Fig_insert2
/// \image html insert2.png "Insertion in a face."

/// \name
/// \cgalAdvancedBegin
/// The following member functions offer more specialized
/// versions of the insertion or removal operations to be used when
/// one knows to be in the corresponding case. The following functions
/// are mainly intended to be used in conjunction with the
/// `find_conflicts()` member functions of Delaunay and constrained
/// Delaunay triangulations to perform insertions.
/// \cgalAdvancedEnd
/// @{

  /*!
  \cgalAdvancedFunction
  \cgalAdvancedBegin
  Inserts the first vertex.
  \cgalAdvancedEnd
  */
  Vertex_handle insert_first(const Point& p);

  /*!
  \cgalAdvancedFunction
  \cgalAdvancedBegin
  Inserts vertex `v` in face
  `f`. Face `f` is modified,
  two new faces are created. If the triangulation contains periodic copies, a point is inserted in all periodic copies.
  \pre The point in vertex `v` lies inside face `f`.
  \cgalAdvancedEnd
  */
  Vertex_handle insert_in_face(const Point& p, Face_handle f);

  /*!
  \cgalAdvancedFunction
  \cgalAdvancedBegin
  Removes a vertex of degree three. Two of the incident faces are
  destroyed, the third one is modified. \pre Vertex `v` is a vertex with degree three.
  \cgalAdvancedEnd
  */
  void remove_degree_3(Vertex_handle v);

  /*!
  \cgalAdvancedFunction
  \cgalAdvancedBegin
  Removes the unique vertex in the
  triangulation.
  \cgalAdvancedEnd
  */
  void
  remove_first(Vertex_handle v);


  /*!
  \cgalAdvancedFunction
  \cgalAdvancedBegin
  creates a new vertex `v` and use it to star the hole
  whose boundary is described by the sequence of edges
  `[edge_begin, edge_end]`. Returns a handle to the new vertex.

  \pre The triangulation is a triangulation of 1 sheet
    \cgalAdvancedEnd
  */
  template<class EdgeIt>
  Vertex_handle star_hole( Point p,
                           EdgeIt edge_begin,
                           EdgeIt edge_end);

  /*!
  \cgalAdvancedFunction
  \cgalAdvancedBegin
  same as above, except that the algorithm
  first recycles faces in the sequence `[face_begin, face_end]`
  and create new ones only when the sequence is exhausted.

  \pre The triangulation is a triangulation of 1 sheet
  \cgalAdvancedEnd
  */
  template<class EdgeIt, class FaceIt>
  Vertex_handle star_hole( Point p,
                           EdgeIt edge_begin,
                           EdgeIt edge_end,
                           FaceIt face_begin,
                           FaceIt face_end);

  /*!
  \cgalAdvancedFunction
  \cgalAdvancedBegin
  Changes the domain. Note that this function calls `clear()`,
  i.e., it erases the existing triangulation.
  \cgalAdvancedEnd
  */
  void set_domain(const Iso_rectangle dom);

/// @}

/// \name Miscellaneous
/// @{

  /*!
  Returns \f$ i+1\f$ modulo 3.\pre \f$0 \leq i \leq 2\f$.
  */
  int ccw(int i) const;

  /*!
  Returns \f$ i+2\f$ modulo 3.\pre \f$0 \leq i \leq 2\f$.
  */
  int cw(int i) const;

  /*
  Returns whether the
  union of the faces `f` and `f->neighbor(i)` form a convex
  quadrilateral.
  */
  void flippable(Face_handle f, int i);

  /*!
  Returns the degree of the
  vertex `v`
  */
  size_t degree(Vertex_handle v);

/// @}

/// \name Checking
/// \cgalAdvancedBegin
/// The responsibility of keeping a valid triangulation
/// belongs to the users if advanced operations are used. Obviously
/// the advanced user, who implements higher levels operations may
/// have to make a triangulation invalid at some times. The following
/// method is provided to help the debugging.
/// \cgalAdvancedEnd
/// @{

  /*!
  \cgalAdvancedFunction
  \cgalAdvancedBegin
  Checks the combinatorial validity of the triangulation and
  also the validity of its geometric embedding.
  This method is mainly a debugging help
  for the users of advanced features.
  \cgalAdvancedEnd
  */
  bool
  is_valid(bool verbose = false, int level = 0) const;

/// @}

}; /* end Periodic_2_triangulation_2 */

/*!
Writes the
triangulation `t` into the stream `os`. \pre The
output operator must be defined for `Point`.
\relates Periodic_2_triangulation_2
*/
ostream& operator<<(ostream& os, const Periodic_2_triangulation_2<Traits, Tds>& t);

/*!
Reads a triangulation from stream
`is` and assigns it to `t`. \pre The input operator
must be defined for `Point`.
\relates Periodic_2_triangulation_2
*/
istream& operator>>(istream& is, Triangulation_2<Traits, Tds>& t);

} /* end namespace CGAL */
