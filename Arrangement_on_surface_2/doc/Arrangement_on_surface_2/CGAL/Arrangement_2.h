namespace CGAL {

/*!
  \ingroup PkgArrangementOnSurface2Ref

  \anchor arr_refarr

  An object `arr` of the class `Arrangement_2` represents the
  planar subdivision induced by a set of \f$ x\f$-monotone curves and isolated
  points into maximally connected cells. The arrangement is represented as
  a doubly-connected edge-list (\dcel) such that each \dcel vertex
  is associated with a point of the plane and each edge is
  associated with an \f$ x\f$-monotone curve whose interior is disjoint from all
  other edges and vertices. Recall that an arrangement
  edge is always comprised of a pair of twin \dcel halfedges.

  The `Arrangement_2` template has two parameters:
  <UL>
  <LI>The `Traits` template-parameter should be instantiated with
  a model of the `ArrangementBasicTraits_2` concept. The traits
  class defines the types of \f$ x\f$-monotone curves and two-dimensional
  points, namely `ArrangementBasicTraits_2::X_monotone_curve_2` and `ArrangementBasicTraits_2::Point_2`,
  respectively, and supports basic geometric predicates on them.
  <LI>The `Dcel` template-parameter should be instantiated with
  a class that is a model of the `ArrangementDcel` concept. The
  value of this parameter is by default
  `Arr_default_dcel<Traits>`.
  </UL>
  The available traits classes and \dcel classes are described below.

  \sa `ArrangementDcel`
  \sa `Arr_default_dcel<Traits>`
  \sa `ArrangementBasicTraits_2`
  \sa `CGAL::overlay()`
  \sa `CGAL::is_valid()`

  Insertion Functions

  \sa `PkgArrangementOnSurface2Insert`
  \sa `CGAL::insert_non_intersecting_curve()`
  \sa `CGAL::insert_non_intersecting_curves()`
  \sa `CGAL::insert_point()`

  Removal functions

  \sa `CGAL::remove_edge()`
  \sa `CGAL::remove_vertex()`

  Input/output functions

  \sa `PkgArrangementOnSurface2Read`
  \sa `PkgArrangementOnSurface2Write`

*/
template< typename Traits, typename Dcel >
class Arrangement_2 {
public:

  /// \name Types

  /// @{

  /*!
    a private type used as an abbreviation of the `Arrangement_2` type hereafter.
  */
  typedef Arrangement_2<Traits_2,Dcel> Self;

  /*!
    the traits class in use.
  */
  typedef Traits Traits_2;

  /*!
    the \dcel representation of the arrangement.
  */
  typedef unspecified_type Dcel;

  /*!
    the point type, as defined by the traits class.
  */
  typedef typename Traits_2::Point_2 Point_2;

  /*!
    the \f$ x\f$-monotone curve type, as defined by the traits class.
  */
  typedef typename Traits_2::X_monotone_curve_2 X_monotone_curve_2;

  /*!
    the size type (equivalent to `size_t`).
  */
  typedef typename Dcel::Size Size;

  /*!
    \ingroup PkgArrangementOnSurface2DCEL
    An object \f$ v\f$ of the class `Vertex` represents an arrangement vertex,
    that is - a \f$ 0\f$-dimensional cell, associated with a point on the plane.
  */
  class Vertex : public Dcel::Vertex {
  public:

    /// \name Creation
    /// @{

    /*!
      default constructor.
    */
    Vertex();

    /// @}

    /// \name Access Functions
    /// All non-const methods listed below also have `const` counterparts
    /// that return constant handles, iterators or circulators:
    /// @{

    /*!
      checks whether the vertex lies at infinity and not associated with
      a point with bounded coordinates.
    */
    bool is_at_open_boundary() const;

    /*!
      checks whether the vertex is isolated (i.e., has no incident edges).
    */
    bool is_isolated() const;

    /*!
      returns the number of edges incident to `v`.
    */
    Size degree() const;

    /*!
      returns a circulator circulator that allows going over the halfedges
      incident to `v` (that have `v` as their target).
      The edges are traversed in a clockwise direction around `v`.
      \pre `v` is <I>not</I> an isolated vertex.
    */
    Halfedge_around_vertex_circulator incident_halfedges();

    /*!
      returns a handle to the face that contains `v` in its interior.
      \pre `v` is an isolated vertex.
    */
    Face_handle face();

    /*!
      returns the point associated with the vertex.
      \pre `v` is not a vertex at infinity.
    */
    const typename Traits::Point_2& point() const;

    /*!
      returns the placement of the \f$ x\f$-coordinate in the parameter space,
      that is, either the left boundary-side, the interior, or the right
      boundary-side.
    */
    Arr_parameter_space parameter_space_in_x () const;

    /*!
      returns the placement of the \f$ y\f$-coordinate in the parameter space,
      that is, either the bottom boundary-side, the interior, or the top
      boundary-side.
    */
    Arr_parameter_space parameter_space_in_y () const;

    /// @}

  }; /* end Vertex */

  /*!
    \ingroup PkgArrangementOnSurface2DCEL
    An object \f$ e\f$ of the class `Halfedge` represents a halfedge in the
    arrangement. A halfedge is directed from its <I>source</I> vertex
    to its <I>target</I> vertex, and has an <I>incident face</I> lying to
    its left. Each halfedge has a <I>twin</I> halfedge directed in the
    opposite direction, where the pair of twin halfedges form together
    an arrangement edge, that is - a \f$ 1\f$-dimensional cell, associated
    with planar \f$ x\f$-monotone curve.

    Halfedges are stored in doubly-connected lists and form chains. These
    chains define the inner and outer boundaries of connected components.

  */
  class Halfedge : public Dcel::Halfedge {
  public:

    /// \name Creation
    /// @{

    /*!
      default constructor.
    */
    Halfedge();

    /// @}

    /// \name Access Functions
    /// All non-const methods listed below also have `const` counterparts
    /// that return constant handles, iterators or circulators:
    /// @{

    /*!
      returns whether the halfedge is fictitious (i.e., connects two vertices at
      infinity and is not associated with a valid curve).
    */
    bool is_fictitious () const;

    /*!
      returns a handle for the source vertex of `e`.
    */
    Vertex_handle source();

    /*!
      returns a handle for the target vertex of `e`.
    */
    Vertex_handle target();

    /*!
      returns the direction of the halfedge: `ARR_LEFT_TO_RIGHT` if
      `e`'s source vertex is lexicographically smaller than it
      target (so the halfedge is directed from left to right), and
      `ARR_RIGHT_TO_LEFT` if it is lexicographically larger than
      the target (so the halfedge is directed from right to left).
    */
    Arr_halfedge_direction direction() const;

    /*!
      returns the face that `e` is incident to (The face lies to
      the left of `e`).
    */
    Face_handle face();

    /*!
      returns the twin halfedge.
    */
    Halfedge_handle twin();

    /*!
      returns `e`'s predecessor in the connected component it belongs to.
    */
    Halfedge_handle prev();

    /*!
      returns `e`'s successor in the connected component it belongs to.
    */
    Halfedge_handle next();

    /*!
      returns a circulator that allows traversing the halfedges of the
      connected component boundary (CCB) that contains `e`.
      The circulator is initialized to point to `e`.
    */
    Ccb_halfedge_circulator ccb();

    /*!
      returns the \f$ x\f$-monotone curve associated with `e`.
      \pre `e` is not a fictitious halfedge.
    */
    const typename Traits::X_monotone_curve_2& curve() const;

    /// @}

  }; /* end Halfedge */

  /*!
    \ingroup PkgArrangementOnSurface2DCEL

    An object of the class `Face` represents an arrangement face,
    namely, a \f$ 2\f$-dimensional arrangement cell. An arrangement that supports
    only bounded curves contains exactly one <I>unbounded</I> face, and a
    number of bounded faces. An arrangement that supports also unbounded
    curves has one or more unbounded faces. Such an arrangement has also
    exactly one fictitious face, which does not correspond to a real
    two-dimensional cell of the arrangement (and thus it is ignored in
    counting the number of faces of the arrangement.)
    Each bounded face has an outer boundary comprising a halfedge chain
    winding in counterclockwise orientation around it. Each unbounded face of
    an arrangement that has a fictitious face also has a boundary comprising
    a counterclockwise halfedge-chain. The edges on the boundary of a face
    incident to the fictitious face are fictitious, as they do not correspond
    to real curves. A face may also contain holes, which are defined by
    clockwise-oriented halfedge chains, and isolated vertices.
  */
  class Face : public Dcel::Face {
  public:

    /// \name Creation
    /// @{

    /*!
      default constructor.
    */
    Face();

    /// @}

    /// \name Access Functions
    /// All non-const methods listed below also have `const` counterparts
    /// that return constant handles, iterators or circulators:
    /// @{

    /*!
      returns a Boolean indicating whether this is the fictitious face,
      which contain the entire arrangement (and does not have an outer CCB).
      An arrangement that supports only bounded curves does not have a
      fictitious face at all.
    */
    bool is_fictitious () const;

    /*!
      returns a Boolean indicating whether the face is unbounded.
    */
    bool is_unbounded() const;

    /*!
      returns a Boolean indicating whether the face has an outer CCB.
      (The fictitious face and the unbounded face of an arrangement that
      does not have a fictitious face do not have outer CCBs.)
    */
    bool has_outer_ccb() const;

    /*!
      returns a circulator that enables traversing the outer boundary of
      `f`. The edges along the CCB are traversed in a counterclockwise
      direction.
      \pre The face `f` has an outer CCB.
    */
    Ccb_halfedge_circulator outer_ccb();

    /*!
      returns an iterator for traversing all the holes (inner CCBs) of
      `f`.
    */
    Hole_iterator holes_begin();

    /*!
      returns a past-the-end iterator for the holes of `f`.
    */
    Hole_iterator holes_end();

    /*!
      returns an iterator for traversing all the isolated vertices
      contained in the interior of `f`.
    */
    Isolated_vertex_iterator isolated_vertices_begin();

    /*!
      returns a past-the-end iterator for the isolated vertices
      contained inside `f`.
    */
    Isolated_vertex_iterator isolated_vertices_end();

    /// @}

  }; /* end Face */


/// @}

/*! \name

The following handles, iterators, and circulators all have respective
constant counterparts (for example, in addition to `Vertex_iterator`
the type `Vertex_const_iterator` is also defined). See \cgalCite{cgal:ms-strg-96}
for a discussion of constant versus mutable iterator
types. The mutable types are assignable to their constant
counterparts. `Vertex_iterator`, `Halfedge_iterator`, and
`Face_iterator` are equivalent to the respective handle types (namely,
`Vertex_handle`, `Halfedge_handle`, and `Face_handle`). Thus, wherever
the handles appear in function parameter lists, the respective
iterators can be passed as well.

All handles are model of `LessThanComparable` and `Hashable`,
that is they can be used as keys in containers such as `std::map`
and `boost::unordered_map`.

*/
/// @{


  /*!
    a handle for an arrangement vertex.
  */
  typedef unspecified_type Vertex_handle;

  /*!
    a handle for a halfedge.
    The halfedge and its twin form together an arrangement edge.
  */
  typedef unspecified_type Halfedge_handle;

  /*!
    a handle for an arrangement face.
  */
  typedef unspecified_type Face_handle;

  /*!
    a bidirectional iterator over the
    vertices of the arrangement. Its value-type is `Vertex`.
  */
  typedef unspecified_type Vertex_iterator;

  /*!
    a bidirectional iterator over the
    halfedges of the arrangement. Its value-type is `Halfedge`.
  */
  typedef unspecified_type Halfedge_iterator;

  /*!
    a bidirectional iterator over the
    edges of the arrangement. (That is, it skips every other halfedge.)
    Its value-type is `Halfedge`.
  */
  typedef unspecified_type Edge_iterator;

  /*!
    a bidirectional iterator over the
    faces of arrangement. Its value-type is `Face`.
  */
  typedef unspecified_type Face_iterator;

  /*!
    a bidirectional iterator over the
    unbounded faces of arrangement. Its value-type is `Face`.
  */
  typedef unspecified_type Unbounded_face_iterator;

  /*!
    a bidirectional circulator
    over the halfedges that have a given vertex as their target.
    Its value-type is `Halfedge`.
  */
  typedef unspecified_type Halfedge_around_vertex_circulator;

  /*!
    a bidirectional circulator over the
    halfedges of a CCB (connected component of the boundary).
    Its value-type is `Halfedge`. Each
    bounded face has a single CCB representing it outer boundary, and may
    have several inner CCBs representing its holes.
  */
  typedef unspecified_type Ccb_halfedge_circulator;

  /*!
    a bidirectional iterator over the holes
    (i.e., inner CCBs) contained inside a given face.
    Its value type is `Ccb_halfedge_circulator`.
  */
  typedef unspecified_type Hole_iterator;

  /*!
    a bidirectional iterator over the
    isolated vertices contained inside a given face.
    Its value type is `Vertex`.
  */
  typedef unspecified_type Isolated_vertex_iterator;

  /// @}

  /// \name Creation
  /// @{

  /*!
    constructs an empty arrangement containing one unbounded face,
    which corresponds to the entire plane.
  */
  Arrangement_2<Traits, Dcel>();

  /*!
    copy constructor.
  */
  Arrangement_2<Traits, Dcel>(const Self& other);

  /*!
    constructs an empty arrangement that uses the given `traits`
    instance for performing the geometric predicates.
  */
  Arrangement_2<Traits, Dcel>(const Traits_2 *traits);

  /// @}

  /// \name Assignment Methods
  /// @{

  /*!
    assignment operator.
  */
  Self& operator= (other);

  /*!
    assigns the contents of another arrangement.
  */
  void assign (const Self& other);

  /*!
    clears the arrangement.
  */
  void clear ();

  /// @}

  /// \name Access Functions

  /// @{

  /*!
    returns the traits object used by the arrangement instance.
    A `const` version is also available.
  */
  Traits_2* get_traits();

  /*!
    determines whether the arrangement is empty (contains only the unbounded
    face, with no vertices or edges).
  */
  bool is_empty() const;


/// @}

/*! \name Accessing the Arrangement Vertices

All `_begin()` and `_end()` methods listed below also have `const` counterparts,
returning constant iterators instead of mutable ones.
*/
/// @{

  /*!
    returns the number of vertices in the arrangement.
  */
  Size number_of_vertices() const;

  /*!
    returns the total number of isolated vertices in the arrangement.
  */
  Size number_of_isolated_vertices() const;

  /*!
    returns the begin-iterator of the vertices in the arrangement.
  */
  Vertex_iterator vertices_begin();

  /*!
    returns the past-the-end iterator of the vertices in the arrangement.
  */
  Vertex_iterator vertices_end();

  /*!
  returns a range over handles of the arrangement vertices .
  */
  unspecified_type vertex_handles();
  /*!
    returns the number of arrangement vertices that lie at infinity and
    are not associated with valid points. Such vertices are not considered
    to be regular arrangement vertices and `arr.number_of_vertices()`
    does not count them.
  */
  Size number_of_vertices_at_infinity() const;


/// @}

/*! \name Accessing the Arrangement Halfedges

All `_begin()` and `_end()` methods listed below also have `const` counterparts,
returning constant iterators instead of mutable ones.
*/
/// @{

  /*!
    returns the number of halfedges in the arrangement.
  */
  Size number_of_halfedges() const;

  /*!
    returns the begin-iterator of the halfedges in the arrangement.
  */
  Halfedge_iterator halfedges_begin();

  /*!
    returns the past-the-end iterator of the halfedges in the arrangement.
  */
  Halfedge_iterator halfedges_end();

  /*!
  returns a range over handles of the arrangement halfedges .
  */
  unspecified_type halfedge_handles();

  /*!
    returns the number of edges in the arrangement (equivalent to
    `arr.number_of_halfedges() / 2`).
  */
  Size number_of_edges() const;

  /*!
    returns the begin-iterator of the edges in the arrangement.
  */
  Edge_iterator edges_begin();

  /*!
    returns the past-the-end iterator of the edges in the arrangement.
  */
  Edge_iterator edges_end();

  /*!
  returns a range over handles of the arrangement edges .
  */
  unspecified_type edge_handles();

/// @}

/*! \name Accessing the Arrangement Faces

All `_begin()` and `_end()` methods listed below also have `const` counterparts,
returning constant iterators instead of mutable ones.
*/
/// @{

  /*!
    returns a handle for an unbounded face of the arrangement.
    In case the arrangement comprises only bounded curves, there is a single
    unbounded face and the function returns a handle to it. Otherwise, a
    handle to an arbitrary unbounded face is returned.
  */
  Face_handle unbounded_face();

  /*!
    returns the number of faces in the arrangement.
  */
  Size number_of_faces() const;

  /*!
    returns the begin-iterator of the faces in the arrangement.
  */
  Face_iterator faces_begin();

  /*!
    returns the past-the-end iterator of the faces in the arrangement.
  */
  Face_iterator faces_end();

  /*!
  returns a range over handles of the arrangement faces .
  */
  unspecified_type face_handles();

  /*!
    returns the number of unbounded faces in the arrangement.
    Note `arr.number_of_faces()` also counts the unbounded faces
    of the arrangement.
  */
  Size number_of_unbounded_faces() const;

  /*!
    returns the begin-iterator of the unbounded faces in the arrangement.
  */
  Unbounded_face_iterator unbounded_faces_begin();

  /*!
    returns the past-the-end iterator of the unbounded faces in the
    arrangement.
  */
  Unbounded_face_iterator unbounded_faces_end();

  /*!
    returns a handle to the fictitious face of the arrangement.
    If the arrangement is not unbounded, there is no fictitious
    face. In this case the result is not deterministic. A const
    version is also available.
  */
  Face_handle fictitious_face();

  /// @}

/*! \name Casting Away Constness

It is sometimes necessary to convert a constant (non-mutable) handle to
a mutable handle. For example, the result of a point-location query is
a non-mutable handle for the arrangement cell containing the query
point. Assume that the query point lies on an edge, so we obtain a
`Halfedge_const_handle`; if we wish to use this handle and remove the
edge, we first need to cast away its "constness".

*/

  /// @{

  /*!
    casts the given constant vertex handle to an equivalent mutable handle.
  */
  Vertex_handle non_const_handle (Vertex_const_handle v);

  /*!
    casts the given constant halfedge handle to an equivalent mutable handle.
  */
  Halfedge_handle non_const_handle (Halfedge_const_handle e);

  /*!
    casts the given constant face handle to an equivalent mutable handle.
  */
  Face_handle non_const_handle (Face_const_handle f);

  /// @}

/// \name Specialized Insertion Methods

  /// @{

  /*!
    inserts the point `p` into the arrangement as an isolated vertex in
    the interior of the face `f` and returns a handle for the newly
    created vertex.
    \pre `p` lies in the interior of the face `f`.
  */
  Vertex_handle insert_in_face_interior (const Point_2& p,
                                         Face_handle f);

  /*!
    inserts the curve `c` that is entirely contained in the interior
    of a given face `f`. If `c` is a bounded curve two new vertices
    that correspond to `c`'s endpoints are created and connected with a
    newly created halfedge pair, which forms a new hole (inner component)
    in the face `f`. If `c` is unbounded, at least one of the two
    vertices that represents its end lies at infinity, and its creation
    modifies the outer boundary of `f`.
    The function returns a handle for one of the new halfedges
    corresponding to the inserted curve, directed in lexicographic
    increasing order (from left to right).
    \pre `c` lies entirely in the interior of the face `f` and is disjoint from all existing arrangement vertices and edges (in particular, both its endpoints are not already associated with existing arrangement vertices).
    \pre In case `c` is an unbounded curve, `f` must be an unbounded face.
  */
  Halfedge_handle insert_in_face_interior
  (const X_monotone_curve_2& c,
   Face_handle f);

  /*!
    inserts the curve `c` into the arrangement, such that its left
    endpoint corresponds to a given arrangement vertex. As a result, a new
    vertex that correspond to `c`'s right endpoint is created and
    connected to `v` with a newly created halfedge pair. If `c` has
    an unbounded right end, the new vertex lies at infinity and the
    unbounded face that contains the interior of the curve is split.
    The function returns a handle for one of the new halfedges corresponding
    to the inserted curve, directed towards the newly created vertex -
    that is, directed in lexicographic increasing order (from left to right).
    \pre The interior of `c` is disjoint from all existing arrangement vertices and edges.
    \pre `v` is associated with the left endpoint of `c`.
    \pre The right endpoint of `c` is not already associated with an existing arrangement vertex.
  */
  Halfedge_handle insert_from_left_vertex
  (const X_monotone_curve_2& c,
   Vertex_handle v);

  /*!
    inserts the curve `c` into the arrangement, such that its right
    endpoint corresponds to a given arrangement vertex. As a result, a new
    vertex that correspond to `c`'s left endpoint is created and
    connected to `v` with a newly created halfedge pair. If `c` has
    an unbounded left end, the new vertex lies at infinity and the
    unbounded face that contains the interior of the curve is split.
    The function returns a handle for one of the new halfedges corresponding
    to the inserted curve, directed to the newly created vertex -
    that is, directed in lexicographic decreasing order (from right to left).
    \pre The interior of `c` is disjoint from all existing arrangement vertices and edges.
    \pre `v` is associated with the right endpoint of `c`.
    \pre The left endpoint of `c` is not already associated with an existing arrangement vertex.
  */
  Halfedge_handle insert_from_right_vertex
  (const X_monotone_curve_2& c,
   Vertex_handle v);

  /*!
    inserts the curve `c` into the arrangement, such that both `c`'s
    endpoints correspond to existing arrangement vertices, given by `v1`
    and `v2`. The function creates a new halfedge pair that connects the
    two vertices, and returns a handle for the halfedge directed from `v1`
    to `v2`.
    \pre The interior of `c` is disjoint from all existing arrangement vertices and edges.
    \pre `c` must not be an unbounded curve.
    \pre `v1` and `v2` are associated with `c`'s endpoints.
    \pre If `v1` and `v2` are already connected by an edge, this edge represents an \f$ x\f$-monotone curve that is interior-disjoint from `c`).
  */
  Halfedge_handle insert_at_vertices (const X_monotone_curve_2& c,
                                      Vertex_handle v1,
                                      Vertex_handle v2);

  /*!
    inserts an unbounded curve `c` into the arrangement, such that `c`
    is entirely contained within a single unbounded face of the arrangement.
    `fict_pred1` specifies the fictitious halfedge that should contain the
    vertex at infinity that corresponds to the unbounded end of `c`. If
    both ends of `c` are unbounded, `fict_pred1` indicated the place
    for its left end and `fict_pred2` indicated a place for its right end.
    The function returns a handle for one of the new halfedges directed
    (lexicographically) from left to right.
    \pre `c` is an unbounded curve disjoint from all existing arrangement vertices and edges.
    \pre `fict_pred1` (and `fict_pred2`) are fictitious halfedges that contains the unbounded end(s) of `c`. If both halfedges are given they must be both incident to the same unbounded face.
  */
  Halfedge_handle insert_in_face_interior
  (const X_monotone_curve_2& c,
   Halfedge_handle fict_pred1,
   Halfedge_handle fict_pred2 = Halfedge_handle());

  /*!
    inserts the curve `c` into the arrangement, such that its left
    endpoint corresponds to a given arrangement vertex. This vertex is the
    target vertex of the halfedge `pred`, such that `c` is inserted
    to the circular list of halfedges around `pred->target()` right
    between `pred` and its successor. The function returns a handle for
    one of the new halfedges directed (lexicographically) from left to right.
    \pre The interior of `c` is disjoint from all existing arrangement vertices and edges.
    \pre `pred->target()` is associated with the left endpoint of `c`, and `c` should be inserted after `pred` in a clockwise order around this vertex.
    \pre The right endpoint of `c` is not already associated with an existing arrangement vertex.
  */
  Halfedge_handle insert_from_left_vertex
  (const X_monotone_curve_2& c,
   Halfedge_handle pred);

  /*!
    inserts an unbounded curve `c` into the arrangement, such that its left
    endpoint is bounded and corresponds to a given arrangement vertex. This
    vertex is the target vertex of the halfedge `pred`, such that `c`
    is inserted to the circular list of halfedges around `pred->target()`
    right between `pred` and its successor. Similarly, `fict_pred`
    specifies the fictitious halfedge that should contain the vertex at infinity
    that corresponds to the unbounded right end of `c`.
    The function returns a handle for one of the new halfedges directed
    (lexicographically) from left to right.
    \pre The interior of `c` is disjoint from all existing arrangement vertices and edges. `c` must have a bounded left endpoint and an unbounded right end.
    \pre `pred->target()` is associated with the left endpoint of `c`, and `c` should be inserted after `pred` in a clockwise order around this vertex.
    \pre `fict_pred` is a fictitious halfedge that contains the unbounded right end of `c`.
  */
  Halfedge_handle insert_from_left_vertex
  (const X_monotone_curve_2& c,
   Halfedge_handle pred,
   Halfedge_handle fict_pred);

  /*!
    inserts the curve `c` into the arrangement, such that its right
    endpoint corresponds to a given arrangement vertex. This vertex is the
    target vertex of the halfedge `pred`, such that `c` is inserted
    to the circular list of halfedges around `pred->target()` right
    between `pred` and its successor. The function returns a handle for
    one of the new halfedges directed (lexicographically) from right to left.
    \pre The interior of `c` is disjoint from all existing arrangement vertices and edges.
    \pre `pred->target()` is associated with the right endpoint of `c`, and `c` should be inserted after `pred` in a clockwise order around this vertex.
    \pre The left endpoint of `c` is not already associated with an existing arrangement vertex.
  */
  Halfedge_handle insert_from_right_vertex
  (const X_monotone_curve_2& c,
   Halfedge_handle pred);

  /*!
    inserts an unbounded curve `c` into the arrangement, such that its right
    endpoint is bounded and corresponds to a given arrangement vertex. This
    vertex is the target vertex of the halfedge `pred`, such that `c`
    is inserted to the circular list of halfedges around `pred->target()`
    right between `pred` and its successor. Similarly, `fict_pred`
    specifies the fictitious halfedge that should contain the vertex at infinity
    that corresponds to the unbounded left end of `c`.
    The function returns a handle for one of the new halfedges directed
    (lexicographically) from right to left.
    \pre The interior of `c` is disjoint from all existing arrangement vertices and edges. `c` must have a bounded right endpoint and an unbounded left end.
    \pre `pred->target()` is associated with the right endpoint of `c`, and `c` should be inserted after `pred` in a clockwise order around this vertex.
    \pre `fict_pred` is a fictitious halfedge that contains the unbounded left end of `c`.
  */
  Halfedge_handle insert_from_right_vertex
  (const X_monotone_curve_2& c,
   Halfedge_handle pred,
   Halfedge_handle fict_pred);

  /*!
    inserts the curve `c` into the arrangement, such that both `c`'s
    endpoints correspond to existing arrangement vertices, given by
    `pred1->target()` and `v2`. The function creates a new halfedge
    pair that connects the two vertices (where the corresponding halfedge is
    inserted right between `pred1` and its successor around `pred1`'s
    target vertex) and returns a handle for the halfedge directed from
    `pred1->target()` to `v2`.
    \pre The interior of `c` is disjoint from all existing arrangement vertices and edges.
    \pre `pred1->target()` and `v2` are associated with `c`'s endpoints.
    \pre If `pred1->target` and `v2` are already connected by an edge, this edge represents an \f$ x\f$-monotone curve that is interior-disjoint from `c`).
  */
  Halfedge_handle insert_at_vertices (const X_monotone_curve_2& c,
                                      Halfedge_handle pred1,
                                      Vertex_handle v2);

  /*!
    inserts the curve `c` into the arrangement, such that both `c`'s
    endpoints correspond to existing arrangement vertices, given by
    `pred1->target()` and `pred2->target()`. The function creates a
    new halfedge pair that connects the two vertices (with `pred1` and
    `pred2` indicate the exact place for these halfedges around the two
    target vertices) and returns a handle for the halfedge directed from
    `pred1->target()` to `pred2->target()`.
    \pre The interior of `c` is disjoint from all existing arrangement vertices and edges.
    \pre `pred1->target()` and `pred2->target()` are associated with `c`'s endpoints.
    \pre If `pred1->target` and `pred2->target()` are already connected by an edge, this edge represents an \f$ x\f$-monotone curve that is interior-disjoint from `c`).
  */
  Halfedge_handle insert_at_vertices (const X_monotone_curve_2& c,
                                      Halfedge_handle pred1,
                                      Halfedge_handle pred2);


/// @}

/// \name Modifying Vertices and Edges

/// @{

  /*!
    sets `p` to be the point associated with the vertex `v`.
    The function returns a handle for the modified vertex (same as `v`).
    \pre `v` is not a vertex at infinity and `p` is geometrically equivalent to the point currently associated with `v`.
  */
  Vertex_handle modify_vertex (Vertex_handle v,
                               const Point_2& p);

  /*!
    removes the isolated vertex `v` from the arrangement. The function
    returns the face `f` that used to contain the isolated vertex.
    \pre `v` is an isolated vertex (has no incident edges).
  */
  Face_handle remove_isolated_vertex (Vertex_handle v);

  /*!
    sets `c` to be the \f$ x\f$-monotone curve associated with the edge `e`.
    The function returns a handle for the modified edge (same as `e`).
    \pre `c` is geometrically equivalent to the curve currently associated with `e`.
  */
  Halfedge_handle modify_edge (Halfedge_handle e,
                               const X_monotone_curve_2& c);

  /*!
    splits the edge `e` into two edges (more precisely, into two halfedge
    pairs), associated with the given subcurves `c1` and `c2`, and
    creates a vertex that corresponds to the split point.
    The function returns a handle for the halfedge, whose source is the same
    as `e->source()` and whose target vertex is the split point.
    \pre Either `c1`'s left endpoint and `c2`'s right endpoint correspond to `e`'s end-vertices such that `c1`'s right endpoint and `c2`'s left endpoint are equal and define the split point - or vice-versa (with change of roles between `c1` and `c2`).
  */
  Halfedge_handle split_edge (Halfedge_handle e,
                              const X_monotone_curve_2& c1,
                              const X_monotone_curve_2& c2);

  /*!
    merges the edges represented by `e1` and `e2` into
    a single edge, associated with the given merged curve `c`.
    Denote `e1`'s end-vertices as \f$ u_1\f$ and \f$ v\f$, while `e2`'s
    end-vertices are denoted \f$ u_2\f$ and \f$ v\f$. The function removes the
    common vertex \f$ v\f$ returns a handle for one of the merged halfedges,
    directed from \f$ u_1\f$ to \f$ u_2\f$.
    \pre `e1` and `e2` share a common end-vertex, such that the two other end-vertices of the two edges are associated with `c`'s endpoints.
  */
  Halfedge_handle merge_edge (Halfedge_handle e1,
                              Halfedge_handle e2,
                              const X_monotone_curve_2& c);

  /*!
    removes the edge `e` from the arrangement. Since the `e` may
    be the only edge incident to its source vertex (or its target vertex),
    this vertex can be removed as well. The flags `remove_source` and
    `remove_target` indicate whether the endpoints of `e` should be
    removed, or whether they should be left as isolated vertices in the
    arrangement.
    If the operation causes two faces to merge, the merged face is returned.
    Otherwise, the face to which the edge was incident is returned.
  */
  Face_handle remove_edge (Halfedge_handle e,
                           bool remove_source = true,
                           bool remove_target = true);

  /// @}

  /// \name Miscellaneous
  /// @{

  /*!
    returns `true` if `arr` represents a valid instance of
    `Arrangement_2`. In particular, the functions checks the topological
    structure of the arrangement and assures that it is valid. In addition,
    the function performs several simple geometric tests to ensure the
    validity of some of the geometric properties of the arrangement. Namely,
    it checks that all arrangement vertices are associated with distinct
    points, and that the halfedges around every vertex are ordered clockwise.
  */
  bool is_valid () const;

  /// @}

}; /* end Arrangement_2 */
} /* end namespace CGAL */
namespace CGAL {

/*!
  \ingroup PkgArrangementOnSurface2Insert insert
  The function `%insert` inserts one or more curves or \f$ x\f$-monotone curves
  into a given arrangement, where no restrictions are imposed on the inserted
  curves. If an inserted curve is not \f$ x\f$-monotone curve, it is subdivided
  into \f$ x\f$-monotone subcurves (and perhaps isolated points), which are
  inserted into the arrangement.

\cgalHeading{Requirements}

  <UL>
  <LI>If the curve is \f$ x\f$-monotone curve then The instantiated
  `Traits` class must model the `ArrangementXMonotoneTraits_2`
  concept. In case that the curve is not \f$ x\f$-monotone then the
  instantiated `Traits` class must model the
  `ArrangementTraits_2` concept. That is, it should define the
  `Curve_2` type, and support its subdivision into \f$ x\f$-monotone
  subcurves (and perhaps isolated points).
  <LI>The point-location object `pl`, must model the
  `ArrangementPointLocation_2` concept.
  </UL>
*/
/// @{

/*!
Inserts the given curve `c` into the arrangement `arr`.
`c` is subdivided into \f$ x\f$-monotone subcurves (and perhaps isolated
points). Each subcurve is in turn inserted into the arrangement by locating
its left endpoint and computing its zone until reaching the right endpoint.

The given point-location object `pl` is used to locate the left
endpoints of the \f$ x\f$-monotone curves. By default, the function uses the
"walk along line" point-location strategy  -  namely an instance of
the class `Arr_walk_along_line_point_location<Arrangement_2<Traits,Dcel> >`.

\pre If provided, `pl` must be attached to the given arrangement `arr`.

*/
template<class Traits, class Dcel,
         class Curve, class PointLocation>
void insert (Arrangement_2<Traits,Dcel>& arr,
             const Curve& c,
             const PointLocation& pl = walk_pl);

/*!
Inserts the<I>\f$ x\f$-monotone (only)</I> curve `xc` into the
arrangement `arr`. The object `obj`, which either
wraps a `Vertex_const_handle`, a `Halfedge_const_handle`, or a
`Face_const_handle`, represents the location of `xc`'s left
endpoint in the arrangement. The zone of `xc` is computed starting
from the feature represented by `obj`. As in the case above, the
zone computation terminates, when the right endpoint is reached.
Thus, point-location is not required.
*/
template<typename Traits, typename Dcel>
void insert (Arrangement_2<Traits,Dcel>& arr,
             const typename Traits::X_monotone_curve_2& xc,
             const Object& obj);


/*!
Aggregately inserts the curves or \f$ x\f$-monotone curves in the range
`[first,last)` into the arrangement `arr` using the sweep-line
framework.
*/
template<class Traits, class Dcel, class InputIterator>
void insert (Arrangement_2<Traits,Dcel>& arr,
             InputIterator first, InputIterator last);

/// @}

/*!
  \ingroup PkgArrangementOnSurface2Funcs

  Checks if a given curve or \f$ x\f$-monotone
  curve intersects an existing arrangement's edges or vertices.

  If the give curve is not an \f$ x\f$-monotone curve then the function
  subdivides the given curve into \f$ x\f$-monotone subcurves and isolated
  vertices . Each subcurve is in turn checked for intersection.
  The function uses the zone algorithm to check if the curve intersects
  the arrangement. First, the curve's left endpoint is located. Then,
  its zone is computed starting from its left endpoint location. The
  zone computation terminates when an intersection with an arrangement's
  edge/vertex is found or when the right endpoint is reached.

  A given point-location object is used for locating the left endpoint
  of the given curve in the existing arrangement. By default, the function
  uses the "walk along line" point-location strategy - namely an
  instance of the class
  `Arr_walk_along_line_point_location<Arrangement_2<Traits,Dcel> >`.

  Checks if the given curve or \f$ x\f$-monotone curve `c` intersects
  edges or vertices of the existing arrangement `arr`.
  \pre If provided, `pl` must be attached to the given arrangement `arr`.

\cgalHeading{Requirements}

  <UL>
  <LI>If `c` is \f$ x\f$-monotone then the instantiated `GeomTraits`
  class must model the `ArrangementXMonotoneTraits_2` concept. If
  `c` is a curve then the instantiated `GeomTraits` class must
  model the `ArrangementTraits_2` concept. That is, it should
  define the `Curve_2` type, and support its subdivision into
  \f$ x\f$-monotone subcurves (and perhaps isolated points).
  <LI>The point-location object `pl`, must model the
  `ArrangementPointLocation_2` concept.
  </UL>

*/
template <class GeomTraits, class TopTraits,
          class Curve, class PointLocation>
bool do_intersect (
  Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
  const Curve& c, const PointLocation& pl);


/*!
  \ingroup PkgArrangementOnSurface2Funcs

  Inserts a given \f$ x\f$-monotone curve into a given
  arrangement, where the interior of the given curve is disjoint from all
  existing arrangement vertices and edges. Under this assumption, it is
  possible to locate the endpoints of the given curve in the arrangement,
  and use one of the specialized insertion member-functions of the
  arrangement according to the results. The insertion operations creates a
  single new edge, that is, two twin halfedges, and the function returns a
  handle for the one directed lexicographically in increasing order (from
  left to right).

  A given point-location object is used for answering the two point-location
  queries on the given curve endpoints. By default, the function uses the
  "walk along line" point-location strategy - namely, an instance of the
  class `Arr_walk_along_line_point_location<Arrangement_2<Traits,Dcel> >`.

  \pre If provided, `pl` must be attached to the given arrangement `arr`.

\cgalHeading{Requirements}

  <UL>
  <LI>The instantiated `Traits` class must model the restricted
  `ArrangementBasicTraits_2` concept, as no intersections are computed.
  <LI>The point-location object `pl` must model the
  `ArrangementPointLocation_2` concept.
  </UL>

*/
template<typename Traits, typename Dcel,
         typename PointLocation>
typename Arrangement_2<Traits,Dcel>::Halfedge_handle
insert_non_intersecting_curve (Arrangement_2<Traits,Dcel>& arr,
                               const typename Traits::X_monotone_curve_2& xc,
                               const PointLocation& pl = walk_pl);


/*!
  \ingroup PkgArrangementOnSurface2Funcs

  Inserts a set of \f$ x\f$-monotone curves in a given
  range into a given arrangement. The insertion is performed in an aggregated
  manner, using the sweep-line algorithm. The input curves should be pairwise
  disjoint in their interior and pairwise interior-disjoint from all existing
  arrangement vertices and edges.

\cgalHeading{Requirements}

  <UL>
  <LI>The instantiated `Traits` class must model the
  `ArrangementBasicTraits_2` concept, as no intersections are computed.
  <LI>`InputIterator::value_type` must be `Traits::X_monotone_curve_2`
  </UL>

*/
template<typename Traits, typename Dcel, InputIterator>
void insert_non_intersecting_curves(Arrangement_2<Traits,Dcel>& arr,
                                    InputIterator first, InputIterator last);


/*!
  \ingroup PkgArrangementOnSurface2Funcs

  Inserts a given point into a given arrangement.
  It uses a given point-location object to locate the given
  point in the given arrangement. If the point conincides with an existing
  vertex, there is nothing left to do; if it lies on an edge, the edge is
  split at the point. Otherwise, the point is contained inside a face, and is
  inserted as an isolated vertex inside this face.
  By default, the function uses the "walk along line" point-location
  strategy - namely, an instance of the class
  `Arr_walk_along_line_point_location<Arrangement_2<Traits,Dcel> >`.
  In either case, the function returns a handle for the vertex associated
  with the point.

  \pre If provided, `pl` must be attached to the given arrangement `arr`.

\cgalHeading{Requirements}

  <UL>
  <LI>The instantiated `Traits` class must model the
  `ArrangementXMonotoneTraits_2` concept. Not all expressions listed
  by this concept are required. In fact the traits class must model the
  `ArrangementBasicTraits_2` concept, and support the splitting
  functionality.
  <LI>The point-location object `pl`, must model the
  `ArrangementPointLocation_2` concept.
  </UL>

*/
template<typename Traits, typename Dcel,
         typename PointLocation>
typename Arrangement_2<Traits,Dcel>::Vertex_handle
insert_point (Arrangement_2<Traits,Dcel>& arr,
              const typename Traits::Point_2& p,
              const PointLocation& pl = walk_pl);


/*!
  \ingroup PkgArrangementOnSurface2Funcs

  Checks the validity of a given arrangement.

  Invokes the member function `arr.is_valid()` to verify the
  topological correctness of the arrangement. Then it performs additional
  validity tests. It checks that all \f$ x\f$-monotone curves associated with
  arrangement edges are pairwise disjoint in their interior. Then it makes
  sure that all holes and all isolated vertices are located within the
  proper arrangement faces. Note that the test carried out by this
  function may take a considerable amount of time; it is recommended to be
  used only for debugging purposes.

\cgalHeading{Requirements}

  The instantiated traits class must model the concept
  `ArranagmentXMonotoneTraits_2`.

*/
template<typename Traits, typename Dcel>
bool is_valid (const Arrangement_2<Traits, Dcel>& arr);


/*!
  \ingroup PkgArrangementOnSurface2Funcs

  Removes an edge given by one of the twin halfedges
  that forms it, from a given arrangement. Once the edge is removed, if the
  vertices associated with its endpoints become isolated, they are removed as
  well. The call `remove_edge(arr, e)` is equivalent to the call
  `arr.remove_edge (e, true, true)`. However, this free function requires
  that `Traits` be a model of the refined concept
  `ArrangementXMonotoneTraits_2`, which requires merge operations
  on \f$ x\f$-monotone curves. If one of the end-vertices of the given edge
  becomes redundant after the edge is removed (see `remove_vertex()`
  for the definition of a redundant vertex), it is removed, and its
  incident edges are merged.
  If the edge-removal operation causes two faces to merge, the merged face
  is returned. Otherwise, the face to which the edge was incident before the
  removal is returned.

\cgalHeading{Requirements}

  <UL>
  <LI>The instantiated traits class must model the concept
  `ArrangementXMonotoneTraits_2`.
  </UL>

*/
template<typename Traits, typename Dcel>
typename Arrangement_2<Traits,Dcel>::Face_handle
remove_edge (Arrangement_2<Traits,Dcel>& arr,
             typename Arrangement_2<Traits,Dcel>::Halfedge_handle e);


/*!
  \ingroup PkgArrangementOnSurface2Funcs

  Attempts to removed a given vertex from a given
  arrangement. The vertex can be removed if it is either an isolated vertex,
  (and has no incident edge,) or if it is a <I>redundant</I> vertex. That
  is, it has exactly two incident edges, whose associated curves can be
  merged to form a single \f$ x\f$-monotone curve.
  The function returns a boolean value that indicates whether it succeeded
  removing the vertex from the arrangement.

\cgalHeading{Requirements}

  <UL>
  <LI>The instantiated `Traits` class must model the
  `ArrangementXMonotoneTraits_2` concept. Not all expressions listed
  by this concept are required. In fact the traits class must model the
  `ArrangementBasicTraits_2` concept and support the merging
  functionality.
  </UL>
*/
template <typename Traits, typename Dcel>
bool remove_vertex (Arrangement_2<Traits,Dcel>& arr,
                    typename Arrangement_2<Traits,Dcel>::Vertex_handle v);

/*!
  \ingroup PkgArrangementOnSurface2Funcs

  Compute the zone of the given \f$ x\f$-monotone
  curve in the existing arrangement. Meaning, it output the
  arrangement's vertices, edges and faces that the \f$ x\f$-monotone curve
  intersects. The order of the objects is the order that they are
  discovered when traversing the \f$ x\f$-monotone curve from left to right.

  A given point-location object is used for answering point-location queries
  during the insertion process. By default, the function uses the
  "walk along line" point-location strategy - namely an instance of the
  class `Arr_walk_along_line_point_location<Arrangement_2<Traits,Dcel> >`.

  Compute the zone of the given \f$ x\f$-monotone curve `c` in the
  arrangement `arr`.
  \pre If provided, `pl` must be attached to the given arrangement `arr`.

\cgalHeading{Requirements}

  <UL>
  <LI>The instantiated `GeomTraits` class must model the
  `ArrangementXMonotoneTraits_2` concept.
  <LI>The point-location object `pl`, must model the
  `ArrangementPointLocation_2` concept.
  </UL>

*/
template <class GeomTraits, class TopTraits,
          class OutputIterator, class PointLocation>
OutputIterator zone (
  Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
  const typename GeomTraits::X_monotone_curve_2& c,
  OutputIterator oi,
  const PointLocation& pl);

} /* namespace CGAL */
