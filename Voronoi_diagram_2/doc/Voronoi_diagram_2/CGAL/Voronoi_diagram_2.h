
namespace CGAL {

/*!
\ingroup PkgVoronoiDiagram2Ref

The class `Voronoi_diagram_2` provides an adaptor that enables us
to view a triangulated Delaunay graph as their dual subdivision, the
Voronoi diagram. The class `Voronoi_diagram_2` is designed to provide an API
that is similar to that of \cgal's arrangements.

\tparam DG corresponds to the triangulated Delaunay graph and must be a model of the `DelaunayGraph_2` concept.

\tparam AT must be a model of the `AdaptationTraits_2` concept.

\tparam AP must be a model of the `AdaptationPolicy_2` concept.  The third template parameter defaults to `CGAL::Identity_policy_2<DG,AT>`.

\cgalRefines{CopyConstructible,Assignable,DefaultConstructible}

\sa `CGAL::Delaunay_triangulation_2<Traits,Tds>`
\sa `CGAL::Regular_triangulation_2<Traits,Tds>`
\sa `CGAL::Triangulation_hierarchy_2<Tr>` provided that `Tr` is a  model of `DelaunayGraph_2`
\sa `CGAL::Segment_Delaunay_graph_2<Gt,DS>`
\sa `CGAL::Segment_Delaunay_graph_hierarchy_2<Gt,STag,DS>`
\sa `CGAL::Apollonius_graph_2<Gt,Agds>`
\sa `CGAL::Apollonius_graph_hierarchy_2<Gt,Agds>`

*/
template< typename DG, typename AT, typename AP >
class Voronoi_diagram_2 {
public:



  /*!
    \ingroup PkgVoronoiDiagram2Ref

    The class `Face` is the class provided by the
    `Voronoi_diagram_2<DG,AT,AP>` class for Voronoi faces. Below we
    present its interface.

    \cgalModels{DefaultConstructible,CopyConstructible,Assignable,EqualityComparable,LessThanComparable}

    \sa `CGAL::Voronoi_diagram_2<DG,AT,AP>`
    \sa \link CGAL::Voronoi_diagram_2::Vertex `CGAL::Voronoi_diagram_2<DG,AT,AP>::Vertex` \endlink
    \sa \link CGAL::Voronoi_diagram_2::Halfedge `CGAL::Voronoi_diagram_2<DG,AT,AP>::Halfedge` \endlink
    \sa `DelaunayGraph_2`

  */

  class Face {
  public:

    /// \name Types
    /// @{

    /*!
      A type for the vertices of the Voronoi diagram.
    */
    typedef unspecified_type Vertex;

    /*!
      A type for the halfedges of the Voronoi diagram.
    */
    typedef unspecified_type Halfedge;

    /*!
      Handle for the vertices of the Voronoi diagram.
    */
    typedef unspecified_type Vertex_handle;

    /*!
      Handle for the faces of the Voronoi diagram.
    */
    typedef unspecified_type Face_handle;

    /*!
      Handle for the halfedges of the Voronoi
      diagram.
    */
    typedef unspecified_type Halfedge_handle;

    /*!
      A type for a bidirectional
      circulator over the halfedges on the boundary of the face. The value
      type of the circulator is
      \ref CGAL::Voronoi_diagram_2::Halfedge "CGAL::Voronoi_diagram_2<DG,AT,AP>::Halfedge", and is convertible
      to `Halfedge_handle`.
    */
    typedef unspecified_type Ccb_halfedge_circulator;

    /*!
      A type for the Delaunay graph. It is a
      model of the `DelaunayGraph_2` concept.
    */
    typedef unspecified_type Delaunay_graph;

    /*!
      A type for the handle of the dual vertex.
    */
    typedef Delaunay_graph::Vertex_handle
    Delaunay_vertex_handle;

    /// @}

    /// \name Access Methods
    /// @{

    /*!
      Returns an incident halfedge on
      the boundary of `f`.
    */
    Halfedge_handle halfedge();

    /*!
      Returns a bidirectional
      circulator for traversing the halfedges on the boundary of
      `f`. The halfedges are traversed in counterclockwise order.
    */
    Ccb_halfedge_circulator ccb();

    /*!
      Returns a handle to the corresponding dual vertex in the Delaunay graph.
    */
    Delaunay_vertex_handle dual();

    /// @}

    /// \name Predicate Methods
    /// @{

    /*!
      Returns `true` iff the face is
      an unbounded face in the Voronoi diagram.
    */
    bool is_unbounded();

    /*!
      Returns
      `true` iff `e` is a halfedge of the boundary of
      `f`.
    */
    bool is_halfedge_on_ccb(Halfedge e);

    /*!
      Returns `true` iff the following
      conditions are met: the face is not rejected by the chosen
      adaptation policy;
      all its adjacent halfedges do not have zero length; all its adjacent
      halfedges return the face as their adjacent face.
    */
    bool is_valid();

    /// @}

  }; /* end Face */


  /*!
    \ingroup PkgVoronoiDiagram2Ref

    The class `Halfedge` is the class provided by the
    `Voronoi_diagram_2<DG,AT,AP>` class for Voronoi halfedges.
    Below we present its interface.

    \cgalModels{DefaultConstructible,CopyConstructible,Assignable,EqualityComparable,LessThanComparable}

    \sa `CGAL::Voronoi_diagram_2<DG,AT,AP>`
    \sa \link CGAL::Voronoi_diagram_2::Vertex `CGAL::Voronoi_diagram_2<DG,AT,AP>::Vertex` \endlink
    \sa \link CGAL::Voronoi_diagram_2::Face `CGAL::Voronoi_diagram_2<DG,AT,AP>::Face` \endlink
    \sa `DelaunayGraph_2`

  */

  class Halfedge {
  public:

    /// \name Types
    /// @{

    /*!
      A type for the vertices of the Voronoi diagram.
    */
    typedef unspecified_type Vertex;

    /*!
      A type for the faces of the Voronoi diagram.
    */
    typedef unspecified_type Face;

    /*!
      Handle for the vertices of the Voronoi diagram.
    */
    typedef unspecified_type Vertex_handle;

    /*!
      Handle for the faces of the Voronoi diagram.
    */
    typedef unspecified_type Face_handle;

    /*!
      Handle for the halfedges of the Voronoi
      diagram.
    */
    typedef unspecified_type Halfedge_handle;

    /*!
      A type for a bidirectional
      circulator over the halfedges of the boundary of a
      Voronoi face. The value type of the circulator is
      \ref CGAL::Voronoi_diagram_2::Halfedge "CGAL::Voronoi_diagram_2<DG,AT,AP>::Halfedge" and is convertible to
      `Halfedge_handle`.
    */
    typedef unspecified_type Ccb_halfedge_circulator;

    /*!
      A type for the Delaunay graph. It is a
      model of the `DelaunayGraph_2` concept.
    */
    typedef unspecified_type Delaunay_graph;

    /*!
      A type for
      the dual edge in the Delaunay graph.
    */
    typedef Delaunay_graph::Edge Delaunay_edge;

    /*!
      A type for vertex handles in the Delaunay graph.
    */
    typedef Delaunay_graph::Vertex_handle Delaunay_vertex_handle;

    /// @}

    /// \name Access Methods
    /// @{

    /*!
      Returns the twin halfedge.
    */
    Halfedge_handle twin();

    /*!
      Same as `e.twin()`.
    */
    Halfedge_handle opposite();

    /*!
      Returns the next halfedge in the
      counterclockwise sense around the boundary of the face that `e`
      is incident to.
    */
    Halfedge_handle next();

    /*!
      Returns the previous halfedge in the
      counterclockwise sense around the boundary of the adjacent face.
    */
    Halfedge_handle previous();

    /*!
      Returns the face that `e` is
      incident to.
    */
    Face_handle face();

    /*!
      Returns the source vertex of
      `e`.
      \pre The source vertex must exist, i.e., `has_source()` must return `true`.
    */
    Vertex_handle source();

    /*!
      Returns the target vertex of
      `e`.
      \pre The target vertex must exist, i.e., `has_target()` must return `true`.
    */
    Vertex_handle target();

    /*!
      Returns a bidirectional
      circulator to traverse the halfedges on the boundary of the Voronoi
      face containing `e`. The circulator is initialized to
      `e`. Applying `operator++` (resp. `operator-`) to this
      circulator returns the next halfedge on the boundary of the face
      containing `e` in the counterclockwise (resp. clockwise) sense.
    */
    Ccb_halfedge_circulator ccb();

    /*!
      Returns the
      corresponding dual edge in the Delaunay graph.
    */
    Delaunay_edge dual();

    /// @}

    /// \name
    /// In the four methods below we consider Voronoi halfedges to be
    /// "parallel" to the \f$ x\f$-axis, oriented from left to right.
    /// @{

    /*!
      Returns a handle to the vertex in
      the Delaunay graph corresponding to the defining site above
      the Voronoi edge.
    */
    Delaunay_vertex_handle up();

    /*!
      Returns a handle to the vertex
      in the Delaunay graph corresponding to the defining site below
      the Voronoi edge.
    */
    Delaunay_vertex_handle down();

    /*!
      Returns a handle to the vertex in
      the Delaunay graph corresponding to the defining site to the left of
      the Voronoi edge.
      \pre `has_source()` must be `true`.
    */
    Delaunay_vertex_handle left();

    /*!
      Returns a handle to the vertex in
      the Delaunay graph corresponding to the defining site to the right of
      the Voronoi edge.
      \pre `has_target()` must be `true`.
    */
    Delaunay_vertex_handle right();

    /// @}

    /// \name Predicate Methods
    /// @{

    /*!
      Returns `true` iff the halfedge
      corresponds to a bisecting segment or a bisecting ray oriented
      appropriately so that its apex is its source.
    */
    bool has_source();

    /*!
      Returns `true` iff the halfedge
      corresponds to a bisecting segment or a bisecting ray oriented
      appropriately so that its apex is its target.
    */
    bool has_target();

    /*!
      Returns `true` iff the source or
      the target of the halfedge does not exist, i.e., if either of
      `has_source()` or `has_target()` return `false`.
    */
    bool is_unbounded();

    /*!
      Returns `true` iff the Voronoi
      edge is an entire bisector.
    */
    bool is_bisector();

    /*!
      Returns `true` iff the Voronoi
      edge has both a source and a target Voronoi vertex.
    */
    bool is_segment();

    /*!
      Returns `true` iff the Voronoi
      edge has either a source or a target Voronoi vertex, but not both;
      in other words it is a bisecting ray.
    */
    bool is_ray();

    /*!
      Returns `true` if the following
      conditions are met: the halfedge is not a rejected
      edge with respect to the chosen adaptation policy;
      the twin edge of its twin edge is itself; its adjacent face is not a
      rejected face with respect to the chosen adaptation policy;
      its source and target vertices are valid (provided
      they exist, of course); the previous of its next halfedge is itself
      and the next of its previous halfedge is itself.
    */
    bool is_valid();

    /// @}

  }; /* end Halfedge */


  /*!
    \ingroup PkgVoronoiDiagram2Ref

    The class `Vertex` is the Voronoi vertex class provided by the class
    `Voronoi_diagram_2<DG,AT,AP>` class. Below we present its interface.

    \cgalModels{DefaultConstructible,CopyConstructible,Assignable,,EqualityComparable,LessThanComparable}

    \sa `CGAL::Voronoi_diagram_2<DG,AT,AP>`
    \sa \link CGAL::Voronoi_diagram_2::Halfedge `CGAL::Voronoi_diagram_2<DG,AT,AP>::Halfedge` \endlink
    \sa \link CGAL::Voronoi_diagram_2::Face `CGAL::Voronoi_diagram_2<DG,AT,AP>::Face` \endlink
    \sa `DelaunayGraph_2`

  */

  class Vertex {
  public:

    /// \name Types
    /// @{

    /*!
      A type for the halfedges of the Voronoi diagram.
    */
    typedef unspecified_type Halfedge;

    /*!
      A type for the faces of the Voronoi diagram.
    */
    typedef unspecified_type Face;

    /*!
      Handle for the vertices of the Voronoi diagram.
    */
    typedef unspecified_type Vertex_handle;

    /*!
      Handle for the faces of the Voronoi diagram.
    */
    typedef unspecified_type Face_handle;

    /*!
      Handle for the halfedges of the Voronoi
      diagram.
    */
    typedef unspecified_type Halfedge_handle;

    /*!
      A type for the point represented by the
      vertex.
    */
    typedef unspecified_type Point_2;

    /*!
      A type for sizes.
    */
    typedef unspecified_type size_type;

    /*!
      A type for a bidirectional
      circulator that allows to traverse all incident halfedges, i.e., all
      halfedges that have the vertex as their target. The value
      type of the circulator is
      \ref CGAL::Voronoi_diagram_2::Halfedge "CGAL::Voronoi_diagram_2<DG,AT,AP>::Halfedge" and is convertible to
      `Halfedge_handle`.
    */
    typedef unspecified_type Halfedge_around_vertex_circulator;

    /*!
      A type for the Delaunay graph. It is a
      model of the `DelaunayGraph_2` concept.
    */
    typedef unspecified_type Delaunay_graph;

    /*!
      A type for the handle of the dual face.
    */
    typedef Delaunay_graph::Face_handle Delaunay_face_handle;

    /*!
      A type for the vertex handles in the Delaunay graph.
    */
    typedef Delaunay_graph::Vertex_handle Delaunay_vertex_handle;

    /// @}

    /// \name Access Methods
    /// @{

    /*!
      Returns an incident halfedge
      that has `v` as its target.
    */
    Halfedge_handle halfedge();

    /*!
      Returns the in-degree of the vertex,
      i.e.\ the number of halfedges that have `v` as their target.
    */
    size_type degree();

    /*!
      Returns the point represented by the vertex.
    */
    Point_2 point();

    /*!
      Returns a handle to the corresponding dual face in the
      Delaunay graph.
    */
    Delaunay_face_handle dual();

    /*!
      Returns a handle to the vertex in the Delaunay graph corresponding to
      the \f$ (i+1)\f$-th generating site of the Voronoi vertex.
      \pre `i` must be smaller or equal to 2.
    */
    Delaunay_vertex_handle site(unsigned int i);

    /*!
      Returns a bidirectional circulator that allows the traversal of the
      halfedges that have `v` as their target. Applying
      `operator++` (resp. `operator-`) to this circulator returns
      the next incident halfedge in the counterclockwise (resp. clockwise) sense.
    */
    Halfedge_around_vertex_circulator incident_halfedges();

    /// @}

    /// \name Predicate Methods
    /// @{

    /*!
      Returns `true`
      if the halfedge `e` is incident to `v`.
    */
    bool is_incident_edge(Halfedge_handle e);

    /*!
      Returns `true`
      if the face `f` is incident to `v`.
    */
    bool is_incident_face(Face_handle e);

    /*!
      Returns `true` if the following
      conditions are met: the dual face is not an infinite face; all
      incident halfedges have the vertex as their target.
    */
    bool is_valid();

    /// @}

  }; /* end Vertex */


/// \name Types
/// @{

/*!
A type for the dual Delaunay graph.
*/
typedef DG Delaunay_graph;

/*!
A type for the adaptation
traits needed by the Voronoi diagram adaptor.
*/
typedef AT Adaptation_traits;

/*!
A type for the adaptation
policy used.
*/
typedef AP Adaptation_policy;

/*!
A type a point.
*/
typedef Adaptation_traits::Point_2 Point_2;

/*!
A type
for the sites of the Voronoi diagram.
*/
typedef Adaptation_traits::Site_2 Site_2;

/*!
A type for sizes.
*/
typedef Delaunay_graph::size_type size_type;

/*!
A type for the geometric traits of the Delaunay graph.
*/
typedef Delaunay_graph::Geom_traits Delaunay_geom_traits;

/*!
A type for the vertex handles of the Delaunay graph.
*/
typedef Delaunay_graph::Vertex_handle Delaunay_vertex_handle;

/*!
A type for the face handles of the Delaunay graph.
*/
typedef Delaunay_graph::Face_handle Delaunay_face_handle;

/*!
A type for the edges of the Delaunay graph.
*/
typedef Delaunay_graph::Edge Delaunay_edge;

/// @}

/// \name
/// The vertices, edges and faces of the Voronoi diagram are accessed
/// through `handles`, `iterators` and `circulators`. The iterators
/// and circulators are all bidirectional and non-mutable. The
/// circulators and iterators are assignable to the corresponding
/// handle types, and they are also convertible to the corresponding
/// handles.
/// @{

/*!
Handle for halfedges.
*/
typedef unspecified_type Halfedge_handle;

/*!
Handle for vertices.
*/
typedef unspecified_type Vertex_handle;

/*!
Handle for faces.
*/
typedef unspecified_type Face_handle;

/*!
A type for an iterator over Voronoi
edges. Edges are considered non-oriented. Its value type is
`Halfedge`.
*/
typedef unspecified_type Edge_iterator;

/*!
A type for an iterator over Voronoi
halfedges. Halfedges are oriented and come in pairs. Its value type
is `Halfedge`.
*/
typedef unspecified_type Halfedge_iterator;

/*!
A type for an iterator over Voronoi
faces. Its value type is `Face`.
*/
typedef unspecified_type Face_iterator;

/*!
A type for an iterator over Voronoi
vertices. Its value type is `Vertex`.
*/
typedef unspecified_type Vertex_iterator;

/*!
A type for a
circulator over the halfedges that have a common vertex as their
target. Its value type is `Halfedge`.
*/
typedef unspecified_type Halfedge_around_vertex_circulator;

/*!
A type for a circulator over
the halfedges on the boundary of a Voronoi face. Its value type of
is `Halfedge`.
*/
typedef unspecified_type Ccb_halfedge_circulator;

/*!
A type for an iterator over
the unbounded faces of the Voronoi diagram. Its value type is
`Face`.
*/
typedef unspecified_type Unbounded_faces_iterator;

/*!
A type for an iterator over
the bounded faces of the Voronoi diagram. Its value type is
`Face`.
*/
typedef unspecified_type Bounded_faces_iterator;

/*!
A type for an iterator over
the unbounded halfedges of the Voronoi diagram. Its value type is
`Halfedge`.
*/
typedef unspecified_type Unbounded_halfedges_iterator;

/*!
A type for an iterator over
the bounded halfedges of the Voronoi diagram. Its value type is
`Halfedge`.
*/
typedef unspecified_type Bounded_halfedges_iterator;

/*!
A type for an iterator over the
sites of the Voronoi diagram. Its value type is `Site_2`.
*/
typedef unspecified_type Site_iterator;

/*!
The result type of the point location queries.
*/
typedef std::variant<Face_handle,Halfedge_handle,Vertex_handle>
Locate_result;

/// @}

/// \name Creation
/// @{

/*!
Creates a Voronoi diagram using `at` as adaptation traits and
`ap` as adaptation policy; the underlying Delaunay graph is
created using `gt` as geometric traits.
*/
Voronoi_diagram_2(Adaptation_traits
at = Adaptation_traits(), Adaptation_policy ap = Adaptation_policy(),
Delaunay_geom_traits gt = Delaunay_geom_traits());

/*!
Creates a Voronoi diagram from the Delaunay graph `dg` and using
`at` as adaptation traits and `ap` as adaptation policy. The
Delaunay graph `dg` is fully copied if `swap_dg` is set to
`false`, or swapped with the one stored internally if
`swap_dg` is set to `true`.
*/
Voronoi_diagram_2(Delaunay_graph dg, bool swap_dg = false,
Adaptation_traits at = Adaptation_traits(), Adaptation_policy ap =
Adaptation_policy());

/*!
Creates a Voronoi diagram using as sites the sites in the iterator
range `[first, beyond)`, `at` as adaptation traits and
`ap` as adaptation policy; the underlying Delaunay graph is
created using `gt` as geometric traits. `Iterator` must be a
model of the `InputIterator` concept and its value type must be
`Site_2`.
*/
template<class Iterator>
Voronoi_diagram_2(Iterator first, Iterator beyond, Adaptation_traits
at = Adaptation_traits(), Adaptation_policy ap = Adaptation_policy(),
Delaunay_geom_traits gt = Delaunay_geom_traits());

/// @}

/// \name Access Methods
/// @{

/*!
Returns a const reference to the dual graph, i.e., the Delaunay graph.
*/
Delaunay_graph dual();

/*!
Returns a handle to the halfedge in the Voronoi diagram that is dual
to the edge `e` in the Delaunay graph.
*/
Halfedge_handle dual(Delaunay_edge e);

/*!
Returns a handle to the face in the Voronoi diagram that is dual to
the vertex corresponding to the vertex handle `v` in the
Delaunay graph.
*/
Face_handle dual(Delaunay_vertex_handle v);

/*!
Returns a handle to the vertex in the Voronoi diagram that is dual to
the face corresponding to the face handle `f` in the Delaunay graph.
*/
Vertex_handle dual(Delaunay_face_handle f);

/*!
Returns a reference to the Voronoi traits.
*/
Adaptation_traits adaptation_traits();

/*!
Returns a reference to the adaptation policy.
*/
Adaptation_policy adaptation_policy();

/*!
Returns the number of Voronoi vertices.
*/
size_type number_of_vertices();

/*!
Returns the number of Voronoi faces (bounded and unbounded).
*/
size_type number_of_faces();

/*!
Returns the number of halfedges (bounded and unbounded) in the
Voronoi diagram. This is always an even number.
*/
size_type number_of_halfedges();

/*!
Returns the number of connected components of the Voronoi skeleton.
*/
size_type number_of_connected_components();

/*!
Returns one of the unbounded
faces of the Voronoi diagram. If no unbounded faces exist (this can
happen if the number of sites is zero) the
default constructed face handle is returned.
*/
Face_handle unbounded_face();

/*!
Returns one of the bounded
faces of the Voronoi diagram. If no bounded faces exist the default
constructed face handle is returned.
*/
Face_handle bounded_face();

/*!
Returns one of the unbounded
halfedges of the Voronoi diagram. If no unbounded halfedges exist the
default constructed halfedge handle is returned.
*/
Halfedge_handle unbounded_halfedge();

/*!
Returns one of the bounded
halfedges of the Voronoi diagram. If no bounded halfedges exist the
default constructed halfedge handle is returned.
*/
Halfedge_handle bounded_halfedge();

/// @}

/// \name Traversal of the Voronoi Diagram
///
/// A Voronoi diagram can be seen as a container of faces, vertices and
/// halfedges. Therefore the Voronoi diagram provides several iterators
/// and circulators that allow to traverse it.

/// \name Iterators
/// The following iterators allow respectively to visit the faces (all
/// or only the unbounded/bounded ones), edges, halfedges (all or only
/// the unbounded/bounded ones) and vertices of the Voronoi
/// diagram. These iterators are non-mutable, bidirectional and their
/// value types are respectively `Face`, `Halfedge`, `Halfedge` and
/// `Vertex`. All iterators are convertible to the corresponding
/// handles and are invalidated by any change in the Voronoi
/// diagram. The following iterator provides access to the sites that
/// define the Voronoi diagram. Its value type is `Site_2`. It is
/// invalidated by any change in the Voronoi diagram.
/// @{

/*!
Starts at an arbitrary Voronoi face.
*/
Face_iterator faces_begin();

/*!
Past-the-end iterator.
*/
Face_iterator faces_end();

/*!
Starts at an arbitrary unbounded Voronoi face.
*/
Unbounded_faces_iterator unbounded_faces_begin();

/*!
Past-the-end iterator.
*/
Unbounded_faces_iterator unbounded_faces_end();

/*!
Starts at an arbitrary bounded Voronoi face.
*/
Bounded_faces_iterator bounded_faces_begin();

/*!
Past-the-end iterator.
*/
Bounded_faces_iterator bounded_faces_end();

/*!
Starts at an arbitrary Voronoi edge.
*/
Edge_iterator edges_begin();

/*!
Past-the-end iterator.
*/
Edge_iterator edges_end();

/*!
Starts at an arbitrary Voronoi halfedge.
*/
Halfedge_iterator halfedges_begin();

/*!
Past-the-end iterator.
*/
Halfedge_iterator halfedges_end();

/*!
Starts at an arbitrary unbounded Voronoi edge.
*/
Unbounded_halfedges_iterator unbounded_halfedges_begin();

/*!
Past-the-end iterator.
*/
Unbounded_halfedges_iterator unbounded_halfedges_end();

/*!
Starts at an arbitrary bounded Voronoi edge.
*/
Bounded_halfedges_iterator bounded_halfedges_begin();

/*!
Past-the-end iterator.
*/
Bounded_halfedges_iterator bounded_halfedges_end();

/*!
Starts at an arbitrary Voronoi vertex.
*/
Vertex_iterator vertices_begin();

/*!
Past-the-end iterator.
*/
Vertex_iterator vertices_end();

/*!
Starts at an arbitrary site.
*/
Site_iterator sites_begin();

/*!
Past-the-end iterator.
*/
Site_iterator sites_end();

/// @}

/// \name Circulators
/// The Voronoi diagram adaptor also provides circulators that allow
/// to visit all halfedges whose target is a given vertex - this is
/// the `Halfedge_around_vertex_circulator`, as well as all halfedges
/// on the boundary of a Voronoi face - this is the
/// `Ccb_halfedge_circulator`. These circulators are non-mutable and
/// bidirectional. The operator `operator++` moves the former
/// circulator counterclockwise around the vertex while the
/// `operator-` moves clockwise. The latter circulator is moved by the
/// operator `operator++` to the next halfedge on the boundary in the
/// counterclockwise sense, while `operator-` moves clockwise. When
/// the `Ccb_halfedge_circulator` is defined over an infinite Voronoi
/// face `f`, then applying `operator++` to a circulator corresponding
/// to a halfedge whose target is not finite moves to the next
/// infinite (or semi-infinite) halfedge of `f` in the
/// counterclockwise sense. Similarly, applying `operator++` to a
/// circulator corresponding to a halfedge whose source is not finite,
/// moves to the previous infinite (or semi-infinite) halfedge of `f`
/// in the clockwise sense. The `Halfedge_around_vertex_circulator`
/// circulator is invalidated by any modification in the faces
/// adjacent to the vertex over which it is defined. The
/// `Ccb_halfedge_circulator` is invalidated by any modification in
/// the face over which it is defined.
/// @{

/*!
Returns a circulator over the halfedges on the boundary of `f`.
The circulator is initialized to an arbitrary halfedge on the
boundary of the Voronoi face `f`.
*/
Ccb_halfedge_circulator ccb_halfedges(Face_handle f);

/*!
Returns a circulator over the halfedges on the boundary of
`f`. The circulator is initialized with the halfedge `h`.
\pre The halfedge `h` must lie on the boundary of `f`.
*/
Ccb_halfedge_circulator ccb_halfedges(Face_handle f,
Halfedge_handle h);

/*!
Returns a circulator over the halfedges whose target is the Voronoi
vertex `v`. The circulator is initialized to an arbitrary
halfedge incident to `v`.
*/
Halfedge_around_vertex_circulator
incident_halfedges(Vertex_handle v);

/*!
Returns a circulator over the halfedges whose target is the Voronoi
vertex `v`. The circulator is initialized with the halfedge
`h`.
\pre The vertex `v` must be the target vertex of the halfedge `h`.
*/
Halfedge_around_vertex_circulator
incident_halfedges(Vertex_handle v, Halfedge_handle h);

/// @}

/// \name Insertion
/// @{

/*!
Inserts the site
`t` in the Voronoi diagram. A handle to the face corresponding
to the Voronoi face of `t` in the Voronoi diagram is
returned. If `t` has an empty Voronoi cell, the default
constructed face handle is returned. This method is supported only
if `Voronoi_traits::Has_inserter` is set to
`CGAL::Tag_true`.
*/
Face_handle insert(Site_2 t);

/*!
Inserts, in the
Voronoi diagram, the sites in the iterator range `[first, beyond)`. The value type of `Iterator` must be
`Site_2`. The number of sites in the iterator range is
returned. This method is supported only if
`Voronoi_traits::Has_inserter` is set to `CGAL::Tag_true`.
*/
template<class Iterator>
size_type insert(Iterator first, Iterator beyond);

/// @}

/// \name Queries
/// @{

/*!
Performs point location for
the query point `p`. In other words, the face, halfedge or
vertex of the Voronoi diagram is found on which the point `p`
lies. This method is supported only if
`Voronoi_traits::Has_nearest_site_2` is set to
`CGAL::Tag_true`.
\pre The Voronoi diagram must contain at least one face.
*/
Locate_result locate(Point_2 p);

/// @}

/// \name I/O
/// @{

/*!
Writes the current state of the Voronoi diagram to the output
stream `os`.

The following operator must be defined:

`std::ostream& operator<<(std::ostream&, Delaunay_graph)`

*/
void file_output(std::ostream& os);

/*!
Reads the current state of the Voronoi diagram from the input
stream `is`.

The following operator must be defined:

`std::istream& operator>>(std::istream&, Delaunay_graph)`

*/
void file_input(std::istream& is);

/*!
Writes the current state of the Voronoi diagram to the output
stream `os`.

The following operator must be defined:

`std::ostream& operator<<(std::ostream&, Delaunay_graph)`
*/
std::ostream& operator<<(std::ostream& os, Voronoi_diagram_2<DG,AT,AP> vd);

/*!
Reads the current state of the Voronoi diagram from the input
stream `is`.

The following operator must be defined:
`std::istream& operator>>(std::istream&, Delaunay_graph)`
*/
std::istream& operator>>(std::istream& is, Voronoi_diagram_2<DG,AT,AP> vd);

/// @}

/// \name Validity Check
/// @{

/*!
Checks the validity of the dual Delaunay
graph and the Voronoi diagram adaptor.
*/
bool is_valid();

/// @}

/// \name Miscellaneous
/// @{

/*!
Clears all contents of the Voronoi diagram.
*/
void clear();

/*!
The Voronoi
diagrams `other` and `vd` are
swapped. `vd`.`swap(other)` should be preferred to
`vd``= other` or to
`vd``(other)` if `other` is deleted afterwards.
*/
void swap(Voronoi_diagram_2<DG,AT,AP> other);

/// @}

}; /* end Voronoi_diagram_2 */
} /* end namespace CGAL */
