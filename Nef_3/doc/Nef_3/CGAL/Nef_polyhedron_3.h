namespace CGAL {

/*!
  \ingroup PkgNef3Ref

  A 3D Nef polyhedron is a subset of the 3-dimensional space that is the
  result of forming complements and intersections starting from a finite
  set `H` of 3-dimensional halfspaces. Nef polyhedra are closed
  under all binary set operations, i.e., intersection, union,
  difference, complement, and under the topological operations boundary,
  closure, and interior.

  A 3D Nef polyhedron can be represented by the local pyramids of the minimal
  elements of its incidence structure. Without going into to much detail, a local
  pyramid essentially reflects the topologic and geometric situation at a certain
  location in a point set. For finite polyhedra the minimal elements
  of the incidence structure are vertices only. This means, that it suffices to
  model the topological and geometric situation of the vertices. For
  3D Nef polyhedra, the local pyramid of a vertex is represented by
  a planar Nef polyhedra embedded on a sphere.

  A `Nef_polyhedron_3` consists of vertices <I>V</I>, a sphere map for each
  vertex in <I>V</I>, edges <I>E</I>, facets <I>F</I>, volumes <I>C</I>, a mark
  for every item, and an incidence relation on them. Each edge and each facet
  is represented by two halfedges or two halffacets, respectively.

\cgalHeading{Template Parameters}

  The first parameter requires one of the following exact kernels:
  `Homogeneous`, `Simple_homogeneous`, `Extended_homogeneous`
  parametrized with `Gmpz`, `leda_integer` or any other number type
  modeling \f$\mathbb{Z}\f$, or `Cartesian`, `Simple_cartesian`,
  `Extended_cartesian` parametrized with `Gmpq`, `leda_rational`,
  `Quotient<Gmpz>` or any other number type modeling \f$\mathbb{Q}\f$.

  The second parameter and the third parameter are for future considerations.
  Neither `Nef_polyhedronItems_3` nor `Nef_polyhedronMarks` is
  specified, yet. Do not use any other than the default types for these two
  template parameters.

  \sa `CGAL::Nef_polyhedron_3::Vertex`
  \sa `CGAL::Nef_polyhedron_3::Halfedge`
  \sa `CGAL::Nef_polyhedron_3::Halffacet`
  \sa `CGAL::Nef_polyhedron_3::Volume`
  \sa `CGAL::Nef_polyhedron_3::SHalfedge`
  \sa `CGAL::Nef_polyhedron_3::SHalfloop`
  \sa `CGAL::Nef_polyhedron_3::SFace`
  \sa `CGAL::Nef_polyhedron_S2<Traits>`
  \sa `CGAL::Polyhedron_3<Traits>`

*/
template< class Nef_polyhedronTraits_3,
          class Nef_polyhedronItems_3 = CGAL::Default_items<Nef_polyhedronTraits_3>
          class Nef_polyhedronMarks = bool
          > class Nef_polyhedron_3 {
public:

/// \name Types
/// @{

/*!
  \ingroup PkgNef3Ref

  A Halfedge has a double meaning. In the global incidence structure of a
  `Nef_polyhedron_3` it is an oriented edge going from one vertex to another.
  A halfedge also coincides with an svertex of the sphere map of its source
  vertex. Because of this, we offer the types `Halfedge` and `SVertex`
  which are the same. Furthermore, the redundant functions `center_vertex()`
  and `source()` are provided. The reason is, that we get the same vertex
  either if we want to have the source vertex of a halfedge, or if we want to
  have the vertex in the center of the sphere map a svertex lies on.
  Figures \ref figureNef3HalfedgeIncidences
  and \ref figureNef3FacetIncidences
  illustrate the incidence of a svertex on a sphere map and of
  a halfedge in the global structure.

  As part of the global incidence structure, the member functions `source`
  and `target` return the source and target vertex of an edge. The member
  function `twin()` returns the opposite halfedge.

  Looking at the incidence structure on a sphere map, the member function
  `out_sedge` returns the first outgoing shalfedge, and `incident_sface`
  returns the incident sface.

\cgalHeading{Creation}

  There is no need for a user to create a `Halfedge` explicitly. The
  class `Nef_polyhedron_3<Traits>` manages the needed halfedges internally.

  \sa `CGAL::Nef_polyhedron_3::Vertex`
  \sa `CGAL::Nef_polyhedron_3::SHalfedge`
  \sa `CGAL::Nef_polyhedron_3::SFace`
  \sa `CGAL::Nef_polyhedron_S2::Sphere_point`

*/
  class Halfedge {
  public:

/// \name Types
/// The following types are the same as in `Nef_polyhedron_3<Traits>`.
/// @{

/*!
  type of mark.
*/
    typedef unspecified_type Mark;

/*!
  sphere point type stored in Halfedge.
*/
    typedef unspecified_type Sphere_point;

/*!
  const handle to vertex.
*/
    typedef unspecified_type Vertex_const_handle;

/*!
  const handle to halfedge.
*/
    typedef unspecified_type Halfedge_const_handle;

/*!
  const handle to SHalfedge.
*/
    typedef unspecified_type SHalfedge_const_handle;

/*!
  const handle to SFace.
*/
    typedef unspecified_type SFace_const_handle;

/// @}

/// \name Operations
/// @{

/*!
  the mark of `e` .
*/
    const Mark& mark() const;

/*!
  the sphere point of `e` .
*/
    const Sphere_point& point() const;

/*!
  returns |true| if `e` has no adjacent sedges.
*/
    bool is_isolated() const;

/*!
  the center vertex of the sphere map `e` belongs to.
*/
    Vertex_const_handle center_vertex() const;

/*!
  the source vertex of `e` .
*/
    Vertex_const_handle source() const;

/*!
  the target vertex `e`.
*/
    Vertex_const_handle target() const;

/*!
  the twin of `e` .
*/
    Halfedge_const_handle twin() const;

/*!
  the first out sedge of `e` .
*/
    SHalfedge_const_handle out_sedge() const;

/*!
  the incident sface of `e` .
*/
    SFace_const_handle incident_sface() const;

/// @}

  }; /* end Halfedge */

/*!
  \ingroup PkgNef3Ref

  The type `Halffacet_cycle_iterator` iterates over a list of
  `Object_handles`. Each item of that list can either be assigned
  to `SHalfedge_handle` or `SHalfloop_handle`. To find out which
  of these assignment works out, the member functions `is_shalfedge()`
  and `is_shalfloop()` are provided.

  \sa `CGAL::Nef_polyhedron_3::SHalfedge`
  \sa `CGAL::Nef_polyhedron_3::SHalfloop`

*/

  class Halffacet_cycle_iterator {
  public:

/// \name Types
/// @{

/*!
  const handle to SHalfedge.
*/
    typedef unspecified_type SHalfedge_handle;

/*!
  const handle to SHalfloop.
*/
    typedef unspecified_type SHalfloop_handle;

/// @}

/// \name Creation
/// @{

/*!
  default constructor.
*/
    Halffacet_cycle_iterator();

/// @}

/// \name Operations
/// @{

/*!
  returns true if `hfc` represents a `SHalfedge_handle`.
*/
    bool is_shalfedge() const;

/*!
  returns true if `hfc` represents a `SHalfloop_handle`.
*/
    bool is_shalfloop() const;

/*!
  casts `hfc` to `SHalfedge_handle`.
*/
    operator SHalfedge_handle() const;

/*!
  casts `hfc` to `SHalfloop_handle`.
*/
    operator SHalfloop_handle() const;

/// @}

  }; /* end Halffacet_cycle_iterator */


/*!
  \ingroup PkgNef3Ref

  A halffacet is an oriented, rectilinear bounded part of a plane. The following
  figure depicts the incidences to halfedges, vertices and the notion of facet
  cycles.


  \anchor figureNef3FacetIncidences
  \image html snc.png
  \image latex snc.png

  The member function `twin()` returns the opposite halffacet, `incident_volume`
  returns the incident volume. A Halffacet cycle either consists of consecutive
  shalfedges along the border (or a hole) of the halffacet, or of a single
  shalfloop on the sphere map of a vertex isolated on the halffacet. The
  iterator range (`halffacet_cycles_begin()`/`halffacet_cycles_end()`)
  provides an entry element for each halffacet cycle of a halffacet.

\cgalHeading{Creation}

  There is no need for a user to create a `Halffacet` explicitly. The
  class `Nef_polyhedron_3<Traits>` manages the needed halffacets internally.

  \sa `CGAL::Nef_polyhedron_3::Volume`
  \sa `CGAL::Nef_polyhedron_3::Halfedge`
  \sa `CGAL::Nef_polyhedron_3::SHalfedge`

*/

  class Halffacet {
  public:

/// \name Types
/// The following types are the same as in `Nef_polyhedron_3<Traits>`.
/// @{

/*!
  type of mark.
*/
    typedef unspecified_type Mark;

/*!
  plane type stored in Halffacet.
*/
    typedef unspecified_type Plane_3;

/*!
  list of Object handles.
*/
    typedef unspecified_type Object_list;

/*!
  const handle to Halffacet.
*/
    typedef unspecified_type Halffacet_const_handle;

/*!
  const handle to volume.
*/
    typedef unspecified_type Volume_const_handle;

/*!
  const iterator over the entries to all halffacet cycles of a halffacet.
*/
    typedef unspecified_type Halffacet_cycle_const_iterator;

/// @}

/// \name Operations
/// @{

/*!
  the mark of `f` .
*/
    const Mark& mark() const;

/*!
  the supporting plane of `f` .
*/
    const Plane_3& plane() const;

/*!
  the twin of `f` .
*/
    Halffacet_const_handle twin() const;

/*!
  the incident volume of `f` .
*/
    Volume_const_handle incident_volume() const;

/*!
  iterator over the entries to all halffacet cycles of `f` .
*/
    Halffacet_cycle_const_iterator facet_cycles_begin() const;

/*!
  past-the-end iterator.
*/
    Halffacet_cycle_const_iterator facet_cycles_end() const;

/// @}

  }; /* end Halffacet */

/*!
  \ingroup PkgNef3Ref

  The type `SFace_cycle_iterator` iterates over a list of
  `Object_handles`. Each item of that list can either be assigned
  to `SVertex_handle`, `SHalfedge_handle` or `SHalfloop_handle`.
  To find out which
  of these assignment works out, the member functions `is_svertex()`,
  `is_shalfedge()` and `is_shalfloop()` are provided.

  \sa `CGAL::Nef_polyhedron_3::Halfedge`
  \sa `CGAL::Nef_polyhedron_3::SHalfedge`
  \sa `CGAL::Nef_polyhedron_3::SHalfloop`

*/

  class SFace_cycle_iterator {
  public:

/// \name Types
/// @{

/*!
  const handle to SVertex.
*/
    typedef unspecified_type SVertex_handle;

/*!
  const handle to SHalfedge.
*/
    typedef unspecified_type SHalfedge_handle;

/*!
  const handle to SHalfloop.
*/
    typedef unspecified_type SHalfloop_handle;

/// @}

/// \name Creation
/// @{

/*!
  default constructor.
*/
    SFace_cycle_iterator();

/// @}

/// \name Operations
/// @{

/*!
  returns true if `sfc` represents a `SVertex_handle`.
*/
    bool is_svertex() const;

/*!
  returns true if `sfc` represents a `SHalfedge_handle`.
*/
    bool is_shalfedge() const;

/*!
  returns true if `sfc` represents a `SHalfloop_handle`.
*/
    bool is_shalfloop() const;

/*!
  casts `sfc` to `SVertex_handle`.
*/
    operator SVertex_handle() const;

/*!
  casts `sfc` to `SHalfedge_handle`.
*/
    operator SHalfedge_handle() const;

/*!
  casts `sfc` to `SHalfloop_handle`.
*/
    operator SHalfloop_handle() const;

/// @}

  }; /* end SFace_cycle_iterator */

/*!
  \ingroup PkgNef3Ref

  An sface is described by its boundaries.
 Figures \ref figureNef3HalfedgeIncidences
  and \ref figureNef3HalfloopIncidences
  illustrate the incidences of an sface.
 An entry item to each boundary cycle can be accessed
  using the iterator range `[sface_cycles_begin(), sface_cycles_end())`.
  Additionally, `Nef_polyhedron_S2` provides the macro
  `CGAL_forall_sface_cylces_of`. The iterators are of type
  `SFace_cycle_const_iterator` and represent either a shalfedge, a shalfloop,
  or a svertex.

\cgalHeading{Creation}

  There is no need for a user to create a `SFace` explicitly. The
  class `Nef_polyhedron_3<Traits>` manages the needed sfaces internally.

  \sa `CGAL::Nef_polyhedron_3::Vertex`
  \sa `CGAL::Nef_polyhedron_3::Volume`

*/

  class SFace {
  public:

/// \name Types
/// The following types are the same as in `Nef_polyhedron_3<Traits>`.
/// @{

/*!
  type of mark.
*/
    typedef unspecified_type Mark;

/*!
  list of Object handles.
*/
    typedef unspecified_type Object_list;

/*!
  const handle to Vertex.
*/
    typedef unspecified_type Vertex_const_handle;

/*!
  const handle to Volume.
*/
    typedef unspecified_type Volume_const_handle;

/*!
  const handle to SFace.
*/
    typedef unspecified_type SFace_const_handle;

/*!
  const iterator over the entries to all sface cycles of a sface.
*/
    typedef unspecified_type SFace_cycle_const_iterator;

/// @}

/// \name Operations
/// @{

/*!
  the mark of `sf` .
*/
    const Mark& mark() const;

/*!
  the center vertex of the sphere map `sf` belongs to.
*/
    Vertex_const_handle center_vertex() const;

/*!
  the volume that corresponds to `sf` in the 3D incidence structure.
*/
    Volume_const_handle volume() const;

/*!
  iterator over the entries to all sface cycles of `sf` .
*/
    SFace_cycle_const_iterator sface_cycle_begin() const;

/*!
  past-the-end iterator.
*/
    SFace_cycle_const_iterator sface_cycle_end() const;

/// @}

  }; /* end SFace */

/*!
  \ingroup PkgNef3Ref

  A shalfedge is a great arc on a sphere map.
  Figure \ref figureNef3HalfedgeIncidences
  depicts the relationship between a shalfedge and its incident
  shalfedges, svertices, and sfaces on a sphere map. A shalfedge is
  an oriented sedge between two svertices. It is always paired with a
  shalfedge pointing in
  the opposite direction. The `twin()` member function returns
  this shalfedge of opposite orientation.

  \anchor figureNef3HalfedgeIncidences
  \image html shalfedge.png
  \image latex shalfedge.png

  The `snext()` member function points
  to the successor shalfedge around this sface while the `sprev()` member
  function points to the preceding shalfedge. An
  successive assignments of the form `se = se->snext()` cycles
  counterclockwise around the sface (or hole).

  Similarly, the successive
  assignments of the form `se = se->snext()->twin()` cycle
  clockwise around the svertex and traverse all halfedges incident to
  this svertex. The assignment `se = se->cyclic_adj_succ()` can be
  used as a shortcut.

  The role of shalfedges in a facet is illustrated in
  Figure \ref figureNef3FacetIncidences.
  The `facet()` member function returns the facet in which
  the shalfedge is part of one of the facet cycles. The successive assignment of
  the form `se = se->next()` cycles counterclockwise around the facet (or a
  hole of the facet).

  A const circulators is provided for each of the three circular orders.
  The circulators are bidirectional and assignable to `SHalfedge_const_handle`.

\cgalHeading{Creation}

  There is no need for a user to create a `SHalfedge` explicitly. The
  class `Nef_polyhedron_3<Traits>` manages the needed shalfedges internally.

  \sa `CGAL::Nef_polyhedron_3::Halfedge`
  \sa `CGAL::Nef_polyhedron_3::Halffacet`
  \sa `CGAL::Nef_polyhedron_3::SFace`
  \sa `CGAL::Nef_polyhedron_S2::Sphere_circle`

*/

  class SHalfedge {
  public:

/// \name Types
/// The following types are the same as in `Nef_polyhedron_3<Traits>`.
/// @{

/*!
  type of mark.
*/
    typedef unspecified_type Mark;

/*!
  sphere circle type stored in SHalfedge.
*/
    typedef unspecified_type Sphere_circle;

/*!
  const handle to Halffacet.
*/
    typedef unspecified_type Halffacet_const_handle;

/*!
  const handle to SVertex.
*/
    typedef unspecified_type SVertex_const_handle;

/*!
  const handle to SHalfedge.
*/
    typedef unspecified_type SHalfedge_const_handle;

/*!
  const handle to SFace.
*/
    typedef unspecified_type SFace_const_handle;

/// @}

/// \name Operations
/// @{

/*!
  the mark of `se` .
*/
    const Mark& mark() const;

/*!
  the sphere circle of `se` .
*/
    const Sphere_circle& circle() const;

/*!
  the twin of `se` .
*/
    SHalfedge_const_handle twin() const;

/*!
  the source svertex of `se` .
*/
    SVertex_const_handle source() const;

/*!
  equals `twin()->source()`.
*/
    SVertex_const_handle target() const;

/*!
  the SHalfedge previous to `se` in a facet cycle.
*/
    SHalfedge_const_handle prev() const;

/*!
  the next SHalfedge of `se` in a facet cycle.
*/
    SHalfedge_const_handle next() const;

/*!
  the SHalfedge previous to `se` in a sface cycle.
*/
    SHalfedge_const_handle sprev() const;

/*!
  the next SHalfedge of `se` in a sface cycle.
*/
    SHalfedge_const_handle snext() const;

/*!
  the edge before `se` in the cyclic ordered adjacency list of source().
*/
    SHalfedge_const_handle cyclic_adj_pred() const;

/*!
  the edge after `se` in the cyclic ordered adjacency list of source().
*/
    SHalfedge_const_handle cyclic_adj_succ() const;

/*!
  the facet that corresponds to `se` in the 3D incidence structure.
*/
    Halffacet_const_handle facet() const;

/*!
  the incident
  sface of `se` .
*/
    SFace_const_handle incident_sface() const;

/*!
  determines whether `se` is
  in an outer sface cycle.
*/
    bool in_outer_sface_cycle() const;

/*!
  determines whether `se` is
  in an inner sface cycle.
*/
    bool in_inner_sface_cycle() const;

/*!
  determines whether `se` is
  in an outer facet cycle.
*/
    bool in_outer_facet_cycle() const;

/*!
  determines whether `se` is
  in an inner facet cycle.
*/
    bool in_inner_facet_cycle() const;

/// @}

  }; /* end SHalfedge */

/*!
  \ingroup PkgNef3Ref

  A shalfloop is a great circle on a sphere map.
  Figure \ref figureNef3HalfloopIncidences
  depicts the relationship between a shalfloop and its incident
  shalfloops, and sfaces on a sphere map. A shalfloop is
  an oriented sloop. It is always paired with a
  shalfloop whose supporting `Sphere_circle` is pointing in
  the opposite direction. The `twin()` member function returns
  this shalfloop of opposite orientation.

  \anchor figureNef3HalfloopIncidences
  \image html shalfloopB.png
  \image latex shalfloopB.png

  A sphere map having a shalfloop models the neighborhood of a vertex which is
  isolated on a facet. That facet is returned by the member function
  `facet`.

\cgalHeading{Creation}

  There is no need for a user to create a `SHalfloop` explicitly. The
  class `Nef_polyhedron_3<Traits>` manages the needed shalfloops internally.

  \sa `CGAL::Nef_polyhedron_3::Halffacet`
  \sa `CGAL::Nef_polyhedron_3::SFace`
  \sa `CGAL::Nef_polyhedron_S2::Sphere_point`

*/

  class SHalfloop {
  public:

/// \name Types
/// The following types are the same as in `Nef_polyhedron_3<Traits>`.
/// @{

/*!
  type of mark.
*/
    typedef unspecified_type Mark;

/*!
  sphere circle type stored in SHalfloop.
*/
    typedef unspecified_type Sphere_circle;

/*!
  const handle to Halffacet.
*/
    typedef unspecified_type Halffacet_const_handle;

/*!
  const handle to SHalfloop.
*/
    typedef unspecified_type SHalfloop_const_handle;

/*!
  const handle to SFace.
*/
    typedef unspecified_type SFace_const_handle;

/// @}

/// \name Operations
/// @{

/*!
  the mark of `se` .
*/
    const Mark& mark() const;

/*!
  the sphere circle of `se` .
*/
    const Sphere_circle& circle() const;

/*!
  the twin of `se` .
*/
    SHalfloop_const_handle twin() const;

/*!
  the facet that corresponds to `se` in the 3D incidence structure.
*/
    Halffacet_const_handle facet() const;

/*!
  the incident sface of `se` .
*/
    SFace_const_handle incident_sface() const;

/// @}

  }; /* end SHalfloop */

/*!
  \ingroup PkgNef3Ref

  A vertex is a point in the 3-dimensional space. Its incidence
  structure can be accessed creating a sphere map of the vertex.
  This is done by the member function `Nef_polyhedron_3::get_sphere_map()`.

\cgalHeading{Creation}

  There is no need for a user to create a `Vertex` explicitly. The
  class `Nef_polyhedron_3<Traits>` manages the needed vertices internally.

  \sa `CGAL::Nef_polyhedron_3<Traits>`
  \sa `CGAL::Nef_polyhedron_S2<Traits>`

*/
  class Vertex {
  public:

/// \name Types
/// The following types are the same as in `Nef_polyhedron_3<Traits>`.
/// @{

/*!
  type of mark.
*/
    typedef unspecified_type Mark;

/*!
  point type stored in Vertex.
*/
    typedef unspecified_type Point_3;

/// @}

/// \name Operations
/// @{

/*!
  the mark of `v` .
*/
    const Mark& mark() const;

/*!
  the point of `v` .
*/
    const Point_3& point() const;

/// @}

  }; /* end Vertex */

/*!
  \ingroup PkgNef3Ref

  A volume is a full-dimensional connected point set in \f$ \mathbb{R}^3\f$. It is
  bounded by several shells, i.e.\ a connected part of the boundary incident
  to a single volume. An entry element to each shell is provided by the
  iterator range (`shells_begin()`/`shells_end()`). A
  `Shell_entry_iterator` is assignable to `SFace_handle`.

\cgalHeading{Creation}

  There is no need for a user to create a `Volume` explicitly. The
  class `Nef_polyhedron_3<Traits>` manages the needed volumes internally.

  \sa `CGAL::Nef_polyhedron_3::SFace`

*/

  class Volume {
  public:

/// \name Types
/// The following types are the same as in `Nef_polyhedron_3<Traits>`.
/// @{

/*!
  type of mark.
*/
    typedef unspecified_type Mark;

/*!
  list of Object handles.
*/
    typedef unspecified_type Object_list;

/*!
  const handle to Volume.
*/
    typedef unspecified_type Volume_const_handle;

/*!
  const iterator over the entries to all shells adjacent to a volume.
*/
    typedef unspecified_type Shell_entry_const_iterator;

/// @}

/// \name Operations
/// @{

/*!
  the mark of `c` .
*/
    const Mark& mark() const;

/*!
  const iterator over the entries to all shells adjacent to `c` .
*/
    Shell_entry_const_iterator shells_begin() const;

/*!
  past-the-end iterator.
*/
    Shell_entry_const_iterator shells_end() const;

/// @}

  }; /* end Volume */

/*!
  traits class selected for `Nef_polyhedronTraits_3`.
*/
  typedef unspecified_type Traits;

/*!
  All object (vertices, edges, etc.) are attributed by a Mark.
  Mark equals bool.
*/
  typedef unspecified_type Mark;

/*!
  size type of `Nef_polyhedron_3`.
*/
  typedef unspecified_type size_type;

/*!
  non-mutable handle to a vertex.
*/
  typedef unspecified_type Vertex_const_handle;

/*!
  non-mutable handle to a halfedge.
*/
  typedef unspecified_type Halfedge_const_handle;

/*!
  non-mutable handle to a halffacet.
*/
  typedef unspecified_type Halffacet_const_handle;

/*!
  non-mutable handle to a volume.
*/
  typedef unspecified_type Volume_const_handle;

/*!
  non-mutable handle to a svertex.
*/
  typedef unspecified_type SVertex_const_handle;

/*!
  non-mutable handle to a shalfedge.
*/
  typedef unspecified_type SHalfedge_const_handle;

/*!
  non-mutable handle to a shalfloop.
*/
  typedef unspecified_type SHalfloop_const_handle;

/*!
  non-mutable handle to a sface.
*/
  typedef unspecified_type SFace_const_handle;

/*!
  non-mutable iterator over all vertices.
*/
  typedef unspecified_type Vertex_const_iterator;

/*!
  non-mutable iterator over all halfeges.
*/
  typedef unspecified_type Halfedge_const_iterator;

/*!
  non-mutable iterator over all halffacets.
*/
  typedef unspecified_type Halffacet_const_iterator;

/*!
  non-mutable iterator over all volumes.
*/
  typedef unspecified_type Volume_const_iterator;

/*!
  non-mutable iterator over all svertices.
*/
  typedef unspecified_type SVertex_const_iterator;

/*!
  non-mutable iterator over all shalfedges.
*/
  typedef unspecified_type SHalfedge_const_iterator;

/*!
  non-mutable iterator over all shalfloops.
*/
  typedef unspecified_type SHalfloop_const_iterator;

/*!
  non-mutable iterator over all sfaces.
*/
  typedef unspecified_type SFace_const_iterator;

/*!
  non-mutable circulator of shalfedges around a svertex (cw).
*/
  typedef unspecified_type SHalfedge_around_svertex_const_circulator;

/*!
  non-mutable circulator of shalfedges around a sface (ccw).
*/
  typedef unspecified_type SHalfedge_around_sface_const_circulator;

/*!
  non-mutable circulator of shalfedges around a halffacet (ccw).
*/
  typedef unspecified_type SHalfedge_around_facet_const_circulator;

/*!
  non-mutable iterator over the cycles of a sface.
*/
  typedef unspecified_type SFace_cycle_const_iterator;

/*!
  non-mutable iterator over the cycles of a halffacet.
*/
  typedef unspecified_type Halffacet_cycle_const_iterator;

/*!
  non-mutable iterator providing an entry to each shell.
*/
  typedef unspecified_type Shell_entry_const_iterator;

/*!
  a generic handle to an object.
  The kind of object `(vertex, halfedge, halffacet, volume, svertex, shalfedge, shalfloop, sface)` can
  be determined and the object can be assigned to a corresponding
  constant handle by one of the following functions:

  `bool assign(Vertex_const_handle& h, Object_handle)`

  `bool assign(Halfedge_const_handle& h, Object_handle)`

  `bool assign(Halffacet_const_handle& h, Object_handle)`

  `bool assign(Volume_const_handle& h, Object_handle)`

  `bool assign(SVertex_const_handle& h, Object_handle)`

  `bool assign(SHalfedge_const_handle& h, Object_handle)`

  `bool assign(SHalfloop_const_handle& h, Object_handle)`

  `bool assign(SFace_const_handle& h, Object_handle)`

  where each function returns `true` iff the assignment to
  `h` could be accomplished.
*/
  typedef unspecified_type Object_handle;

/*!
  location of vertices.
*/
  typedef unspecified_type Point_3;

/*!
  segment represented by a halfedge.
*/
  typedef unspecified_type Segment_3;

/*!
  direction of a halfedge.
*/
  typedef unspecified_type Vector_3;

/*!
  plane of a halffacet lies in.
*/
  typedef unspecified_type Plane_3;

/*!
  affine transformation.
*/
  typedef unspecified_type Aff_transformation_3;

/*!
  tag for calling polyline constructor.
*/
  typedef unspecified_type Polylines_tag;

/*!
  tag for calling point constructor.
*/
  typedef unspecified_type Points_tag;

/*!
  construction selection.
*/
  enum Boundary { EXCLUDED, INCLUDED };

/*!
  construction selection.
*/
  enum Content { EMPTY, COMPLETE };

/*!
  intersection selection.
*/
  enum Intersection_mode { CLOSED_HALFSPACE, OPEN_HALFSPACE,
                           PLANE_ONLY };

/*!
  a sphere map.
*/
  typedef unspecified_type Nef_polyhedron_S2;

/*!
  a polyhedral surface.
*/
  typedef unspecified_type Polyhedron;

/// @}

/// \name Creation
/// @{

/*!

  creates a Nef polyhedron and initializes it to the empty set if `plane
  == EMPTY` and to the whole space if `space == COMPLETE`.
*/
  Nef_polyhedron_3(Content space = EMPTY);

/*!
  creates a
  Nef polyhedron containing the halfspace on the negative side of
  `p` including `p` if `b==INCLUDED`, excluding `p` if
  `b==EXCLUDED`.
*/
  Nef_polyhedron_3(const Plane_3& p,
                   Boundary b = INCLUDED);

/*!

  creates a Nef polyhedron, which represents the same point set as
  the polyhedral surface `P` does.
*/
  Nef_polyhedron_3(Polyhedron& P);

/*!
  creates a Nef polyhedron, which represents the same point set as
  the polyhedral surface `pm` does. `him` and `fim` must be both initialized
  so that halfedges and faces are indexed in `[0, num_halfedges(pm)[`
  and `[0, num_faces(pm)[` respectively. If `PolygonMesh` has an internal
  halfedge index map and an internal face index map, the last two parameters
  can be omitted.
  \tparam PolygonMesh a model of `FaceListGraph` and `VertexListGraph`.
  \tparam HalfedgeIndexMap a class model of `ReadablePropertyMap` with
          `boost::graph_traits<PolygonMesh>::%halfedge_descriptor` as key type
           a value type convertible to `std::size_t`
  \tparam FaceIndexMap a class model of `ReadablePropertyMap` with
          `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type
           a value type convertible to `std::size_t`
*/
 template <class PolygonMesh, class HalfedgeIndexMap, class FaceIndexMap>
 explicit Nef_polyhedron_3(const PolygonMesh& pm,
                           const HalfedgeIndexMap& him,
                           const FaceIndexMap& fim);

/*!
  creates a Nef polyhedron consisting of a single polygon
  spanned by the list of points in the iterator range
  `[begin,end)`. If the points do not lie on a common
  supporting plane, the constructor tries to triangulate
  the polygon into multiple facets.If the construction does
  not succeed, the empty set is created.
*/
  Nef_polyhedron_3(Input_iterator begin, Input_iterator end);

/*!
  creates a Nef polyhedron consisting of polylines.
  The iterator range [it, end) defines a range of polylines.
  Each polyline is defined as a range of points, first and
  past-the-end iterators being provided as a `std::pair` of iterators.
*/
  template <class Forward_iterator>
  Nef_polyhedron_3(Forward_iterator it, Forward_iterator end,
                   Polylines_tag);

/*!
  creates a Nef polyhedron that consists only of points.
  The iterator range [it, end) defines a range of points.
*/
  template <class Forward_iterator>
  Nef_polyhedron_3(Forward_iterator it, Forward_iterator end,
                   Points_tag);

/*!
  creates a Nef polyhedron that consists of point p.
*/
  explicit
  Nef_polyhedron_3(const Point_3& p);

/*!
  creates a Nef polyhedron that consists of segment s.
*/
  explicit
  Nef_polyhedron_3(const Segment_3& s);

/// @}

/// \name Access Member Functions
/// The following macros are provided: `CGAL_forall_vertices(v,N)`,
/// `CGAL_forall_halfedges(e,N)`, `CGAL_forall_edges(e,N)`,
/// `CGAL_forall_halffacets(f,N)`, `CGAL_forall_facets(f,N)`,
/// `CGAL_forall_volumes(c,N)` where `N` is a `Nef_polyhedron_3`.
/// @{

/*!
  returns true, if `N` is a 2-manifold.
*/
  bool is_simple() const;

/*!
  checks the integrity of `N` .
*/
  bool is_valid() const;

/*!
  returns the number of vertices.
*/
  Size_type number_of_vertices() const;

/*!
  return the number of halfedges.
*/
  Size_type number_of_halfedges() const;

/*!
  returns the number of halfedge pairs.
*/
  Size_type number_of_edges() const;

/*!
  returns the number of halffacets.
*/
  Size_type number_of_halffacets() const;

/*!
  returns the number of halffacet pairs.
*/
  Size_type number_of_facets() const;

/*!
  returns the number of volumes.
*/
  Size_type number_of_volumes() const;

/*!
  iterator over all vertices.
*/
  Vertex_const_iterator vertices_begin() const;

/*!
  past-the-end iterator.
*/
  Vertex_const_iterator vertices_end() const;

/*!
  iterator over all halfedges.
*/
  Halfedge_const_iterator halfedges_begin()const;

/*!
  past-the-end iterator.
*/
  Halfedge_const_iterator halfedges_end() const;

/*!
  iterator over all halffacets.
*/
  Halffacet_const_iterator halffacets_begin() const;

/*!
  past-the-end iterator.
*/
  Halffacet_const_iterator halffacets_end() const;

/*!
  iterator over all volumes.
*/
  Volume_const_iterator volumes_begin() const;

/*!
  past-the-end iterator.
*/
  Volume_const_iterator volumes_end() const;

/*!
  returns a generic handle to the object (vertex, edge, facet or volume) which contains the point p in its relative interior.
*/
  Object_handle locate(const Point_3& p) const;

/*!
  returns the neighborhood of a vertex modeled by a `Nef_polyhedron_S2`.
*/
  Nef_polyhedron_S2 get_sphere_map(Vertex_const_iterator v) const;

/// @}

/// \name Point Set Predicates
/// @{

/*!
  returns true, if `N` is the
  empty point set.
*/
  bool is_empty() const;

/*!
  returns true, if `N` is the
  complete 3D space.
*/
  bool is_space() const;

/*!
  returns true, if `N` and N1 comprise the same point sets.
*/
  bool operator==(const Nef_polyhedron_3<Traits>& N1) const;

/*!
  returns true, if `N` and N1 comprise different point sets.
*/
  bool operator!=(const Nef_polyhedron_3<Traits>& N1) const;

/*!
  returns true, if `N` is a proper subset of N1.
*/
  bool operator<(const Nef_polyhedron_3<Traits>& N1) const;

/*!
  returns true, if `N` is a proper superset of N1.
*/
  bool operator>(const Nef_polyhedron_3<Traits>& N1) const;

/*!
  returns true, if `N` is a subset of N1.
*/
  bool operator<=(const Nef_polyhedron_3<Traits>& N1) const;

/*!
  returns true, if `N` is a superset of N1.
*/
  bool operator>=(const Nef_polyhedron_3<Traits>& N1) const;

/// @}

/// \name Unary Set Operations
/// @{

/*!
  returns the complement of `N` .
*/
  Nef_polyhedron_3<Traits> complement() const;

/*!
  returns the interior of `N` .
*/
  Nef_polyhedron_3<Traits> interior() const;

/*!
  returns the boundary of `N` .
*/
  Nef_polyhedron_3<Traits> boundary() const;

/*!
  returns the closure of `N` .
*/
  Nef_polyhedron_3<Traits> closure() const;

/*!
  returns the regularization, i.e.\ the closure of the interior, of `N` .
*/
  Nef_polyhedron_3<Traits> regularization() const;

/*!
  returns the complement of `N` .
*/
  Nef_polyhedron_3<Traits> operator!() const;

/// @}

/// \name Binary Set Operations
/// @{

/*!
  return the intersection of `N` and N1.
*/
  Nef_polyhedron_3<Traits> intersection(const Nef_polyhedron_3<Traits>& N1) const;

/*!
  return the union of `N` and N1. (Note that ''union'' is a C++ keyword and cannot be used for this operation.)
*/
  Nef_polyhedron_3<Traits> join(const Nef_polyhedron_3<Traits>& N1) const;

/*!
  return the difference between `N` and N1.
*/
  Nef_polyhedron_3<Traits> difference(const Nef_polyhedron_3<Traits>& N1) const;

/*!
  return the symmetric difference of `N` and N1.
*/
  Nef_polyhedron_3<Traits> symmetric_difference(const Nef_polyhedron_3<Traits>& N1) const;

/*!
  returns intersection of `N` with
  plane (`im=PLANE_ONLY`), open halfspace (`im=OPEN_HALFSPACE`), or closed
  halfspace (`im=CLOSED_HALFSPACE`). In the latter two cases, the
  halfspaces are on the negative side of the plane `p`. The function is
  written for the use with standard kernels, where halfspaces are not
  part of the domain. The function does not work in combination with
  an extended kernels or with an unbounded polyhedron.
*/
  Nef_polyhedron_3<Traits> intersection(const Plane_3& p, Intersection_mode im) const;

/*!
  return the intersection of `N` and N1.
*/
  Nef_polyhedron_3<Traits> operator*(const Nef_polyhedron_3<Traits>& N1) const;

/*!
  return the union of `N` and N1.
*/
  Nef_polyhedron_3<Traits> operator+(const Nef_polyhedron_3<Traits>& N1) const;

/*!
  return the difference between `N` and N1.
*/
  Nef_polyhedron_3<Traits> operator-(const Nef_polyhedron_3<Traits>& N1) const;

/*!
  return the symmetric difference of `N` and N1.
*/
  Nef_polyhedron_3<Traits> operator^(const Nef_polyhedron_3<Traits>& N1) const;

/*!
  intersects `N` and N1.
*/
  void operator*=(const Nef_polyhedron_3<Traits>& N1);

/*!
  unites `N` with N1.
*/
  void operator+=(const Nef_polyhedron_3<Traits>& N1);

/*!
  subtracts N1 from `N` .
*/
  void operator-=(const Nef_polyhedron_3<Traits>& N1);

/*!
  performs a symmetric intersection of `N` and N1.
*/
  void operator^=(const Nef_polyhedron_3<Traits>& N1);

/// @}

/// \name Operations
/// @{

/*!
  make `N` the empty set if `space == EMPTY` and
  the complete 3D space if `space == COMPLETE`.
*/
  void clear(Content space = EMPTY);

/*!
  applies an affine transformation to `N` .
*/
  void transform(const Aff_transformation_3& aff);

/*!
  converts `N` into a triangulated Polyhedron.
  For conversion to other types, see `convert_nef_polyhedron_to_polygon_mesh()`.
  \pre `N` is simple.
*/
  void convert_to_polyhedron(Polyhedron& P) const;

/*!
  calls the visit function of `V` for every item which belongs to the same shell as `sf`.
*/
  void visit_shell_objects(SFace_const_handle f, Visitor& V) const;

/// @}

}; /* end Nef_polyhedron_3 */
} /* end namespace CGAL */
