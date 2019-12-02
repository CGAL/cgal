namespace CGAL {

/*! \ingroup PkgArrangementOnSurface2Ref
 *
 * \anchor arr_refarr
 *
 * An object `arr` of the class `Arrangement_2` represents the planar
 * subdivision induced by a set of \f$ x\f$-monotone curves and isolated points
 * into maximally connected cells. The arrangement is represented as a
 * doubly-connected edge-list (<span class="textsc">Dcel</span>) such that each
 * <span class="textsc">Dcel</span> vertex is associated with a point of the
 * plane and each edge is associated with an \f$ x\f$-monotone curve whose
 * interior is disjoint from all other edges and vertices. Recall that an
 * arrangement edge is always comprised of a pair of twin <span
 * class="textsc">Dcel</span> halfedges.
 *
 * The `Arrangement_2` template has two parameters:
 * <UL>
 * <LI>The `Traits` template-parameter should be instantiated with
 * a model of the `ArrangementBasicTraits_2` concept. The traits
 * class defines the types of \f$ x\f$-monotone curves and two-dimensional
 * points, namely `ArrangementBasicTraits_2::X_monotone_curve_2` and
 * `ArrangementBasicTraits_2::Point_2`,
 *   respectively, and supports basic geometric predicates on them.
 * <LI>The `Dcel` template-parameter should be instantiated with
 * a class that is a model of the `ArrangementDcel` concept. The
 * value of this parameter is by default
 * `Arr_default_dcel<Traits>`.
 * </UL>
 *
 * The available traits classes and <span class="textsc">Dcel</span> classes are
 * described below.
 *
 * \sa `ArrangementDcel`
 * \sa `Arr_default_dcel<Traits>`
 * \sa `ArrangementBasicTraits_2`
 * \sa `CGAL::overlay()`
 * \sa `CGAL::is_valid()`
 *
 * Insertion Functions
 *
 * \sa `PkgArrangementOnSurface2Insert`
 * \sa `CGAL::insert_non_intersecting_curve()`
 * \sa `CGAL::insert_non_intersecting_curves()`
 * \sa `CGAL::insert_point()`
 *
 * Removal functions
 *
 * \sa `CGAL::remove_edge()`
 * \sa `CGAL::remove_vertex()`
 *
 * Input/output functions
 *
 * \sa `PkgArrangementOnSurface2Read`
 * \sa `PkgArrangementOnSurface2Write`
 */
template <typename Traits, typename Dcel>
class Arrangement_2 : public Arrangement_on_surface_2<Traits, typename Default_planar_topology<Traits, Dcel>::Traits> {
public:
  /// \name Types
  /// @{

  /*! a private type used as an abbreviation of the `Arrangement_2` type
   * hereafter.
   */
  typedef Arrangement_2<Traits_2, Dcel> Self;

  /*! the traits class in use. */
  typedef Traits Traits_2;

  /*! the <span class="textsc">Dcel</span> representation of the arrangement. */
  typedef unspecified_type Dcel;

  /// @}

  /// \name Types inherited from the base Arrangement_on_surface_2
  /// @{

  typedef typename Arrangement_on_surface_2::Point_2            Point_2;
  typedef typename Arrangement_on_surface_2::X_monotone_curve_2 X_monotone_curve_2;
  typedef typename Arrangement_on_surface_2::Size               Size;
  typedef typename Arrangement_on_surface_2::Vertex              Vertex;
  typedef typename Arrangement_on_surface_2::Halfedge            Halfedge;
  typedef typename Arrangement_on_surface_2::Face                Face;
  typedef typename Arrangement_on_surface_2::Vertex_handle       Vertex_handle;
  typedef typename Arrangement_on_surface_2::Halfedge_handle     Halfedge_handle;
  typedef typename Arrangement_on_surface_2::Face_handle         Face_handle;
  typedef typename Arrangement_on_surface_2::Vertex_iterator     Vertex_iterator;
  typedef typename Arrangement_on_surface_2::Halfedge_iterator   Halfedge_iterator;
  typedef typename Arrangement_on_surface_2::Edge_iterator       Edge_iterator;
  typedef typename Arrangement_on_surface_2::Face_iterator       Face_iterator;
  typedef typename Arrangement_on_surface_2::Unbounded_face_iterator Unbounded_face_iterator;
  typedef typename Arrangement_on_surface_2::Halfedge_around_vertex_circulator Halfedge_around_vertex_circulator;
  typedef typename Arrangement_on_surface_2::Ccb_halfedge_circulator Ccb_halfedge_circulator;
  typedef typename Arrangement_on_surface_2::Hole_iterator       Hole_iterator;
  typedef typename Arrangement_on_surface_2::Isolated_vertex_iterator Isolated_vertex_iterator;

  /// @}

  /// \name Creation
  /// @{

  /*! constructs an empty arrangement containing one unbounded face,
   * which corresponds to the entire plane.
   */
  Arrangement_2<Traits, Dcel>();

  /*! copy constructor. */
  Arrangement_2<Traits, Dcel>(const Self& other);

  /*! constructs an empty arrangement that uses the given `traits`
   * instance for performing the geometric predicates.
   */
  Arrangement_2<Traits, Dcel>(const Traits_2 *traits);

  /// @}

  /// \name Assignment Methods
  /// @{

  /*! assignment operator. */
  Self& operator=(other);

  /*! assigns the contents of another arrangement. */
  void assign(const Self& other);

  /*! clears the arrangement. */
  void clear();

  /// @}

  /// \name Access Functions
  /// @{

  /*! obtains the traits object used by the arrangement instance.
   * A `const` version is also available.
   */
  Traits_2* traits();

  /// @}
}; /* end Arrangement_2 */
} /* end namespace CGAL */

namespace CGAL {

/*! \ingroup PkgArrangementOnSurface2Insert insert The function `%insert`
 * inserts one or more curves or \f$ x\f$-monotone curves into a given
 * arrangement, where no restrictions are imposed on the inserted curves. If an
 * inserted curve is not \f$ x\f$-monotone curve, it is subdivided into \f$
 * x\f$-monotone subcurves (and perhaps isolated points), which are inserted
 * into the arrangement.
 *
 * \cgalHeading{Requirements}
 *
 * <UL>
 * <LI>If the curve is \f$ x\f$-monotone curve then The instantiated
 *   `Traits` class must model the `ArrangementXMonotoneTraits_2`
 *   concept. In case that the curve is not \f$ x\f$-monotone then the
 *   instantiated `Traits` class must model the
 *   `ArrangementTraits_2` concept. That is, it should define the
 *   `Curve_2` type, and support its subdivision into \f$ x\f$-monotone
 *   subcurves (and perhaps isolated points).
 * <LI>The point-location object `pl`, must model the
 *   `ArrangementPointLocation_2` concept.
 * </UL>
 */
/// @{

/*! Inserts the given curve `c` into the arrangement `arr`.
 * `c` is subdivided into \f$ x\f$-monotone subcurves (and perhaps isolated
 * points). Each subcurve is in turn inserted into the arrangement by locating
 * its left endpoint and computing its zone until reaching the right endpoint.
 *
 * The given point-location object `pl` is used to locate the left
 * endpoints of the \f$ x\f$-monotone curves. By default, the function uses the
 * "walk along line" point-location strategy  -  namely an instance of
 * the class `Arr_walk_along_line_point_location<Arrangement_2<Traits,Dcel> >`.
 *
 * \pre If provided, `pl` must be attached to the given arrangement `arr`.
 */
template <typename Traits, typename Dcel,
          typename Curve, typename PointLocation>
void insert(Arrangement_2<Traits,Dcel>& arr, const Curve& c,
            const PointLocation& pl = walk_pl);

/*! Inserts the<I>\f$ x\f$-monotone (only)</I> curve `xc` into the arrangement
 * `arr`. The object `obj`, which either wraps a `Vertex_const_handle`, a
 * `Halfedge_const_handle`, or a `Face_const_handle`, represents the location of
 * `xc`'s left endpoint in the arrangement. The zone of `xc` is computed
 * starting from the feature represented by `obj`. As in the case above, the
 * zone computation terminates, when the right endpoint is reached.  Thus,
 * point-location is not required.
 */
template <typename Traits, typename Dcel>
void insert(Arrangement_2<Traits,Dcel>& arr,
            const typename Traits::X_monotone_curve_2& xc, const Object& obj);

/*! Aggregately inserts the curves or \f$ x\f$-monotone curves in the range
 * `[first,last)` into the arrangement `arr` using the sweep-line framework.
 */
template <typename Traits, typename Dcel, class InputIterator>
void insert(Arrangement_2<Traits,Dcel>& arr,
            InputIterator first, InputIterator last);

/// @}

/*! \ingroup PkgArrangementOnSurface2Funcs
 *
 * Checks if a given curve or \f$ x\f$-monotone curve intersects an existing
 * arrangement's edges or vertices.
 *
 * If the give curve is not an \f$ x\f$-monotone curve then the function
 * subdivides the given curve into \f$ x\f$-monotone subcurves and isolated
 * vertices . Each subcurve is in turn checked for intersection.  The function
 * uses the zone algorithm to check if the curve intersects the
 * arrangement. First, the curve's left endpoint is located. Then, its zone is
 * computed starting from its left endpoint location. The zone computation
 * terminates when an intersection with an arrangement's edge/vertex is found or
 * when the right endpoint is reached.
 *
 * A given point-location object is used for locating the left endpoint
 * of the given curve in the existing arrangement. By default, the function
 * uses the "walk along line" point-location strategy - namely an
 * instance of the class
 * `Arr_walk_along_line_point_location<Arrangement_2<Traits,Dcel> >`.
 *
 * Checks if the given curve or \f$ x\f$-monotone curve `c` intersects
 * edges or vertices of the existing arrangement `arr`.
 * \pre If provided, `pl` must be attached to the given arrangement `arr`.
 *
 * \cgalHeading{Requirements}
 *
 * <UL>
 * <LI>If `c` is \f$ x\f$-monotone then the instantiated `GeomTraits`
 * class must model the `ArrangementXMonotoneTraits_2` concept. If
 * `c` is a curve then the instantiated `GeomTraits` class must
 * model the `ArrangementTraits_2` concept. That is, it should
 * define the `Curve_2` type, and support its subdivision into
 * \f$ x\f$-monotone subcurves (and perhaps isolated points).
 * <LI>The point-location object `pl`, must model the
 * `ArrangementPointLocation_2` concept.
 * </UL>
 */
template <typename GeomTraits, typename TopTraits,
          typename Curve, typename PointLocation>
bool do_intersect(Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
                  const Curve& c, const PointLocation& pl);

/*! \ingroup PkgArrangementOnSurface2Funcs
 *
 * Inserts a given \f$ x\f$-monotone curve into a given arrangement, where the
 * interior of the given curve is disjoint from all existing arrangement
 * vertices and edges. Under this assumption, it is possible to locate the
 * endpoints of the given curve in the arrangement, and use one of the
 * specialized insertion member-functions of the arrangement according to the
 * results. The insertion operations creates a single new edge, that is, two
 * twin halfedges, and the function returns a handle for the one directed
 * lexicographically in increasing order (from left to right).
 *
 * A given point-location object is used for answering the two point-location
 * queries on the given curve endpoints. By default, the function uses the "walk
 * along line" point-location strategy - namely, an instance of the class
 * `Arr_walk_along_line_point_location<Arrangement_2<Traits,Dcel> >`.
 *
 * \pre If provided, `pl` must be attached to the given arrangement `arr`.
 *
 * \cgalHeading{Requirements}
 *
 * <UL>
 * <LI>The instantiated `Traits` class must model the restricted
 * `ArrangementBasicTraits_2` concept, as no intersections are computed.
 * <LI>The point-location object `pl` must model the
 * `ArrangementPointLocation_2` concept.
 * </UL>
 */
template <typename Traits, typename Dcel,
          typename PointLocation>
typename Arrangement_2<Traits,Dcel>::Halfedge_handle
insert_non_intersecting_curve(Arrangement_2<Traits,Dcel>& arr,
                              const typename Traits::X_monotone_curve_2& xc,
                              const PointLocation& pl = walk_pl);

/*! \ingroup PkgArrangementOnSurface2Funcs
 *
 * Inserts a set of \f$ x\f$-monotone curves in a given range into a given
 * arrangement. The insertion is performed in an aggregated manner, using the
 * sweep-line algorithm. The input curves should be pairwise disjoint in their
 * interior and pairwise interior-disjoint from all existing arrangement
 * vertices and edges.
 *
 * \cgalHeading{Requirements}
 *
 * <UL>
 * <LI>The instantiated `Traits` class must model the
 * `ArrangementBasicTraits_2` concept, as no intersections are computed.
 * <LI>`InputIterator::value_type` must be `Traits::X_monotone_curve_2`
 * </UL>
 */
template <typename Traits, typename Dcel, InputIterator>
void insert_non_intersecting_curves(Arrangement_2<Traits,Dcel>& arr,
                                    InputIterator first, InputIterator last);

/*! \ingroup PkgArrangementOnSurface2Funcs
 *
 * Inserts a given point into a given arrangement.  It uses a given
 * point-location object to locate the given point in the given arrangement. If
 * the point conincides with an existing vertex, there is nothing left to do; if
 * it lies on an edge, the edge is split at the point. Otherwise, the point is
 * contained inside a face, and is inserted as an isolated vertex inside this
 * face.  By default, the function uses the "walk along line" point-location
 * strategy - namely, an instance of the class
 * `Arr_walk_along_line_point_location<Arrangement_2<Traits,Dcel> >`.  In either
 * case, the function returns a handle for the vertex associated with the point.
 *
 * \pre If provided, `pl` must be attached to the given arrangement `arr`.
 *
 * \cgalHeading{Requirements}
 *
 * <UL>
 * <LI>The instantiated `Traits` class must model the
 * `ArrangementXMonotoneTraits_2` concept. Not all expressions listed
 * by this concept are required. In fact the traits class must model the
 * `ArrangementBasicTraits_2` concept, and support the splitting functionality.
 * <LI>The point-location object `pl`, must model the
 * `ArrangementPointLocation_2` concept.
 * </UL>
 */
template<typename Traits, typename Dcel, typename PointLocation>
typename Arrangement_2<Traits,Dcel>::Vertex_handle
insert_point(Arrangement_2<Traits,Dcel>& arr,
             const typename Traits::Point_2& p,
             const PointLocation& pl = walk_pl);

/*! \ingroup PkgArrangementOnSurface2Funcs
 *
 * Checks the validity of a given arrangement.
 *
 * Invokes the member function `arr.is_valid()` to verify the topological
 * correctness of the arrangement. Then it performs additional validity
 * tests. It checks that all \f$ x\f$-monotone curves associated with
 * arrangement edges are pairwise disjoint in their interior. Then it makes sure
 * that all holes and all isolated vertices are located within the proper
 * arrangement faces. Note that the test carried out by this function may take a
 * considerable amount of time; it is recommended to be used only for debugging
 * purposes.
 *
 * \cgalHeading{Requirements}
 *
 * The instantiated traits class must model the concept
 * `ArranagmentXMonotoneTraits_2`.
 */
template<typename Traits, typename Dcel>
bool is_valid(const Arrangement_2<Traits, Dcel>& arr);

/*! \ingroup PkgArrangementOnSurface2Funcs
 *
 * Removes an edge given by one of the twin halfedges that forms it, from a
 * given arrangement. Once the edge is removed, if the vertices associated with
 * its endpoints become isolated, they are removed as well. The call
 * `remove_edge(arr, e)` is equivalent to the call `arr.remove_edge (e, true,
 * true)`. However, this free function requires that `Traits` be a model of the
 * refined concept `ArrangementXMonotoneTraits_2`, which requires merge
 * operations on \f$ x\f$-monotone curves. If one of the end-vertices of the
 * given edge becomes redundant after the edge is removed (see `remove_vertex()`
 * for the definition of a redundant vertex), it is removed, and its incident
 * edges are merged.  If the edge-removal operation causes two faces to merge,
 * the merged face is returned. Otherwise, the face to which the edge was
 * incident before the removal is returned.
 *
 * \cgalHeading{Requirements}
 *
 * <UL>
 * <LI>The instantiated traits class must model the concept
 * `ArrangementXMonotoneTraits_2`.
 * </UL>
 */
template <typename Traits, typename Dcel>
typename Arrangement_2<Traits,Dcel>::Face_handle
remove_edge(Arrangement_2<Traits,Dcel>& arr,
            typename Arrangement_2<Traits,Dcel>::Halfedge_handle e);

/*! \ingroup PkgArrangementOnSurface2Funcs
 *
 * Attempts to removed a given vertex from a given arrangement. The vertex can
 * be removed if it is either an isolated vertex, (and has no incident edge,) or
 * if it is a <I>redundant</I> vertex. That is, it has exactly two incident
 * edges, whose associated curves can be merged to form a single \f$
 * x\f$-monotone curve.  The function returns a boolean value that indicates
 * whether it succeeded removing the vertex from the arrangement.
 *
 * \cgalHeading{Requirements}
 *
 * <UL>
 * <LI>The instantiated `Traits` class must model the
 * `ArrangementXMonotoneTraits_2` concept. Not all expressions listed
 * by this concept are required. In fact the traits class must model the
 * `ArrangementBasicTraits_2` concept and support the merging functionality.
 * </UL>
 */
template <typename Traits, typename Dcel>
bool remove_vertex(Arrangement_2<Traits,Dcel>& arr,
                   typename Arrangement_2<Traits,Dcel>::Vertex_handle v);

/*! \ingroup PkgArrangementOnSurface2Funcs
 *
 * Compute the zone of the given \f$ x\f$-monotone curve in the existing
 * arrangement. Meaning, it output the arrangement's vertices, edges and faces
 * that the \f$ x\f$-monotone curve intersects. The order of the objects is the
 * order that they are discovered when traversing the \f$ x\f$-monotone curve
 * from left to right.
 *
 * A given point-location object is used for answering point-location queries
 * during the insertion process. By default, the function uses the "walk along
 * line" point-location strategy - namely an instance of the class
 * `Arr_walk_along_line_point_location<Arrangement_2<Traits,Dcel> >`.
 *
 * Compute the zone of the given \f$ x\f$-monotone curve `c` in the arrangement
 * `arr`.
 *
 * \pre If provided, `pl` must be attached to the given arrangement `arr`.
 *
 * \cgalHeading{Requirements}
 *
 * <UL>
 * <LI>The instantiated `GeomTraits` class must model the
 * `ArrangementXMonotoneTraits_2` concept.
 * <LI>The point-location object `pl`, must model the
 * `ArrangementPointLocation_2` concept.
 * </UL>
 */
template <typename Traits_2, typename Dcel,
          typename OutputIterator, typename PointLocation>
OutputIterator zone(Arrangement_2<Traits_2, Dcel>& arr,
                    const typename Traits_2::X_monotone_curve_2& c,
                    OutputIterator oi, const PointLocation& pl);

} /* namespace CGAL */
