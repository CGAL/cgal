namespace CGAL {

/*! \ingroup PkgArrangementOnSurface2Ref
 *
 * \anchor arr_refarr_with_hist
 *
 * An object `arr` of the class `Arrangement_with_history_2` represents the
 * planar subdivision induced by a set of input curves \f$ \cal C\f$.  The
 * arrangement is represented as a doubly-connected edge-list (\dcel).  As is
 * the case for the `Arrangement_2<Traits,Dcel>`, each \dcel vertex is
 * associated with a point and each edge is associated with an \f$ x\f$-monotone
 * curve whose interior is disjoint from all other curves and points. Each such
 * \f$ x\f$-monotone curve is a subcurve of some \f$ C \in \cal C\f$, or may
 * represent an overlap among several curves in \f$ \cal C\f$.
 *
 * The `Arrangement_with_history_2` class-template extends the `Arrangement_2`
 * class-template by keeping an additional container of input curves
 * representing \f$ \cal C\f$, and by maintaining a cross-mapping between these
 * curves and the arrangement edges they induce. This way it is possible to
 * determine the inducing curve(s) of each arrangement edge. This mapping also
 * allows the traversal of input curves, and the traversal of edges induced by
 * each curve.
 *
 * The `Arrangement_with_history_2` template has two parameters:
 * <UL>
 * <LI>The `Traits` template-parameter should be substituted by a model of
 * the `ArrangementTraits_2` concept. The traits class defines the `Curve_2`
 * type, which represents an input curve.  It also defines the types of \f$
 * x\f$-monotone curves and two-dimensional points, namely
 * `ArrangementTraits_2::X_monotone_curve_2` and `ArrangementTraits_2::Point_2`,
 * respectively, and supports basic geometric predicates on them.
 * <LI>The `Dcel` template-parameter should be substituted by a class that is
 * a model of the `ArrangementDcelWithRebind` concept. The value of this
 * parameter is by default `Arr_default_dcel<Traits>`.
 * </UL>
 *
 * \sa `ArrangementDcel`
 * \sa `Arr_default_dcel<Traits>`
 * \sa `ArrangementTraits_2`
 * \sa `Arrangement_2<Traits,Dcel>`
 * \sa `insertion functions`
 * \sa `removal functions`
 * \sa `overlaying arrangements`
 */
template <typename Traits, typename Dcel>
class Arrangement_with_history_2 : public Arrangement_on_surface_with_history_2<Traits, typename Default_planar_topology<Traits, Dcel>::Traits> {
public:

  /// \name Types
  /// @{

  //! the geometry traits class.
  typedef Traits                                  Geometry_traits;

  //! The topology traits.
  typedef typename Default_planar_topology<Geometry_traits, Dcel>::Traits
    Topology_traits;

  //! The base arrangement on surface type.
  typedef Arrangement_on_surface_with_history_2<Geometry_traits, Topology_traits>
    Base;

  /// @}

  /// \name Types inherited from the base Arrangement_on_surface_2_with_history_2.
  /// @{

  typedef typename Base::Point_2                    Point_2;
  typedef typename Base::X_monotone_curve_2         X_monotone_curve_2;
  typedef typename Base::Curve_2                    Curve_2;
  typedef typename Base::Size                       Size;

  typedef typename Base::Curve_handle               Curve_handle;
  typedef typename Base::Curve_iterator             Curve_iterator;
  typedef typename Base::Induced_edge_iterator      Induced_edge_iterator;
  typedef typename Base::Originating_curve_iterator Originating_curve_iterator;

  /// @}

  /// \name Creation
  /// @{

  /*! constructs an empty arrangement containing one unbounded face, which
   * corresponds to the whole plane.
   */
  Arrangement_with_history_2<Traits, Dcel>();

  /*! copy constructor. */
  Arrangement_with_history_2<Traits, Dcel>(const Arrangement_with_history_2<Traits, Dcel>& other);

  /*! constructs an empty arrangement that uses the given `traits` instance for
   * performing the geometric predicates.
   */
  Arrangement_with_history_2<Traits, Dcel>(const Traits* traits);

  /// @}

  /// \name Assignment Methods
  /// @{

  /*! assignment operator. */
  Arrangement_with_history_2<Traits, Dcel>& operator= (other);

  /*! assigns the contents of another arrangement. */
  void assign(const Arrangement_with_history_2<Traits, Dcel>& other);

  /*! clears the arrangement. */
  void clear();

  /// @}

  /// \name Access Functions
  /// @{

  /*! obtains the traits object used by the arrangement instance.
   * A `const` version is also available.
   */
  Traits* traits();

  /// @}
}; /* end Arrangement_with_history_2 */

/*! \ingroup PkgArrangementOnSurface2Insert
 *
 * Inserts the given curve `c` into the arrangement with history `arr`, and
 * returns a handle to the inserted curve. `c` is subdivided into \f$
 * x\f$-monotone subcurves (and perhaps isolated points). Each subcurve is in
 * turn inserted into the arrangement by locating its left endpoint and
 * computing its zone until reaching the right endpoint.
 *
 * The given point-location object `pl` is used to locate the left endpoints of
 * the \f$ x\f$-monotone curves. By default, the function uses the "walk along
 * line" point-location strategy - namely an instance of the class
 * `Arr_walk_along_line_point_location<Arrangement_2<Traits,Dcel> >`.
 *
 * \pre If provided, `pl` is attached to the given arrangement `arr`.
 */
template<typename Traits, typename Dcel, typename PointLocation>
typename Arrangement_with_history_2<Traits, Dcel>::Curve_handle
insert(Arrangement_with_history_2<Traits, Dcel>& arr,
       const typename Traits::Curve_2& c,
       const PointLocation& pl = walk_pl);

/*! \ingroup PkgArrangementOnSurface2Insert
 * Aggregately inserts the curves in the range `[first,last)` into the
 * arrangement with history `arr` using the sweep-line framework.
 * \param arr the target arrangement with history.
 * \param first the iterator to the first element in the range of curves.
 * \param last the past-the-end iterator of the range of curves.
 */
template <typename Traits, typename Dcel, typename InputIterator>
void insert(Arrangement_with_history_2<Traits,Dcel>& arr,
            InputIterator first, InputIterator last);

/*! \ingroup PkgArrangementOnSurface2Funcs
 * Removes a given curvespecified by its handle `ch`, from a given arrangement
 * `arr`, deleting all the edges it induces. The function returns the number
 * of deleted edges.
 */
template <typename Traits, typename Dcel>
Size remove_curve(Arrangement_with_history_2<Traits,Dcel>& arr,
                  typename Arrangement_with_history_2<Traits,Dcel>::Curve_handle ch);


/*! \addtogroup PkgArrangementOnSurface2Overlay
 * Computes the overlay of two arrangements with history `arr1` and `arr2`, and
 * sets the output arrangement with history `res` to represent the overlaid
 * arrangement. The function also constructs a consolidated set of curves that
 * induce `res`.
 *
 * \pre `res` does not refer to either `arr1` or `arr2`.
 */
template<typename Traits, typename Dcel1, typename Dcel2,
         typename ResDcel, typename OverlayTraits>
void overlay(const Arrangement_with_history_2<Traits,Dcel1>& arr1,
             const Arrangement_with_history_2<Traits,Dcel2>& arr2,
             Arrangement_with_history_2<Traits,ResDcel>& res,
             OverlayTraits& ovl_tr);


/*! \addtogroup PkgArrangementOnSurface2Overlay
 *
 * Computes the (simple) overlay of two arrangements with history `arr1` and
 * `arr2`, and sets the output arrangement with history `res` to represent the
 * overlaid arrangement. The function also constructs a consolidated set of
 * curves that induce `res`. It employs the default overlay-traits, which
 * practically does nothing.
 *
 * \pre `res` does not refer to either `arr1` or `arr2` (that is, "self overlay"
 * is not supported).
 */
template<typename Traits, typename Dcel1, typename Dcel2,
         typename ResDcel>
void overlay(const Arrangement_with_history_2<Traits,Dcel1>& arr1,
             const Arrangement_with_history_2<Traits,Dcel2>& arr2,
             Arrangement_with_history_2<Traits,ResDcel>& res);

} /* namespace CGAL */
