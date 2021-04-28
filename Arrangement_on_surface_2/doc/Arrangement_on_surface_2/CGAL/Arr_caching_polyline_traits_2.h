/// \ingroup PkgArrangementOnSurface2Macros
/// @{
/*!
 * If the macro is set to one, then \f$x\f$-monotone curves are always
 * directed from left-to-right.
 */
#define CGAL_ALWAYS_LEFT_TO_RIGHT
/// @}

namespace CGAL {

/*! \ingroup PkgArrangementOnSurface2TraitsClasses
 *
 * The traits class `Arr_caching_polyline_traits_2` handles piecewise linear
 * curves, commonly referred to as polylines. Each polyline is a chain of
 * segments, where each two neighboring segments in the chain share a common
 * endpoint; that is, the polyline is continuous. Furthermore, the target of the
 * \f$i\f$th segement of a polyline has to coincide with the source of the
 * \f$i+1\f$st segment; that is, the polyline has to be \a well-oriented.
 *
 * The traits class template uses an internal type to represent polylines, which
 * in turn uses an internal type to represent segments that compose polyline.
 *
 * \tparam Kernel a type that represents a geometric kernel.
 * The number type used by the injected kernel should support exact
 * rational arithmetic (that is, the number type should support the arithmetic
 * operations \f$ +\f$, \f$ -\f$, \f$ \times\f$ and \f$ \div\f$ carried out
 * without loss of precision), in order to avoid robustness problems, although
 * other inexact number types could be used at the user's own risk.
 *
 * \tparam Range a type that represents a valid range of points.
 *
 * A polyline that comprises \f$n > 0\f$ segments has \f$ n+1 \f$ points, and
 * they are represented as objects of type `Kernel::Point_2`. Since the
 * notion of a \a vertex is reserved to 0-dimensional elements of an
 * arrangement, we use, in this context, the notion of \a points in order to
 * refer to the vertices of a polyline. For example, an arrangement induced by a
 * single non-self intersecting polyline has exactly two vertices regardless of
 * the number of points.
 *
 * \cgalModels `ArrangementTraits_2`
 * \cgalModels `ArrangementDirectionalXMonotoneTraits_2`
 * \cgalModels `ArrangementConstructXMonotoneCurveTraits_2`
 * \cgalModels `ArrangementApproximateTraits_2`
 *
 * \sa `Arr_polyline_traits_2<SubcurveTraits_2>`
 */

template <typename Kernel, typename Range>
class Arr_caching_polyline_traits_2 :
    public Arr_polycurve_traits_2<SegmentTraits_2> {
public:

  /// \name Types
  /// @{
  /*!
   */
  typedef typename Kernel::Point_2                                Point_2;

  /*!
   */
  typedef internal::X_monotone_caching_polyline_2<Kernel, Range>
    X_monotone_curve_2;

  /*!
   */
  typedef internal::X_monotone_caching_polyline_2<Kernel, Range>  Curve_2;
  typedef
  /// @}

  /*! Construction functor of a polyline. Its `operator()` is oveloaded to
   * support various input types.
   */
  class Construct_curve_2 {
  public:
    /// \name Operations
    /// @{

    /*! Obtain a polyline connecting the two given endpoints.
     * \param p the first point.
     * \param q the second point.
     * \pre `p` and `q` are distinct.
     * \return a polyline that consists of a single segment connecting `p` and
     *         `q`.
     */
    Curve_2 operator()(const Point_2& p, const Point_2& q) const;

    /*! Construct a well-oriented polyline from a range of points.
     *
     * \param range a range of points
     * \param duplicate_first indicates whether the first point should be
     *        duplicated to construct a closed polyline.
     * \pre contains no duplicated points.
     * \return a polyline that consists of the given range of points.
     */
    template <typename Range>
    Curve_2 operator()(const Range& range, bool duplicate_first = false) const;

    /// @} /* end of operations */
  }; /* end of Arr_caching_polyline_traits_2::Construct_curve_2 */

  /*! Construction functor of a \f$x\f$-monotone polyline. Its `operator()` is
   * oveloaded to support various input types.
   */
  class Construct_x_monotone_curve_2 {
  };
  public:
    /// \name Operations
    /// @{

    /*! Obtain a polyline connecting the two given endpoints.
     * \param p the first point.
     * \param q the second point.
     * \pre `p` and `q` are distinct.
     * \return a polyline that consists of a single segment connecting `p` and
     *         `q`.
     */
    X_monotone_curve_2 operator()(const Point_2& p, const Point_2& q) const;

    /*! Construct a well-oriented polyline from a range of points.
     *
     * \param range a range of points.
     * \param duplicate_first indicates whether the first point should be
     *        duplicated to construct a closed polyline.
     * \pre contains no duplicated points.
     * \return a polyline that consists of the given range of points.
     */
    template <typename Range>
    X_monotone_curve_2 operator()(const Range& range,
                                  bool duplicate_first = false) const;

    /// @} /* end of operations */

  /// \name Accessing Functor Objects
  /// @{

  /*!
   */
  Construct_curve_2 construct_curve_2_object() const;

  /*!
   */
  Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object() const;

  /// @} /* End Accessing Functor Objects */

}; /* end Arr_caching_polyline_traits_2 */

} /* end namespace CGAL */
