namespace CGAL {

/*! \ingroup PkgSurfaceSweep2Ref
 *
 * Given a range of curves, compute all intersection points between two (or
 * more) input curves. When the flag `report_endpoints` is `true`, this function
 * reports all the curve endpoints as well. If a curve endpoint is also an
 * intersection point, it is reported once (regardless of the value of the
 * `report_endpoints` flag). The value-type of `InputIterator` is a curve type
 * and the value-type of `OutputIterator` is a point type. The output points are
 * reported in an increasing \f$xy\f$-lexicographical order.
 */
template <typename InputIterator, typename OutputIterator>
OutputIterator compute_intersection_points(InputIterator curves_begin,
                                           InputIterator curves_end,
                                           OutputIterator points,
                                           bool report_endpoints = false);

/*! \ingroup PkgSurfaceSweep2Ref
 *
 * Given a range of curves, compute all intersection points between two (or
 * more) input curves. When the flag `report_endpoints` is `true`, this function
 * reports all the curve endpoints as well. If a curve endpoint is also an
 * intersection point, it is reported once (regardless of the value of the
 * `report_endpoints` flag). The `Traits` type must be a model of the
 * `AosTraits_2` concept, such that the value-type of `InputIterator` is
 * `Traits::Curve_2`, and the value-type of `OutputIterator` is
 * `Traits::Point_2`.  The output points are reported in an increasing
 * \f$xy\f$-lexicographical order.
*/
template <typename InputIterator, typename OutputIterator, class Traits>
OutputIterator compute_intersection_points(InputIterator curves_begin,
                                           InputIterator curves_end,
                                           OutputIterator points,
                                           bool report_endpoints = false,
                                           Traits traits);

/*! \ingroup PkgSurfaceSweep2Ref
 *
 * Given a range of curves, compute all \f$ x\f$-monotone subcurves that are
 * pairwise disjoint in their interior, as induced by the input curves.  If the
 * flag `multiple_overlaps` is `true`, then a subcurve that represents an
 * overlap of \f$ k\f$ input curves is reported \f$ k\f$ times; otherwise, each
 * subcurve is reported only once.  The value-type of `InputIterator` is a curve
 * type, and the value-type of `OutputIterator` is an \f$x\f$-monotone curve
 * type.
 */
template <typename InputIterator, typename OutputIterator>
OutputIterator compute_subcurves(InputIterator curves_begin,
                                 InputIterator curves_end,
                                 OutputIterator subcurves,
                                 bool multiple_overlaps = false);

/*! \ingroup PkgSurfaceSweep2Ref
 *
 * Given a range of curves, compute all \f$ x\f$-monotone subcurves that are
 * pairwise disjoint in their interior, as induced by the input curves.  If the
 * flag `multiple_overlaps` is `true`, then a subcurve that represents an
 * overlap of \f$ k\f$ input curves is reported \f$ k\f$ times; otherwise, each
 * subcurve is reported only once. The `Traits` type must be a model of the
 * `AosTraits_2` concept, such that the value-type of `InputIterator` is
 * `Traits::Curve_2`, and the value-type of `OutputIterator` is
 * `Traits::X_monotone_curve_2`.
*/
template <typename InputIterator, typename OutputIterator, typename Traits>
OutputIterator compute_subcurves(InputIterator curves_begin,
                                 InputIterator curves_end,
                                 OutputIterator subcurves,
                                 bool multiple_overlaps = false,
                                 Traits traits = Default_traits());

/*! \deprecated Call Surface_sweep_2::do_intersect(curves_begin, curves_end, false) instead
 * \ingroup PkgSurfaceSweep2Ref
 *
 * Given a range of curves, check whether there is at least one pair of curves
 * that intersect in their interior. The function returns `true` if such a pair
 * is found, and `false` if all curves are pairwise disjoint in their
 * interior. The value-type of `InputIterator` is a curve type.
 *
 * The call `CGAL::do_curves_intersect(begin, end, traits)` should be replaced by the call
 * `CGAL::Surface_sweep_2::do_intersect(begin, end, false, traits)`
 */
template <typename InputIterator>
bool do_curves_intersect(InputIterator curves_begin,
                         InputIterator curves_end);

/*! \deprecated Call Surface_sweep_2::do_intersect(curves_begin, curves_end, false, traits) instead
 * \ingroup PkgSurfaceSweep2Ref
 *
 * Given a range of curves, check whether there is at least one pair of curves
 * that intersect in their interior. The function returns `true` if such a pair
 * is found, and `false` if all curves are pairwise disjoint in their
 * interior. The `Traits` type must be a model of the `AosTraits_2` concept,
 * such that the value-type of `InputIterator` is `Traits::Curve_2`.
 *
 * The call `CGAL::do_curves_intersect(begin, end)` should be replaced by the call
 * `CGAL::Surface_sweep_2::do_intersect(begin, end, false)`
 */
template <typename InputIterator, typename Traits>
bool do_curves_intersect(InputIterator curves_begin,
                         InputIterator curves_end,
                         Traits traits = Default_traits());

namespace Surface_sweep_2 {

/*! \ingroup PkgSurfaceSweep2Ref
 *
 * Given a range of curves, check whether there is at least one pair of curves
 * that intersect. The value-type of `InputIterator` is a curve type. If closed
 * is `true`, the curves are considered close. Otherwise, they are considered
 * open. If the curve are considered close, the function returns `true` if the
 * curves in the range pairwise intersect, and false otherwise. If the curves
 * are considered open, the function returns `true` if the curves in the range
 * pairwise intersect in their interiros, and false otherwise.
 */
template <typename InputIterator>
bool do_intersect(InputIterator curves_begin,
                  InputIterator curves_end,
                  bool closed = true);

/*! \ingroup PkgSurfaceSweep2Ref
 *
 * Given a range of curves, check whether there is at least one pair of curves
 * that intersect. The value-type of `InputIterator` is a curve type. If closed
 * is `true`, the curves are considered close. Otherwise, they are considered
 * open. If the curve are considered close, the function returns `true` if the
 * curves in the range pairwise intersect, and false otherwise. If the curves
 * are considered open, the function returns `true` if the curves in the range
 * pairwise intersect in their interiros, and false otherwise. The `Traits` type
 * must be a model of the `AosTraits_2` concept, such that the value-type of
 * `InputIterator` is either `Traits::Curve_2` or `Traits::X_monotone_curve_2`.
 */
template <typename InputIterator, typename Traits>
bool do_intersect(InputIterator curves_begin,
                  InputIterator curves_end,
                  bool closed = true,
                  Traits traits = Default_traits());

} // namespace Surface_sweep_2
} // namespace CGAL
