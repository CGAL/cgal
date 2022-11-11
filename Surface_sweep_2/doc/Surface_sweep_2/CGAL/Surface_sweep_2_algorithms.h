namespace CGAL {


/*!
\ingroup PkgSurfaceSweep2Ref

Given a range of curves, compute all intersection points between two (or more)
input curves. When the flag `report_endpoints` is `true`,
this function reports all the curve endpoints as well. If a curve
endpoint is also an intersection point, it is reported once (regardless
of the value of the `report_endpoints` flag). The value-type of
`InputIterator` is a curve type and the value-type of `OutputIterator`
is a point type. The output points are reported in an increasing
\f$xy\f$-lexicographical order.
*/
template <class InputIterator, class OutputIterator>
OutputIterator compute_intersection_points (InputIterator curves_begin,
                                            InputIterator curves_end,
                                            OutputIterator points,
                                            bool report_endpoints = false);

/*!
\ingroup PkgSurfaceSweep2Ref

Given a range of curves, compute all intersection points between two (or more)
input curves. When the flag `report_endpoints` is `true`,
this function reports all the curve endpoints as well. If a curve
endpoint is also an intersection point, it is reported once (regardless
of the value of the `report_endpoints` flag). The `Traits` type
must be a model of the `ArrangementTraits_2` concept, such that the
value-type of `InputIterator` is `Traits::Curve_2`, and the
value-type of `OutputIterator` is `Traits::Point_2`.
The output points are reported in an increasing \f$ xy\f$-lexicographical order.
*/
template <class InputIterator, class OutputIterator, class Traits>
OutputIterator compute_intersection_points (InputIterator curves_begin,
                                            InputIterator curves_end,
                                            OutputIterator points,
                                            bool report_endpoints = false,
                                            Traits traits);

/*!
\ingroup PkgSurfaceSweep2Ref

Given a range of curves, compute all \f$ x\f$-monotone subcurves that are pairwise
disjoint in their interior, as induced by the input curves.
If the flag `multiple_overlaps` is `true`, then a subcurve that
represents an overlap of \f$ k\f$ input curves is reported \f$ k\f$ times; otherwise,
each subcurve is reported only once.
The value-type of `InputIterator` is a curve type, and the value-type of
`OutputIterator` is an \f$x\f$-monotone curve type.
*/
template <class InputIterator, class OutputIterator>
OutputIterator compute_subcurves (InputIterator curves_begin,
                                  InputIterator curves_end,
                                  OutputIterator subcurves,
                                  bool multiple_overlaps = false);

/*!
\ingroup PkgSurfaceSweep2Ref

Given a range of curves, compute all \f$ x\f$-monotone subcurves that are pairwise
disjoint in their interior, as induced by the input curves.
If the flag `multiple_overlaps` is `true`, then a subcurve that
represents an overlap of \f$ k\f$ input curves is reported \f$ k\f$ times; otherwise,
each subcurve is reported only once. The `Traits` type must be a model
of the `ArrangementTraits_2` concept, such that the value-type of
`InputIterator` is `Traits::Curve_2`, and the value-type of
`OutputIterator` is `Traits::X_monotone_curve_2`.
*/
template <class InputIterator, class OutputIterator, class Traits>
OutputIterator compute_subcurves (InputIterator curves_begin,
                                  InputIterator curves_end,
                                  OutputIterator subcurves,
                                  bool multiple_overlaps = false,
                                  Traits traits = Default_traits());

/*!
\ingroup PkgSurfaceSweep2Ref

Given a range of curves, check whether there is at least one pair of curves
that intersect in their interior. The function returns `true` if such
a pair is found, and `false` if all curves are pairwise disjoint in
their interior. The value-type of `InputIterator` is a curve type.
*/
template <class InputIterator>
bool do_curves_intersect (InputIterator curves_begin,
                          InputIterator curves_end);

/*!
\ingroup PkgSurfaceSweep2Ref

Given a range of curves, check whether there is at least one pair of curves
that intersect in their interior. The function returns `true` if such
a pair is found, and `false` if all curves are pairwise disjoint in
their interior. The `Traits` type must be a model
of the `ArrangementTraits_2` concept, such that the value-type of
`InputIterator` is `Traits::Curve_2`.
*/
template <class InputIterator, class Traits>
bool do_curves_intersect (InputIterator curves_begin,
                          InputIterator curves_end,
                          Traits traits = Default_traits());

} /* namespace CGAL */

