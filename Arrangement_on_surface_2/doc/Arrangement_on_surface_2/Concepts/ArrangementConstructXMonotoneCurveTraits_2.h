/*!
\ingroup PkgArrangementOnSurface2ConceptsTraits
\cgalConcept

The concept `ArrangementConstructXMonotoneCurveTraits_2` refines the basic
traits concept `ArrangementBasicTraits_2`. A model of this concept is able to
construct an \f$ x\f$-monotone curve from two points.

\cgalRefines `ArrangementBasicTraits_2`

\cgalHasModel `CGAL::Arr_conic_traits_2<RatKernel,AlgKernel,NtTraits>`
\cond \cgalHasModel `CGAL::Arr_geodesic_arc_on_sphere_traits_2` \endcond
\cgalHasModel `CGAL::Arr_linear_traits_2<Kernel>`
\cgalHasModel `CGAL::Arr_non_caching_segment_traits_2<Kernel>`
\cgalHasModel `CGAL::Arr_segment_traits_2<Kernel>`
\cgalHasModel `CGAL::Arr_polyline_traits_2<SegmentTraits_2>`
\cgalHasModel `CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1>`

\sa `ArrangementApproximateTraits_2`,
    `ArrangementXMonotoneTraits_2`, and
    `ArrangementTraits_2`
*/
class ArrangementConstructXMonotoneCurveTraits_2 {
public:
  typedef ArrTraits::X_monotone_curve_2                 X_monotone_curve_2;
  typedef ArrTraits::Point_2                            Point_2;

  /*!
   * returns an \f$ x\f$-monotone curve connecting `p1` and `p2` (i.e., the
   * two input points are its endpoints).
   */
  X_monotone_curve_2 operator() (Point_2 p1, Point_2 p2);
}
