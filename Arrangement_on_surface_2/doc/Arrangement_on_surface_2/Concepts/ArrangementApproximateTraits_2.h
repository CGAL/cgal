/*!
\ingroup PkgArrangement2ConceptsTraits
\cgalConcept

The concept `ArrangementApproximateTraits_2` refines the basic traits concept
`ArrangementBasicTraits_2`. A model of this concept is able to approximate a
point.

\cgalRefines `ArrangementBasicTraits_2`

\cgalHasModel `CGAL::Arr_conic_traits_2<RatKernel,AlgKernel,NtTraits>`
\cond \cgalHasModel `CGAL::Arr_geodesic_arc_on_sphere_traits_2` \endcond
\cgalHasModel `CGAL::Arr_linear_traits_2<Kernel>`
\cgalHasModel `CGAL::Arr_non_caching_segment_traits_2<Kernel>`
\cgalHasModel `CGAL::Arr_segment_traits_2<Kernel>`
\cgalHasModel `CGAL::Arr_polycurve_traits_2<GeometryTraits_2>`
\cgalHasModel `CGAL::Arr_polyline_traits_2<SegmentTraits_2>`
\cgalHasModel `CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1>`

\sa ArrangementConstructXMonotoneCurveTraits_2`,
    ArrangementXMonotoneTraits_2`, and
    `ArrangementTraits_2`
*/
class ArrangementApproximateTraits_2 {
public:
  /// \name Types
  /// @{

  /*!
   * the number type used to approximate point coordinates, e.g., double.
   */
  typedef unspecified_type Approximate_number_type;

  /// @}

  /// \name Functor Types
  /// @{

  /*!
   * models the concept `ArrTraits::Approximate_2`.
   */
  typedef unspecified_type Approximate_2;

  /// \name Accessing Functor Objects
  /// @{

  /*!
   *
   */
  Approximate_2 approximate_2_object() const;

  /// @}
}
