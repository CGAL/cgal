
/*!
\ingroup PkgArrangementOnSurface2ConceptsTraits
\cgalConcept

The concept `ArrangementXMonotoneTraits_2` refines the basic arrangement-traits concept.
A model of this concept is able to handle \f$ x\f$-monotone curves that
intersect in their interior (and points that coincide with curve
interiors). This is necessary for constructing arrangements of sets of
intersecting \f$ x\f$-monotone curves.

As the resulting structure, represented by the `Arrangement_2` class,
stores pairwise interior-disjoint curves, the input curves are split at
the intersection points before being inserted into the arrangement.
A model of this refined concept therefore needs to compute the intersections
(and possibly overlaps) between two \f$ x\f$-monotone curves and to support
curve splitting.

\cgalRefines{ArrangementBasicTraits_2}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Arr_segment_traits_2<Kernel>}
\cgalHasModels{CGAL::Arr_non_caching_segment_traits_2<Kernel>}
\cgalHasModels{CGAL::Arr_linear_traits_2<Kernel>}
\cgalHasModels{CGAL::Arr_polyline_traits_2<SegmentTraits>}
\cgalHasModels{CGAL::Arr_circle_segment_traits_2<Kernel>}
\cgalHasModels{CGAL::Arr_line_arc_traits_2<CircularKernel>}
\cgalHasModels{CGAL::Arr_circular_arc_traits_2<CircularKernel>}
\cgalHasModels{CGAL::Arr_circular_line_arc_traits_2<CircularKernel>}
\cgalHasModels{CGAL::Arr_conic_traits_2<RatKernel,AlgKernel,NtTraits>}
\cgalHasModels{CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1>}
\cgalHasModels{CGAL::Arr_Bezier_curve_traits_2<RatKernel,AlgKernel,NtTraits>}
\cgalHasModels{CGAL::Arr_algebraic_segment_traits_2<Coefficient>}
\cgalHasModels{CGAL::Arr_curve_data_traits_2<Tr,XData,Mrg,CData,Cnv>}
\cgalHasModels{CGAL::Arr_consolidated_curve_data_traits_2<Traits,Data>}
\cgalHasModelsEnd

\sa `ArrangementBasicTraits_2`

*/

class ArrangementXMonotoneTraits_2 {
public:

/// \name Types
/// @{

/*!
the multiplicity type.
*/
typedef unspecified_type Multiplicity;

/// @}

/// \name Tags
/// @{

/*!
indicates whether the nested functors `Are_mergeable_2` and
`Merge_2` are provided.
*/
typedef unspecified_type Has_merge_category;

/// @}

/// \name Functor Types
/// @{

/*!
models the concept `ArrTraits::Intersect_2`.
*/
typedef unspecified_type Intersect_2;

/*!
models the concept `ArrTraits::Split_2`.
*/
typedef unspecified_type Split_2;

/// @}

/// \name
/// \attention The two following function-object types are
/// optional. If they are supported, the `Has_merge_category` tag
/// should be defined as `Tag_true` and otherwise as `Tag_false`.
/// @{

/*!
models the concept `ArrTraits::AreMergeable_2`.



*/
typedef unspecified_type Are_mergeable_2;

/*!
models the concept `ArrTraits::Merge_2`.
*/
typedef unspecified_type Merge_2;

/// @}

/// \name Accessing Functor Objects
/// @{

/*!

*/
Intersect_2 intersect_2_object() const;

/*!

*/
Split_2 split_2_object() const;

/// @}

/// \name
/// The two following methods are optional. If they are supported, the
/// `Has_merge_category` tag should be defined as `Tag_true` and otherwise
/// as `Tag_false`.
/// @{

/*!

*/
Are_mergeable_2 are_mergeable_2_object() const;

/*!

*/
Merge_2 merge_2_object() const;

/// @}

}; /* end ArrangementXMonotoneTraits_2 */

