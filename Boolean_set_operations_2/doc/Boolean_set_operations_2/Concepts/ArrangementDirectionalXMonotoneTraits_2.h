
/*!
\ingroup PkgBooleanSetOperations2Concepts
\cgalConcept

This concept refines the basic arrangement \f$ x\f$-monotone traits concept.
A model of this concept is able to handle <I>directed</I> \f$ x\f$-monotone curves
that intersect in their interior. Namely, an instance of the
`X_monotone_curve_2` type defined by a model of the concept
`ArrangementXMonotoneTraits_2` is only required to have a <I>left</I>
(lexicographically smaller) endpoint and a <I>right</I> endpoint.
If the traits class is also a model of
`ArrangementDirectionalXMonotoneTraits_2`, the \f$ x\f$-monotone curve is
also required to have a direction, namely one of these two endpoint serves
as its <I>source</I> and the other as its <I>target</I>.

\cgalRefines{ArrangementXMonotoneTraits_2}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Arr_segment_traits_2<Kernel>}
\cgalHasModels{CGAL::Arr_non_caching_segment_traits_2<Kernel>}
\cgalHasModels{CGAL::Arr_circle_segment_traits_2<Kernel>}
\cgalHasModels{CGAL::Arr_conic_traits_2<RatKernel,AlgKernel,NtTraits>}
\cgalHasModels{CGAL::Arr_rational_function_traits_2<AlgKernel,NtTraits>}
\cgalHasModels{CGAL::Arr_Bezier_curve_traits_2<RatKernel,AlgKernel,NtTraits>}
\cgalHasModels{CGAL::Arr_algebraic_segment_traits_2<Coefficient>}
\cgalHasModelsEnd

\sa `ArrangementXMonotoneTraits_2`

*/

class ArrangementDirectionalXMonotoneTraits_2 {
public:

/// \name Functor Types
/// @{

/*!
models the concept `ArrDirectionalTraits::CompareEndpointsXy_2`.
*/
typedef unspecified_type Compare_endpoints_xy_2;

/*!
models the concept `ArrDirectionalTraits::ConstructOpposite_2`.
*/
typedef unspecified_type Construct_opposite_2;

/*!
models the concept `ArrDirectionalTraits::Intersect_2`.
*/
typedef unspecified_type Intersect_2;

/*!
models the concept `ArrDirectionalTraits::Split_2`.
*/
typedef unspecified_type Split_2;

/*!
models the concept `ArrDirectionalTraits::AreMergeable_2`.
*/
typedef unspecified_type Are_mergeable_2;

/*!
models the concept `ArrDirectionalTraits::Merge_2`.
*/
typedef unspecified_type Merge_2;

/// @}

/// \name Creation
/// @{

/*!
default constructor.
*/
ArrangementDirectionalXMonotoneTraits_2();

/*!
copy constructor
*/
ArrangementDirectionalXMonotoneTraits_2(ArrangementDirectionalXMonotoneTraits_2 other);

/*!
assignment operator.
*/
ArrangementDirectionalXMonotoneTraits_2 operator=(other);

/// @}

/// \name Accessing Functor Objects
/// @{

/*!

*/
Compare_endpoints_xy_2 compare_endpoints_xy_2_object();

/*!

*/
Construct_opposite_2 construct_opposite_2_object();

/// @}

}; /* end ArrangementDirectionalXMonotoneTraits_2 */

