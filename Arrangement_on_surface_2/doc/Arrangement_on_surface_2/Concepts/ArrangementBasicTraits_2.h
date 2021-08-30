
/*!
\ingroup PkgArrangementOnSurface2ConceptsTraits
\cgalConcept

The concept `ArrangementBasicTraits_2` defines the minimal set of geometric
predicates needed for the construction and maintenance of objects of the
class `Arrangement_2`, as well as performing simple queries (such as
point-location queries) on such arrangements.

A model of this concept must define nested `Point_2` and
`X_monotone_curve_2` types, which represent planar points and
continuous \f$ x\f$-monotone curves (a vertical segment is also considered to be
<I>weakly</I> \f$ x\f$-monotone), respectively. The \f$ x\f$-monotone curves are assumed
to be pairwise disjoint in their interiors, so they do not intersect
except at their endpoints.

The `X_monotone_curve_2` curves of an arrangement are confined to an
iso-rectangular area called the parameter space. The iso-rectangule can
be unbounded, open, or closed. The set of predicates provided by a model
the concept `ArrangementBasicTraits_2` is sufficient for constructing arrangements of
\f$ x\f$-monotone curves that do not reach or approach the boundary of the
parameter space. The nature of the input curves, whether they are
expected to reach or approach the left, right, bottom, or top side of the
boundary of the parameter space, are conveyed through the definition of
four additional nested types, namely `Left_side_category`,
`Right_side_category`, `Bottom_side_category`, and
`Top_side_category`. Each such type must be convertible to the type
`Arr_oblivious_side_tag`.

\cgalRefines DefaultConstructible
\cgalRefines CopyConstructible
\cgalRefines Assignable

\cgalHasModel `CGAL::Arr_segment_traits_2<Kernel>`
\cgalHasModel `CGAL::Arr_non_caching_segment_basic_traits_2<Kernel>`
\cgalHasModel `CGAL::Arr_non_caching_segment_traits_2<Kernel>`
\cgalHasModel `CGAL::Arr_linear_traits_2<Kernel>`
\cgalHasModel `CGAL::Arr_polyline_traits_2<SegmentTraits>`
\cgalHasModel `CGAL::Arr_circle_segment_traits_2<Kernel>`
\cgalHasModel `CGAL::Arr_line_arc_traits_2<CircularKernel>`
\cgalHasModel `CGAL::Arr_circular_arc_traits_2<CircularKernel>`
\cgalHasModel `CGAL::Arr_circular_line_arc_traits_2<CircularKernel>`
\cgalHasModel `CGAL::Arr_conic_traits_2<RatKernel,AlgKernel,NtTraits>`
\cgalHasModel `CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1>`
\cgalHasModel `CGAL::Arr_Bezier_curve_traits_2<RatKernel,AlgKernel,NtTraits>`
\cgalHasModel `CGAL::Arr_algebraic_segment_traits_2<Coefficient>`
\cgalHasModel `CGAL::Arr_curve_data_traits_2<Tr,XData,Mrg,CData,Cnv>`
\cgalHasModel `CGAL::Arr_consolidated_curve_data_traits_2<Traits,Data>`

*/

class ArrangementBasicTraits_2 {
public:

/// \name Types
/// @{

/*!
models the concept `ArrTraits::Point_2`.
*/
typedef unspecified_type Point_2;

/*!
models the concept `ArrTraits::XMonotoneCurve_2`.
*/
typedef unspecified_type X_monotone_curve_2;

/// @}

/// \name Categories
/// @{

/*!
indicates whether the nested functor `Compare_at_x_left_2` is
provided.
*/
typedef unspecified_type Has_left_category;

/*!
Must be convertible to `Arr_oblivious_side_tag`.
*/
typedef unspecified_type Left_side_category;

/*!
Must be convertible to `Arr_oblivious_side_tag`.
*/
typedef unspecified_type Bottom_side_category;

/*!
Must be convertible to `Arr_oblivious_side_tag`.
*/
typedef unspecified_type Top_side_category;

/*!
Must be convertible to `Arr_oblivious_side_tag`.
*/
typedef unspecified_type Right_side_category;

/// @}

/// \name Functor Types
/// @{

/*!
models the concept `ArrTraits::CompareX_2`.
*/
typedef unspecified_type Compare_x_2;

/*!
models the concept `ArrTraits::CompareXy_2`.
*/
typedef unspecified_type Compare_xy_2;

/*!
models the concept `ArrTraits::ConstructMinVertex_2`.
*/
typedef unspecified_type Construct_min_vertex_2;

/*!
models the concept `ArrTraits::ConstructMaxVertex_2`.
*/
typedef unspecified_type Construct_max_vertex_2;

/*!
models the concept `ArrTraits::IsVertical_2`.
*/
typedef unspecified_type Is_vertical_2;

/*!
models the concept `ArrTraits::CompareYAtX_2`.
*/
typedef unspecified_type Compare_y_at_x_2;

/*!
models the concept `ArrTraits::CompareYAtXLeft_2`.
Required only if the `Has_left_category` category is convertible to
`Tag_true`.
*/
typedef unspecified_type Compare_y_at_x_left_2;

/*!
models the concept `ArrTraits::CompareYAtXRight_2`.
*/
typedef unspecified_type Compare_y_at_x_right_2;

/*!
models the concept `ArrTraits::Equal_2`.
*/
typedef unspecified_type Equal_2;

/// @}

/// \name Accessing Functor Objects
/// @{

/*!

*/
Compare_x_2 compare_x_2_object() const;

/*!

*/
Compare_xy_2 compare_xy_2_object() const;

/*!

*/
Construct_min_vertex_2 construct_min_vertex_2_object() const;

/*!

*/
Construct_max_vertex_2 construct_max_vertex_2_object() const;

/*!

*/
Is_vertical_2 is_vertical_2_object() const;

/*!

*/
Compare_y_at_x_2 compare_y_at_x_2_object() const;

/*!

*/
Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const;

/*!

*/
Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const;

/*!

*/
Equal_2 equal_2_object() const;

/// @}

}; /* end ArrangementBasicTraits_2 */

