
/*!
\ingroup PkgArrangement2ConceptsTraits
\cgalConcept

The concept `ArrangementTraits_2` allows the construction of arrangement 
of <I>general</I> planar curves. Models of this concept are used 
by the free \link PkgArrangement2Insert `CGAL::insert()` \endlink functions of the arrangement package and 
by the `CGAL::Arrangement_with_history_2` class. 

A model of this concept must define the nested `Curve_2` type, which 
represents a general planar curve that is not necessarily \f$ x\f$-monotone 
and is not necessarily connected. Such curves are eventually subdivided into 
\f$ x\f$-monotone subcurves and isolated points (represented by the `Point_2` 
and `X_monotone_curve_2` types, defined in the basic traits concept). 

A model of the concept `ArrangementTraits_2` that handles arbitrary curves, 
which are always \f$ x\f$-monotone, such as a traits class that handles 
linear curves may define the nested types `Curve_2` and 
`X_monotone_curve_2` to be of equivalent types. Moreover, defining 
them as of equivalent types is advantageous, as it enables a generic 
simple implementation of the nested Functor 
`Make_x_monotone_2`. 

On the other hand, a model of the `ArrangementTraits_2` concept that handles 
arbitrary curves, which may be not \f$ x\f$-monotone must define the 
`Curve_2` and `X_monotone_curve_2` nested types to be of 
different types to allow proper dispatching of the free functions 
that accept such curves, such as `intsert()`. 

\cgalRefines `ArrangementXMonotoneTraits_2` 

\cgalHasModel `CGAL::Arr_segment_traits_2<Kernel>`
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

\sa `ArrangementBasicTraits_2` 
\sa `ArrangementXMonotoneTraits_2` 
\sa `ArrangementLandmarkTraits_2` 

*/

class ArrangementTraits_2 {
public:

/// \name Types 
/// @{

/*!
models the concept `ArrTraits::Curve_2`. 
*/ 
typedef unspecified_type Curve_2; 

/// @} 

/// \name Functor Types 
/// @{

/*!
models the concept `ArrTraits::MakeXMonotone_2`. 
*/ 
typedef unspecified_type Make_x_monotone_2; 

/// @} 

/// \name Accessing Functor Objects 
/// @{

/*!

*/ 
Make_x_monotone_2 make_x_monotone_2_object() const; 

/// @}

}; /* end ArrangementTraits_2 */

