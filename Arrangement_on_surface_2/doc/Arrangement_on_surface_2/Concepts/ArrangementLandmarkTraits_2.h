
/*!
\ingroup PkgArrangement2ConceptsTraits
\cgalConcept

The concept `ArrangementLandmarkTraits_2` refines the general traits concept by adding 
operations needed for the landmarks point-location strategy, namely - 
approximating points and connecting points with a simple \f$ x\f$-monotone 
curve. 

A model of this concept must define the `Approximate_number_type`, which 
is used to approximate the coordinates of `Point_2` instances. It is 
recommended to define the approximated number type as the built-in 
`double` type. 

\cgalRefines `ArrangementTraits_2` 

\cgalHasModel `CGAL::Arr_non_caching_segment_traits_2<Kernel>` 
\cgalHasModel `CGAL::Arr_segment_traits_2<Kernel>` 
\cgalHasModel `CGAL::Arr_linear_traits_2<Kernel>` 
\cgalHasModel `CGAL::Arr_polyline_traits_2<SegmentTraits>` 
\cgalHasModel `CGAL::Arr_conic_traits_2<RatKernel,AlgKernel,NtTraits>` 

\sa `ArrangementTraits_2` 

*/

class ArrangementLandmarkTraits_2 {
public:

/// \name Types 
/// @{

/*! 
the number type used to approximate point coordinates. 
*/ 
typedef Hidden_type Approximate_number_type; 

/// @} 

/// \name Functor Types 
/// @{

/*! 
models the concept `ArrTraits::Approximate_2`. 
*/ 
typedef Hidden_type Approximate_2; 

/*! 
models the concept `ArrTraits::ConstructXMonotoneCurve_2`. 
*/ 
typedef Hidden_type Construct_x_monotone_curve_2; 

/// @} 

/// \name Accessing Functor Objects 
/// @{

/*! 

*/ 
Approximate_2 approximate_2_object() const; 

/*! 

*/ 
Construct_x_monotone_curve_2 
construct_x_monotone_curve_2_object() const; 

/// @}

}; /* end ArrangementLandmarkTraits_2 */

