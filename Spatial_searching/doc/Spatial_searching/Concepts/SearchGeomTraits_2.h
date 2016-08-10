/*!
\ingroup PkgSpatialSearchingDConcepts
\cgalConcept

The concept `SearchGeomTraits_2` defines the requirements for the template 
parameter of the search traits classes.

\cgalHasModel Any \cgal Kernel
*/

class SearchGeomTraits_2 {
public:

/// \name Types 
/// @{

/*!
Point type. `CGAL::Kernel_traits` has to be 
specialized for this type. 
*/ 
typedef unspecified_type Point_2; 

/*!
The number type of the %Cartesian coordinates of types `Point_2` 
*/ 
typedef unspecified_type FT; 
/*!
Iso box type.
*/ 
typedef unspecified_type Iso_rectangle_2;
/*!
Circle type.
*/
typedef unspecified_type Circle_2;

/*!
Functor model of `Kernel::ConstructMinVertex_2`
*/
typedef unspecified_type Construct_min_vertex_2;

/*!
Functor model of `Kernel::ConstructMaxVertex_2`
*/
typedef unspecified_type Construct_max_vertex_2;

/*!
Functor model of `Kernel::ConstructCenter_2`
*/
typedef unspecified_type Construct_center_2;

/*!
Functor model of `Kernel::ComputeSquaredRadius_2 `
*/
typedef unspecified_type Compute_squared_radius_2;

/*!
Functor model of `Kernel::ConstructIsoRectangle_2 `
*/
typedef unspecified_type Construct_iso_rectangle_2;

/*!
A random access iterator type to enumerate the 
%Cartesian coordinates of a point. 
*/ 
typedef unspecified_type Cartesian_const_iterator_2;

/*!
Functor with operators to construct iterators on the 
first and the past-the-end iterator for the %Cartesian coordinates of a point. This functor must 
provide the type `result_type` that must be the same a `Cartesian_const_iterator_2`. 
*/ 
typedef unspecified_type Construct_cartesian_const_iterator_2; 

/// @} 



}; /* end SearchTraits */
