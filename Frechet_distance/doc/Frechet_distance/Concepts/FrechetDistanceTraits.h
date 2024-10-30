/*!
\ingroup PkgFrechetDistanceConcepts
\cgalConcept

The concept `FrechetDistanceTraits` defines the requirements of the
first template parameter of the functions `CGAL::is_Frechet_distance_larger()`
and `CGAL::approximate_Frechet_distance()`.


\cgalHasModelsBegin
\cgalHasModels{CGAL::Frechet_distance_traits_2}
\cgalHasModels{CGAL::Frechet_distance_traits_3}
\cgalHasModels{CGAL::Frechet_distance_traits_d}
\cgalHasModelsEnd
*/

class FrechetDistanceTraits {

/// \name Types
/// @{

/*!
Dimension type. Either `CGAL::Dimension_tag`
or `CGAL::Dynamic_dimension_tag`.
*/
typedef unspecified_type Dimension;

/*!
Point type.
*/
typedef unspecified_type Point_d;

/*!
The number type of the %Cartesian coordinates of types `Point_d`
*/
typedef unspecified_type FT;

/*!
A random access iterator type to enumerate the
%Cartesian coordinates of a point.
*/
typedef unspecified_type Cartesian_const_iterator_d;

/*!
Functor with operators to construct iterators on the
first and the past-the-end iterator for the %Cartesian coordinates of a point. This functor must
provide the type `result_type` that must be the same a `Cartesian_const_iterator_d`.
*/
typedef unspecified_type Construct_cartesian_const_iterator_d;

/*!
Functor with operator to construct the bounding box of an object of type `Point_d`.
Its result_type is either `Bbox_2`, `Bbox_3` or `Bbox` depending on `Dimension`.
*/
typedef unspecified_type Construct_bbox_d;

/*!
Functor with operator taking two `Point_d` objects and returning the squared distance between then.
Its result_type is either `FT`.
*/
typedef unspecified_type Compute_squared_distance_d;

/// @}

/// \name Operations
/// @{

/*!
Function used to construct an object of type `Construct_cartesian_const_iterator_d`.
*/
Construct_cartesian_const_iterator_d construct_construct_cartesian_const_iterator_d_object(const Point_d& p) const;

/*!
Function used to construct an object of type `Construct_bbox_d`.
*/
Construct_bbox_d construct_construct_construct_bbox_d_object() const;

/*!
Function used to construct an object of type `Compute_squared_distance_d`.
*/
Compute_squared_distance_d construct_compute_squared_distance_d_object() const;

/// @}
};
