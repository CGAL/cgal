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
public:
/// \name Types
/// @{

/*!
Dimension type, being an instance of `CGAL::Dimension_tag`.
*/
typedef unspecified_type Dimension;

/*!
Point type.
*/
typedef unspecified_type Point_d;

/*!
The number type of the %Cartesian coordinates of types `Point_d`.
For a given `FT n`, `to_interval(n)` must be a valid expression and it must
return an interval containing `n`, represented by a `std::pair<double, double>`.
*/
typedef unspecified_type FT;

/*!
A random access iterator type to enumerate the
%Cartesian coordinates of a point, with `FT` as value type.
*/
typedef unspecified_type Cartesian_const_iterator_d;

/*!
Functor model of `Kernel_d::ConstructCartesianConstIterator_d` to get
iterators over %Cartesian coordinates of a point.
*/
typedef unspecified_type Construct_cartesian_const_iterator_d;

/*!
Functor with operator to construct the bounding box of an object of type `Point_d`,
result type is either `CGAL::Bbox_2`, `CGAL::Bbox_3` or `CGAL::Bbox` depending on `Dimension`.
*/
typedef unspecified_type Construct_bbox_d;

/// @}

/// \name Operations
/// @{

/*!
Function used to construct an object of type `Construct_cartesian_const_iterator_d`.
*/
Construct_cartesian_const_iterator_d construct_cartesian_const_iterator_d_object() const;

/*!
Function used to construct an object of type `Construct_bbox_d`.
*/
Construct_bbox_d construct_bbox_d_object() const;

/// @}
};
