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
Dimension type with the ambient dimension `d`.
*/
using Dimension = CGAL::Dimension_tag<d>;

/*!
Point type.
*/
using Point_d = unspecified_type;

/*!
The number type of the %Cartesian coordinates of type `Point_d`.  It must be a model of `FieldNumberType`.
For a given `FT n`, `to_interval(n)` must be a valid expression and it must
return an interval containing `n`, represented by a `std::pair<double, double>`.
*/
using FT = unspecified_type;

/*!
A random access iterator type to enumerate the
%Cartesian coordinates of a point, with `FT` as value type.
*/
using Cartesian_const_iterator_d = unspecified_type;

/*!
Functor with `operator()(const Point_d&)` and `operator()(const Point_d&, int)` for constructing
a begin and past-the-end `Cartesian_const_iterator_d`, respectively.
*/
using Construct_cartesian_const_iterator_d = unspecified_type;

/*!
Functor with operator to construct the bounding box of an object of type `Point_d`,
result type is either `CGAL::Bbox_2`, `CGAL::Bbox_3` or `CGAL::Bbox_d` depending on `Dimension`.
*/
using Construct_bbox_d = unspecified_type;

/*!
Functor with operator to compare the squared distance of two objects of type `Point_d` with a bound of type `FT`,
returning a `CGAL::Comparison_result`.
\note Only needed by class template `CGAL::Frechet_distance_Neighor_search`.
*/
using Compare_squared_distance_d = unspecified_type;
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

/*!
Function used to construct an object of type `Compare_squared_distance_d`.
*/
Compare_squared_distance_d compare_squared_distance_d_object() const;
/// @}
};
