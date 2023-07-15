/*!
\ingroup PkgSpatialSearchingDConcepts
\cgalConcept

The concept `SearchGeomTraits_3` defines the requirements for the template
parameter of the search traits classes.

\cgalHasModelsBegin
\cgalHasModelsBare{All models of the \cgal concept `Kernel`}
\cgalHasModelsEnd

*/

class SearchGeomTraits_3 {
public:

/// \name Types
/// @{

/*!
Point type.
*/
typedef unspecified_type Point_3;

/*!
The number type of the %Cartesian coordinates of types `Point_3`
*/
typedef unspecified_type FT;
/*!
Iso box type.
*/
typedef unspecified_type Iso_cuboid_3;
/*!
Sphere type.
*/
typedef unspecified_type Sphere_3;

/*!
Functor model of `Kernel::ConstructIsoCuboid_3`
*/
typedef unspecified_type Construct_iso_cuboid_3;

/*!
Functor model of `Kernel::ConstructMinVertex_3`
*/
typedef unspecified_type Construct_min_vertex_3;

/*!
Functor model of `Kernel::ConstructMaxVertex_3`
*/
typedef unspecified_type Construct_max_vertex_3;

/*!
Functor model of `Kernel::ConstructCenter_3`
*/
typedef unspecified_type Construct_center_3;

/*!
Functor model of `Kernel::ComputeSquaredRadius_3 `
*/
typedef unspecified_type Compute_squared_radius_3;

/*!
A random access iterator type to enumerate the
%Cartesian coordinates of a point.
*/
typedef unspecified_type Cartesian_const_iterator_3;

/*!
Functor with operators to construct iterators on the
first and the past-the-end iterator for the %Cartesian coordinates of a point. This functor must
provide the type `result_type` that must be the same a `Cartesian_const_iterator_3`.
*/
typedef unspecified_type Construct_cartesian_const_iterator_3;

/// @}



}; /* end SearchTraits */
