
/*!
\ingroup PkgLinearCellComplexConcepts
\cgalconcept

Required types and functors for the `LinearCellComplexTraits` concept. This 
geometric traits concept is used in the `Linear_cell_complex` 
class. 

\hasModel `CGAL::Linear_cell_complex_traits<d,K>`

\sa `CGAL::Linear_cell_complex<d,d2,LCCTraits,Items,Alloc>`

*/

class LinearCellComplexTraits {
public:

/// \name Constants 
/// @{

/*! 
The ambient dimension, must be \f$ >\f$1. 
*/ 
static unsigned int ambient_dimension; 

/// @} 

/// \name Types 
/// @{

/*! 
a number type that is a model of FieldNumberType. 
*/ 
typedef Hidden_type FT; 

/*! 
point type. 
*/ 
typedef Hidden_type Point; 

/*! 
vector type. 
*/ 
typedef Hidden_type Vector; 

/// @} 

/// \name Constructions 
/// @{

/*! 
Functor that provides `Point operator() (const Point& p, const Vector& v)`, 
which constructs the translation of point `p` by vector `v`, and 
`Point operator() (const CGAL::Origin&, const Vector& v)`, 
which constructs the translation of a point at the origin by vector `v` 
(used in `Linear_cell_complex::barycenter`). 
*/ 
typedef Hidden_type Construct_translated_point; 

/*! 
Functor that provides `Vector operator() (const Point& p1, const Point& p2)` 
which constructs a vector as the difference of points `p2-p1`, and 
`Vector operator() (const CGAL::Origin&, const Point& p)` 
which constructs a vector as the difference of point `p` and a point at the origin 
(used in `Linear_cell_complex::barycenter` and `CGAL::import_from_plane_graph`). 
*/ 
typedef Hidden_type Construct_vector; 

/*! 
Functor that provides `Vector operator() (const Vector& v1, const Vector& v2)` 
which constructs a vector as the sum of vectors `v1+v2` 
(used in `Linear_cell_complex::barycenter`, `CGAL::compute_normal_of_cell_0` 
and `CGAL::compute_normal_of_cell_2`). 
*/ 
typedef Hidden_type Construct_sum_of_vectors; 

/*! 
Functor that provides `Vector operator() (const Vector& v, FT scale)` 
which constructs a vector equal to vector `v` scaled by `scale` factor 
(used in `Linear_cell_complex::barycenter` , `CGAL::compute_normal_of_cell_0` 
and `CGAL::compute_normal_of_cell_2`). 
*/ 
typedef Hidden_type Construct_scaled_vector; 

/*! 
Functor that provides `Point operator() (const Point& p1, const Point& p2)` 
which constructs the midpoint of points `p1` and `p2` 
(used in `Linear_cell_complex::barycenter`). 
*/ 
typedef Hidden_type Construct_midpoint; 

/// @}

/// \name
/// If `ambient_dimension==2`
/// @{

/*! 
a model of `Direction_2`. 
*/ 
typedef Hidden_type Direction_2; 

/*! 
a model of `ConstructDirection_2` (used in `CGAL::import_from_plane_graph`). 
*/ 
typedef Hidden_type Construct_direction_2; 

/// @}

/// \name
/// If `ambient_dimension==3`
/// @{

/*! 
a model of `ConstructNormal_3` (used in `CGAL::compute_normal_of_cell_2`). 
*/ 
typedef Hidden_type Construct_normal_3; 

/*! 
a model of `Collinear_3` (used in `CGAL::compute_normal_of_cell_2`). 
*/ 
typedef Hidden_type Collinear_3; 

/// @}

}; /* end LinearCellComplexTraits */

