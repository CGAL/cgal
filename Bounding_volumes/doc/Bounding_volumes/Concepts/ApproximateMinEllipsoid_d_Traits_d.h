
/*!
\ingroup PkgBoundingVolumesConcepts
\cgalConcept

This concept defines the requirements for traits classes of 
`CGAL::Approximate_min_ellipsoid_d<Traits>`. 

\cgalRefines `DefaultConstructible` 
\cgalRefines `CopyConstructible` 
\cgalRefines `Assignable` 

\cgalHasModel `CGAL::Approximate_min_ellipsoid_d_traits_2<K,ET>` 
\cgalHasModel `CGAL::Approximate_min_ellipsoid_d_traits_3<K,ET>` 
\cgalHasModel `CGAL::Approximate_min_ellipsoid_d_traits_d<K,ET>` 

\sa `CGAL::Min_ellipse_2<Traits>` 

*/

class ApproximateMinEllipsoid_d_Traits_d {
public:

/// \name Types 
/// @{

/*! 
`typedef double FT` 
*/ 
typedef Hidden_type FT; 

/*! 
Some model of concept `RingNumberType` that 
provides exact arithmetic, meaning that 
`CGAL::Number_type_traits<ET>::Has_exact_ring_operations` 
must be `CGAL::Tag_true`. In addition, `ET` must be 
able to exactly represent any finite `double` value. (An 
example for such a type is `CGAL::MP_Float`.). 
The type 
`ET` is to be used by the 
`Approximate_min_ellipsoid_d<Traits>` class for internal, 
exact computations. 
*/ 
typedef Hidden_type ET; 

/*! 
Type of the input points. `Point` must 
provide the default and copy constructor, and must be a model 
of `DefaultConstructible`, `CopyConstructible`, and 
`Assignable`. 
*/ 
typedef Hidden_type Point; 

/*! 
Model for the STL concept 
`RandomAccessIterator` whose value type must be convertible to 
`double`. This type is used to iterate over the %Cartesian 
coordinates of an instance of type `Point`, see 
`cartesian_begin()` below. 
*/ 
typedef Hidden_type Cartesian_const_iterator; 

/// @} 

/// \name Access Functions 
/// @{

/*! 
returns the 
dimension of a point `p`. 
*/ 
int dimension(const Point& p); 

/*! 
returns an input iterator over the Euclidean coordinates 
of the point \f$ p\f$. The range of the iterator must have size 
`dimension(p)`. 
*/ 
Cartesian_const_iterator cartesian_begin(const 
Point& p); 

/// @}

}; /* end ApproximateMinEllipsoid_d_Traits_d */

