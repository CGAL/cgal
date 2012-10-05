/*!
\ingroup PkgSpatialSearchingDConcepts
\cgalconcept

The concept `SearchTraits` defines the requirements for the template 
parameter of the search classes. 

\hasModel `CGAL::Cartesian_d<FT>` 
\hasModel `CGAL::Homogeneous_d<RT>` 
\hasModel `CGAL::Search_traits_2<Kernel>` 
\hasModel `CGAL::Search_traits_3<Kernel>`
\hasModel `CGAL::Search_traits<NT,Point,CartesianCoordinateIterator,ConstructCartesianCoordinateIterator,ConstructMinVertex,ConstructMaxVertex>`

\sa `RangeSearchTraits` 
\sa `CGAL::Search_traits_adapter<Key,PointPropertyMap,BaseTraits>` 

*/

class SearchTraits {
public:

/// \name Types 
/// @{

/*! 
Point type. `CGAL::Kernel_traits` has to be 
specialized for this type. 
*/ 
typedef Hidden_type Point_d; 

/*! 
The number type of the %Cartesian coordinates of types `Point_d` 
*/ 
typedef Hidden_type FT; 

/*! 
A random access iterator type to enumerate the 
%Cartesian coordinates of a point. 
*/ 
typedef Hidden_type Cartesian_const_iterator_d; 

/*! 
Functor with operators to construct iterators on the 
first and the past-the-end iterator for the %Cartesian coordinates of a point. This functor must 
provides the type `result_type` that must be the same a `Cartesian_const_iterator_d`. 
*/ 
typedef Hidden_type Construct_cartesian_const_iterator_d; 

/// @} 

/// \name Operations 
/// @{

/*! 
Function used to construct an object of type `Construct_cartesian_const_iterator_d`. 
*/ 
Construct_cartesian_const_iterator_d construct_construct_cartesian_const_iterator_d_object(const Point_d& p) const; 

/// @}

}; /* end SearchTraits */
