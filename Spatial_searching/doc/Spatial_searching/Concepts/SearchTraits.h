/*!
\ingroup PkgSpatialSearchingDConcepts
\cgalConcept

The concept `SearchTraits` defines the requirements for the template 
parameter of the search classes. 

\cgalHasModel `CGAL::Cartesian_d<FT>` 
\cgalHasModel `CGAL::Homogeneous_d<RT>` 
\cgalHasModel `CGAL::Search_traits_2<Kernel>` 
\cgalHasModel `CGAL::Search_traits_3<Kernel>`
\cgalHasModel `CGAL::Search_traits_d<Kernel,Dim>`
\cgalHasModel `CGAL::Search_traits<NT,Point,CartesianCoordinateIterator,ConstructCartesianCoordinateIterator,ConstructMinVertex,ConstructMaxVertex>`

\sa `RangeSearchTraits` 
\sa `CGAL::Search_traits_adapter<Key,PointPropertyMap,BaseTraits>` 

*/

class SearchTraits {
public:

/// \name Types 
/// @{

/*!
Dimension type. Either `CGAL::Dimension_tag`
or `CGAL::Dynamic_dimension_tag`.
*/
typedef unspecified_type Dimension;

/*!
Point type. `CGAL::Kernel_traits` has to be 
specialized for this type. 
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

/// @} 

/// \name Operations 
/// @{

/*!
Function used to construct an object of type `Construct_cartesian_const_iterator_d`. 
*/ 
Construct_cartesian_const_iterator_d construct_construct_cartesian_const_iterator_d_object(const Point_d& p) const; 

/// @}

}; /* end SearchTraits */
