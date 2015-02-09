namespace CGAL {

/*!
\ingroup SearchTraitsClasses

The class `Search_traits` can be used as a template parameter of the kd tree 
and the search classes. It is a mere wrapper for the geometric types needed 
by these classes. 
\cgalModels `SearchTraits`

\sa `Search_traits_2<Kernel>` 
\sa `Search_traits_3<Kernel>` 
\sa `Search_traits_d<Kernel>` 

*/
template< typename NT, typename Point, typename CartesianIterator, typename ConstructCartesianIterator, typename ConstructMinVertex, typename ConstructMaxVertex, typename Dim >
class Search_traits {
public:

/// \name Types 
/// @{

/*!
Dimension type. Either `Dimension_tag<int dim>`
or `Dynamic_dimension_tag`.
*/
typedef unspecified_type Dimension;

/*!
The number type of the coordinates. 
*/ 
typedef NT FT; 

/*!
Point type. 
*/ 
typedef Point Point_d; 

/*!
An iterator over the coordinates. 
*/ 
typedef CartesianIterator Cartesian_const_iterator_d; 

/*!
A functor with 
two function operators, which return the begin and past the end iterator for the %Cartesian coordinates. 
The functor for begin has as argument a `Point_d`. The functor for the past the end iterator, 
has as argument a `Point_d` and an `int`. 
*/ 
typedef ConstructCartesianIterator Construct_Cartesian_const_iterator_d; 

/*!
Functor with operator to construct 
the vertex with lexicographically smallest coordinates of an object of type `Iso_box_d`. 
*/ 
typedef ConstructMinVertex Construct_min_vertex_d; 

/*!
Functor with operator to construct 
the vertex with lexicographically largest coordinates of an object of type `Iso_box_d`. 
*/ 
typedef Kernel::ConstructMaxVertex Construct_max_vertex_d; 

/// @}

}; /* end Search_traits */
} /* end namespace CGAL */
