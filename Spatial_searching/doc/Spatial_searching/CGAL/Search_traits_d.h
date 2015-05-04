namespace CGAL {

/*!
\ingroup SearchTraitsClasses

The class `Search_traits_d` can be used as a template parameter of the kd tree 
and the search classes.

\tparam Kernel must be a model of the concept `Kernel_d` (for example `Cartesian_d<double>`) 

\tparam Dim must be a `Dimension_tag`
(default value is `Dynamic_dimension_tag`).

\cgalModels `SearchTraits`
\cgalModels `RangeSearchTraits`

\sa `Search_traits_2<Kernel>` 
\sa `Search_traits_3<Kernel>` 
\sa `Search_traits<Point,CartesianConstIterator,ConstructCartesianConstIterator>` 

*/
template< typename Kernel, typename Dim >
class Search_traits_d {
public:

/// \name Types 
/// @{

/*!
Dimension type. Either `Dimension_tag<int dim>`
or `Dynamic_dimension_tag`.
*/
typedef Dim Dimension;

/*!
Number type. 
*/ 
typedef Kernel::FT NT; 

/*!
Point type. 
*/ 
typedef Kernel::Point_d Point_d; 

/*!
Iso box type. 
*/ 
typedef Kernel::Iso_box_d Iso_box_d; 

/*!
Sphere type. 
*/ 
typedef Kernel::Sphere_d Sphere_d; 

/*!
An iterator over the %Cartesian coordinates. 
*/ 
typedef Kernel::Cartesian_const_iterator_d Cartesian_const_iterator; 

/*!
A functor with 
two function operators, which return the begin and past the end iterator for the %Cartesian coordinates. 
The functor for begin has as argument a \link Search_traits_d::Point_d `Point_d`\endlink. The functor for the past the end iterator, 
has as argument a \link Search_traits_d::Point_d `Point_d`\endlink and an `int`. 
*/ 
typedef Kernel::Construct_cartesian_const_iterator_d Construct_cartesian_const_iterator; 

/*!
Functor with operator to construct 
the vertex with lexicographically smallest coordinates of an object of type \link Search_traits_d::Iso_box_d `Iso_box_d`\endlink. 
*/ 
typedef Kernel::Construct_min_vertex_d Construct_min_vertex_d; 

/*!
Functor with operator to construct 
the vertex with lexicographically largest coordinates of an object of type \link Search_traits_d::Iso_box_d `Iso_box_d`\endlink. 
*/ 
typedef Kernel::Construct_max_vertex_d Construct_max_vertex_d; 

/// @}

}; /* end Search_traits_d */
} /* end namespace CGAL */
