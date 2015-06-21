namespace CGAL {

/*!
\ingroup SearchTraitsClasses

The class `Search_traits_2` can be used as a template parameter of the kd tree 
and the search classes. 

\tparam Kernel must be a model of the concept `Kernel`, 
for example `Simple_cartesian<double>` or `Simple_cartesian<Gmpq>`. 

\cgalModels `SearchTraits`
\cgalModels `RangeSearchTraits`

\sa `Search_traits_3<Kernel>` 
\sa `Search_traits<NT,Point,CartesianConstIterator,ConstructCartesianConstIterator,Dim>` 

*/
template< typename Kernel >
class Search_traits_2 {
public:

/// \name Types 
/// @{

/*!
Dimension type.
*/
typedef Dimension_tag<2> Dimension;

/*!
Number type. 
*/ 
typedef Kernel::FT FT; 

/*!
Point type. 
*/ 
typedef Kernel::Point_2 Point_d; 

/*!
Iso box type. 
*/ 
typedef Kernel::Iso_rectangle_2 Iso_box_d; 

/*!
Sphere type. 
*/ 
typedef Kernel::Circle_2 Sphere_d; 

/*!
An iterator over the %Cartesian coordinates. 
*/ 
typedef Kernel::Cartesian_const_iterator_2 Cartesian_const_iterator_d; 

/*!
A functor with 
two function operators, which return the begin and past the end iterator for the %Cartesian coordinates. 
The functor for begin has as argument a \link Search_traits_2::Point_d `Point_d`\endlink. The functor for the past the end iterator, 
has as argument a \link Search_traits_2::Point_d `Point_d`\endlink and an `int`. 
*/ 
typedef Kernel::Construct_cartesian_const_iterator_2 Construct_cartesian_const_iterator_d; 

/*!
Functor with operator to construct 
the iso box from two points. 
*/ 
typedef Kernel::Construct_iso_rectangle_2 Construct_iso_box_d; 

/*!
Functor with operator to construct 
the center of an object of type \link Search_traits_2::Sphere_d `Sphere_d`\endlink. 
*/ 
typedef Kernel::Construct_center_2 Construct_center_d; 

/*!
Functor with operator to compute 
the squared radius of a an object of type \link Search_traits_2::Sphere_d `Sphere_d`\endlink. 
*/ 
typedef Kernel::Compute_squared_radius_2 Construct_squared_radius_d; 

/*!
Functor with operator to construct 
the vertex with lexicographically smallest coordinates of an object of type \link Search_traits_2::Iso_box_d `Iso_box_d`\endlink. 
*/ 
typedef Kernel::Construct_min_vertex_2 Construct_min_vertex_d; 

/*!
Functor with operator to construct 
the vertex with lexicographically largest coordinates of an object of type \link Search_traits_2::Iso_box_d `Iso_box_d`\endlink. 
*/ 
typedef Kernel::Construct_max_vertex_2 Construct_max_vertex_d; 

/// @}

}; /* end Search_traits_2 */
} /* end namespace CGAL */
