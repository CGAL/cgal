namespace CGAL {

/*!
\ingroup SearchTraitsClasses

The class `Search_traits_d` can be used as a template parameter of the kd tree 
and the search classes. `Kernel` must be a \cgal kernel. 

`Kernel` must be a d-dimensional \cgal kernel. 

### Parameters ###

Expects for the template argument a model of the concept `Kernel_d`, 
for example `CGAL::Cartesian_d<double>` or `CGAL::Homogeneous_d<CGAL::Gmpz>`. 

\cgalModels `SearchTraits`
\cgalModels `RangeSearchTraits`

\sa `Search_traits_2<Kernel>` 
\sa `Search_traits_3<Kernel>` 
\sa `Search_traits<Point,CartesianConstIterator,ConstructCartesianConstIterator>` 

*/
template< typename Kernel >
class Search_traits_d {
public:

/// \name Types 
/// @{

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
An iterator over the Cartesian 
coordinates. 
*/ 
typedef Kernel::Cartesian_const_iterator_d Cartesian_const_iterator; 

/*! 
A functor with 
two function operators, which return the begin and past the end iterator for the Cartesian coordinates. 
The functor for begin has as argument a `Point_d`. The functor for the past the end iterator, 
has as argument a `Point_d` and an `int`. 
*/ 
typedef Kernel::Construct_cartesian_const_iterator_d Construct_cartesian_const_iterator; 

/*! 
Functor with operator to construct 
the vertex with lexicographically smallest coordinates of an object of type `Iso_box_d`. 
*/ 
typedef Kernel::Construct_min_vertex_d Construct_min_vertex_d; 

/*! 
Functor with operator to construct 
the vertex with lexicographically largest coordinates of an object of type `Iso_box_d`. 
*/ 
typedef Kernel::Construct_max_vertex_d Construct_max_vertex_d; 

/// @}

}; /* end Search_traits_d */
} /* end namespace CGAL */
