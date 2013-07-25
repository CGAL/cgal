
namespace CGAL {

/// \addtogroup kernel_dimension

/// @{

/*!
The class `Ambient_dimension` allows to retrieve the dimension of the ambient space of 
a type `T` in a kernel `K`. 

\cgalHeading{Parameters}

The parameter `K` has the default value `Kernel_traits<T>::Kernel>`. 

\cgalHeading{Example}

The following retrieves the dimension of a point type. 

\code
typedef K::Point_2 Point; 
int dimension = Ambient_dimension<Point, K>::value; 
assert(dimension == 2); 
\endcode

\sa `CGAL::Dimension_tag<int dim>` 
\sa `CGAL::Dynamic_dimension_tag` 
\sa `CGAL::Feature_dimension<T, K>` 

*/
template< typename T, typename K = typename Kernel_traits<T>::Kernel >
class Ambient_dimension {
public:

/// \name Constants 
/// @{

/*!
The dimension value as a compile-time 
integral constant. It is implemented as `K::Ambient_dimension<T>::%type::%value`.
It exists only when the dimension is a compile-time constant. 
*/ 
static const int value; 

/// @} 

/// \name Types 
/// @{

/*!
Either `Dimension_tag<dim>` if the dimension is a 
compile-time constant of value `dim`, or `Dynamic_dimension_tag` 
otherwise. It is implemented as `K::Ambient_dimension<T>::%type`.
*/ 
typedef unspecified_type type; 

/// @}

}; /* end Ambient_dimension */

/*!


An object of the class `Dimension_tag` is an empty object which can be used 
for dispatching functions based on the dimension of an object, as provided 
by the `dim` parameter. It is useful in cases where it is not more 
practical to pass the dimension as a template parameter directly. 

\cgalHeading{Example}

The following code declares two functions constructing two points at the origin, 
either in 2D or in 3D. 

\code
Point_2<K> get_origin(Dimension_tag<2>) { return Point_2<K>(ORIGIN); } 
Point_3<K> get_origin(Dimension_tag<3>) { return Point_3<K>(ORIGIN); } 

std::cout << get_origin(Dimension_tag<2>())) << std::endl; 
\endcode

\sa `CGAL::Ambient_dimension<T, K>` 
\sa `CGAL::Feature_dimension<T, K>` 
\sa `CGAL::Dynamic_dimension_tag` 

*/
template< typename int dim >
class Dimension_tag {
public:

/// \name Constants 
/// @{

/*!
The value of the `dim` parameter. 
*/ 
static const int value; 

/// @}

}; /* end Dimension_tag */

/*!


An object of the class `Dynamic_dimension_tag` is an empty object which can be used 
for dispatching functions based on the dimension of an object. 
`Dynamic_dimension_tag` indicates that the dimension is not known at compile-time. 
`Dimension_tag` is the tag class dealing with compile-time dimensions. 

\cgalHeading{Example}

The following code declares two functions constructing two points at the origin, 
either in 2D or in 3D. 

\code
Point_2<K> get_origin(Dimension_tag<2>) { return Point_2<K>(ORIGIN); } 
Point_3<K> get_origin(Dimension_tag<3>) { return Point_3<K>(ORIGIN); } 
Point_d<K> get_origin(Dynamic_dimension_tag) { return Point_d<K>(ORIGIN); } 

std::cout << get_origin(Dynamic_dimension_tag())) << std::endl; 
\endcode

\sa `CGAL::Dimension_tag<int dim>` 
\sa `CGAL::Ambient_dimension<T, K>` 
\sa `CGAL::Feature_dimension<T, K>` 

*/

class Dynamic_dimension_tag {
public:

/// @}

}; /* end Dynamic_dimension_tag */

/*!

The class `Feature_dimension` allows to retrieve the geometric dimension of a type `T` 
in a kernel `K`. 

\cgalHeading{Parameters}

The parameter `K` has the default value `Kernel_traits<T>::Kernel`. 

\cgalHeading{Example}

The following retrieves the dimension of a point type. 

\code
typedef K::Point_2 Point; 
int dimension = Feature_dimension<Point, K>::value; 
assert(dimension == 0); 
\endcode

\sa `CGAL::Dimension_tag<int dim>` 
\sa `CGAL::Dynamic_dimension_tag` 
\sa `CGAL::Ambient_dimension<T, K>` 

*/
template< typename T, typename K = typename Kernel_traits<T>::Kernel >
class Feature_dimension {
public:

/// \name Constants 
/// @{

/*!
The dimension value as a compile-time 
integral constant. It is implemented as `K::Feature_dimension<T>::%type::%value`. 
It exists only when the dimension is a compile-time constant. 
*/ 
static const int value; 

/// @} 

/// \name Types 
/// @{

/*!
Either `Dimension_tag<dim>` if the dimension is a 
compile-time constant of value `dim`, or `Dynamic_dimension_tag` 
otherwise. It is implemented as `K::Feature_dimension<T>::%type`. 
*/ 
typedef unspecified_type type; 

/// @}

}; /* end Feature_dimension */

/// @}

} /* end namespace CGAL */

