namespace CGAL {
namespace cpp11 {

/*!
\ingroup PkgStlExtension

An object of the class `cpp11::array` represents an array of elements 
of type `T`, the number of which is specified by the second template argument. 

There is actually no class in namespace `CGAL::cpp11` with this name, but a using declaration which 
imports a class from another namespace. By order of priority: the one in namespace 
`std` is used (provided by C++0x), if not found, then the one in namespace 
`std::tr1` is used (provided by TR1), and finally, the fallback solution 
is taken from Boost. 

\cgalHeading{Parameters}

The parameter `T` is the value type. The second parameter is the
dimension of the array. 

\cgalHeading{Extensions}

\cgal provides a `make_array` function for this purpose, up to a
certain number of arguments.
*/
template< typename T, int >
class array {
}; /* end cpp11::array */

/*!
\relates cpp11::array 

\returns `array<T, N>` where `N` is the number of arguments given to
the function. The position of each argument in the array is the same
as its position in the argument list.

The maximal number of arguments is `6`.
*/ 
template <class T> array<T, N> make_array(const T&...); 

} /* end namespace cpp11 */
} /* end namespace CGAL */
