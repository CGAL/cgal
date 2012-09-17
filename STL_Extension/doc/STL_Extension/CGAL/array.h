
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



### Parameters ###

The parameter `T` is the value type. The second parameter is the 
dimension of the array. 

### Construction Functions ###

The array class does not provide a constructor which can be used to initialize 
data members. \cgal therefore provides a `make_array` function for 
this purpose, up to a certain number of arguments. 

*/
template< typename T, int >
class array {
public:


}; /* end cpp11::array */

/*! 
\relates cpp11::array 
returns an array of dimension 1 whose first element is `a`. 
*/ 
template <class T> array<T, 1> make_array(const T& a); 

/*! 
\relates cpp11::array 
returns an array of dimension 2 whose first element is `a1` 
and second element is `a2`. 
*/ 
template <class T> array<T, 2> make_array(const T& a1, const T& a2); 

} /* end namespace cpp11 */
} /* end namespace CGAL */
