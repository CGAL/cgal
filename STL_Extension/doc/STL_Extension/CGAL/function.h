namespace CGAL {
namespace cpp11 {

/*!
\ingroup PkgStlExtensionRef

An object of the class `cpp11::function` is a general-purpose polymorphic function wrapper.

There is actually no class in namespace `CGAL::cpp11` with this name, but a using declaration which 
imports a class from another namespace. By order of priority: the one in namespace 
`std` is used (provided by C++0x), if not found, then the one in namespace 
`boost` is used.

Please refer to the [C++ manual](http://en.cppreference.com/w/cpp/utility/functional/function) for the documentation of this class.

\cgalHeading{Parameters}

The parameter `R` is the return type. The parameters `Args` are the
parameters of the function.

*/
  
template< typename R, typename... Args >
class function {
}; /* end cpp11::function */

} /* end namespace cpp11 */


} /* end namespace CGAL */
