
namespace CGAL {
namespace cpp11 {

/*!
\ingroup PkgStlExtension



An object of the class `cpp11::tuple` represents a heterogeneous tuple of elements 
of the types specified in parameters, which are in variadic number. 

There is actually no class in namespace `CGAL::cpp11` with this name, but a using declaration which 
imports a class from another namespace. By order of priority: the one in namespace 
`std` is used (provided by C++0x), if not found, then the one in namespace 
`std::tr1` is used (provided by TR1), and finally, the fallback solution 
is taken from Boost. 



\cgalHeading{Parameters}

The parameters `...` are the value types. 

\cgalHeading{Free functions and helper classes}

Some free functions part of the standard interface of `tuple` are also 
brought in namespace `CGAL::cpp11` with using declarations, these are `make_tuple`, 
`get`, `tie`. Like in C++0x, the `get` function template is 
specialized so that it can take `std::pair` as argument. 
Two standard helper classes are also provided for convenience (`tuple_size` and `tuple_element`). 


*/
template< typename ... >
class tuple {
public:


}; /* end cpp11::tuple */
} /* end namespace cpp11 */
} /* end namespace CGAL */

