
namespace CGAL {

/*!
\ingroup PkgKdsFrameworkOtherClasses

The class `Ref_counted` implements a base class for objects which are 
reference counted. To use it simply inherit from `Ref_counted` (passing 
the type to be reference counted as the template argument) and then 
access the object through `Handle` objects rather than bare C++ 
pointers. 

\cgalHeading{Operations}

There are no methods which should be called by users of this class. 

*/
template< typename T >
class Ref_counted {
public:

/// \name Types 
/// @{

/*!
A reference counted pointer to an Object. 
*/ 
typedef unspecified_type Handle; 

/*!
A const reference counted pointer to an Object. 
*/ 
typedef unspecified_type Const_handle; 

/// @} 

/// \name Creation 
/// @{

/*!
default constructor. 
*/ 
Ref_counted(); 

/// @}

}; /* end Ref_counted */
} /* end namespace CGAL */
