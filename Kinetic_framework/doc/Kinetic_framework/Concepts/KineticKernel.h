namespace Kinetic {
/*!
\ingroup PkgKdsFrameworkConcepts
\cgalConcept

The concept `Kinetic::Kernel` acts as the kinetic analog of a \cgal kernel. 
It provides some set of primitives and predicats acting on them. The 
predicates are instances of `Kinetic::CertificateGenerator` and 
can be used to either create `Certificate`s or to evaluate 
instantaneous predicates. 

\cgalHasModel `CGAL::Kinetic::Cartesian<FunctionKernel>`

*/

class Kernel {
public:

/// \name Types 
/// @{

/*!
The type which is used to represent 
coordinates of moving primitives. It is a model of the concept 
`FunctionKernel::Function`. This is the analog of the CGAL 
kernel `RT`. 
*/ 
typedef unspecified_type Motion_function; 

/*!
The type representing the results of 
predicates. See `Kinetic::Certificate`. 
*/ 
typedef unspecified_type Certificate; 

/*!
The type of the function kernel used. 
See `Kinetic::FunctionKernel`. 
*/ 
typedef unspecified_type Function_kernel; 

/// @} 

/// \name Operations 
/// @{

/*!
Gets a copy 
of the function kernel. 
*/ 
Function_kernel function_kernel_object() const; 

/// @}

}; /* end Kinetic::Kernel */

} /* end namespace KineticConcepts */
