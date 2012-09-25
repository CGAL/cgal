
/*!
\ingroup PkgKernelDKernelConcept
\cgalconcept

A model for this must provide: 

*/

class Kernel_d::Component_accessor_d {
public:

/// \name See Also 
/// @{

/*! 
returns 
the dimension of \f$ p\f$. 
*/ 
int dimension(const Kernel_d::Point_d& p); 

/*! 
returns the ith homogeneous coordinate of \f$ p\f$.

\pre `0 <= i <= dimension(p)`. 
*/ 
Kernel_d::RT homogeneous(const Kernel_d::Point_d& p, 
int i); 

/*! 
returns the ith %Cartesian coordinate of \f$ p\f$.

\pre `0 <= i < dimension(p)`. 
*/ 
Kernel_d::FT cartesian(const Kernel_d::Point_d& p, int 
i); 

/// @}

}; /* end Kernel_d::Component_accessor_d */

