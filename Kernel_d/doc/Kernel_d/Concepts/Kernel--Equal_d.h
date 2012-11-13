
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

A model for this must provide: 

*/

class Kernel_d::Equal_d {
public:

/// \name See Also 
/// @{

/*! 
returns true iff \f$ p\f$ and \f$ q\f$ are equal (as 
\f$ d\f$-dimensional points).

\pre `p` and `q` have the same dimension. 
*/ 
bool operator()(const Kernel_d::Point_d&p, const 
Kernel_d::Point_d&q); 

/// @}

}; /* end Kernel_d::Equal_d */

