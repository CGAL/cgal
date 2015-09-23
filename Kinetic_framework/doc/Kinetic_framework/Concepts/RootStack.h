namespace Kinetic {
/*!
\ingroup PkgKdsFrameworkOtherConcepts
\cgalConcept

The concept `Kinetic::RootStack` enumerates through roots of a function 
contained in a half open interval [lb \f$ \dots\f$ ub). 

\sa `CGAL::Kinetic::FunctionKernel`
\sa `CGAL::Kinetic::Certificate`

*/

class RootStack {
public:

/// \name Types 
/// @{

/*!
The root of a function. 
*/ 
typedef unspecified_type Root; 

/*!
The traits class for this concept. 
*/ 
typedef unspecified_type Traits; 

/// @} 

/// \name Creation 
/// @{

/*!
default constructor. 
*/ 
RootStack(); 

/*!
Construct a `Kinetic::RootStack` over the roots of `f` in the half open interval [`lb` to `ub`). 
*/ 
RootStack(Function f, Root lb, Root ub, Traits tr); 

/// @} 

/// \name Operations 
/// @{

/*!
Advance to the next root. As a precondition, empty() must be false. 
*/ 
void pop(); 

/*!
Return the current root. As a precondition, empty() must be false. Note that the `Root` returned might not actually be in the interval (since the solver has not yet proved that there are no more roots). 
*/ 
Root top(); 

/*!
Return true if there are known to be no more roots left. There might not actually be any roots of the polynomial left in the interval, but the work necessary to prove this has been delayed. 
*/ 
bool empty(); 

/// @}

}; /* end Kinetic::RootStack */

} /* end namespace KineticConcepts */
