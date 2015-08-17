namespace Kinetic {
/*!
\ingroup PkgKdsFrameworkOtherConcepts
\cgalConcept

The concept `Kinetic::Certificate` represents certificate. Its main purpose is 
to provide a way of creating `Time` objects corresponding to when 
the certificate fails and to cache any useful work done in find the 
`Time` for later. 

\sa `Kinetic::Kernel`

*/

class Certificate {
public:

/// \name Operations 
/// @{

/*!
Returns the next failure time. 
*/ 
Time failure_time(); 

/*!
Returns true if the certificate will ever fail. 
*/ 
bool will_fail(); 

/*!
Advances to the next failure time (the next root of the certificate functions). 
*/ 
void pop_failure_time(); 

/// @}

}; /* end Kinetic::Certificate */

} /* end namespace KineticConcepts */
