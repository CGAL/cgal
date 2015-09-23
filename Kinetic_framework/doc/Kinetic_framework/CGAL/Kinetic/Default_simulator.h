
namespace CGAL {
namespace Kinetic {

/*!
\ingroup PkgKdsFrameworkClasses

The class `Kinetic::Default_simulator` controls kinetic data structures by maintaining 
a concept of time and ensuring that events are processed when 
necessary. 

\cgalModels `Kinetic::Simulator`

*/
template< typename FunctionKernel, typename EventQueue >
class Default_simulator {
public:

/// \name Creation 
/// @{

/*!
Construct a `Kinetic::Default_simulator` which will process events between times start and end (events outside this window will be discarded). 
*/ 
Default_simulator(const Time start=Time(0), const Time end= Time::infinity()); 

/// @}

}; /* end Kinetic::Default_simulator */
} /* end namespace Kinetic */
} /* end namespace CGAL */
