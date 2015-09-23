namespace Kinetic {
/*!
\ingroup PkgKdsFrameworkConcepts
\cgalConcept

This concept ties together the parts needed in order to run a kinetic 
data structure. All support trajectories defined by polynomial coordinates. The 
`Exact` vs `Inexect` picks whether the roots of the 
certificate functions are compared exactly or approximated 
numerically. The regular triangulation models have weighted points of 
the appropriate dimension as the primitive used in the 
`Kinetic::InstantaneousKernel` and the 
`Kinetic::ActiveObjectsTable`. 

\cgalHasModel `CGAL::Kinetic::Exact_simulation_traits`
\cgalHasModel `CGAL::Kinetic::Inexact_simulation_traits`
\cgalHasModel `CGAL::Kinetic::Regular_triangulation_exact_simulation_traits`
\cgalHasModel `CGAL::Kinetic::Regular_triangulation_inexact_simulation_traits`

*/

class SimulationTraits {
public:

/// \name Types 
/// @{

/*!
The number type used for representation. 
*/ 
typedef unspecified_type NT; 

/*!
A model of 
`Kinetic::InstantaneousKernel` which can be used to apply static CGAL 
data structures to snapshots of moving data. 
*/ 
typedef unspecified_type Instantaneous_kernel; 

/*!
A model of `Kinetic::Kernel`. 
*/ 
typedef unspecified_type Kinetic_kernel; 

/*!
A model of `Kinetic::FunctionKernel`. 
*/ 
typedef unspecified_type Function_kernel; 

/*!
A model of 
`Kinetic::ActiveObjectsTable` which holds the relevant kinetic 
primitives. 
*/ 
typedef unspecified_type Active_points_@omtable; 

/*!
A model of `Kinetic::Simulator` which will be 
used by all the kinetic data structures. 
*/ 
typedef unspecified_type Simulator; 

/// @} 

/// \name Operations 
/// @{

/*!
Get a new instantaneous kernel. 
*/ 
Instantaneous_kernel instantaneous_kernel_object(); 

/*!
Get a new kinetic kernel. 
*/ 
Kinetic_kernel kinetic_kernel_object(); 

/*!
Get a new function kernel. 
*/ 
Function_kernel function_kernel_object(); 

/*!
Return a pointer to the `Kinetic::Simulator` which is to be used in the simulation. 
*/ 
Simulator::Handle simulator_handle(); 

/*!
Return a pointer to the 
`Kinetic::ActiveObjectsTable` which is to be used in the 
simulation. 
*/ 
Active_points_[123]_table::Handle 
active_points_[123]_table_handle(); 

/// @}

}; /* end Kinetic::SimulationTraits */

} /* end namespace KineticConcepts */
