namespace KineticConcepts {
/*!
\ingroup PkgKdsFrameworkConcepts
\cgalconcept

This concept ties together the parts needed in order to run a kinetic 
data structure. We provide several models of this concept: 
<UL> 
<LI>`Kinetic::Exact_simulation_traits` 
<LI>`Kinetic::Inexact_simulation_traits` 
<LI>`Kinetic::Regular_triangulation_exact_simulation_traits` 
<LI>`Kinetic::Regular_triangulation_inexact_simulation_traits` 
</UL> 

All support trajectories defined by polynomial coordinates. The 
`Exact` vs `Inexect` picks whether the roots of the 
certificate functions are compared exactly or approximated 
numerically. The regular triangulation models have weighted points of 
the appropriate dimension as the primitive used in the 
`Kinetic::InstantaneousKernel` and the 
`Kinetic::ActiveObjectsTable`. 

\hasModel `Kinetic::Exact_simulation_traits`
\hasModel `Kinetic::Inexact_simulation_traits`
\hasModel `Kinetic::Regular_triangulation_exact_simulation_traits`
\hasModel `Kinetic::Regular_triangulation_inexact_simulation_traits`

*/

class SimulationTraits {
public:

/// \name Types 
/// @{

/*! 
The number type used for representation. 
*/ 
typedef Hidden_type NT; 

/*! 
A model of 
`Kinetic::InstantaneousKernel` which can be used to apply static CGAL 
data structures to snapshots of moving data. 
*/ 
typedef Hidden_type Instantaneous_kernel; 

/*! 
A model of `Kinetic::Kernel`. 
*/ 
typedef Hidden_type Kinetic_kernel; 

/*! 
A model of `Kinetic::FunctionKernel`. 
*/ 
typedef Hidden_type Function_kernel; 

/*! 
A model of 
`Kinetic::ActiveObjectsTable` which holds the relevant kinetic 
primitives. 
*/ 
typedef Hidden_type Active_points_@omtable; 

/*! 
A model of `Kinetic::Simulator` which will be 
used by all the kinetic data structures. 
*/ 
typedef Hidden_type Simulator; 

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
