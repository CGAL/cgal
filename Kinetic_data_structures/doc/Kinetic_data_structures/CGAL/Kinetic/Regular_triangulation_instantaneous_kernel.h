
namespace CGAL {
namespace Kinetic {

/*!
\ingroup PkgKdsTri

The class `Kinetic::Regular_triangulation_instantaneous_kernel` is an instantaneous kernel for use with a 
regular triangulation data structure. There is not currently a reason 
for the user to call this directly unless the user wants created their 
own simulation traits as it is included as part of the 
`Kinetic::Regular_triangulation_exact_simulation_traits`. 

\cgalModels `Kinetic::RegularTriangulationTraits_3`
\cgalModels `Kinetic::InstantaneousKernel`

\sa `Kinetic::Regular_triangulation_3<Traits, Visitor, Triangulation>`

*/
template< typename ActiveObjectsTable, typename StaticKernel >
class Regular_triangulation_instantaneous_kernel {
public:

/// @}

}; /* end Kinetic::Regular_triangulation_instantaneous_kernel */
} /* end namespace Kinetic */
} /* end namespace CGAL */
