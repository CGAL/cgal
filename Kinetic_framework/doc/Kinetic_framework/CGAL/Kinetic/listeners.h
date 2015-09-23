
namespace CGAL {
namespace Kinetic {

/*!
\ingroup PkgKdsFrameworkOtherClasses

The class `Kinetic::Active_objects_listener` acts as an intermediate between a moving object 
table and a KDS. It translates the 
`ActiveObjectsTable::Listener::IS_EDITING` notification events into 
appropriate calls to `KDS::insert(Key)`, `KDS::set(Key)`, 
`KDS::erase(Key)`. 

Kinetic data structures can still take advantage of the batch editing 
if they are careful. The methods (such as `KDS::set(Key)` are 
called in lexicographical order in the `Key`s. So, when a KDS is 
preparing to update some certificate involving a recently set object, 
it can first check if the certificate involves another changed object 
which is lexicographically prior. If so, then the certificate has 
already been updated and can be skipped. 

\sa `Kinetic::Active_objects_vector<MovingObject>`
\sa `Kinetic::ActiveObjectsTable`

*/
template< typename ActiveObjectsTable, typename KDS >
class Active_objects_listener {
}; /* end Kinetic::Active_objects_listener */

/*!
\ingroup PkgKdsFrameworkOtherClasses

The class `Kinetic::Simulator_listener` acts as a helper class for kinetic data 
structures which want to respond to 
`Simulator::Listener::HAS_AUDIT_TIME` notifications. When kinetic 
data structures can audit themselves, the `Kinetic::Simulator_listener` calls the 
`audit()` method on the kinetic data structure. 

\sa `Kinetic::Simulator`
\sa `Listener<Interface>`

*/
template< typename Listener, typename KDS >
class Simulator_listener {
public:

/// \name Creation 
/// @{

/*!
default constructor. 
*/ 
Simulator_listener(Simulator::Handle, KDS *kds); 

/// @}

}; /* end Kinetic::Simulator_listener */
} /* end namespace Kinetic */
} /* end namespace CGAL */
