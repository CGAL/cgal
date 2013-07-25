
namespace CGAL {
namespace Kinetic {

/*!
\ingroup PkgKdsFrameworkOtherClasses

The class `Kinetic::Simulator_objects_listener` is a helper for classes which wish to react to 
`Simulator::Listener::DIRECTION_OF_TIME` notifications. The helper 
object translates such notifications `reverse_time` function calls 
on the responder. See `Kinetic::Qt_moving_points_2` for a 
simple example of using this helper function. 

\sa `Kinetic::Listener`

*/
template< typename Simulator_listener, typename KDS >
class Simulator_objects_listener {
public:

/// \name Creation 
/// @{

/*!
default constructor. 
*/ 
Simulator_objects_listener(Simulator::Handle, KDS*); 

/// @}

}; /* end Kinetic::Simulator_objects_listener */
} /* end namespace Kinetic */
} /* end namespace CGAL */
