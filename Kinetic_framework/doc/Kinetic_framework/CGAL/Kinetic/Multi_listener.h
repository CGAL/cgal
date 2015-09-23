
namespace CGAL {

/*!
\ingroup PkgKdsFrameworkOtherClasses

The class `Multi_listener` implements a base class for listeners where more 
than one listener is allowed to subscribe to a notifier. See 
`Listener` for full documentation. This uses the function calls 
`new_listener()` and `delete_listener()` to register and 
unrester the listener (instead of `set_listener()`). 

\sa `Listener<Interface>`

*/
template< typename Interface >
class Multi_listener {
public:

/// @}

}; /* end Multi_listener */
} /* end namespace CGAL */
