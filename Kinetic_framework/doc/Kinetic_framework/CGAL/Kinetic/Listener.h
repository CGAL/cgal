
namespace CGAL {

/*!
\ingroup PkgKdsFrameworkOtherClasses

The `Listener` class provides the core of the run time notification 
system used by the kinetic data structures package. In short, 
notifications are handled through proxy objects called listeners. In 
order to listen for notifications from an object, called the notifier, 
you make define a small class called a listener proxy, which inherits 
from the Listener interface defined by the notifier. When constructing 
your listner poxy, you pass a reference counted pointer to the 
notifier, which is used to register the proxy for notifications. When 
a notification occurs, the notifier calls the `new_notification` 
method on the proxy, passing the type of the notification. The proxy 
stores a reference counted pointer to the notifier, ensuring that 
there are never any dangling pointers in the system. 

The class `Listener` provides base class for listener proxy objects. 
A notifier should provide a class which inherits from this base. To 
use this base class, implement a class, here called `Interface`, 
which defines a type `Interface::Notification_type` and a type 
`Interface::Notifier_handle`. 

The `Notification_type` is generally an enum with one value for 
each type of notification which can be used. 

The `Notifier_handle` is the type of a (ref counted) pointer to 
the object providing the notifications. The ref counter pointer must 
provide a nested type `Pointer` which is the type of a raw 
pointer. 

The `Listener` maintains a ref counted pointer to the object 
performing notifications. It is registered for notifications on 
construction and unregistered on destruction using the function 
`set_listener(Listener<Interface>*)` on the object providing the 
notifications. The use of ref counted pointers means that as long as 
the notification object exists, the object providing the notifications 
must exist, ensuring that the object providing the notifications is 
not prematurely destroyed. 

These objects cannot be copied since the notifier only support one 
listener. If copying and more than one listener are desired, the 
`Multi_listener<Interface>` base class should be used instead. 

As a side note, Boost provides a similar functionality in the Boost.Signal 
package. However, it is quite a bit more complex (and flexible). This 
complexity add significantly to compile time and (although I did not 
test this directly), I suspect it is much slower at runtime due to the 
overhead of worrying about signal orders and not supporting single 
signals. In addition, it does not get on well with Qt due to 
collisions with the Qt moc keywords. 

There is also the TinyTL library which implements signals. As of 
writing it did not have any easy support for making sure all pointers 
are valid, so it did not seem to offer significant code saving over 
writing my own. 

\sa `Multi_listener<Interface>`

\cgalHeading{Example}

Here is a simplier class that provides notifications: 

\code{.cpp} 

struct Notifier: public CGAL::Kinetic::Ref_counted<Notifier> 
{ 
public: 
Notifier(): data_(0), listener_(NULL){} 

struct Listener_interface 
{ 
public: 
typedef enum Notification_type {DATA_CHANGED} Notification_type; 
typedef Notifier::Handle Notifier_handle; 
}; 

typedef CGAL::Kinetic::Listener<Listener_interface> Listener; 
friend class CGAL::Kinetic::Listener<Listener_interface>; 

void set_data(int d) { 
data_=d; 
if (listener_ != NULL) listener_->new_notification(Listener_interface::DATA_CHANGED); 
} 

protected: 
void set_listener(Listener *l) { 
listener_= l; 
} 
Listener* listener() const {return listener_;} 

int data_; 
Listener *listener_; 
}; 

\endcode 

Now the listener: 

\code{.cpp} 

struct My_listener: public Notifier::Listener, public CGAL::Kinetic::Ref_counted<My_listener> 
{ 
typedef Notifier::Listener::Notifier_handle PP; 
My_listener(PP p): P(p){} 

void new_notification(P::Notification_type nt) { 
... 
} 
}; 

\endcode 

*/
template< typename Interface >
class Listener {
public:

/// \name Types 
/// @{

/*!
This type is inherited from the `Interface` template argument. It is a reference counted pointer type for the object providing notifications. 
*/ 
typedef unspecified_type Notifier_handle; 

/*!
The type (usually an enum) used to distinguish different types of notifications. This is inherited from the `Interface` template argument. 
*/ 
typedef unspecified_type Notification_type; 

/// @} 

/// \name Creation 
/// @{

/*!
The `Listener` subscribes to events coming from the notifier and stores a pointer to the notifier. 
*/ 
Listener(Notifier_handle np); 

/// @} 

/// \name Operations 
/// @{

/*!
Return a pointer to the notifier. 
*/ 
Notifier_handle notifier(); 

/*!
This method is pure virtual. A class which wishes to receive events must inherit from this class and implement this method. The method will then be called whenever there is a notification. 
*/ 
virtual void new_notification(Notification_type); 

/// @}

}; /* end Listener */
} /* end namespace CGAL */
