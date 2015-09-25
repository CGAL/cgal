namespace Kinetic {
/*!
\ingroup PkgKdsFrameworkOtherConcepts
\cgalConcept

The concept for priority queues used by the `Simulator`. The concept
basically defines a priority queue which supports deletions and
changes of items in the queue (but not their priorities). Items in the
queue must implement the `Event` concept.

\cgalHasModel `CGAL::Kinetic::Two_list_pointer_event_queue<FunctionKernel>`
\cgalHasModel `CGAL::Kinetic::Heap_pointer_event_queue<FunctionKernel>`

*/

class EventQueue {
public:

/// \name Types 
/// @{

/*!
The type used to access items in the queue in order 
to change or delete them. 
*/ 
typedef unspecified_type Key; 

/*!
The priority type for items in the queue. This 
is typically the same as `Kinetic::Simulator::Time` 
*/ 
typedef unspecified_type Priority; 

/// @} 

/// \name Creation 
/// @{

/*!
Construct a queue which will start at time start and 
run until time end. 
*/ 
EventQueue(Priority start, Priority end, int 
size_hint); 

/// @} 

/// \name Operations 
/// @{

/*!
Insert 
an event into the event queue. A `Key` which can be used to 
manipulated the event is returned. 
*/ 
template <class Event> Key insert(Priority, Event); 

/*!
Erase an event from the queue. 
*/ 
void erase(Key); 

/*!
Change the data in the event referred to by the key. 
*/ 
template <class Event> void set(Key, Event); 

/*!
Access the event referred to by the passed key. 
*/ 
template <class Event> Event& get(Key) const; 

/*!
Return the priority of the event. 
*/ 
Priority priority(Key) const; 

/*!
Return true if the queue is empty. 
*/ 
bool empty(); 

/*!
Return the priority of the next event in the queue. 
*/ 
Priority next_priority() const; 

/*!
Process the next `Event` by calling its process method with its `Priority`. 
*/ 
void process_next(); 

/*!
Set the priority beyond which to ignore events. 
*/ 
void set_end_priority(); 

/*!
Return true if the queue contains 
the event and false if it does not. This is used for auditing 
events and can be slow if needed. 
*/ 
bool contains(Key) const; 

/// @}

}; /* end Kinetic::EventQueue */

} /* end namespace KineticConcepts */
