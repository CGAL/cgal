namespace Kinetic {
/*!
\ingroup PkgKdsFrameworkConcepts
\cgalConcept

This container holds a set of objects of a particular type. It creates 
notifications using the standard `Multi_listener<Interface>` 
interface when a primitive changes or is added or deleted. Objects 
which are listening for events can then ask which primitives changed. 

For speed, modifications to the `Kinetic::ActiveObjectsTable` can be grouped into 
editing sessions. A session is begun by calling 
`set_is_editing(true)` and ended by calling `set_is_editing(false)`. 
There is one type of notification, namely, `Listener::IS_EDITING` 
which occurs when the editing mode is set to false, signaling that a 
batch of changes is completed. 

As an convenience, the change methods can be called without setting 
the editing state to true, this acts as if it were set to true for 
that one function call. 

\cgalHasModel `CGAL::Kinetic::Active_objects_vector<MovingObject>`

\sa `Multi_listener<Interface>`
\sa `CGAL::Kinetic::Active_objects_listener_helper<ActiveObjectsTable, KDS>` 

*/

class ActiveObjectsTable {
public:

/// \name Types 
/// @{

/*!
A key identifying an object in the table. 
*/ 
typedef unspecified_type Key; 

/*!
The type being stored in the table. 
*/ 
typedef unspecified_type Data; 

/*!
The base class to derive from for listening for runtime events. 
*/ 
typedef unspecified_type Listener; 

/// @}

/// \name Iterators
/// The following types are iterators. Each type, `Foo_iterator` has
/// two corresponding methods `foo_begin` and `foo_end` which allow
/// you to iterate through the objects in the set `Foo`.
/// @{

/*!
An iterator through all the valid keys in the table. 
*/ 
typedef unspecified_type Key_iterator; 

/*!
An iterator through all the objects which have been changed in the last editing session. The iterator iterates through the objects in lexicographical order. 
*/ 
typedef unspecified_type Changed_iterator; 

/*!
An iterator through all the objects which were added in the last editing session. 
*/ 
typedef unspecified_type Inserted_iterator; 

/*!
An iterator through all the objects which were deleted in the last editing session. 
*/ 
typedef unspecified_type Erased_iterator; 

/// @} 

/// \name Operations 
/// @{

/*!
Access the object referenced by the key. 
*/ 
Data operator[](Key key) const; 

/*!
Access the object referenced by the key. 
*/ 
Data at(Key key) const; 

/*!
Set the editing state of 
the object. A notification is sent when the editing state is set to 
false after it has been true, i.e.\ the editing session is finished. 
This allows changes to be batched together. 
*/ 
void set_is_editing(bool is_editing); 

/*!
Access the editing state. 
*/ 
bool is_editing() const; 

/*!
This method changes 
the motion of one moving object. The position at the current time 
should not be different from the previous current position. However, 
at the moment I do not check this as there is no reference to time 
in the `Kinetic::ActiveObjectsTable`. If `is_editing()` is not true, then it is as 
if the calls `set_is_editing(true)`, `set(key, value)` and 
finally `set_is_editing(false)` were made. If it is true, then 
no notifications are created. 
*/ 
void set(Key key, Data object); 

/*!
Insert a new object into the 
table and return a `Key` which can be used to refer to it. See 
`set(Key, Data)` for a description of editing modes. 
*/ 
Key insert_object(Data ob); 

/*!
Delete an object from the table. The 
object with Key key must already be in the table. This does not 
necessarily decrease the amount of storage used at all. In fact, it 
is unlikely to do so. See `set(Key,Data)` for an explainating 
of how the editing modes are used. 
*/ 
void erase(Key key); 

/*!
Remove all objects from the table and free all storage. 
*/ 
void clear(); 

/*!
Returns true if the object has been 
set in the currently notified editing session. Note that this method 
can only be saftely called when processing an `IS_EDITING` 
notification. 
*/ 
void is_set(Key) const; 

/*!
Returns true if the object has been 
set in the currently notified editing session. 
*/ 
void is_new(Key) const; 

/*!
Returns true if the object has been 
set in the currently notified editing session. 
*/ 
void clear(Key) const; 

/// @}

}; /* end Kinetic::ActiveObjectsTable */

} /* end namespace KineticConcepts */
