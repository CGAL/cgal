namespace Kinetic {

/*!
\ingroup PkgKdsFrameworkOtherConcepts
\cgalConcept

The concept `Key` is a unique identifier for something in some sort of
table. In general, they can be only created by the table and are
returned when a appropriate `new_foo()` method is called on the
table. There are two classes of values for a `Key`, valid and
invalid. The latter cannot refer to something in a table. Use the
method `is_valid()` to differentiate.

\cgalHasModel `CGAL::Kinetic::Simulator::Event_key`
\cgalHasModel `CGAL::Kinetic::Active_objects_vector::Key`

*/

class Key {
public:

/*!
The default constructor is guaranteed to construct an invalid key (i.e.\ one which is false when cast to a bool. 
*/ 
Key(); 

/// \name Operations 
/// @{

/*!
This method returns false if the key 
was created using the default constructor or was otherwise created 
to be invalid. 
*/ 
bool is_valid() const; 

/*!
Write a text description of the key to a standard stream.
*/ 
std::ostream& operator<<(std::ostream&, Key);

/// @}

}; /* end Key */

} /* end namespace Kinetic */
