
namespace CGAL {

/*!
\ingroup PkgIOstreams

The class `Ostream_iterator` is an output iterator adaptor for the 
output stream class `Stream` and value type `T`. 

### Operations ###

`o` fulfills the requirements for an output iterator. 

### Implementation ###

The `operator*()` in class `Ostream_iterator` uses a proxy class. 

*/
template< typename T, typename Stream >
class Ostream_iterator {
public:

/// \name Creation 
/// @{

/*! 
creates an output 
iterator `o` writing to `s`. 
*/ 
Ostream_iterator( Stream& s); 

/// @}

}; /* end Ostream_iterator */
} /* end namespace CGAL */
