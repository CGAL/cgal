
namespace CGAL {

/*!
\ingroup PkgIOstreams

The class `Verbose_ostream` can be used as an output stream. The stream 
output operator `<<` is defined for any type. The class 
`Verbose_ostream` stores in an internal state a stream and whether the 
output is active or not. If the state is active, the stream output 
operator `<<` uses the internal stream to output its argument. If 
the state is inactive, nothing happens. 

\cgalHeading{Example}

The class `Verbose_ostream` can be conveniently used to implement for 
example the `is_valid()` member function for triangulations or 
other complex data structures. 

\code{.cpp}
bool is_valid( bool verbose = false, int level = 0) { 
Verbose_ostream verr( verbose); 
verr << "Triangulation::is_valid( level = " << level << ')' << endl; 
verr << " Number of vertices = " << size_of_vertices() << endl; 
// ... 
} 
\endcode
*/

class Verbose_ostream {
public:

/// \name Creation 
/// @{

/*!
creates an output stream with state set to 
`active` that writes to the stream `out`. 
*/ 
Verbose_ostream( bool active = false, 
std::ostream& out = std::cerr); 

/// @} 

/// \name Operations 
/// @{

/*!

*/ 
template < class T > 
Verbose_ostream& operator<<( const T& t); 

/// @}

}; /* end Verbose_ostream */
} /* end namespace CGAL */
