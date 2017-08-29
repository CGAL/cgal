namespace CGAL {


namespace IO{
/*!
  \ingroup PkgIOstreams

All classes in the \cgal `Kernel` provide input and output operators for
IOStreams.  The basic task of such an operator is to produce a
representation of an object that can be written as a sequence of
characters on devices as a console, a file, or a pipe. The enum `Mode` distinguish between three different printing formats.

In  `ASCII` mode, numbers 
e.g. the coordinates of a point or
the coefficients of a line, are written
in a machine independent format. 
In <span class="textsc">BINARY</span> mode, data are written
in a binary format, e.g. a double is represented
as a sequence of four byte. The format depends on the machine.
 The mode  <span class="textsc">PRETTY</span>
serves mainly for debugging as the type of the geometric
object is written, as well as the data defining the object. For example
for a point at the origin with %Cartesian double coordinates, the output
would be `PointC2(0.0, 0.0)`.  At the moment \cgal does not
provide input operations for pretty printed data. By default a stream
is in <span class="textsc">Ascii</span> mode.

\sa `CGAL::set_mode()`
\sa `CGAL::set_ascii_mode()`
\sa `CGAL::set_binary_mode()`
\sa `CGAL::set_pretty_mode()`
\sa `CGAL::get_mode()`
\sa `CGAL::is_ascii()`
\sa `CGAL::is_binary()`
\sa `CGAL::is_pretty()`
*/
enum Mode { ASCII = 0, BINARY, PRETTY };
}

/*!
\ingroup PkgIOstreams

returns the printing mode of the %IO stream `s`.

\sa `CGAL::IO::Mode`
\sa `CGAL::set_mode()`
\sa `CGAL::set_ascii_mode()`
\sa `CGAL::set_binary_mode()`
\sa `CGAL::set_pretty_mode()`
\sa `CGAL::is_ascii()`
\sa `CGAL::is_binary()`
\sa `CGAL::is_pretty()`


*/
IO::Mode get_mode(std::ios& s);

/*!
\ingroup PkgIOstreams

sets the mode of the %IO stream `s` to be the `IO::ASCII` mode.
Returns the previous mode of `s`.


\sa `CGAL::IO::Mode`
\sa `CGAL::set_mode()`
\sa `CGAL::set_binary_mode()`
\sa `CGAL::set_pretty_mode()`
\sa `CGAL::get_mode()`
\sa `CGAL::is_ascii()`
\sa `CGAL::is_binary()`
\sa `CGAL::is_pretty()`
*/
IO::Mode set_ascii_mode(std::ios& s);

/*!
\ingroup PkgIOstreams

\sa `CGAL::IO::Mode`
\sa `CGAL::set_mode()`
\sa `CGAL::set_ascii_mode()`
\sa `CGAL::set_pretty_mode()`
\sa `CGAL::get_mode()`
\sa `CGAL::is_ascii()`
\sa `CGAL::is_binary()`
\sa `CGAL::is_pretty()`

sets the mode of the %IO stream `s` to be the `IO::BINARY` mode.
Returns the previous mode of `s`.
*/
IO::Mode set_binary_mode(std::ios& s);

/*!
\ingroup PkgIOstreams

sets the printing mode of the %IO stream `s`.

\sa `CGAL::IO::Mode`
\sa `CGAL::set_ascii_mode()`
\sa `CGAL::set_binary_mode()`
\sa `CGAL::set_pretty_mode()`
\sa `CGAL::get_mode()`
\sa `CGAL::is_ascii()`
\sa `CGAL::is_binary()`
\sa `CGAL::is_pretty()`
*/
IO::Mode set_mode(std::ios& s, IO::Mode m);

/*!
\ingroup PkgIOstreams

sets the mode of the %IO stream `s` to be the `IO::PRETTY` mode.
Returns the previous mode of `s`.

\sa `CGAL::IO::Mode`
\sa `CGAL::set_mode()`
\sa `CGAL::set_ascii_mode()`
\sa `CGAL::set_binary_mode()`
\sa `CGAL::get_mode()`
\sa `CGAL::is_ascii()`
\sa `CGAL::is_binary()`
\sa `CGAL::is_pretty()`

*/
IO::Mode set_pretty_mode(std::ios& s);


/*!
\ingroup PkgIOstreams

The definition of `Input_rep` is completely symmetric to `Output_rep`. 

*/
template< typename T, typename F >
class Input_rep {

}; /* end Input_rep */

/*!
\ingroup PkgIOstreams

The purpose of `Output_rep` is to provide a way to control output formatting that works independently of the object's stream output operator. 

If you dont specialize `Output_rep` for `T`, `T`'s stream output operator is called from within `Output_rep`, by default. If you want another behaviour for your type `T`, you have to provide a specialization for that type. Furthermore, you can provide specializations with a second template parameter (a formatting tag). The second template parameter defaults to `Null_tag` and means *default behaviour*. 

Specializations of `Output_rep` should provide the following features: 

\code{.cpp} 

template< class F > 
struct Output_rep< Some_type, F > { 
  static const bool is_specialized = true;
  Output_rep( const Some_type& t );
  std::ostream& operator()( std::ostream& out ) const;
}; 

\endcode 

You can also specialize for a formatting tag `F`. 

The constant `is_specialized` can be tested by meta-programming tools to
verify that a given type can be used with `oformat()`. Its value has to be
`true` in a specialization of `Output_rep`. When there is no specialization
for a type, the class template `Output_rep` defines `is_specialized` to the
default value `false`.

*/
template< typename T, typename F >
class Output_rep {
}; /* end Output_rep */


/*!
\ingroup PkgIOstreams

checks if the %IO stream `s` is in `IO::ASCII` mode.

\sa `CGAL::IO::Mode`
\sa `CGAL::set_mode()`
\sa `CGAL::set_ascii_mode()`
\sa `CGAL::set_binary_mode()`
\sa `CGAL::set_pretty_mode()`
\sa `CGAL::get_mode()`
\sa `CGAL::is_binary()`
\sa `CGAL::is_pretty()`
*/
bool is_ascii(std::ios& s);


/*!
\ingroup PkgIOstreams

checks if the %IO stream `s` is in `IO::BINARY` mode.

\sa `CGAL::IO::Mode`
\sa `CGAL::set_mode()`
\sa `CGAL::set_ascii_mode()`
\sa `CGAL::set_binary_mode()`
\sa `CGAL::set_pretty_mode()`
\sa `CGAL::get_mode()`
\sa `CGAL::is_ascii()`
\sa `CGAL::is_pretty()`
*/
bool is_binary(std::ios& s);

/*!
\ingroup PkgIOstreams

checks if the %IO stream `s` is in `IO::PRETTY` mode.

\sa `CGAL::IO::Mode`
\sa `CGAL::set_mode()`
\sa `CGAL::set_ascii_mode()`
\sa `CGAL::set_binary_mode()`
\sa `CGAL::set_pretty_mode()`
\sa `CGAL::get_mode()`
\sa `CGAL::is_ascii()`
\sa `CGAL::is_binary()`

*/
bool is_pretty(std::ios& s);

/*!
\ingroup PkgIOstreams

Convenience function to construct an output representation (`Output_rep`) for type `T`. 

Generic IO for type `T`.
*/
template <class T> Output_rep<T> oformat( const T& t);

/*!
\ingroup PkgIOstreams

The definition of this function is completely symmetric to `oformat()`.
*/
template <class T> Input_rep<T> iformat( const T& t);


/*!
\ingroup PkgIOstreams

Convenience function to construct an output representation (`Output_rep`) for type `T`. 

Generic IO for type `T` with formatting tag.
*/
template <class T, typename F> Output_rep<T,F> oformat( const T& t, F );

/*!
\ingroup IOstreamOperators

\brief Inserts object `c` in the stream `os`. Returns `os`.

\cgal defines output operators for classes that are derived 
from the class `ostream`. This allows to write to ostreams 
as `cout` or `cerr`, as well as to `std::ostringstream` 
and `std::ofstream`. 
The output operator is defined for all classes in the \cgal `Kernel` and for the class `Color` as well. 

\sa `CGAL::set_mode()`
\sa `CGAL::set_ascii_mode()`
\sa `CGAL::set_binary_mode()`
\sa `CGAL::set_pretty_mode()`
\sa `CGAL::get_mode()`
\sa `CGAL::is_ascii()`
\sa `CGAL::is_binary()`
\sa `CGAL::is_pretty()`
*/
ostream& operator<<(ostream& os, Class c);

/*!
\ingroup IOstreamOperators

\brief \cgal defines input operators for classes that are derived
from the class `istream`. This allows to read from istreams
as `std::cin`, as well as from `std::istringstream` and `std::ifstream`.
The input operator is defined for all classes in the \cgal `Kernel`.



\sa `CGAL::set_mode()`
\sa `CGAL::set_ascii_mode()`
\sa `CGAL::set_binary_mode()`
\sa `CGAL::set_pretty_mode()`
\sa `CGAL::get_mode()`
\sa `CGAL::is_ascii()`
\sa `CGAL::is_binary()`
\sa `CGAL::is_pretty()`
*/
istream& operator>>(istream& is, Class c);

}
