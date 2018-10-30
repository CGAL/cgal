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


//! \ingroup PkgIOstreams
//! \brief `read_point_WKT()` fills a `Point` from a  WKT stream. The first line starting with POINT 
//! in the stream will be used. 
//! 
//! \tparam Point can be a `CGAL::Point_2` or `CGAL::Point_3`.
//! \attention Only Cartesian Kernels are supported.
//! 
//! \see `CGAL::Point_2`
//! \see `CGAL::Point_3`
template<typename Point>
std::istream&
read_point_WKT( std::istream& in,
                Point& point );

//! \ingroup PkgIOstreams
//! \brief `read_multi_point_WKT()` overwrites the content of a `MultiPoint` 
//! with the first line starting with MULTIPOINT in the stream.
//! 
//! \tparam MultiPoint must be a model of `RandomAccessRange` of `CGAL::Point_2` or `CGAL::Point_3`, 
//! and have:
//! - a function `push_back()` that takes the same point type,
//! - a function `clear()`,
//! - a function `resize()` that takes an `size_type`
//! - an `operator[]()` that takes a `size_type`.
//! 
//! \attention Only Cartesian Kernels are supported.
//! 
//! \see `CGAL::Point_2`
//! \see `CGAL::Point_3`
template<typename MultiPoint>
std::istream&
read_multi_point_WKT( std::istream& in,
                      MultiPoint& mp );


//! \ingroup PkgIOstreams
//! \brief `read_linestring_WKT()` fills a `Linestring` from a WKT stream.
//! The first line starting with LINESTRING in the stream will be used.
//!
//! \tparam Linestring must be a model of `RandomAccessRange` of `CGAL::Point_2`, 
//! and have:
//! - a function `push_back()` that takes a `CGAL::Point_2`.
//! - a function `clear()`,
//! - a function `resize()` that takes an `size_type`
//! - an `operator[]()` that takes a `size_type`.
//! \attention Only Cartesian Kernels are supported.
//! \see `CGAL::Point_2`
template<typename LineString>
std::istream&
read_linestring_WKT( std::istream& in,
                     LineString& polyline );

//! \ingroup PkgIOstreams
//! \brief `read_multi_linestring_WKT()` overwrites the content of a `MultiLineString` 
//! with the first line starting with MULTILINESTRING in the stream.
//!
//! \tparam MultiLineString must be a model of `RandomAccessRange` of `Linestring`, 
//! and have:
//! - a function `push_back()` that takes a `Linestring`,
//! - a function `clear()`,
//! - a function `resize()` that takes an `size_type`
//! - an `operator[]()` that takes a `size_type`.
//! \attention Only Cartesian Kernels are supported.
//! 
//! \see `CGAL::Point_2`
template<typename MultiLineString>
std::istream&
read_multi_linestring_WKT( std::istream& in,
                           MultiLineString& mls );

//! \ingroup PkgIOstreams
//! \brief `read_polygon_WKT()` fills `Polygon` from a WKT stream.
//! The first line starting with POLYGON in the stream will be used.
//! 
//! \tparam Polygon is a `CGAL::General_polygon_with_holes_2`.
//! \attention Only Cartesian Kernels are supported.
//! 
//! \see `CGAL::General_polygon_with_holes_2`
template<typename Polygon>
std::istream&
read_polygon_WKT( std::istream& in,
                  Polygon& polygon );

//! \ingroup PkgIOstreams
//! \brief `read_multi_polygon_WKT()` overwrites the content of a `MultiPolygon` 
//! with the first line starting with MULTIPOLYGON in the stream.
//!
//! \tparam Multipolygon must be a model of `RandomAccessRange` of `CGAL::General_polygon_with_holes_2`, 
//! and have:
//! - a function `push_back()` that takes a `CGAL::General_polygon_with_holes_2`,
//! - a function `clear()`,
//! - a function `resize()` that takes an `size_type`
//! - an `operator[]()` that takes a `size_type`.
//! \attention Only Cartesian Kernels are supported.
//! \see `CGAL::General_polygon_with_holes_2`

template<typename MultiPolygon>
std::istream&
read_multi_polygon_WKT( std::istream& in,
                        MultiPolygon& polygons );

//! \ingroup PkgIOstreams
//! \brief `write_point_WKT()` writes `point` into a WKT stream.
//! \tparam Point is a `CGAL::Point_2`
//! \attention Only Cartesian Kernels are supported.
//! \see `CGAL::Point_2`
template<typename Point>
std::ostream&
write_point_WKT( std::ostream& out,
                 const Point& point );

//! \ingroup PkgIOstreams
//! \brief `write_polygon_WKT()` writes `poly` into a WKT stream.
//! \tparam Polygon  must be a `CGAL::General_polygon_with_holes_2`
//! \attention Only Cartesian Kernels are supported.
//! \see `CGAL::General_polygon_with_holes_2`
template<typename Polygon>
std::ostream&
write_polygon_WKT( std::ostream& out,
                   const Polygon& poly );

//! \ingroup PkgIOstreams
//! \brief `write_linestring_WKT()` writes the content of `ls` 
//! into a WKT stream.
//! \tparam LineString must be a `RandomAccessRange` of `CGAL::Point_2`,
//! and have:
//! - a function `push_back()` that takes a `CGAL::Point_2`,
//! - a function `clear()`,
//! - a function `resize()` that takes an `size_type`
//! - an `operator[]()` that takes a `size_type`.
//! \attention Only Cartesian Kernels are supported.
//!\see `CGAL::Point_2`
template<typename LineString>
std::ostream&
write_linestring_WKT( std::ostream& out,
                      LineString ls );

//! \ingroup PkgIOstreams
//! \brief `write_multi_point_WKT()` writes the content of `mp` 
//! into a WKT stream.
//! \tparam MultiPoint must be a `RandomAccessRange` of `CGAL::Point_2`,
//! and have:
//! - a function `push_back()` that takes a `CGAL::Point_2`,
//! - a function `clear()`,
//! - a function `resize()` that takes an `size_type`
//! - an `operator[]()` that takes a `size_type`.
//! \attention Only Cartesian Kernels are supported.
//!\see `CGAL::Point_2`
template<typename MultiPoint>
std::ostream&
write_multi_point_WKT( std::ostream& out,
                       MultiPoint& mp );

//! \ingroup PkgIOstreams
//! \brief `write_multi_polygon_WKT()` writes the content of `polygons` 
//! into a WKT stream.
//! \tparam MultiPolygon must be a `RandomAccessRange` of `CGAL::General_polygon_with_holes_2`,
//! and have:
//! - a function `push_back()` that takes a `CGAL::General_polygon_with_holes_2`,
//! - a function `clear()`,
//! - a function `resize()` that takes an `size_type`
//! - an `operator[]()` that takes a `size_type`.
//! \attention Only Cartesian Kernels are supported.
//!\see `CGAL::General_polygon_with_holes_2`
template<typename MultiPolygon>
std::ostream&
write_multi_polygon_WKT( std::ostream& out,
                         MultiPolygon& polygons );

//! \ingroup PkgIOstreams
//! \brief `write_multi_linestring_WKT()` writes the content of `mls` 
//! into a WKT stream. 
//! \tparam MultiLineString must be a `RandomAccessRange` of `LineString`,
//! and have:
//! - a function `push_back()` that takes a `Linestring`,
//! - a function `clear()`,
//! - a function `resize()` that takes an `size_type`
//! - an `operator[]()` that takes a `size_type`.
//! \attention Only Cartesian Kernels are supported.
//! \see `CGAL::write_linestring_WKT()`
template<typename MultiLineString>
std::ostream&
write_multi_linestring_WKT( std::ostream& out,
                            MultiLineString& mls );

//! \ingroup PkgIOstreams
//!  reads the content of a WKT stream and fills 
//! `points`, `polylines` and `polygons` with all the POINT, MULTIPOINT, 
//! LINESTRING, MULTILINESTRING, POLYGON and MULTIPOLYGON it finds in `input`.
//! \tparam MultiPoint must be a model of `RandomAccessRange` of `CGAL::Point_2` or `CGAL::Point_3`,
//! and have:
//! - a function `push_back()` that takes the same point type,
//! - a function `clear()`,
//! - a function `resize()` that takes an `size_type`
//! - an `operator[]()` that takes a `size_type`.
//! \tparam MultiLineString must be a `RandomAccessRange` of `Linestring`,
//! and have:
//! - a function `push_back()` that takes a `Linestring`,
//! - a function `clear()`,
//! - a function `resize()` that takes an `size_type`
//! - an `operator[]()` that takes a `size_type`.
//! \tparam MultiPolygon must be a model of `RandomAccessRange` of `CGAL::General_polygon_with_holes_2`,
//! and have:
//! - a function `push_back()` that takes a `CGAL::General_polygon_with_holes_2`,
//! - a function `clear()`,
//! - a function `resize()` that takes an `size_type`
//! - an `operator[]()` that takes a `size_type`.
//! \attention Only Cartesian Kernels are supported.
//! \see `CGAL::read_linestring_WKT()`
template<typename MultiPoint,
         typename MultiLineString,
         typename MultiPolygon>
std::istream&
read_WKT( std::istream& input,
                        MultiPoint& points,   
                        MultiLineString& polylines,
                        MultiPolygon& polygons);
}
