// Copyright (c) 1997-2020
//
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_IO_H
#define CGAL_IO_H

#include <CGAL/disable_warnings.h>

#include <CGAL/IO/io_tags.h>
#include <CGAL/IO/Color.h>

#include <CGAL/assertions.h>
#include <CGAL/Fraction_traits.h>
#include <CGAL/tags.h>

#include <cstdio>
#include <cctype>
#include <string>
#include <locale>
#include <iostream>
#include <optional>
#include <variant>

namespace CGAL {

namespace IO {

class Static
{
public:
  static int get_mode()
  {
    static const int mode = std::ios::xalloc();
    return mode;
  }
};

/*!
\ingroup PkgStreamSupportEnumRef

All classes in the \cgal `Kernel` provide input and output operators for
IOStreams.  The basic task of such an operator is to produce a
representation of an object that can be written as a sequence of
characters on devices as a console, a file, or a pipe. The enum `Mode` distinguish between three different printing formats.

In `ASCII` mode, numbers
e.g. the coordinates of a point or
the coefficients of a line, are written
in a machine independent format.
In `BINARY` mode, data are written
in a binary format, e.g. a double is represented
as a sequence of four byte. The format depends on the machine.
 The mode `PRETTY` serves mainly for debugging as the type of the geometric
object is written, as well as the data defining the object. For example
for a point at the origin with %Cartesian double coordinates, the output
would be `PointC2(0.0, 0.0)`.  At the moment \cgal does not
provide input operations for pretty printed data. By default a stream
is in `ASCII` mode.

\sa `CGAL::IO::set_mode()`
\sa `CGAL::IO::set_ascii_mode()`
\sa `CGAL::IO::set_binary_mode()`
\sa `CGAL::IO::set_pretty_mode()`
\sa `CGAL::IO::get_mode()`
\sa `CGAL::IO::is_ascii()`
\sa `CGAL::IO::is_binary()`
\sa `CGAL::IO::is_pretty()`
*/
enum Mode {ASCII = 0, PRETTY, BINARY};

} // namespace IO

#ifdef DOXYGEN_RUNNING
/*!
\ingroup IOstreamOperators

\brief Inserts object `c` in the stream `os`. Returns `os`.
\cgal defines output operators for classes that are derived
from the class `ostream`. This allows to write to ostreams
as `cout` or `cerr`, as well as to `std::ostringstream`
and `std::ofstream`.
The output operator is defined for all classes in the \cgal `Kernel` and for the class `Color` as well.

\sa `CGAL::IO::set_mode()`
\sa `CGAL::IO::set_ascii_mode()`
\sa `CGAL::IO::set_binary_mode()`
\sa `CGAL::IO::set_pretty_mode()`
\sa `CGAL::IO::get_mode()`
\sa `CGAL::IO::is_ascii()`
\sa `CGAL::IO::is_binary()`
\sa `CGAL::IO::is_pretty()`
*/
ostream& operator<<(ostream& os, Class c);

/*!
\ingroup IOstreamOperators

\brief \cgal defines input operators for classes that are derived
from the class `istream`. This allows to read from istreams
as `std::cin`, as well as from `std::istringstream` and `std::ifstream`.
The input operator is defined for all classes in the \cgal `Kernel`.

\sa `CGAL::IO::set_mode()`
\sa `CGAL::IO::set_ascii_mode()`
\sa `CGAL::IO::set_binary_mode()`
\sa `CGAL::IO::set_pretty_mode()`
\sa `CGAL::IO::get_mode()`
\sa `CGAL::IO::is_ascii()`
\sa `CGAL::IO::is_binary()`
\sa `CGAL::IO::is_pretty()`
*/
istream& operator>>(istream& is, Class c);
#endif

template <typename Dummy>
struct IO_rep_is_specialized_aux
{
  static const bool is_specialized = true;
};

template< class Dummy >
const bool IO_rep_is_specialized_aux<Dummy>::is_specialized;

template <typename Dummy>
struct IO_rep_is_not_specialized_aux
{
  static const bool is_specialized = false;
};

template< class Dummy >
const bool IO_rep_is_not_specialized_aux<Dummy>::is_specialized;

typedef IO_rep_is_specialized_aux<void> IO_rep_is_specialized;
typedef IO_rep_is_not_specialized_aux<void> IO_rep_is_not_specialized;

/*!
\ingroup PkgStreamSupportRef

The purpose of `Output_rep` is to provide a way to control output formatting that works independently of the object's stream output operator.

If you dont specialize `Output_rep` for `T`, `T`'s stream output operator is called from within `Output_rep`, by default. If you want another behavior for your type `T`, you have to provide a specialization for that type. Furthermore, you can provide specializations with a second template parameter (a formatting tag). The second template parameter defaults to `Null_tag` and means *default behavior*.

Specializations of `Output_rep` should provide the following features:

\code{.cpp}

template< class F >
struct Output_rep< Some_type, F > {
  static const bool is_specialized = true;
  Output_rep( const Some_type& t );
  std::ostream& operator()( std::ostream& os ) const;
};

\endcode

You can also specialize for a formatting tag `F`.

The constant `is_specialized` can be tested by meta-programming tools to
verify that a given type can be used with `oformat()`. Its value has to be
`true` in a specialization of `Output_rep`. When there is no specialization
for a type, the class template `Output_rep` defines `is_specialized` to the
default value `false`.
*/
template <class T, class F = ::CGAL::Null_tag >
class Output_rep
  : public IO_rep_is_not_specialized
{
  const T& t;

public:
  //! initialize with a const reference to \a t.
  Output_rep( const T& tt, F = {}) : t(tt) {}
  //! perform the output, calls \c operator\<\< by default.
  std::ostream& operator()( std::ostream& os) const { return (os << t); }
};

template <class T, class F>
class Output_rep<std::optional<T>, F>
{
  const std::optional<T>& t;

public:
  Output_rep( const std::optional<T>& tt) : t(tt) {}
  std::ostream& operator()( std::ostream& os) const
  {
    if (t==std::nullopt) return (os << "--");
    return (os << t.value());
  }
};

template <class ... T, class F>
class Output_rep<std::variant<T...>, F>
{
   const std::variant<T...>& t;

public:
  Output_rep( const std::variant<T...>& tt) : t(tt) {}
  std::ostream& operator()( std::ostream& os) const
  {
    std::visit([&os](auto&& v) { os << v; }, t);
    return os;
  }
};

/*!
  \relates Output_rep
  \brief stream output of the \c Output_rep calls its \c operator().

  \cgal defines output operators for classes that are derived from the class `ostream`.
  This enables to write to ostreams as `cout` or `cerr`, as well as to `std::ostringstream`
  and `std::ofstream`.
  The output operator is defined for all classes in the \cgal `Kernel` and for the class `Color` as well.
*/
template <class T, class F>
std::ostream& operator<<( std::ostream& os, Output_rep<T,F> rep) { return rep(os); }

namespace IO {

/*!
\ingroup PkgStreamSupportRef

Convenience function to construct an output representation (`Output_rep`) for type `T`.

Generic IO for type `T`.
*/
template <class T>
Output_rep<T> oformat(const T& t) { return Output_rep<T>(t); }

/*!
\ingroup PkgStreamSupportRef

Convenience function to construct an output representation (`Output_rep`) for type `T`.

Generic IO for type `T` with formatting tag.
*/
template <class T, class F>
Output_rep<T,F> oformat( const T& t, F format) {
  if constexpr (std::is_constructible_v<Output_rep<T,F>, const T&, F>)
    return Output_rep<T,F>(t, format);
  else
    return Output_rep<T,F>(t);
}

} // namespace IO

/*!
\ingroup PkgStreamSupportRef

The definition of `Input_rep` is completely symmetric to `Output_rep`.
*/
template <class T>
class Input_rep
  : public IO_rep_is_not_specialized
{
  T& t;

public:
  //! initialize with a reference to \a t.
  Input_rep( T& tt) : t(tt) {}

  //! perform the input, calls \c operator\>\> by default.
  std::istream& operator()( std::istream& is) const { return (is >> t); }
};

template <class T>
class Input_rep<std::optional<T>>
{
  std::optional<T>& t;

public:
  //! initialize with a reference to \a t.
  Input_rep( std::optional<T>& tt) : t(tt) {}

  //! perform the input, calls \c operator\>\> by default.
  std::istream& operator()( std::istream& is) const {
    T v;
    if(is >> v) t = v;
    return is;
  }
};

#if CGAL_FORCE_IFORMAT_DOUBLE || \
  ( ( _MSC_VER > 1600 ) && ( _MSC_VER < 1910 ) && (! defined( CGAL_NO_IFORMAT_DOUBLE )) )

template <>
class Input_rep<double>
  : public IO_rep_is_specialized
{
  double& t;

public:
  //! initialize with a reference to \a t.
  Input_rep( double& tt) : t(tt) {}

  std::istream& operator()( std::istream& is) const
  {
    typedef std::istream istream;
    typedef istream::char_type char_type;
    typedef istream::int_type int_type;
    typedef istream::traits_type traits_type;

    std::string buffer;
    buffer.reserve(32);

    char_type c;
    do
    {
      const int_type i = is.get();
      if(i == traits_type::eof())
        return is;

      c = static_cast<char_type>(i);
    }
    while (std::isspace(c));

    if(c == '-')
    {
      buffer += '-';
    }
    else if(c != '+')
    {
      is.unget();
    }

    for(;;)
    {
      const int_type i = is.get();
      if(i == traits_type::eof())
      {
        is.clear(is.rdstate() & ~std::ios_base::failbit);
        break;
      }

      c = static_cast<char_type>(i);
      if(std::isdigit(c) || (c =='.') || (c =='E') || (c =='e') || (c =='+') || (c =='-'))
      {
        buffer += c;
      }
      else
      {
        is.unget();
        break;
      }
    }

    if(sscanf_s(buffer.c_str(), "%lf", &t) != 1)
    {
      // if a 'buffer' does not contain a double, set the fail bit.
      is.setstate(std::ios_base::failbit);
    }

    return is;
  }
};

template <>
class Input_rep<float>
{
  float& t;

public:
  //! initialize with a reference to \a t.
  Input_rep( float& tt) : t(tt) {}

  std::istream& operator()( std::istream& is) const
  {
    typedef std::istream istream;
    typedef istream::char_type char_type;
    typedef istream::int_type int_type;
    typedef istream::traits_type traits_type;

    std::string buffer;
    buffer.reserve(32);

    char_type c;
    do
    {
      const int_type i = is.get();
      if(i == traits_type::eof())
        return is;

      c = static_cast<char_type>(i);
    }
    while (std::isspace(c));

    if(c == '-')
    {
      buffer += '-';
    }
    else if(c != '+')
    {
      is.unget();
    }

    for(;;)
    {
      const int_type i = is.get();
      if(i == traits_type::eof())
      {
        is.clear(is.rdstate() & ~std::ios_base::failbit);
        break;
      }

      c = static_cast<char_type>(i);
      if(std::isdigit(c) || (c =='.') || (c =='E') || (c =='e') || (c =='+') || (c =='-'))
      {
        buffer += c;
      }
      else
      {
        is.unget();
        break;
      }
    }

    if(sscanf_s(buffer.c_str(), "%f", &t) != 1)
    {
      // if a 'buffer' does not contain a double, set the fail bit.
      is.setstate(std::ios_base::failbit);
    }

    return is;
  }
};
#endif

/*! \relates Input_rep
    \brief stream input to the \c Input_rep calls its \c operator().

\brief \cgal defines input operators for classes that are derived
from the class `istream`. This allows to read from istreams
as `std::cin`, as well as from `std::istringstream` and `std::ifstream`.
The input operator is defined for all classes in the \cgal `Kernel`.
*/
template <class T>
std::istream& operator>>( std::istream& is, Input_rep<T> rep) { return rep(is); }

namespace IO {

/*!
\ingroup PkgStreamSupportRef

The definition of this function is completely symmetric to `oformat()`.
*/
template <class T>
Input_rep<T> iformat( T& t) { return Input_rep<T>(t); }

} // namespace IO

template <class T, class F = Null_tag >
class Benchmark_rep
{
  const T& t;

public:
  //! initialize with a const reference to \a t.
  Benchmark_rep( const T& tt) : t(tt) {}
  //! perform the output, calls \c operator\<\< by default.
  std::ostream& operator()( std::ostream& os) const { return os << t; }

  // static function to get the benchmark name
  static std::string get_benchmark_name() { return ""; }
};

template <class T, class F>
std::ostream& operator<<( std::ostream& os, Benchmark_rep<T,F> rep) { return rep(os); }

namespace IO {

template <class T>
Benchmark_rep<T> bmformat( const T& t) { return Benchmark_rep<T>(t); }

template <class T, class F>
Benchmark_rep<T,F> bmformat( const T& t, F) { return Benchmark_rep<T,F>(t); }

/*!
\ingroup PkgStreamSupportRef

returns the printing mode of the %IO stream `s`.

\link PkgStreamSupportEnumRef `CGAL::IO::Mode`\endlink
\sa `CGAL::IO::set_mode()`
\sa `CGAL::IO::set_ascii_mode()`
\sa `CGAL::IO::set_binary_mode()`
\sa `CGAL::IO::set_pretty_mode()`
\sa `CGAL::IO::is_ascii()`
\sa `CGAL::IO::is_binary()`
\sa `CGAL::IO::is_pretty()`
*/
inline Mode get_mode(std::ios& s)
{
  return static_cast<Mode>(s.iword(Static::get_mode()));
}

/*!
\ingroup PkgStreamSupportRef

sets the mode of the %IO stream `s` to be the `ASCII` mode.
Returns the previous mode of `s`.

\link PkgStreamSupportEnumRef `CGAL::IO::Mode`\endlink
\sa `CGAL::IO::set_mode()`
\sa `CGAL::IO::set_binary_mode()`
\sa `CGAL::IO::set_pretty_mode()`
\sa `CGAL::IO::get_mode()`
\sa `CGAL::IO::is_ascii()`
\sa `CGAL::IO::is_binary()`
\sa `CGAL::IO::is_pretty()`
*/
inline Mode set_ascii_mode(std::ios& s)
{
  Mode m = get_mode(s);
  s.iword(Static::get_mode()) = ASCII;
  return m;
}

/*!
\ingroup PkgStreamSupportRef

sets the mode of the %IO stream `s` to be the `BINARY` mode.
Returns the previous mode of `s`.

\link PkgStreamSupportEnumRef `CGAL::IO::Mode`\endlink
\sa `CGAL::IO::set_mode()`
\sa `CGAL::IO::set_ascii_mode()`
\sa `CGAL::IO::set_pretty_mode()`
\sa `CGAL::IO::get_mode()`
\sa `CGAL::IO::is_ascii()`
\sa `CGAL::IO::is_binary()`
\sa `CGAL::IO::is_pretty()`
*/
inline Mode set_binary_mode(std::ios& s)
{
  Mode m = get_mode(s);
  s.iword(Static::get_mode()) = BINARY;
  return m;
}

/*!
\ingroup PkgStreamSupportRef

sets the mode of the %IO stream `s` to be the `PRETTY` mode.
Returns the previous mode of `s`.

\link PkgStreamSupportEnumRef `CGAL::IO::Mode`\endlink
\sa `CGAL::IO::set_mode()`
\sa `CGAL::IO::set_ascii_mode()`
\sa `CGAL::IO::set_binary_mode()`
\sa `CGAL::IO::get_mode()`
\sa `CGAL::IO::is_ascii()`
\sa `CGAL::IO::is_binary()`
\sa `CGAL::IO::is_pretty()`
*/
inline Mode set_pretty_mode(std::ios& s)
{
  Mode m = get_mode(s);
  s.iword(Static::get_mode()) = PRETTY;
  return m;
}

/*!
\ingroup PkgStreamSupportRef

sets the printing mode of the %IO stream `s`.

\link PkgStreamSupportEnumRef `CGAL::IO::Mode`\endlink
\sa `CGAL::IO::set_ascii_mode()`
\sa `CGAL::IO::set_binary_mode()`
\sa `CGAL::IO::set_pretty_mode()`
\sa `CGAL::IO::get_mode()`
\sa `CGAL::IO::is_ascii()`
\sa `CGAL::IO::is_binary()`
\sa `CGAL::IO::is_pretty()`
*/
inline Mode set_mode(std::ios& s, Mode m)
{
  Mode old = get_mode(s);
  s.iword(Static::get_mode()) = m;
  return old;
}

/*!
\ingroup PkgStreamSupportRef

checks if the %IO stream `s` is in `PRETTY` mode.

\link PkgStreamSupportEnumRef `CGAL::IO::Mode`\endlink
\sa `CGAL::IO::set_mode()`
\sa `CGAL::IO::set_ascii_mode()`
\sa `CGAL::IO::set_binary_mode()`
\sa `CGAL::IO::set_pretty_mode()`
\sa `CGAL::IO::get_mode()`
\sa `CGAL::IO::is_ascii()`
\sa `CGAL::IO::is_binary()`
*/
inline bool is_pretty(std::ios& s) { return s.iword(Static::get_mode()) == PRETTY; }

/*!
\ingroup PkgStreamSupportRef

checks if the %IO stream `s` is in `ASCII` mode.

\link PkgStreamSupportEnumRef `CGAL::IO::Mode`\endlink
\sa `CGAL::IO::set_mode()`
\sa `CGAL::IO::set_ascii_mode()`
\sa `CGAL::IO::set_binary_mode()`
\sa `CGAL::IO::set_pretty_mode()`
\sa `CGAL::IO::get_mode()`
\sa `CGAL::IO::is_binary()`
\sa `CGAL::IO::is_pretty()`
*/
inline bool is_ascii(std::ios& s) { return s.iword(Static::get_mode()) == ASCII; }

/*!
\ingroup PkgStreamSupportRef

checks if the %IO stream `s` is in `BINARY` mode.

\link PkgStreamSupportEnumRef `CGAL::IO::Mode`\endlink
\sa `CGAL::IO::set_mode()`
\sa `CGAL::IO::set_ascii_mode()`
\sa `CGAL::IO::set_binary_mode()`
\sa `CGAL::IO::set_pretty_mode()`
\sa `CGAL::IO::get_mode()`
\sa `CGAL::IO::is_ascii()`
\sa `CGAL::IO::is_pretty()`
*/
inline bool is_binary(std::ios& s) { return s.iword(Static::get_mode()) == BINARY; }

} // namespace IO

template < class T >
inline void write(std::ostream& os, const T& t, const io_Read_write&)
{
  os.write(reinterpret_cast<const char*>(&t), sizeof(t));
}

template < class T >
inline void write(std::ostream& os, const T& t, const io_Operator&)
{
  os << IO::oformat(t);
}

template < class T >
inline void write(std::ostream& os, const T& t, const io_Extract_insert&)
{
  insert(os, t);
}

template < class T >
inline void write(std::ostream& os, const T& t)
{
  write(os, t, typename Io_traits<T>::Io_tag());
}

template < class T >
inline void read(std::istream& is, T& t, const io_Read_write&)
{
  is.read(reinterpret_cast<char*>(&t), sizeof(t));
}

template < class T >
inline void read(std::istream& is, T& t, const io_Operator&)
{
  is >> IO::iformat(t);
}

template < class T >
inline void read(std::istream& is, T& t, const io_Extract_insert&)
{
  extract(is, t);
}

template < class T >
inline void read(std::istream& is, T& t)
{
  read(is, t, typename Io_traits<T>::Io_tag());
}

namespace IO {

inline std::ostream& operator<<( std::ostream& os, const Color& col)
{
  switch(get_mode(os))
  {
    case ASCII :
      return os << static_cast<int>(col.red())   << ' '
                << static_cast<int>(col.green()) << ' '
                << static_cast<int>(col.blue()) << ' '
                << static_cast<int>(col.alpha());
    case BINARY :
      os.write(reinterpret_cast<const char*>(col.to_rgba().data()), 4);
      return os;
    default:
      return os << "Color(" << static_cast<int>(col.red()) << ", "
                << static_cast<int>(col.green()) << ", "
                << static_cast<int>(col.blue()) << ", "
                << static_cast<int>(col.alpha()) << ")";
  }
}

inline std::istream &operator>>(std::istream &is, Color& col)
{
  unsigned char r = 0, g = 0, b = 0, a = 0;
  int ir = 0, ig = 0, ib = 0, ia = 0;

  switch(get_mode(is))
  {
    case ASCII :
      is >> ir >> ig >> ib >> ia;
      r = (unsigned char)ir;
      g = (unsigned char)ig;
      b = (unsigned char)ib;
      a = (unsigned char)ia;
      break;
    case BINARY :
      read(is, r);
      read(is, g);
      read(is, b);
      read(is, a);
      break;
    default:
      std::cerr << "" << std::endl;
      std::cerr << "Stream must be in ASCII or binary mode" << std::endl;
      break;
  }

  col = Color(r,g,b,a);
  return is;
}

inline const char* mode_name( IO::Mode m )
{
  static const char* const names[] = {"ASCII", "PRETTY", "BINARY" };
  CGAL_assertion( IO::ASCII <= m && m <= IO::BINARY );
  return names[m];
}

namespace internal {

template <class P> constexpr auto has_exact(int) -> decltype(exact(P()), bool()) { return true; }
template <class P> constexpr bool has_exact(...) { return false; }

}

template <class P>
auto
serialize(const P& p) {
  if constexpr (internal::has_exact<P>(0)) {
    return exact(p);
  } else {
    return p;
  }
}



} // IO namespace

#ifndef CGAL_NO_DEPRECATED_CODE
using IO::oformat;
using IO::iformat;
using IO::bmformat;
using IO::get_mode;
using IO::set_ascii_mode;
using IO::set_binary_mode;
using IO::set_pretty_mode;
using IO::set_mode;
using IO::is_pretty;
using IO::is_ascii;
using IO::is_binary;
using IO::mode_name;
#endif

// From polynomial.h TODO: Where to put this?
inline void swallow(std::istream &is, char d)
{
  char c;
  do { is.get(c); } while(isspace(c));
  if(c != d)
  {
    std::stringstream msg;
    msg << "input error: expected '" << d << "' but got '" << c << "'";
    CGAL_error_msg( msg.str().c_str() );
  }
}

inline void swallow(std::istream &is, const std::string& s)
{
  std::string t;
  is >> t;
  if(s != t)
  {
    std::stringstream msg;
    msg << "input error: expected '" << s << "' but got '" << t << "'";
    CGAL_error_msg( msg.str().c_str() );
  }
}

namespace internal {

inline void eat_white_space(std::istream &is)
{
  std::istream::int_type c;
  do
  {
    c = is.peek();
    if(c == std::istream::traits_type::eof())
    {
      return;
    }
    else
    {
      std::istream::char_type cc= static_cast<std::istream::char_type>(c);
      if(std::isspace(cc, std::locale::classic()))
      {
        is.get();
        // since peek succeeded, this should too
        CGAL_assertion(!is.fail());
      }
      else
      {
        return;
      }
    }
  }
  while (true);
}

inline bool is_space(const std::istream& /*is*/, std::istream::int_type c)
{
  return (c == std::istream::traits_type::eof()) ||
          std::isspace(static_cast<std::istream::char_type>(c),
                       std::locale::classic() );
}

inline bool is_eof(const std::istream& /*is*/, std::istream::int_type c)
{
  return c == std::istream::traits_type::eof();
}

inline bool is_digit(const std::istream& /*is*/, std::istream::int_type c)
{
  CGAL_assertion(c != std::istream::traits_type::eof());
  return std::isdigit(static_cast<std::istream::char_type>(c),
                      std::locale::classic());
}

inline std::istream::int_type peek(std::istream& is)
{
  // Workaround for a bug in the version of libc++ that is shipped with
  // Apple-clang-3.2. See the long comment in the function
  // gmpz_new_read() in <CGAL/GMP/Gmpz_type.h>.

  if(is.eof())
    return std::istream::traits_type::eof();
  else
    return is.peek();
}

template <typename ET>
inline void read_float_or_quotient(std::istream & is, ET& et)
{
  is >> et;
}

template <typename Int, typename Rat>
inline void read_float_or_quotient(std::istream& is, Rat &z)
{
  // To build a rational from numerator and denominator. Hope that `Int`
  // and `Fraction_traits::(Numerator|Denominator)_type` are consistent...
  typename Fraction_traits<Rat>::Compose compose;

  // reads rational and floating point literals.
  const std::istream::char_type zero = '0';
  std::istream::int_type c;
  std::ios::fmtflags old_flags = is.flags();

  is.unsetf(std::ios::skipws);
  internal::eat_white_space(is);

  Int n(0); // unsigned number before '/' or '.'
  Int d(1); // number after '/', or denominator (fp-case)
  bool negative = false; // do we have a leading '-'?
  bool digits = false;   // for fp-case: are there any digits at all?

  c = internal::peek(is);
  if(c != '.')
  {
    // is there a sign?
    if(c == '-' || c == '+')
    {
      is.get();
      negative = (c == '-');
      internal::eat_white_space(is);
      c = internal::peek(is);
    }

    // read n (could be empty)
    while (!internal::is_eof(is, c) && internal::is_digit(is, c))
    {
      digits = true;
      n = n*10 + (c-zero);
      is.get();
      c = internal::peek(is);
    }

    // are we done?
    if(internal::is_eof(is, c) || internal::is_space(is, c))
    {
      is.flags(old_flags);
      if(digits && !is.fail())
        z = negative? compose(-n,1): compose(n,1);
      return;
    }
  }
  else
  {
    n = 0;
  }

  // now we have read n, we are not done, and c is the next character
  // in the stream
  if(c == '/' || c == '.')
  {
    is.get();
    if(c == '/')
    {
      // rational case
      is >> d;
      is.flags(old_flags);
      if(!is.fail())
        z = negative? compose(-n,d) : compose(n,d);
      return;
    }

    // floating point case; read number after '.' (may be empty)
    for(;;)
    {
      c = internal::peek(is);
      if(internal::is_eof(is, c) || !internal::is_digit(is, c))
        break;
      // now we have a digit
      is.get();
      digits = true;
      d *= 10;
      n = n*10 + (c-zero);
    }
  }

  // now we have read all digits after '.', and c is the next character;
  // read the exponential part (optional)
  int e = 0;
  if(c == 'e' || c == 'E')
  {
    is.get();
    is >> e;
  }

  // now construct the Gmpq
  if(!digits)
  {
    // illegal floating-point number
    is.setstate(std::ios_base::failbit);
    is.flags(old_flags);
    return;
  }

  // handle e
  if(e > 0)
    while (e--) n *= 10;
  else
    while (e++) d *= 10;

  is.flags(old_flags);
  if(!is.fail())
    z = (negative ? compose(-n,d) : compose(n,d));
}

} // namespace internal

} // namespace CGAL

#if CGAL_CAN_USE_CXX20_FORMAT
#  include <format>
#  include <sstream>

namespace std {

template <typename T, typename F, typename CharT>
struct formatter<CGAL::Output_rep<T, F>, CharT> : public std::formatter<std::basic_string<CharT>>
{
  constexpr auto parse(std::basic_format_parse_context<CharT>& ctx)
  {
    auto it = ctx.begin();
    const auto end = ctx.end();
    if(it == end)
      return it;
    if(*it != CharT('.')) {
      if(*it == CharT('}')) return it;
      throw std::format_error("formatter for CGAL::Output_rep only support precision, like `{:.6}`");
    }
    if(++it == end)
      throw std::format_error("Missing precision");
    if(*it < CharT('0') || *it > CharT('9'))
      throw std::format_error("Invalid value for precision");
    precision = *it - CharT('0');
    while(++it != end) {
      if(*it < CharT('0') || *it > CharT('9'))
        return it;
      precision = precision * 10 + (*it - CharT('0'));
    }
    return it;
  }

  template <typename FormatContext>
  auto format(const CGAL::Output_rep<T, F> &rep, FormatContext& ctx) const
  {
    std::basic_stringstream<CharT> ss;
    ss.precision(precision);
    ss << rep;
    return std::formatter<std::basic_string<CharT>>::format(ss.str(), ctx);
  }

  int precision = 17;
};

} // namespace std
#endif // CGAL_CAN_USE_CXX20_FORMAT

#include <CGAL/enable_warnings.h>

#endif // CGAL_IO_H
