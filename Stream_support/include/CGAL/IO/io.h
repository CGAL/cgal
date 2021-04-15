// Copyright (c) 1997
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

#include <cstdio>
#include <cctype>
#include <string>
#include <locale>
#include <iostream>
#include <CGAL/tags.h>
#include <CGAL/IO/io_tags.h>
#include <CGAL/IO/Color.h>
#include <CGAL/assertions.h>
#include <CGAL/Fraction_traits.h>


namespace CGAL {



namespace IO {

class Static {
public:

  static int get_mode()
  {
    static const int mode = std::ios::xalloc();
    return mode;
  }

};

  enum Mode {ASCII = 0, PRETTY, BINARY};

}



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

typedef IO_rep_is_specialized_aux<void>     IO_rep_is_specialized;
typedef IO_rep_is_not_specialized_aux<void> IO_rep_is_not_specialized;

template <class T, class F = ::CGAL::Null_tag >
class Output_rep : public IO_rep_is_not_specialized {
    const T& t;
public:
    //! initialize with a const reference to \a t.
    Output_rep( const T& tt) : t(tt) {}
    //! perform the output, calls \c operator\<\< by default.
    std::ostream& operator()( std::ostream& out) const { return (out << t); }
};

/*! \relates Output_rep
    \brief stream output of the \c Output_rep calls its \c operator().
*/
template <class T, class F>
std::ostream&
operator<<( std::ostream& out, Output_rep<T,F> rep) {
    return rep( out);
}

//! generic IO output format manipulator.
template <class T>
Output_rep<T>
oformat( const T& t) { return Output_rep<T>(t); }

//! generic IO output format manipulator with formatting tag.
template <class T, class F>
Output_rep<T,F>
oformat( const T& t, F) { return Output_rep<T,F>(t); }



/*!\brief
 * input functor class created by the generic IO input manipulator.
 *
 * It holds a reference to the input object. Default implementation
 * calls the stream input operator. Specializations can be written
 * for external types not supporting our stream IO format.
 */
template <class T>
class Input_rep : public IO_rep_is_not_specialized {
    T& t;
public:
    //! initialize with a reference to \a t.
    Input_rep( T& tt) : t(tt) {}
    //! perform the input, calls \c operator\>\> by default.
    std::istream& operator()( std::istream& in) const { return (in >> t); }
};

#if CGAL_FORCE_IFORMAT_DOUBLE || \
  ( ( _MSC_VER > 1600 ) && ( _MSC_VER < 1910 ) && (! defined( CGAL_NO_IFORMAT_DOUBLE )) )
template <>
class Input_rep<double> : public IO_rep_is_specialized {
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
    do {
      const int_type i = is.get();
      if(i == traits_type::eof()) {
        return is;
      }
      c = static_cast<char_type>(i);
    }while (std::isspace(c));
    if(c == '-'){
      buffer += '-';
    } else if(c != '+'){
      is.unget();
    }
    do {
      const int_type i = is.get();
      if(i == traits_type::eof()) {
        is.clear(is.rdstate() & ~std::ios_base::failbit);
        break;
      }
      c = static_cast<char_type>(i);
      if(std::isdigit(c) || (c =='.') || (c =='E') || (c =='e') || (c =='+') || (c =='-')){
        buffer += c;
      }else{
        is.unget();
        break;
      }
    }while(true);
    if(sscanf_s(buffer.c_str(), "%lf", &t) != 1) {
      // if a 'buffer' does not contain a double, set the fail bit.
      is.setstate(std::ios_base::failbit);
    }
    return is;
  }
};

template <>
class Input_rep<float> {
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
    do {
      const int_type i = is.get();
      if(i == traits_type::eof()) {
        return is;
      }
      c = static_cast<char_type>(i);
    }while (std::isspace(c));
    if(c == '-'){
      buffer += '-';
    } else if(c != '+'){
      is.unget();
    }
    do {
      const int_type i = is.get();
      if(i == traits_type::eof()) {
        is.clear(is.rdstate() & ~std::ios_base::failbit);
        break;
      }
      c = static_cast<char_type>(i);
      if(std::isdigit(c) || (c =='.') || (c =='E') || (c =='e') || (c =='+') || (c =='-')){
        buffer += c;
      }else{
        is.unget();
        break;
      }
    }while(true);
    if(sscanf_s(buffer.c_str(), "%f", &t) != 1) {
      // if a 'buffer' does not contain a double, set the fail bit.
      is.setstate(std::ios_base::failbit);
    }
    return is;
  }
};
#endif

/*! \relates Input_rep
    \brief stream input to the \c Input_rep calls its \c operator().
*/
template <class T>
std::istream&
operator>>( std::istream& in, Input_rep<T> rep) {
    return rep( in);
}

//! generic IO input format manipulator.
template <class T>
Input_rep<T>
iformat( T& t) { return Input_rep<T>(t); }


template <class T, class F = Null_tag >
class Benchmark_rep {
    const T& t;
public:
    //! initialize with a const reference to \a t.
    Benchmark_rep( const T& tt) : t(tt) {}
    //! perform the output, calls \c operator\<\< by default.
    std::ostream& operator()( std::ostream& out) const {
        return out << t;
    }

    // static function to get the benchmark name
    static std::string get_benchmark_name() {
        return "";
    }
};

template <class T, class F>
std::ostream& operator<<( std::ostream& out, Benchmark_rep<T,F> rep) {
    return rep( out);
}

template <class T>
Benchmark_rep<T> bmformat( const T& t) { return Benchmark_rep<T>(t); }

template <class T, class F>
Benchmark_rep<T,F> bmformat( const T& t, F) { return Benchmark_rep<T,F>(t); }


CGAL_EXPORT
IO::Mode
get_mode(std::ios& i);

CGAL_EXPORT
IO::Mode
set_ascii_mode(std::ios& i);

CGAL_EXPORT
IO::Mode
set_binary_mode(std::ios& i);

CGAL_EXPORT
IO::Mode
set_pretty_mode(std::ios& i);

CGAL_EXPORT
IO::Mode
set_mode(std::ios& i, IO::Mode m);

CGAL_EXPORT
bool
is_pretty(std::ios& i);

CGAL_EXPORT
bool
is_ascii(std::ios& i);

CGAL_EXPORT
bool
is_binary(std::ios& i);


template < class T >
inline
void
write(std::ostream& os, const T& t, const io_Read_write&)
{
    os.write(reinterpret_cast<const char*>(&t), sizeof(t));
}


template < class T >
inline
void
write(std::ostream& os, const T& t, const io_Operator&)
{
    os << oformat(t);
}


template < class T >
inline
void
write(std::ostream& os, const T& t, const io_Extract_insert&)
{
    insert(os, t);
}


template < class T >
inline
void
write(std::ostream& os, const T& t)
{
    write(os, t, typename Io_traits<T>::Io_tag());
}


template < class T >
inline
void
read(std::istream& is, T& t, const io_Read_write&)
{
    is.read(reinterpret_cast<char*>(&t), sizeof(t));
}


template < class T >
inline
void
read(std::istream& is, T& t, const io_Operator&)
{
    is >> iformat(t);
}


template < class T >
inline
void
read(std::istream& is, T& t, const io_Extract_insert&)
{
    extract(is, t);
}


template < class T >
inline
void
read(std::istream& is, T& t)
{
    read(is, t, typename Io_traits<T>::Io_tag());
}


inline
std::ostream& operator<<( std::ostream& out, const Color& col)
{
    switch(get_mode(out)) {
    case IO::ASCII :
        return out << static_cast<int>(col.red())   << ' '
                  << static_cast<int>(col.green()) << ' '
                  << static_cast<int>(col.blue()) << ' '
                  << static_cast<int>(col.alpha());
    case IO::BINARY :
        out.write(reinterpret_cast<const char*>(col.to_rgba().data()), 4);
        return out;
    default:
        return out << "Color(" << static_cast<int>(col.red()) << ", "
                  << static_cast<int>(col.green()) << ", "
                  << static_cast<int>(col.blue()) << ", "
                  << static_cast<int>(col.alpha()) << ")";
    }
}

inline
std::istream &operator>>(std::istream &is, Color& col)
{
    unsigned char r = 0, g = 0, b = 0, a = 0;
    int ir = 0, ig = 0, ib = 0, ia = 0;
    switch(get_mode(is)) {
    case IO::ASCII :
        is >> ir >> ig >> ib >> ia;
        r = (unsigned char)ir;
        g = (unsigned char)ig;
        b = (unsigned char)ib;
        a = (unsigned char)ia;
        break;
    case IO::BINARY :
        read(is, r);
        read(is, g);
        read(is, b);
        read(is, a);
        break;
    default:
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
    }
    col = Color(r,g,b,a);
    return is;
}

CGAL_EXPORT
const char* mode_name( IO::Mode m );

// From polynomial.h TODO: Where to put this?
CGAL_EXPORT
void swallow(std::istream &is, char d);

CGAL_EXPORT
void swallow(std::istream &is, const std::string& s );


  namespace internal {
inline
void eat_white_space(std::istream &is)
{
  std::istream::int_type c;
  do {
    c= is.peek();
    if (c== std::istream::traits_type::eof())
      return;
    else {
      std::istream::char_type cc= static_cast<std::istream::char_type>(c);
      if ( std::isspace(cc, std::locale::classic()) ) {
        is.get();
        // since peek succeeded, this should too
        CGAL_assertion(!is.fail());
      } else {
        return;
      }
    }
  } while (true);
}


  inline
  bool is_space (const std::istream& /*is*/, std::istream::int_type c)
  {
    return (c == std::istream::traits_type::eof()) ||
      std::isspace(static_cast<std::istream::char_type>(c),
                   std::locale::classic() );
  }

  inline
  bool is_eof (const std::istream& /*is*/, std::istream::int_type c)
  {
    return c == std::istream::traits_type::eof();
  }

  inline
  bool is_digit (const std::istream& /*is*/, std::istream::int_type c)
  {
    CGAL_assertion(c != std::istream::traits_type::eof());
    return std::isdigit(static_cast<std::istream::char_type>(c),
                        std::locale::classic() );
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

  Int n(0);             // unsigned number before '/' or '.'
  Int d(1);             // number after '/', or denominator (fp-case)
  bool negative = false; // do we have a leading '-'?
  bool digits = false;   // for fp-case: are there any digits at all?

  c = internal::peek(is);
  if (c != '.') {
    // is there a sign?
    if (c == '-' || c == '+') {
      is.get();
      negative = (c == '-');
      internal::eat_white_space(is);
      c=internal::peek(is);
    }
    // read n (could be empty)
    while (!internal::is_eof(is, c) && internal::is_digit(is, c)) {
      digits = true;
      n = n*10 + (c-zero);
      is.get();
      c = internal::peek(is);
    }
    // are we done?
    if (internal::is_eof(is, c) || internal::is_space(is, c)) {
      is.flags(old_flags);
      if (digits && !is.fail())
        z = negative? compose(-n,1): compose(n,1);
      return;
    }
  } else
    n = 0;

  // now we have read n, we are not done, and c is the next character
  // in the stream
  if (c == '/' || c == '.') {
    is.get();
    if (c == '/') {
      // rational case
      is >> d;
      is.flags(old_flags);
      if (!is.fail())
        z = negative? compose(-n,d): compose(n,d);
      return;
    }

    // floating point case; read number after '.' (may be empty)
    while (true) {
      c = internal::peek(is);
      if (internal::is_eof(is, c) || !internal::is_digit(is, c))
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
  if (c == 'e' || c == 'E') {
    is.get();
    is >> e;
  }

  // now construct the Gmpq
  if (!digits) {
    // illegal floating-point number
    is.setstate(std::ios_base::failbit);
    is.flags(old_flags);
    return;
  }

  // handle e
  if (e > 0)
    while (e--) n *= 10;
  else
    while (e++) d *= 10;
  is.flags(old_flags);
  if (!is.fail())
    z = (negative ? compose(-n,d) : compose(n,d));

}


  } // namespace internal

} //namespace CGAL


#ifdef CGAL_HEADER_ONLY
#include <CGAL/IO/io_impl.h>
#endif // CGAL_HEADER_ONLY

#include <CGAL/enable_warnings.h>

#endif // CGAL_IO_H
