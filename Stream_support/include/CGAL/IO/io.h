// Copyright (c) 1997  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Andreas Fabri


#ifndef CGAL_IO_H
#define CGAL_IO_H


#include <cstdio>
#include <cctype>
#include <string>
#include <iostream>
#include <CGAL/tags.h>
#include <CGAL/IO/io_tags.h>
#include <CGAL/IO/Color.h>


namespace CGAL {

class IO {
public:
    CGAL_EXPORT static int mode;
    enum Mode {ASCII = 0, PRETTY, BINARY};
};

template <class T, class F = ::CGAL::Null_tag >
class Output_rep {
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
class Input_rep {
    T& t;
public:
    //! initialize with a reference to \a t.
    Input_rep( T& tt) : t(tt) {}
    //! perform the input, calls \c operator\>\> by default.
    std::istream& operator()( std::istream& in) const { return (in >> t); }
};


#if CGAL_FORCE_IFORMAT_DOUBLE || \
  ( ( _MSC_VER > 1600 ) && (! defined( CGAL_NO_IFORMAT_DOUBLE )) )
template <>
class Input_rep<double> {
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

    if(sscanf(buffer.c_str(), "%lf", &t) != 1) {
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
    switch(out.iword(IO::mode)) {
    case IO::ASCII :
        return out << static_cast<int>(col.red())   << ' '
		   << static_cast<int>(col.green()) << ' '
		   << static_cast<int>(col.blue());
    case IO::BINARY :
        write(out, static_cast<int>(col.red()));
        write(out, static_cast<int>(col.green()));
        write(out, static_cast<int>(col.blue()));
        return out;
    default:
        return out << "Color(" << static_cast<int>(col.red()) << ", "
		   << static_cast<int>(col.green()) << ", "
                   << static_cast<int>(col.blue()) << ')';
    }
}

inline
std::istream &operator>>(std::istream &is, Color& col)
{
    int r = 0, g = 0, b = 0;
    switch(is.iword(IO::mode)) {
    case IO::ASCII :
        is >> r >> g >> b;
        break;
    case IO::BINARY :
        read(is, r);
        read(is, g);
        read(is, b);
        break;
    default:
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
    }
    col = Color((unsigned char)r,(unsigned char)g,(unsigned char)b);
    return is;
}

CGAL_EXPORT
const char* mode_name( IO::Mode m );

// From polynomial.h TODO: Where to put this?
CGAL_EXPORT
void swallow(std::istream &is, char d);

CGAL_EXPORT
void swallow(std::istream &is, const std::string& s );



} //namespace CGAL

#endif // CGAL_IO_H
