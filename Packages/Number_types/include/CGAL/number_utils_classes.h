// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/number_utils_classes.h
// package       : Number_types
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
//           
//                 slight modification of Michael's 1.1 version
// ======================================================================

// to be included by number_utils.h

#ifndef CGAL_NUMBER_UTILS_CLASSES_H
#define CGAL_NUMBER_UTILS_CLASSES_H 1

#ifndef CGAL_CONFIG_H
#include <CGAL/config.h>
#endif // CGAL_CONFIG_H

#include <algorithm>
#include <CGAL/functional_base.h>

CGAL_BEGIN_NAMESPACE

template < class NT >
struct Is_zero :public CGAL_STD::unary_function< NT, bool > {
  typedef Arity_tag< 1 > Arity;
  bool operator()( const NT& x) const
  { return is_zero( x); }
};

template < class NT >
struct Is_one :public CGAL_STD::unary_function< NT, bool > {
  typedef Arity_tag< 1 > Arity;
  bool operator()( const NT& x) const
  { return is_one( x); }
};

template < class NT >
struct Is_negative :public CGAL_STD::unary_function< NT, bool > {
  typedef Arity_tag< 1 > Arity;
  bool operator()( const NT& x) const
  { return is_negative( x); }
};

template < class NT >
struct Is_positive :public CGAL_STD::unary_function< NT, bool > {
  typedef Arity_tag< 1 > Arity;
  bool operator()( const NT& x) const
  { return is_positive( x); }
};

// Sign would result in a name clash with enum.h
template < class NT >
struct Sgn :public CGAL_STD::unary_function< NT, int > {
  typedef Arity_tag< 1 > Arity;
  Sign operator()( const NT& x) const
  { return sign( x); }
};

template < class NT >
struct Lexicographical_sign
  :public CGAL_STD::binary_function< NT, NT, int > {
  typedef Arity_tag< 2 > Arity;

  Sign operator()( const NT& x, const NT& y) const
  { return lexicographical_sign( x, y); }
};

template < class NT >
struct Abs :public CGAL_STD::unary_function< NT, NT > {
  typedef Arity_tag< 1 > Arity;
  NT operator()( const NT& x) const
  { return abs( x); }
};

template < class NT, class Compare = std::less< NT > >
struct Min :public CGAL_STD::binary_function< NT, NT, NT > {
  typedef Arity_tag< 2 > Arity;
  Min() {}
  Min(const Compare& c_) : c(c_) {}
  NT operator()( const NT& x, const NT& y) const
  { return std::min( x, y, c); }
protected:
  Compare c;
};

template < class NT, class Compare = std::less< NT > >
struct Max :public CGAL_STD::binary_function< NT, NT, NT > {
  typedef Arity_tag< 2 > Arity;
  Max() {}
  Max(const Compare& c_) : c(c_) {}
  NT operator()( const NT& x, const NT& y) const
  { return std::max( x, y, c); }
protected:
  Compare c;
};

template < class NT >
struct Compare
  :public CGAL_STD::binary_function< NT, NT, Comparison_result > {
  typedef Arity_tag< 2 > Arity;

  Comparison_result
  operator()( const NT& x, const NT& y) const
  { return compare( x, y); }
};

template < class NT >
struct Square :public CGAL_STD::unary_function< NT, NT > {
  typedef Arity_tag< 2 > Arity;

  NT
  operator()( const NT& x) const
  { return square( x ); }
};

CGAL_END_NAMESPACE

#endif // CGAL_NUMBER_UTILS_CLASSES_H
