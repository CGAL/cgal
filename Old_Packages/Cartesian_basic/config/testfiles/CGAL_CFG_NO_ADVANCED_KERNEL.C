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
// release       : $CGAL_Revision: CGAL-2.1-I-4 $
// release_date  : $CGAL_Date: 1999/07/06 $
//
// file          : config/testfiles/CGAL_CFG_NO_ADVANCED_KERNEL.C
// package       : Cartesian_basic
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : Herve.Bronnimann@sophia.inria.fr
//
// coordinator   : INRIA Sophia-Antipolis
//
// ======================================================================

// CGAL_CFG_NO_ADVANCED_KERNEL.C

// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| We evaluate the deferred instantiation scheme of Michael Hoffmann
//| versus the more conservative approach of Stefan Schirra. Michael uses
//! partial specialisation, passing of Self as template parameter of a base
//! class, and other features which are not available for some
//! compilers. This program is used to detect this problem.

class Cartesian_tag {};

// Basic kernel scheme

template < class R, class T = typename R::Tag > class Point_2;
template < class R, class T = typename R::Tag > class Line_2;

template < class R >
class Point_2 < R, Cartesian_tag >
{
public:
  typename R::FT  x, y;
};

template < class R >
class Line_2 < R, Cartesian_tag >
{
public:
  typename R::FT  a, b, c;
};

template < class R, class _FT >
class Cartesian_base
{
public:
  typedef _FT                                 FT;
  typedef Cartesian_tag                       Tag;
  // Strange error:  no type named `Tag' in `class Cartesian<double>'
  // typedef ::Point_2<R>                        Point_2;
  // typedef ::Line_2<R>                         Line_2;
  typedef ::Point_2<R,Tag>                     Point_2;
  typedef ::Line_2<R,Tag>                      Line_2;
};

template < class _FT > 
class Cartesian
: public Cartesian_base< Cartesian< _FT >, _FT >
{};

// Kernel extensibility

template < class R, class T = typename R::Tag > class Line_2p_2;

template < class R >
class Line_2p_2 < R, Cartesian_tag >
{
public:
  typename R::Point_2  p, q;
};

template < class R, class _FT >
class My_cartesian_base : public Cartesian_base<R,_FT>
{
  typedef Line_2p_2<R,Tag>                     Line_2;
};

template < class _FT > 
class My_cartesian
: public My_cartesian_base< My_cartesian< _FT >, _FT >
{};

// Test

int
main()
{
  typedef double             FT;
  typedef Cartesian<FT>      Rep1;
  typedef Point_2<Rep1>      Pt1;
  typedef Rep1::Point_2      Qt1;
  typedef My_cartesian<FT>   Rep2;
  typedef Rep2::Point_2      Pt2;
  typedef Rep2::Line_2       Ln2;

  Pt1 p;
  p.x=1.0;
  p.y=2.0;
  Qt1 q;
  q = p;    // verifies equivalence between R::Point_2 and Point_2<R>
  Ln2 l;
  // l.p = p; // Error: l.p is of type Point_2<Rep2> and p is of type Point_2<Rep1>
  l.p = Pt2();
  l.q = Pt2();
  return 0;
}
