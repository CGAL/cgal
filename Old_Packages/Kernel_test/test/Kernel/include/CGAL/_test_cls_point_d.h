// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// file          : _test_cls_point_d.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL__TEST_CLS_POINT_D_H
#define CGAL__TEST_CLS_POINT_D_H


template <class R>
bool
_test_cls_point_d( const R& )
{
  double coord1[6] = {0.0, 1.0, 2.0, 3.0, 4.0, 1.0};
  double coord2[6] = {0.0, 2.0, 4.0, 6.0, 8.0, 2.0};
  CGAL::Point_d<R> p (5, coord1, coord1+6);
  CGAL::Point_d<R> pp (5, coord1, coord1+5);
  CGAL::Point_d<R> q (5, coord2, coord2+6);
  CGAL::Point_d<R> s (p);   // constructors
  CGAL::Point_d<R> t;
  t = s;                    // assignment

  assert ( p == pp);
  assert ( p == q);
  assert ( p == s);         // equality test
  assert ( p == t);

  for (int i=0; i<5; ++i)
  {
    assert (q.cartesian(i) == (double)i);   // cartesian method
    assert (q[i] == (double)i);             // operator[]
    assert (p.homogeneous(i) == (double)i);
  }
  assert (s.dimension() == 5);              // dimension

  // I/O

  const int size = 100;                 // suffices to store a point
  char buffer[size];
  std::ostrstream ost (buffer, size);
  CGAL::set_ascii_mode (ost);
  ost << p;                    // write p

  std::istrstream ist (buffer, size);
  CGAL::Point_d<R> u;
  CGAL::set_ascii_mode (ist);
  ist >> u;                    // read p back in as u

  assert (p == u);

  return true;
}

#endif // CGAL__TEST_CLS_POINT_D_H
