// Copyright (c) 1999
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
// SPDX-License-Identifier: LGPL-3.0-or-later
//
//
// Author(s)     : Stefan Schirra


#ifndef CGAL__TEST_CLS_OBJECT_H
#define CGAL__TEST_CLS_OBJECT_H

#include <CGAL/Testsuite/use.h>
#include <cassert>

using CGAL::internal::use;

// Test that we can derive from Object.

class Object_handle
  : public CGAL::Object
{
  typedef CGAL::Object Base;
public:
  Object_handle() : Base() {}
  Object_handle(const CGAL::Object& o) : Base(o) {}
  Object_handle(const Object_handle& h) : Base(h) {}

  Object_handle&
  operator=(const Object_handle& v)=default;
};

Object_handle return_obj()
{
  Object_handle o;
  return o;
}


template <class R>
bool
_test_cls_object(const R&)
{
  typedef typename  R::RT   RT;
  std::cout << "testing class object" ;

  CGAL::Object o1;
  assert( o1.is_empty() );
  CGAL::Object o2;
  assert( o2.is_empty() );
  CGAL::Object o3;
  assert( o3.is_empty() );
  CGAL::Point_2<R>   p21( RT(1), RT(1) );
  CGAL::Point_2<R>   p22;
  CGAL::Point_2<R>   p23;
  CGAL::Point_2<R>   p24( RT(4), RT(4) );
  CGAL::Point_3<R>   p31;
  CGAL::Point_3<R>   p32;
  CGAL::Line_2<R>    l21( p21, p24);
  CGAL::Line_2<R>    l22;
  CGAL::Line_3<R>    l31;
  CGAL::Line_3<R>    l32;

  std::cout << '.';

  o1 = o2;
  o2 = CGAL::make_object(p21);
  assert( o1.is_empty() );
  o1 = CGAL::make_object(p31);
  assert(   CGAL::assign( p22, o2 ) );
  assert( p22 == p21 );
  assert( ! CGAL::assign( p32, o2 ) );
  assert( ! CGAL::assign( l22, o2 ) );
  o3 = o2;
  assert(   CGAL::assign( p23, o3 ) );
  assert( p23 == p21 );
  assert( ! CGAL::assign( p32, o2 ) );
  assert( ! CGAL::assign( l32, o2 ) );
  assert(   CGAL::assign( p32, o1 ) );
  assert( ! CGAL::assign( p21, o1 ) );
  assert( p21 == p23 );

  std::cout << '.';

  o2 = CGAL::make_object(l21);
  assert(   CGAL::assign( l22, o2 ) );
  assert( l22 == l21 );
  assert( ! CGAL::assign( l32, o2 ) );
  assert( ! CGAL::assign( p22, o2 ) );
  o3 = CGAL::make_object(l31);
  assert(   CGAL::assign( l32, o3 ) );
  assert( ! CGAL::assign( p32, o3 ) );

  // Test that deriving from Object works.
  Object_handle o4;
  o4 = return_obj();

  // Test .type().
  CGAL::Object o5;
  std::cout << o5.type().name() << std::endl;
  assert( o5.type() == typeid(void) );

  CGAL::Object o6 = CGAL::make_object(2);
  std::cout << o6.type().name() << std::endl;
  assert( o6.type() == typeid(int) );

  // Test object_cast<>().
  const int *i = CGAL::object_cast<int>(&o6);
  assert( i != nullptr );
  int j = CGAL::object_cast<int>(o6);
  use(j);

  const double *d = CGAL::object_cast<double>(&o6);
  assert( d == nullptr );
  try {
    // This case must throw.
    double k = CGAL::object_cast<double>(o6);
    use(k);
    assert(false);
  }
  catch (CGAL::Bad_object_cast&) {}


  std::cout << "done" << std::endl;
  return true;
}

#endif // CGAL__TEST_CLS_OBJECT_H
