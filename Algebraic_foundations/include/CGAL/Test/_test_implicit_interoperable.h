// Copyright (c) 2006-2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id:$
//
//
// Author(s)     : Michael Hemmer    <hemmer@mpi-inf.mpg.de>
//
// =============================================================================


#ifndef CGAL_TEST_IMPLICIT_INTEROPERABLE_H
#define CGAL_TEST_IMPLICIT_INTEROPERABLE_H


#include <CGAL/basic.h>

#include <cstddef>
#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>

#include <cassert>
#include <CGAL/tags.h>
#include <CGAL/Real_embeddable_traits.h>
#include <CGAL/Algebraic_structure_traits.h>



CGAL_BEGIN_NAMESPACE

template <typename A, typename B>
void test_implicit_interoperable_for_real_embeddable (CGAL::Tag_false){};

template <typename A, typename B>
void test_implicit_interoperable_for_real_embeddable (CGAL::Tag_true){
  // two sided test for interoperability with int  
  A a;
  B b;
  a = A(-5);
  b = B(-2);
  // a < b 
  assert (!(a == b));
  assert ( (a != b));
  assert ( (a <  b));
  assert ( (a <= b));
  assert (!(a >  b));
  assert (!(a >= b));
    
  assert (!(b == a));
  assert ( (b != a));
  assert (!(b <  a));
  assert (!(b <= a));
  assert ( (b >  a));
  assert ( (b >= a)); 

  // a > b 
  a = A(5);
  b = B(2);
  assert (!(a == b));
  assert ( (a != b));
  assert (!(a <  b));
  assert (!(a <= b));
  assert ( (a >  b));
  assert ( (a >= b));
    
  assert (!(b == a));
  assert ( (b != a));
  assert ( (b <  a));
  assert ( (b <= a));
  assert (!(b >  a));
  assert (!(b >= a));

  // a == b
  a = A(3);
  b = B(3);
  assert ( (a == b));
  assert (!(a != b));
  assert (!(a <  b));
  assert ( (a <= b));
  assert (!(a >  b));
  assert ( (a >= b));
    
  assert ( (b == a));
  assert (!(b != a));
  assert (!(b <  a));
  assert ( (b <= a));
  assert (!(b >  a));
  assert ( (b >= a));
};


template <typename A, typename B>
void test_implicit_interoperable_for_algebraic_structure 
(CGAL::Null_tag){};

template <typename A, typename B>
void test_implicit_interoperable_for_algebraic_structure 
(CGAL::Integral_domain_without_division_tag){
  typedef CGAL::Coercion_traits<A,B> CT; 
  typedef typename CT::Type C; 
  A a(6); B b(2);
  assert(a + b == C(8));
  assert(a - b == C(4));
  assert(a * b == C(12));
  
  assert(b + a == C(8));
  assert(b - a == C(-4));
  assert(b * a == C(12));

  assert((C(4)+= A(3)) == C(7));
  assert((C(4)-= A(3)) == C(1));
  assert((C(4)*= A(3)) == C(12));
  
  assert((C(4)+= B(3)) == C(7));
  assert((C(4)-= B(3)) == C(1));
  assert((C(4)*= B(3)) == C(12));  
};

template <typename A, typename B>
void test_implicit_interoperable_for_algebraic_structure 
(CGAL::Field_tag){
  test_implicit_interoperable_for_algebraic_structure<A,B>
    (CGAL::Integral_domain_without_division_tag());
  
  typedef CGAL::Coercion_traits<A,B> CT; 
  typedef typename CT::Type C; 
  A a(6); B b(2);
  assert(a / b == C(3));
  assert(b / a == C(2)/C(6));
  assert((C(4)/= A(2)) == C(2));
  assert((C(4)/= B(2)) == C(2));
};

  template <typename A, typename B>
void test_implicit_interoperable(){
  typedef CGAL::Coercion_traits<A,B> CT; 
  
  //typedef typename CT::Are_implicit_interoperable Are_implicit_interoperable;
  //assert((::boost::is_same<Are_implicit_interoperable, CGAL::Tag_true>::value));  
  
  typedef typename CT::Type C; 
  
  typename CGAL::Real_embeddable_traits<C>::Is_real_embeddable is_real_embeddable;
  test_implicit_interoperable_for_real_embeddable<A,B>(is_real_embeddable);
  
  typename CGAL::Algebraic_structure_traits<C>::Algebraic_category category;
  test_implicit_interoperable_for_algebraic_structure<A,B>(category);
}

CGAL_END_NAMESPACE

#endif // CGAL_TEST_IMPLICIT_INTEROPERABLE_H
