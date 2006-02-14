// Copyright (c) 2005  Stanford University (USA).
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
// $URL$
// $Id$ $Date$
// 
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KINETIC_IO_KERNEL_CERTIFICATE_H
#define CGAL_KINETIC_IO_KERNEL_CERTIFICATE_H
#include <CGAL/Kinetic/basic.h>
#include <limits>

CGAL_KINETIC_BEGIN_INTERNAL_NAMESPACE

template <class Function_kernel_t>
class Certificate {
public:
  typedef typename Function_kernel_t::Function Function;
  typedef typename Function_kernel_t::Root Time;
  Certificate(const Function &f, const Function_kernel_t& fk, const Time &b, const Time &e): rs_(fk.root_stack_object(f, b, e)) {
  }
  Certificate(){}

  const Time &failure_time() const {
    static Time inf= std::numeric_limits<Time>::infinity();
    if (rs_.empty()) return inf;
    else return rs_.top();
  }
  void pop_failure_time() {
    if (!rs_.empty()) rs_.pop();
  }
private:
  typename Function_kernel_t::Root_stack rs_;
};

template <class KK_t, class Generator>
struct Certificate_function_generator: public Generator {
  typedef Certificate<typename KK_t::Function_kernel> result_type;
  typedef typename KK_t::Function_kernel::Root Time;
  Certificate_function_generator(typename KK_t::Function_kernel fk): fk_(fk){}
  

  template <class A, class B, class C, class D, class E>
  result_type operator()(const A &a, const B &b, const C &c, const D &d, const E &e, const Time &begin, const Time &end) const {
    return result_type(Generator::operator()(a, b, c, d, e), fk_, begin, end);
  }
  template <class A, class B, class C, class D>
  result_type operator()(const A &a, const B &b, const C &c, const D &d, const Time &begin, const Time &end) const {
    return result_type(Generator::operator()(a, b, c, d), fk_, begin, end);
  }
  template <class A, class B, class C>
  result_type operator()(const A &a, const B &b, const C &c, const Time &begin, const Time &end) const {
    return result_type(Generator::operator()(a, b, c), fk_, begin, end);
  }
  template <class A, class B>
  result_type operator()(const A &a, const B &b, const Time &begin, const Time &end) const {
    return result_type(Generator::operator()(a, b), fk_, begin, end);
  }
  template <class A>
  result_type operator()(const A &a, const Time &begin, const Time &end) const {
    return result_type(Generator::operator()(a), fk_,begin,end);
  }

  


  

  typename KK_t::Function_kernel fk_;
};

CGAL_KINETIC_END_INTERNAL_NAMESPACE


#endif
