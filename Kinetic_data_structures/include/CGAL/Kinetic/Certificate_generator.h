// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
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
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KINETIC_IO_KERNEL_CERTIFICATE_GENERATOR_H
#define CGAL_KINETIC_IO_KERNEL_CERTIFICATE_GENERATOR_H
#include <CGAL/Kinetic/basic.h>
//#include <CGAL/Kinetic/internal/Kernel/Certificate.h>
namespace CGAL { namespace Kinetic {
template <class KK_t, class Generator>
struct Certificate_generator {
  typedef typename KK_t::Certificate result_type;
  typedef typename KK_t::Function_kernel::Root Time;
  Certificate_generator(typename KK_t::Function_kernel fk): fk_(fk){}
  Certificate_generator(){}
  
  enum When {AT, AFTER};

  template <class A, class B, class C, class D, class E, class Time>
  CGAL::Sign sign_at(const A &a, const B &b, const C &c, const D &d, const E &e,
                     const Time &begin) const {
    return eval_sign_at(gen_(a,b,c,d,e), begin);
  }


  template <class A, class B, class C, class D, class E, class Time>
  CGAL::Sign sign_after(const A &a, const B &b, const C &c, const D &d, const E &e,
                        const Time &begin) const {
    //if (when==AFTER) {
      return eval_sign_after(gen_(a, b, c, d, e), begin);
      /*} elsee {
      return sign_at(gen_(a,b,c,d,e), begin);
      }*/
  }

  template <class A, class B, class C, class D>
  CGAL::Sign sign_at(const A &a, const B &b, const C &c, const D &d,
                     const Time &begin) const {
      return eval_sign_at(gen_(a, b, c, d), begin);
  }

  template <class A, class B, class C, class D>
  CGAL::Sign sign_after(const A &a, const B &b, const C &c, const D &d,
                        const Time &begin) const {
      return eval_sign_after(gen_(a, b, c, d), begin);
  }


  template <class A, class B, class C>
  CGAL::Sign sign_at(const A &a, const B &b, const C &c, 
                     const Time &begin) const {
    return eval_sign_at(gen_(a,b,c), begin);
  }

  template <class A, class B, class C>
  CGAL::Sign sign_after(const A &a, const B &b, const C &c, 
                        const Time &begin) const {
    return eval_sign_after(gen_(a,b,c), begin);
  }


  template <class A, class B>
  CGAL::Sign sign_at(const A &a, const B &b,
                     const Time &begin) const {
    return eval_sign_at(gen_(a, b), begin);
  }
  template <class A, class B>
  CGAL::Sign sign_after(const A &a, const B &b,
                        const Time &begin) const {
    return eval_sign_after(gen_(a, b), begin);
  }

  template <class A>
  CGAL::Sign sign_at(const A &a, const Time &begin) const {
    return eval_sign_at(gen_(a), begin);
  }

  template <class A>
  CGAL::Sign sign_after(const A &a, const Time &begin) const {
    return eval_sign_after(gen_(a), begin);
  }

  template <class A, class B, class C, class D, class E>
  result_type operator()(const A &a, const B &b, const C &c, const D &d, const E &e,
                         const Time &begin, const Time &end) const {
    return result_type(gen_(a, b, c, d, e), fk_, begin, end);
  }

  template <class A, class B, class C, class D>
  result_type operator()(const A &a, const B &b, const C &c, const D &d, 
                         const Time &begin, const Time &end) const {
    return result_type(gen_(a, b, c, d), fk_, begin, end);
  }

  template <class A, class B, class C>
  result_type operator()(const A &a, const B &b, const C &c, const Time &begin, 
                         const Time &end) const {
    return result_type(gen_(a, b, c), fk_, begin, end);
  }

  template <class A, class B>
  result_type operator()(const A &a, const B &b, const Time &begin, const Time &end) const {
    return result_type(gen_(a, b), fk_, begin, end);
  }

  template <class A>
  result_type operator()(const A &a, const Time &begin, const Time &end) const {
    return result_type(gen_(a), fk_,begin,end);
  }


protected:
  CGAL::Sign eval_sign_after(const typename KK_t::Function_kernel::Function &f,
			const Time &t) const {
    return fk_.sign_after_object()(f,t);
  }
  CGAL::Sign eval_sign_at(const typename KK_t::Function_kernel::Function &f,
			const Time &t) const {
    return fk_.sign_at_object()(f,t);
  }
  
  Generator gen_;
  typename KK_t::Function_kernel fk_;
};

} } //namespace CGAL::Kinetic

#endif
