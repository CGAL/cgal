// Copyright (c) 2017 Inria.
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
// Author: Marc Glisse <marc.glisse@inria.fr>

#ifndef CGAL_GMPXX_ARITHMETIC_KERNEL_H
#define CGAL_GMPXX_ARITHMETIC_KERNEL_H

#include <CGAL/Arithmetic_kernel/Arithmetic_kernel_base.h>
#include <CGAL/Get_arithmetic_kernel.h>

#include <CGAL/boost_mp.h>
#include <boost/multiprecision/cpp_int.hpp>
#ifdef CGAL_USE_GMP
#include <boost/multiprecision/gmp.hpp>
#endif

// FIXME: the could be several kernels based on Boost.Multiprecision.

namespace CGAL {
/** \ingroup CGAL_Arithmetic_kernel
 *  \brief The Boost.Multiprecision set of exact number types
 */
struct BOOST_MP_arithmetic_kernel : internal::Arithmetic_kernel_base {
  typedef boost::multiprecision::cpp_int Integer;
  typedef boost::multiprecision::cpp_rational Rational;
};

template <class B1, boost::multiprecision::expression_template_option E1>
struct Get_arithmetic_kernel<boost::multiprecision::number<B1, E1> > {
  typedef BOOST_MP_arithmetic_kernel Arithmetic_kernel;
};
template <class T1, class T2, class T3, class T4, class T5>
struct Get_arithmetic_kernel<boost::multiprecision::detail::expression<T1,T2,T3,T4,T5> > {
  typedef BOOST_MP_arithmetic_kernel Arithmetic_kernel;
};
} //namespace CGAL
#endif
