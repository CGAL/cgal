// Copyright (c) 2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
// 
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_KERNEL_PROFILER_H
#define CGAL_KERNEL_PROFILER_H

// This file contains the definition of a kernel traits profiler.

#include <CGAL/basic.h>
#include <typeinfo>

namespace CGAL {

// Primitive wrapper which handles the profiling.
template < typename P >
struct Primitive_profiler
  : public P
{
    typedef typename P::result_type  result_type;

// #define CGAL_KERNEL_PROFILER CGAL_PROFILER(CGAL_PRETTY_FUNCTION);
#define CGAL_KERNEL_PROFILER \
        CGAL_PROFILER(typeid(static_cast<const P&>(*this)).name())

    Primitive_profiler(const P& p = P())
      : P(p) {}

    template <class A1>
    result_type
    operator()(const A1 &a1) const
    {
	CGAL_KERNEL_PROFILER;
	return P::operator()(a1);
    }

    template <class A1, class A2>
    result_type
    operator()(const A1 &a1, const A2 &a2) const
    {
	CGAL_KERNEL_PROFILER;
	return P::operator()(a1, a2);
    }

    template <class A1, class A2, class A3>
    result_type
    operator()(const A1 &a1, const A2 &a2, const A3 &a3) const
    {
	CGAL_KERNEL_PROFILER;
	return P::operator()(a1, a2, a3);
    }

    template <class A1, class A2, class A3, class A4>
    result_type
    operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4) const
    {
	CGAL_KERNEL_PROFILER;
	return P::operator()(a1, a2, a3, a4);
    }

    template <class A1, class A2, class A3, class A4, class A5>
    result_type
    operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	       const A5 &a5) const
    {
	CGAL_KERNEL_PROFILER;
	return P::operator()(a1, a2, a3, a4, a5);
    }

    template <class A1, class A2, class A3, class A4, class A5, class A6>
    result_type
    operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	       const A5 &a5, const A6 &a6) const
    {
	CGAL_KERNEL_PROFILER;
	return P::operator()(a1, a2, a3, a4, a5, a6);
    }

    template <class A1, class A2, class A3, class A4,
              class A5, class A6, class A7>
    result_type
    operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	       const A5 &a5, const A6 &a6, const A7 &a7) const
    {
	CGAL_KERNEL_PROFILER;
	return P::operator()(a1, a2, a3, a4, a5, a6, a7);
    }

    template <class A1, class A2, class A3, class A4,
              class A5, class A6, class A7, class A8>
    result_type
    operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	       const A5 &a5, const A6 &a6, const A7 &a7, const A8 &a8) const
    {
	CGAL_KERNEL_PROFILER;
	return P::operator()(a1, a2, a3, a4, a5, a6, a7, a8);
    }

    // ...
};

// We inherit all geometric objects from K, and just replace the primitives.
template < typename K >
struct Kernel_profiler
  : public K
{
#define CGAL_prof_prim(X, Y) \
    typedef Primitive_profiler<typename K::X> X; \
    X Y() const { return X(static_cast<const K&>(*this).Y()); }

#define CGAL_Kernel_pred(X, Y)  CGAL_prof_prim(X, Y)
#define CGAL_Kernel_cons(X, Y)  CGAL_prof_prim(X, Y)

#include <CGAL/Kernel/interface_macros.h>
};

} //namespace CGAL

#endif // CGAL_KERNEL_PROFILER_H
