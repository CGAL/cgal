// Copyright (c) 2001  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_KERNEL_CHECKER_H
#define CGAL_KERNEL_CHECKER_H

// This file contains the definition of a kernel traits checker.
//
// TODO:
// - At the moment, only predicates are checked.  To handle constructions as
//   well, the best approach is probably to have objects be pairs, and do
//   everything in parallel.
//   So the template parameter will be a comparator, not a converter.

#include <CGAL/basic.h>
#include <utility>

CGAL_BEGIN_NAMESPACE

// Class used by Kernel_checker.
template <class O1, class O2, class Conv>
class Predicate_checker
{
    O1 o1;
    O2 o2;
    Conv c;

public:

    Predicate_checker(const O1 &oo1 = O1(), const O2 &oo2 = O2())
	: o1(oo1), o2(oo2) {}

    typedef typename O1::result_type result_type;
    typedef typename O1::Arity       Arity;

    template <class A1>
    result_type
    operator()(const A1 &a1) const
    {
	typename O1::result_type res1 = o1(a1);
	typename O2::result_type res2 = o2(c(a1));
	if (res1 != res2)
	{
	    std::cerr << "Kernel_checker error : " << res1 << " != " << res2
		      << " for the inputs : " << std::endl;
	    std::cerr << a1 << std::endl;
#ifdef __GNUG__
	    std::cerr << __PRETTY_FUNCTION__ << std::endl;
#endif
	    CGAL_kernel_assertion(false);
	}
	return res1;
    }

    template <class A1, class A2>
    result_type
    operator()(const A1 &a1, const A2 &a2) const
    {
	typename O1::result_type res1 = o1(a1, a2);
	typename O2::result_type res2 = o2(c(a1), c(a2));
	if (res1 != res2)
	{
	    std::cerr << "Kernel_checker error : " << res1 << " != " << res2
		      << " for the inputs : " << std::endl;
	    std::cerr << a1 << ", " << a2 << std::endl;
#ifdef __GNUG__
	    std::cerr << __PRETTY_FUNCTION__ << std::endl;
#endif
	    CGAL_kernel_assertion(false);
	}
	return res1;
    }

    template <class A1, class A2, class A3>
    result_type
    operator()(const A1 &a1, const A2 &a2, const A3 &a3) const
    {
	typename O1::result_type res1 = o1(a1, a2, a3);
	typename O2::result_type res2 = o2(c(a1), c(a2), c(a3));
	if (res1 != res2)
	{
	    std::cerr << "Kernel_checker error : " << res1 << " != " << res2
		      << " for the inputs : " << std::endl;
	    std::cerr << a1 << ", " << a2 << ", " << a3 << std::endl;
#ifdef __GNUG__
	    std::cerr << __PRETTY_FUNCTION__ << std::endl;
#endif
	    CGAL_kernel_assertion(false);
	}
	return res1;
    }

    template <class A1, class A2, class A3, class A4>
    result_type
    operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4) const
    {
	typename O1::result_type res1 = o1(a1, a2, a3, a4);
	typename O2::result_type res2 = o2(c(a1), c(a2), c(a3), c(a4));
	if (res1 != res2)
	{
	    std::cerr << "Kernel_checker error : " << res1 << " != " << res2
		      << " for the inputs : " << std::endl;
	    std::cerr << a1 << ", " << a2 << ", " << a3 << ", " << a4
		      << std::endl;
#ifdef __GNUG__
	    std::cerr << __PRETTY_FUNCTION__ << std::endl;
#endif
	    CGAL_kernel_assertion(false);
	}
	return res1;
    }

    template <class A1, class A2, class A3, class A4, class A5>
    result_type
    operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	       const A5 &a5) const
    {
	typename O1::result_type res1 = o1(a1, a2, a3, a4, a5);
	typename O2::result_type res2 = o2(c(a1), c(a2), c(a3), c(a4), c(a5));
	if (res1 != res2)
	{
	    std::cerr << "Kernel_checker error : " << res1 << " != " << res2
		      << " for the inputs : " << std::endl;
	    std::cerr << a1 << ", " << a2 << ", " << a3 << ", " << a4
		      << ", " << a5 << std::endl;
#ifdef __GNUG__
	    std::cerr << __PRETTY_FUNCTION__ << std::endl;
#endif
	    CGAL_kernel_assertion(false);
	}
	return res1;
    }

    // Same thing with more arguments...
};

// For now, we inherit all geometric objects and constructions from K1, and
// just overload the predicates.
template <class K1, class K2, class Conv>
class Kernel_checker
  : public K1
{
    typedef K1     Kernel1;
    typedef K2     Kernel2;

    Kernel2 k2;

    typedef Conv   c;

    // typedef std::pair<K1::Point_2, K2::Point_2>  Point_2;
    // ...  Same thing for all objects.

#define CGAL_check_pred(X, Y) \
    typedef Predicate_checker<typename K1::X, typename K2::X, Conv> X; \
    X Y() const { return X(K1::Y(), k2.Y()); }

#define CGAL_Kernel_pred(Y,Z) CGAL_check_pred(Y, Z)
#define CGAL_Kernel_cons(Y,Z)

public:

#include <CGAL/Kernel/interface_macros.h>
};

CGAL_END_NAMESPACE

#endif // CGAL_KERNEL_CHECKER_H
