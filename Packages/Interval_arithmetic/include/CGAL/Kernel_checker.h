// ============================================================================
//
// Copyright (c) 2001 The CGAL Consortium
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
// file          : include/CGAL/Kernel_checker.h
// revision      : $Revision$
// revision_date : $Date$
// package       : ???
// author(s)     : Sylvain Pion
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_KERNEL_CHECKER_H
#define CGAL_KERNEL_CHECKER_H

// This file contains the definition of a kernel traits checker.
//
// TODO:
// - have a look at the PM_checker from Tel-Aviv, and Stefan's NT checker.


#include <CGAL/basic.h>
#include <pair>

CGAL_BEGIN_NAMESPACE

class Default_O_Comp
{
public:
    template <class C1, class C2>
    bool
    operator()(const C1 &c1, const C2 &c2)
    {
	return c1 == c2;
    }
};

// sub class used by Kernel_checker.
// predicate vs construction could be selected by a tag
// predicate/construction as additional parameter.
template <class O1, class O2, class O_Comp>
class Object_checker
{
    O1 o1;
    O2 o2;

public:

    // for predicate
    typedef O1::result_type result_type;

    // for construction
    typedef std::pair<O1::result_type, O2::result_type> result_type;

    template <class A1, class A2>
    result_type
    operator()(const A1 &a1, const A2 &a2)
    {
	O1::result_type res1 = o1(a1.first, a2.first);
	O2::result_type res2 = o2(a1.second, a2.second);
	CGAL_assertion(O_Comp()(res1, res2));
	return res1; // predicate
	return std::make_pair(res1, res2); // construction.
    }
    // Same thing with more arguments...
    
};

template <class K1, class K2, class Comp = Default_Comparator>
class Kernel_checker
{
    typedef K1     Kernel1;
    typedef K2     Kernel2;
    typedef Comp   Comparator;

    typedef std::pair<K1::Point_2, K2::Point_2>  Point_2;
    // ...  Same thing for all objects.

    typedef Object_checker<K1::Orientation_2, K2::Orientation_2, Comp>
	    Orientation_2;
};

CGAL_END_NAMESPACE

#endif // CGAL_KERNEL_CHECKER_H
