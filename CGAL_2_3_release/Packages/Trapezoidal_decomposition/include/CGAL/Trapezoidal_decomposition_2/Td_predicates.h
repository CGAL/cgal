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
// release		 : 
// release_date  : 1999, October 13
//
// file 		 : include/CGAL/Trapezoidal_decomposition_2/Td_predicates.h
// package		 : Trapezoidal decomposition 2
// source		 : 
// revision 	 : 
// revision_date : 
// author(s)	 : Oren Nechushtan <theoren@math.tau.ac.il>
//		   Iddo Hanniel <hanniel@math.tau.ac.il>
//
// maintainer(s) : Oren Nechushtan <theoren@math.tau.ac.il>
//
//
// coordinator	 : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter		 : 
// ======================================================================
#ifndef CGAL_TD_PREDICATES_H
#define CGAL_TD_PREDICATES_H

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif

#include <functional>

CGAL_BEGIN_NAMESPACE

template < class Td_traits> class Trapezoidal_decomposition_2;

template <class X_trapezoid>
struct Td_active_trapezoid : public std::unary_function<X_trapezoid,bool>
{
	bool operator()(const X_trapezoid& tr) const
	{
		return tr.is_active();
	}
};
template <class X_trapezoid,class Traits>
struct Td_active_non_degenerate_trapezoid : 
public std::unary_function<X_trapezoid,bool>
{
	Td_active_non_degenerate_trapezoid(Traits& t) : traits(t) {}
	bool operator()(const X_trapezoid& tr) const
	{
		return tr.is_active() && !traits.is_degenerate(tr);
	}
	protected:
		const Traits& traits;
};
template <class X_trapezoid,class Traits>
struct Td_active_right_degenerate_curve_trapezoid:
public std::unary_function<X_trapezoid,bool>
{
	Td_active_right_degenerate_curve_trapezoid(Traits& t) : traits(t) {}
	bool operator()(const X_trapezoid& tr) const
	{
		return tr.is_active() && traits.is_degenerate_curve(tr) && 
			!tr.right_bottom_neighbour();
	}
	protected:
		const Traits& traits;
};

CGAL_END_NAMESPACE

#endif //CGAL_TD_PREDICATES_H












