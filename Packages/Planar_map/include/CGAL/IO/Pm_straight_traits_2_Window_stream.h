// Copyright (c) 1997  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Oren Nechushtan <theoren@math.tau.ac.il>
//                 Shai Hirsch     <shaihi@post.tau.ac.il>
#ifndef CGAL_IO_PM_STRAIGHT_TRAITS_2_WINDOW_STREAM_H
#define CGAL_IO_PM_STRAIGHT_TRAITS_2_WINDOW_STREAM_H

#ifndef CGAL_PM_STRAIGHT_TRAITS_2_H
#include <CGAL/Pm_straight_traits_2.h>
#endif

#ifndef CGAL_LEDA_WINDOW_H
#include <CGAL/IO/leda_window.h>
#endif

CGAL_BEGIN_NAMESPACE

template <class R> 
leda_window& operator<<(leda_window& os,const Straight_2_<R>& cv)
{
	typedef Straight_2_<R> Straight;
	switch(cv.current_state())
	{
	case Straight::EMPTY:
	  return os;
	case Straight::POINT:
		{
			Point_2<R> p;
			cv.current(p);
			return os << p;
		}
	case Straight::SEGMENT:
		{
			Segment_2<R> seg;
			cv.current(seg);
			return os << seg;
		}
	case Straight::RAY:
		{
			Ray_2<R> ray;
			cv.current(ray);
			return os << ray;
		}
	case Straight::LINE:
		{
			Line_2<R> line;
			cv.current(line);
			return os << line;
		}
	}
	CGAL_assertion_msg(
		cv.current_state()==Straight::EMPTY||
		cv.current_state()==Straight::POINT||
		cv.current_state()==Straight::SEGMENT||
		cv.current_state()==Straight::RAY||
		cv.current_state()==Straight::LINE,
		"\nUnknown type in  leda_window& operator<<(leda_window& os,\
const Straight& cv)");
	return os;
}
template <class R> 
leda_window& operator>>(leda_window& os,const Straight_2_<R>& cv)
{
	typedef Straight_2_<R> Straight;
	switch(cv.current_state())
	{
	case Straight::EMPTY:
	  return os;
	case Straight::POINT:
		{
			Point_2<R> p;
			cv.current(p);
			return os >> p;
		}
	case Straight::SEGMENT:
		{
			Segment_2<R> seg;
			cv.current(seg);
			return os >> seg;
		}
	case Straight::RAY:
		{
			Ray_2<R> ray;
			cv.current(ray);
			return os >> ray;
		}
	case Straight::LINE:
		{
			Line_2<R> line;
			cv.current(line);
			return os >> line;
		}
	}
	CGAL_assertion_msg(
		cv.get_type()==Straight::EMPTY||
		cv.get_type()==Straight::POINT||
		cv.get_type()==Straight::SEGMENT||
		cv.get_type()==Straight::RAY||
		cv.get_type()==Straight::LINE,
		"\nUnknown type in  leda_window& operator>>(leda_window& os,"
		<< "const Straight& cv)");
	return os;
}
template <class R>
Window_stream& write(
		     Window_stream& os, 
		     const typename Pm_straight_traits_2<R>::X_curve& cv,
		     const Pm_straight_traits_2<R>& traits)
{
	typedef Pm_straight_traits_2<R> Traits; 
	typedef typename Traits::X_bounded_curve X_bounded_curve;
	return os << X_bounded_curve(
		traits.curve_source(cv),
		traits.curve_target(cv));
}

CGAL_END_NAMESPACE

#endif // CGAL_IO_PM_STRAIGHT_TRAITS_2_WINDOW_STREAM_H
// EOF
