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
// release       : 
// release_date  : 1999, October 13
//
// file          : include/CGAL/Pm_default_point_location.C
// package       : pm (4.08)
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Oren Nechushtan <theoren@math.tau.ac.il>
//
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
#ifndef CGAL_PM_DEFAULT_POINT_LOCATION_C
#define CGAL_PM_DEFAULT_POINT_LOCATION_C

#ifndef CGAL_PM_DEFAULT_POINT_LOCATION_H
#include <CGAL/Pm_default_point_location.h>
#endif // CGAL_PM_DEFAULT_POINT_LOCATION_H

CGAL_BEGIN_NAMESPACE

//IMPLEMENTATION
//if unbounded face - returns NULL or some edge on unbounded face 
//if its a vertex returns a halfedge pointing _at_ it
	/* postconditions:
	The output Halfedge_handle represents a 
	planar map subdivision region that contains the 
	input Point in its interior or equal to it.
	The input Locate_type is equal to the type 
	of this region.
	*/
template <class Planar_map>
Pm_default_point_location<Planar_map>::Halfedge_handle
Pm_default_point_location<Planar_map>::locate(const Point& p, Locate_type& lt)
  const
{
		//there are different internal compiler errors if we
		// typedef the Locate_type
		typename TD::Locate_type td_lt; 
		
		const X_curve_plus& cv = td.locate(p,td_lt).top();
		// treat special case, where trapezoid is unbounded.
		//	for then get_parent() is not defined
		if (td_lt==TD::UNBOUNDED_TRAPEZOID)
		{
			lt=PM::UNBOUNDED_FACE;
			return halfedge_representing_unbounded_face();
		}
		Halfedge_handle h = cv.get_parent();
		lt=convert(p,td_lt,h);
		
		return h;
    }

template <class Planar_map>
Pm_default_point_location<Planar_map>::Halfedge_handle
Pm_default_point_location<Planar_map>::locate(const Point& p, Locate_type& lt){
	((Bounding_box*)get_bounding_box())->insert(p);
	Halfedge_handle h=((cPLp)this)->locate(p,lt);
	if (!((Bounding_box*)get_bounding_box())->locate(p,lt,h))
		h=((cPLp)this)->locate(p,lt);
	return h;
}

	/* postconditions:
	The output Halfedge_handle represents a planar map 
	subdivision region that contains the first Point 
	on the closed vertical ray eminating from the input 
	Point in upward or downward direction depending on 
	the input bool in its interior or equal to it.
	The input Locate_type is equal to the type 
	of this region.
	*/
template <class Planar_map>
Pm_default_point_location<Planar_map>::Halfedge_handle
Pm_default_point_location<Planar_map>::vertical_ray_shoot(
	const Point& p, Locate_type& lt, bool up) const{

		//trying to workaround internal compiler error
		typename TD::Locate_type td_lt;
		
		X_curve_plus cv = td.vertical_ray_shoot(p,td_lt,up);
		// treat special case, where trapezoid is unbounded.
		//	for then get_parent() is not defined
		if (td_lt==TD::UNBOUNDED_TRAPEZOID)
		{
			lt=PM::UNBOUNDED_FACE;
			return halfedge_representing_unbounded_face();
		}
		Halfedge_handle h=cv.get_parent();
		lt=convert(p,td_lt,h,up);
		
		return h;
    }

template <class Planar_map>
Pm_default_point_location<Planar_map>::Halfedge_handle
Pm_default_point_location<Planar_map>::vertical_ray_shoot(
	const Point& p, Locate_type& lt, bool up){
/* Make sure the source point is in the bounding box on the output */
	((Bounding_box*)get_bounding_box())->insert(p);
	Halfedge_handle h=((cPLp)this)->vertical_ray_shoot(p,lt,up);
/* Apply the bounding box on the output */
	if (!((Bounding_box*)get_bounding_box())->vertical_ray_shoot(p,lt,
								     up,h))
	{
		h=((cPLp)this)->vertical_ray_shoot(p,lt,up);
		CGAL_assertion(lt!=Planar_map::UNBOUNDED_FACE);
	}
	return h;
}

CGAL_END_NAMESPACE

#endif // CGAL_PM_DEFAULT_POINT_LOCATION_C
