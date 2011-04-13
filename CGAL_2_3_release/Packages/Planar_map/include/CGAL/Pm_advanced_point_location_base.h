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
// file          : include/CGAL/Pm_advanced_point_location_base.h
// package       : pm (4.08)
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Iddo Hanniel <hanniel@math.tau.ac.il>
//                 Oren Nechushtan <theoren@math.tau.ac.il>
//
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
#ifndef CGAL_PM_ADVANCED_POINT_LOCATION_BASE_H
#define CGAL_PM_ADVANCED_POINT_LOCATION_BASE_H

#ifndef CGAL_PM_BOUNDING_BOX_BASE_H
#include <CGAL/Pm_bounding_box_base.h>
#endif

CGAL_BEGIN_NAMESPACE


////////////////////////////////////////////////////////////////////
//               ABSTRACT BASE CLASS OF STRATEGY
//////////////////////////////////////////////////////////////////

template <class _Planar_map>
class Pm_advanced_point_location_base : 
	public Pm_advanced_point_location_base<_Planar_map>{
public:
	typedef _Planar_map Planar_map;
	typedef typename Planar_map::Traits Traits;
	typedef typename Planar_map::Locate_type Locate_type;
	typedef typename Planar_map::Ccb_halfedge_circulator 
		Ccb_halfedge_circulator;
	typedef typename Planar_map::Halfedge_handle Halfedge_handle;
	typedef typename Planar_map::Halfedge Halfedge;
	typedef typename Traits::X_curve X_curve;
	typedef typename Traits::Point Point;
	typedef Pm_bounding_box_base<Planar_map> Bounding_box;
	
	Pm_advanced_point_location_base():Pm_point_location_base {}
	
	virtual Halfedge_handle ray_shoot(const Point& p, 
					  Locate_type& lt, 
					  const X_curve&) = 0;
	virtual Halfedge_handle x_curve_shoot(const Point& p, 
					      Locate_type& lt, 
					      const Ray&) = 0;
	
};


CGAL_END_NAMESPACE

#endif //CGAL_PM_ADVANCED_POINT_LOCATION_BASE_H














