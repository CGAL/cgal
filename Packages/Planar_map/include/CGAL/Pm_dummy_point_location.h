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
// Author(s)     : Efi Fogel <efif@post.math.tau.ac.il>

#ifndef CGAL_PM_DUMMY_POINT_LOCATION_H
#define CGAL_PM_DUMMY_POINT_LOCATION_H

#include <CGAL/Pm_point_location_base.h>
#include <CGAL/Planar_map_2/Pm_traits_wrap_2.h>

CGAL_BEGIN_NAMESPACE

template <class Planar_map_>
class Pm_dummy_point_location : public Pm_point_location_base<Planar_map_> {

public:
  typedef Planar_map_                                   Planar_map;
  typedef typename Planar_map::Traits                   Traits;
  typedef Pm_point_location_base<Planar_map>            Base;
  typedef typename Planar_map::Traits_wrap              Traits_wrap;
  typedef typename Planar_map::Locate_type              Locate_type;
  typedef typename Planar_map::Halfedge_handle          Halfedge_handle;
  typedef typename Planar_map::Halfedge_const_handle    Halfedge_const_handle;
  typedef typename Planar_map::Halfedge_iterator        Halfedge_iterator;
  typedef typename Planar_map::Halfedge                 Halfedge;
  typedef typename Planar_map::Vertex_handle            Vertex_handle;
  typedef typename Traits::Point_2                      Point_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef Pm_bounding_box_base<Planar_map>              Bounding_box;
  typedef typename Base::Halfedge_handle_iterator Halfedge_handle_iterator;
  typedef typename Base::Token                          Token;
	
public:	
  Pm_dummy_point_location() :
    Pm_point_location_base<Planar_map>(), m_traits(0) {}
  Pm_dummy_point_location(Planar_map *, Traits_wrap * traits) : 
    Pm_point_location_base<Planar_map>(), m_traits(traits) {}
	
  void init(Planar_map &, const Traits & traits)
  { m_traits = (const Traits_wrap*)(&traits); }
  void insert(Halfedge_handle, const X_monotone_curve_2 &) {}
	
  inline const Traits * get_traits() const {return m_traits;}

  Halfedge_const_handle locate(const Point_2 &, Locate_type &) const
  {
    CGAL_assertion_msg(false, "Dummy point location - locate not allowed");
    Halfedge_const_handle h; return h;
  }
	
  Halfedge_handle locate(const Point_2 &, Locate_type &)
  {
    CGAL_assertion_msg(false, "Dummy point location - locate not allowed");
    Halfedge_handle h; return h;
  }

  Halfedge_const_handle vertical_ray_shoot(const Point_2 &, Locate_type &, bool)
    const
  {
    CGAL_assertion_msg(false, "Dummy point location - locate not allowed");
    Halfedge_const_handle h; return h;
  }
	
  Halfedge_handle vertical_ray_shoot(const Point_2 &, Locate_type &, bool)
  {
    CGAL_assertion_msg(false, "Dummy point location - locate not allowed");
    Halfedge_handle h; return h;
  }

  void split_edge(const X_monotone_curve_2 &,
		  Halfedge_handle, Halfedge_handle,
		  const X_monotone_curve_2 &, const X_monotone_curve_2 &) {}

  void merge_edge(const X_monotone_curve_2 &, const X_monotone_curve_2 &,
                  Halfedge_handle,
                  const X_monotone_curve_2 &) {}

  void remove_edge(Halfedge_handle) {}

  void remove_edge(const Halfedge_handle_iterator &,
                   const Halfedge_handle_iterator &) {}

  void clear() {}

  void update(const Halfedge_handle_iterator &,
	      const Halfedge_handle_iterator &,
	      const Token & token) {}

#ifdef CGAL_PM_DEBUG
  void debug(){}
#endif
protected:
  const Traits_wrap * m_traits;
};

CGAL_END_NAMESPACE

#endif
