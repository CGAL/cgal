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
// Author(s)     : Iddo Hanniel <hanniel@post.tau.ac.il>
//                 Oren Nechushtan <theoren@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>
#ifndef CGAL_PM_NAIVE_POINT_LOCATION_H
#define CGAL_PM_NAIVE_POINT_LOCATION_H

#include <CGAL/Pm_point_location_base.h>
#include <CGAL/Planar_map_2/Pm_traits_wrap_2.h>

CGAL_BEGIN_NAMESPACE

template <class Planar_map_>
class Pm_naive_point_location : public Pm_point_location_base<Planar_map_> {
public:
  typedef Planar_map_ Planar_map;
  typedef typename Planar_map::Traits Traits;
  typedef Pm_point_location_base<Planar_map> Base;
  typedef Pm_naive_point_location<Planar_map> Self;
  typedef typename Planar_map::Traits_wrap Traits_wrap;
  typedef typename Planar_map::Locate_type Locate_type;
  typedef typename Planar_map::Face_handle Face_handle;
  typedef typename Planar_map::Face_const_handle Face_const_handle;
  typedef typename Planar_map::Ccb_halfedge_circulator 
    Ccb_halfedge_circulator;
  typedef typename Planar_map::Ccb_halfedge_const_circulator 
    Ccb_halfedge_const_circulator;
  typedef typename Planar_map::Halfedge_handle Halfedge_handle;
  typedef typename Planar_map::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Planar_map::Halfedge_iterator Halfedge_iterator;
  typedef typename Planar_map::Halfedge_const_iterator Halfedge_const_iterator;
  typedef typename Planar_map::Halfedge Halfedge;
  typedef typename Planar_map::Vertex_handle Vertex_handle;
  typedef typename Planar_map::Vertex_const_handle Vertex_const_handle;
  typedef typename Traits::Point Point;
  typedef typename Traits::X_curve X_curve;
  typedef Pm_bounding_box_base<Planar_map> Bounding_box;
  typedef typename Base::Halfedge_handle_iterator Halfedge_handle_iterator;
  typedef typename Base::Token Token;
  
public:
  // Constructors
  Pm_naive_point_location() : 
    Pm_point_location_base<Planar_map>(),
    pm(0),
    traits(0) 
  {}

  Pm_naive_point_location(Planar_map * _pm,Traits_wrap * _traits) : 
    Pm_point_location_base<Planar_map>(), traits(_traits), pm(_pm) {}

  inline void init(Planar_map & pmp, const Traits & tr) 
  {
    CGAL_precondition_msg(pm == NULL,
    "Point location instance should be uninitialized "
    "(Do not use the same instance for more than one map).");

    pm = &pmp;
    traits = (Traits_wrap*)(&tr);
  }
  
  inline void insert(Halfedge_handle, const X_curve &) {}
  
  Halfedge_const_handle locate(const Point & p, Locate_type & lt) const;
  Halfedge_handle locate(const Point & p, Locate_type & lt);
  
  Halfedge_const_handle vertical_ray_shoot(const Point & p,
                                     Locate_type & lt, bool up) const;
  Halfedge_handle vertical_ray_shoot(const Point & p,
                                     Locate_type & lt, bool up);
  
  inline void split_edge(const X_curve &,
                         Halfedge_handle, Halfedge_handle,
                         const X_curve &, const X_curve &) {}
  
  inline void merge_edge(const X_curve &, const X_curve &,
                         Halfedge_handle, const X_curve &) {}
  
  inline void remove_edge(Halfedge_handle) {}
  inline void remove_edge(const Halfedge_handle_iterator &,
                          const Halfedge_handle_iterator &) {};
  inline void clear() {}
  inline void update(const Halfedge_handle_iterator &,
                     const Halfedge_handle_iterator &,
                     const Token & token)
  { token.rebuild_bounding_box(this); }

public:
  inline const Bounding_box * get_bounding_box() const 
  {return pm->get_bounding_box();}
  inline const Traits * get_traits() const {return traits;}
  
protected:
  Halfedge_const_handle find_lowest(Vertex_const_handle v, bool highest) const;
  
#ifdef CGAL_PM_DEBUG
  void debug(){}
#endif

protected:
  typedef const Self * const_Self_ptr;
  
protected:
  Planar_map * pm;
  const Traits_wrap * traits;
};

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/Pm_naive_point_location.C>
#endif

#endif
