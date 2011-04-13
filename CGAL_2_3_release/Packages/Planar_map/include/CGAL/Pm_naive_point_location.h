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
// file          : include/CGAL/Pm_naive_point_location.h
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
#ifndef CGAL_PM_NAIVE_POINT_LOCATION_H
#define CGAL_PM_NAIVE_POINT_LOCATION_H

#ifndef CGAL_PM_POINT_LOCATION_BASE_H
#include <CGAL/Pm_point_location_base.h>
#endif

#ifndef CGAL_PLANAR_MAP_MISC_H
#include <CGAL/Planar_map_2/Planar_map_misc.h>
#endif


CGAL_BEGIN_NAMESPACE

////////////////////////////////////////////////////////
//NAIVE STRATEGY
////////////////////////////////////////////////////////

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
  typedef typename Planar_map::Ccb_halfedge_circulator 
    Ccb_halfedge_circulator;
  typedef typename Planar_map::Halfedge_handle Halfedge_handle;
  typedef typename Planar_map::Halfedge_iterator Halfedge_iterator;
  typedef typename Planar_map::Halfedge Halfedge;
  typedef typename Planar_map::Vertex_handle Vertex_handle;
  typedef typename Traits::Point Point;
  typedef typename Traits::X_curve X_curve;
  typedef Pm_bounding_box_base<Planar_map> Bounding_box;
  typedef typename Base::Halfedge_handle_iterator Halfedge_handle_iterator;
  typedef typename Base::Token Token;

public:	
  Pm_naive_point_location() : Pm_point_location_base<Planar_map>(),traits(0) {}
  Pm_naive_point_location(Planar_map* _pm,Traits_wrap* _traits) : 
    Pm_point_location_base<Planar_map>(),traits(_traits),pm(_pm) {}
  
  inline void init(Planar_map& pmp, Traits& tr) {
    pm = &pmp;
    traits = (Traits_wrap*)(&tr);
  }
  
  inline void insert(Halfedge_handle h
                     ,const X_curve& cv
                     ) {}
  
  Halfedge_handle locate(const Point& p, Locate_type& lt) const;
  Halfedge_handle locate(const Point& p, Locate_type& lt);
  
  Halfedge_handle vertical_ray_shoot(const Point& p, Locate_type& lt, bool up)
    const;
  Halfedge_handle vertical_ray_shoot(const Point& p, Locate_type& lt, bool up);
  
  inline void split_edge(const X_curve &cv,
                         Halfedge_handle e1,
                         Halfedge_handle e2
                         ,const X_curve& cv1, const X_curve& cv2
                         ) {}
  
  inline void merge_edge(const X_curve &cv1,
                         const X_curve &cv2,
                         Halfedge_handle e
                         ,const X_curve& cv
                         ) {}
  
  inline void remove_edge(Halfedge_handle e) {}
  inline void remove_edge(const Halfedge_handle_iterator& begin,
		const Halfedge_handle_iterator& end) {};
  inline void clear() {}
  inline void update(const Halfedge_handle_iterator&,
                     const Halfedge_handle_iterator&,
                     const Token& token)
  { token.rebuild_bounding_box(this); }

public:
  inline const Bounding_box* get_bounding_box() const 
  {return pm->get_bounding_box();}	
  inline const Traits* get_traits() const {return traits;}
  
protected:
  Halfedge_handle find_lowest(Vertex_handle v,Traits_wrap *traits, 
			      bool highest) const;
  
#ifdef CGAL_PM_DEBUG
  void debug(){}
#endif

protected:
  typedef const Self* cPLp;
  
protected:
  Planar_map* pm;
  Traits_wrap* traits;
};

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/Pm_naive_point_location.C>
#endif

#endif //CGAL_PM_NAIVE_POINT_LOCATION_H








