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
//                 Iddo Hanniel <hanniel@math.tau.ac.il>
#ifndef CGAL_PM_WALK_ALONG_LINE_POINT_LOCATION_H
#define CGAL_PM_WALK_ALONG_LINE_POINT_LOCATION_H

#include <CGAL/Pm_point_location_base.h>
#include <CGAL/Planar_map_2/Pm_traits_wrap_2.h>

//#define CGAL_PM_WALK_DEBUG
//#define CGAL_PM_DEBUG

#ifdef CGAL_PM_DEBUG
#ifndef CGAL_NAIVE_POINT_LOCATION_H
#include <CGAL/Pm_naive_point_location.h>
#endif
#endif

////////////////////////////////////////////////////////
//    WALK_ALONG_LINE STRATEGY
////////////////////////////////////////////////////////

CGAL_BEGIN_NAMESPACE

template <class Planar_map_> class Pm_walk_along_line_point_location : 
  public Pm_point_location_base<Planar_map_> {
public:
  typedef Planar_map_ Planar_map;
  typedef Pm_point_location_base<Planar_map> Base;
  typedef Pm_walk_along_line_point_location<Planar_map> Self;
  typedef typename Planar_map::Traits Traits;
  typedef typename Planar_map::Traits_wrap Traits_wrap;
  typedef typename Traits_wrap::Point Point;
  typedef typename Traits_wrap::X_curve X_curve;
  typedef typename Traits::Point_2                Point_2;
  typedef typename Traits::X_monotone_curve_2     X_monotone_curve_2;
  typedef typename Planar_map::Locate_type Locate_type;
  typedef typename Planar_map::Vertex_handle Vertex_handle;
  typedef typename Planar_map::Vertex_const_handle Vertex_const_handle;
  typedef typename Planar_map::Halfedge_handle Halfedge_handle;
  typedef typename Planar_map::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Planar_map::Face_handle Face_handle;
  typedef typename Planar_map::Face_const_handle Face_const_handle;
  typedef typename Planar_map::Ccb_halfedge_circulator Ccb_halfedge_circulator;
  typedef typename Planar_map::Ccb_halfedge_const_circulator
    Ccb_halfedge_const_circulator;
  typedef typename Planar_map::Holes_iterator Holes_iterator;
  typedef typename Planar_map::Holes_const_iterator Holes_const_iterator;
  typedef typename Planar_map::Halfedge_around_vertex_const_circulator Avcc;
  typedef typename Planar_map::Halfedge_iterator Halfedge_iterator;
  typedef Pm_bounding_box_base<Planar_map> Bounding_box;
  typedef typename Base::Halfedge_handle_iterator Halfedge_handle_iterator;
  typedef typename Base::Token Token;

protected:
  typedef const Self* const_Self_ptr;

public:
  // Constructor
  Pm_walk_along_line_point_location() : 
    Pm_point_location_base<Planar_map>(),
    pm(0),
    traits(0) {}
  
  void init(Planar_map & pmp, const Traits & tr) 
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

  Halfedge_const_handle vertical_ray_shoot(const Point& p, Locate_type& lt,
                                           bool up) const;
  Halfedge_handle vertical_ray_shoot(const Point& p, Locate_type& lt, bool up);

  inline void split_edge(const X_curve &, Halfedge_handle, Halfedge_handle,
                         //additions by iddo for arrangement
                         const X_curve &, const X_curve &) {}

  inline void merge_edge(const X_curve &, const X_curve &, Halfedge_handle, 
                         //additions by iddo for arrangement
                         const X_curve &) {}

  inline void remove_edge(Halfedge_handle) {}
  inline void remove_edge(const Halfedge_handle_iterator &,
                          const Halfedge_handle_iterator &) {};
  inline void clear() {}
  inline void update(const Halfedge_handle_iterator &,
                     const Halfedge_handle_iterator &,
                     const Token& token) 
  { token.rebuild_bounding_box(this); }

private:

  void walk_along_line(const Point & p, bool up, bool including,
                       Halfedge_const_handle & e, Locate_type & lt) const ;
  /* Simulates a walk along a vertical ray shoot whose shape is determined by 
     'up' and 'including'. e is the returned edge. */

  Halfedge_const_handle find_vertex_representation(Halfedge_const_handle e,
						   const Point & p,
						   bool up) const
  /* find the first halfedge pointing to p, when going clockwise
    if up==true - start from 6 oclock, else start from 12 oclock
    precondition:    e points to p.
    postcondition:   returned value points to p. */
{

#ifdef CGAL_PM_DEBUG

  CGAL_precondition(traits->point_equal(e->target()->point(),p));

#endif

  Avcc first = e,curr=first;
  ++curr;
  if (up)
    while(curr!=first)
    {
      if (traits->curves_compare_y_at_x_from_bottom(curr->curve(),e->curve(),p)
          ==SMALLER) 
      {
	e=curr;
        break;// this can't be improved
      }
      ++curr;
    }
  else
    while(curr!=first)
      {
        if (traits->curves_compare_y_at_x_from_top(curr->curve(),e->curve(),p)
	    ==SMALLER) 
        {
	  e=curr;
          break;// this can't be improved
        }
        ++curr;
      }

#ifdef CGAL_PM_DEBUG

  CGAL_postcondition(e!=pm->halfedges_end());
  CGAL_postcondition(traits->point_equal(e->target()->point(),p));

#endif

  return e;

}

  bool find_closest(const Point & p, const Ccb_halfedge_const_circulator & c,
                    bool up, bool including,
                    Halfedge_const_handle & e, Locate_type & lt) const;
  /* Finds the closest halfedge on a ccb along a vertical ray shoot.
     The bools 'up' and 'including' set the vertical ray shoot's shape.
     The return value is true iff such an halfedge exists.
     Locate type is 
     UNBOUNDED_FACE if point is outside ccb
     FACE if inside
     EDGE of on edge-boundary
     VERTEX if on vertex-boundary.
  */

#ifdef CGAL_PM_DEBUG

  void debug() {}

  void debug(const Halfedge_handle& e) const
    {
      {
        if (e!=pm->halfedges_end()) 
          std::cerr << "(" << e->source()->point() << "," 
		    << e->target()->point() << ")" << std::flush;
        else std::cerr << "(oo)";
      }
    }

#endif

public:
  inline const Traits * get_traits() const {return traits;}

protected:
  inline const Bounding_box * get_bounding_box() const 
  {
    return pm->get_bounding_box();
  }

  Planar_map * pm;
  const Traits_wrap * traits;
};
  
CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/Pm_walk_along_line_point_location.C>
#endif

#endif //PM_WALK_ALONG_LINE_POINT_LOCATION_H
