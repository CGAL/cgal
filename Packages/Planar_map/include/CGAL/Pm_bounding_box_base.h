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
#ifndef CGAL_PM_BOUNDING_BOX_BASE_H
#define CGAL_PM_BOUNDING_BOX_BASE_H

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif
#include <vector>

CGAL_BEGIN_NAMESPACE


////////////////////////////////////////////////////////////////////
//               ABSTRACT BASE CLASS OF BOUNDING BOX
//////////////////////////////////////////////////////////////////

template <class Planar_map_>
class Pm_point_location_base;

template <class Planar_map_>
struct Pm_token {
  virtual void rebuild_bounding_box(const Pm_point_location_base<Planar_map_>*) const {};
};

template <class Planar_map_>
class Pm_bounding_box_base {
public:
  typedef Planar_map_ Planar_map;
  typedef typename Planar_map::Traits Traits;
  typedef typename Planar_map::Locate_type Locate_type;
  typedef typename Planar_map::Halfedge_handle Halfedge_handle;
  
  typedef typename Traits::Point Point;
  typedef typename Traits::X_curve X_curve;
  typedef typename std::vector<Point>::iterator Point_iterator;
  typedef typename std::vector<X_curve>::iterator X_curve_iterator;

  typedef Pm_point_location_base<Planar_map> Point_location_base;
  
  typedef Pm_token<Planar_map> Token;

  Pm_bounding_box_base() {}
  
  virtual void init(Planar_map& pmp, Traits& tr) = 0;
  virtual ~Pm_bounding_box_base() {}
  
  /* Returns true if bounding box remained unchanged */
  virtual bool insert(const Point& p) = 0;
  virtual bool insert(const Point_iterator& begin,const Point_iterator& end
#ifdef _MSC_VER
                      ,Point* dummy=0
#endif
                      )=0;

  virtual bool insert(const X_curve& cv) = 0;
  virtual bool insert(const X_curve_iterator& begin,
                      const X_curve_iterator& end
#ifdef _MSC_VER
                      ,X_curve* dummy=0
#endif
                      )=0;
  // workaround for MSVC6.0

  /* The point location query function may updates the resulting 
     halfedge handle and locate type as expected from the bounding box */
  /* Returns true if bounding box remained unchanged */
  virtual bool locate(const Point& p, Locate_type& lt,Halfedge_handle& h) = 0;
  virtual bool vertical_ray_shoot(const Point& p, Locate_type& lt, bool up,
                                  Halfedge_handle& h) = 0;
  
  //the function is called after the combinatoric split
  //cv is the original curve , e1 e2 are the new halfedges returned 
  virtual void split_edge(const X_curve &cv,
                          Halfedge_handle e1,
                          Halfedge_handle e2,
                          const X_curve& cv1, 
                          const X_curve& cv2
                          ) = 0;
  
  virtual void split_boundary_edge(const Halfedge_handle &h,
                                   Halfedge_handle h1,
                                   Halfedge_handle h2,
                                   const Point& p) =0;
  
  //called after combinatoric merge
  //e is the new edge cv1,cv2 are the original curves
  virtual void merge_edge(const X_curve &cv1,
                          const X_curve &cv2,
                          Halfedge_handle e 
                          //additions by iddo for arrangement
                          ,const X_curve& cv
                          //end additions
                          ) = 0; 
  
  //called before combinatoric deletion
  virtual void remove_edge(Halfedge_handle e) = 0;
  virtual void clear(){}
  virtual bool is_empty() const = 0;
  
#ifdef CGAL_PM_DEBUG
  virtual void debug() const = 0;
#endif
  
};

CGAL_END_NAMESPACE

#endif //CGAL_PM_BOUNDING_BOX_BASE_H














