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

#ifndef CGAL_PM_UNBOUNDING_BOX_H
#define CGAL_PM_UNBOUNDING_BOX_H

#include <CGAL/Planar_map_2.h>
#include <CGAL/Pm_bounding_box_base.h>

CGAL_BEGIN_NAMESPACE

template <class Planar_map_>
class Pm_unbounding_box : public Pm_bounding_box_base<Planar_map_> {
public:
  typedef Planar_map_ Planar_map;
  typedef Pm_unbounding_box<Planar_map> Self;
  //  typedef Planar_map_2<Dcel,Traits,Self> Planar_map;
  /*
  typedef Planar_map_2<Dcel,Traits> Base;
  typedef Planar_map_Bbox_2<Dcel,Traits> Self;
  typedef Pm_traits_wrap_2<Traits> Traits_wrap;
  */
  typedef typename Planar_map::Traits                   Traits;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits::Point_2                      Point_2;

  typedef typename std::vector<Point_2>::iterator       Point_iterator;
  typedef typename std::vector<X_monotone_curve_2>::iterator
                                                        X_curve_iterator;

  /*  typedef typename Traits::Bounding_box Bounding_box;
  typedef typename Traits::Boundary_type Boundary_type;
  typedef typename Traits::Point_boundary_container 
  Point_boundary_container;
  typedef typename Traits::X_curve_boundary_container 
  X_curve_boundary_container;
  */
  typedef typename Planar_map::Halfedge_handle          Halfedge_handle;
  typedef typename Planar_map::Face_handle              Face_handle;
  typedef typename Planar_map::Vertex_handle            Vertex_handle;
  typedef typename Planar_map::Vertex_const_handle      Vertex_const_handle;
  typedef typename Planar_map::Halfedge_const_handle    Halfedge_const_handle;
  typedef typename Planar_map::Face_const_handle        Face_const_handle;
  typedef typename Planar_map::Vertex_iterator          Vertex_iterator;
  typedef typename Planar_map::Halfedge_iterator        Halfedge_iterator;
  typedef typename Planar_map::Face_iterator            Face_iterator;
  typedef typename Planar_map::Vertex_const_iterator    Vertex_const_iterator;
  typedef typename Planar_map::Halfedge_const_iterator Halfedge_const_iterator;
  typedef typename Planar_map::Face_const_iterator      Face_const_iterator;
  typedef typename Planar_map::Locate_type              Locate_type;
  
  /*
  typedef typename Base::Halfedge_handle Halfedge_handle;
  typedef typename Base::Face_handle Face_handle;
  typedef typename Base::Vertex_handle Vertex_handle;
  typedef typename Base::Vertex_const_handle Vertex_const_handle;
  typedef typename Base::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Base::Face_const_handle Face_const_handle;
  typedef typename Base::Vertex_iterator Vertex_iterator;
  typedef typename Base::Halfedge_iterator Halfedge_iterator;
  typedef typename Base::Face_iterator Face_iterator;
  typedef typename Base::Vertex_const_iterator Vertex_const_iterator;
  typedef typename Base::Halfedge_const_iterator Halfedge_const_iterator;
  typedef typename Base::Face_const_iterator Face_const_iterator;
  typedef typename Base::Locate_type Locate_type;
  typedef Pm_point_location_base<Base> Point_location_base;
  //  typedef std::list<X_curve> X_curve_container;
  //  typedef Topological_map<_Dcel> TPM;
  typedef typename Base::Halfedge_around_vertex_circulator 
  Halfedge_around_vertex_circulator;
  typedef typename Base::Holes_iterator Holes_iterator;
  typedef typename Base::Holes_const_iterator Holes_const_iterator;
  typedef typename Base::Ccb_halfedge_const_circulator 
  Ccb_halfedge_const_circulator;
  typedef typename Base::Ccb_halfedge_circulator Ccb_halfedge_circulator;
  typedef typename Base::Size Size;
  */

  Pm_unbounding_box(){}
  ~Pm_unbounding_box(){}

  void init(Planar_map &, Traits &) {
    /*
    pm = &pmp;
    traits = (Traits_wrap*)(&tr);
    */
  }

  bool insert(const Point_2 &) {return true;}
#ifndef _MSC_VER
  bool insert(const Point_iterator &, const Point_iterator &)
#else
  // workaround for MSVC6.0
  bool insert(const Point_iterator &, const Point_iterator &,
              Point_2 * dummy = 0)
#endif
  {return true;}
  bool insert(const X_monotone_curve_2 &) {return true;}
#ifndef _MSC_VER
  bool insert(const X_curve_iterator &, const X_curve_iterator &)
#else
  // workaround for MSVC6.0
  bool insert(const X_curve_iterator &, const X_curve_iterator &,
              X_monotone_curve_2 * dummy = 0)
#endif          
  {return true;}

  /* The point location query function may updates the resulting 
     halfedge handle and locate type as expected from the bounding box */
  bool locate(const Point_2 &, Locate_type &, Halfedge_handle &) {return true;}

  bool vertical_ray_shoot(const Point_2 &, Locate_type &, bool,
                          Halfedge_handle &) {return true;}

  void split_edge(const X_monotone_curve_2 &, Halfedge_handle, Halfedge_handle,
                  const X_monotone_curve_2 &, const X_monotone_curve_2 &) {}

  void split_boundary_edge(const Halfedge_handle &,
                           Halfedge_handle, Halfedge_handle,
                           const Point_2 &) {}

  void merge_edge(const X_monotone_curve_2 &, const X_monotone_curve_2 &,
                  Halfedge_handle,
                  //additions by iddo for arrangement
                  const X_monotone_curve_2 &) {}

  void remove_edge(Halfedge_handle) {}
  inline bool is_empty() const {return false;}

#ifdef CGAL_PM_DEBUG
  void debug() const {}
#endif
};

CGAL_END_NAMESPACE

#endif //CGAL_PM_UNBOUNDING_BOX_H
