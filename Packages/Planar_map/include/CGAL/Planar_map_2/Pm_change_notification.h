// Copyright (c) 2000  Tel-Aviv University (Israel).
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
// Author(s)     : Eyal flato <flato@math.tau.ac.il>
#ifndef CGAL_PM_CHANGE_NOTIFICATION_H
#define CGAL_PM_CHANGE_NOTIFICATION_H

CGAL_BEGIN_NAMESPACE

template<class Planar_map_>
class Pm_change_notification
{
public:
  typedef Planar_map_ Planar_map;
  typedef typename Planar_map::Traits Traits;
  
  virtual void add_edge(const typename Traits::X_curve &, 
                        typename Planar_map::Halfedge_handle, 
                        bool /* original_direction */, bool overlap = false)
  {
      (void) overlap;
  }

  virtual void split_edge(typename Planar_map::Halfedge_handle /* org */,
                          typename Planar_map::Halfedge_handle /* new */,
                          const typename Traits::X_curve &,
                          const typename Traits::X_curve &)
  {
  }

  //  virtual void merge_edge(typename Planar_map::Halfedge_handle orig_edge, 
  //                          typename Planar_map::Halfedge_handle new_edge,
  //                          const typename Traits::X_curve & c)
  //    {
  //    }

  //   virtual void remove_edge(typename Planar_map::Halfedge_handle orig_edge)
  //    {
  //    }

  virtual void split_face(typename Planar_map::Face_handle /* org */, 
                          typename Planar_map::Face_handle /* new */)
  {
  }

  virtual void add_hole(typename Planar_map::Face_handle /* in_face */, 
                        typename Planar_map::Halfedge_handle /* new_hole */)
  {
  }

  virtual const typename Traits::X_curve &
  edge_support_curve(typename Planar_map::Halfedge_handle edge)
  {
    return edge->curve();
  }

  virtual bool have_support_curve()
  {
    return false;
  }

};

CGAL_END_NAMESPACE

#endif  // PM_CHANGE_NOTIFICATION
