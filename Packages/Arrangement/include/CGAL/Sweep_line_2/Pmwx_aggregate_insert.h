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
// Author(s)     : Tali Zvi <talizvi@post.tau.ac.il>

#ifndef CGAL_PMWX_AGGREGATE_INSERT_H
#define CGAL_PMWX_AGGREGATE_INSERT_H

#include <CGAL/Sweep_line_2/Sweep_line_2_impl.h>
#include <CGAL/Sweep_line_2/Pmwx_sweep_line_event.h>
#include <CGAL/Sweep_line_2/Pmwx_sweep_line_curve.h>
#include <CGAL/Sweep_line_2/Pmwx_sweep_line_notification.h>
#include <CGAL/assertions.h>
#include <list>
#include <vector>

CGAL_BEGIN_NAMESPACE

template <class SweepLineTraits_2, 
          class Arr,
          class Change_notification>
class Pmwx_aggregate_insert 
{
  typedef typename Arr::Halfedge_handle                             Halfedge_handle;
  typedef typename Arr::Edge_iterator                               Edge_iterator;
  typedef SweepLineTraits_2                                         Traits;
  typedef typename Traits::X_monotone_curve_2                       X_monotone_curve_2;
  typedef Pmwx_sweep_line_notification<Traits,
                                       Arr,
                                       Change_notification>         Visitor;

  typedef Pmwx_sweep_line_curve<Traits,Visitor, Halfedge_handle>    Subcurve; 
  typedef Pmwx_sweep_line_event<Traits, Subcurve, Visitor>          Event;

 
  typedef Sweep_line_2_impl<Traits,
                            Event,
                            Subcurve,
                            Visitor,
                            CGAL_ALLOCATOR(int)>       Sweep_line;
 


public:

  Pmwx_aggregate_insert(Arr *arr, Change_notification* ntf = NULL):
      m_traits(new Traits()),
      m_traits_owner(true),
      m_arr(arr),
      m_notif(ntf),
      m_visitor(arr, ntf, m_traits),
      m_sweep_line(m_traits, &m_visitor, false)
      {}


  Pmwx_aggregate_insert(Traits *traits, Arr *arr, Change_notification* ntf = NULL):
      m_traits(traits),
      m_traits_owner(false),
      m_arr(arr),
      m_notif(ntf),
      m_visitor(arr, ntf, m_traits),
      m_sweep_line(m_traits, &m_visitor, false)
      {}

  template<class CurveInputIterator>
  void insert_curves(CurveInputIterator begin, 
                     CurveInputIterator end)
  {
    if(m_arr->edges_begin() == m_arr->edges_end())
    {
       m_sweep_line.init(begin, end);
       m_sweep_line.sweep();
    }
    else
    {
      std::vector<X_monotone_curve_2>      curves_vec(begin,end);
      for (Edge_iterator eit = m_arr->edges_begin(); eit != m_arr->edges_end(); ++eit) 
      {
        curves_vec.push_back(eit->curve());
      }
      m_arr->clear();
      m_sweep_line.init(curves_vec.begin(), curves_vec.end());
      m_sweep_line.sweep();
    }
  }


  ~Pmwx_aggregate_insert()
  {
    if(m_traits_owner)
      delete m_traits;
  }


              
protected:

  Traits*              m_traits;
  bool                 m_traits_owner;
  Arr*                 m_arr;
  Change_notification* m_notif;
  Visitor              m_visitor;
  Sweep_line           m_sweep_line;
};

CGAL_END_NAMESPACE

#endif
