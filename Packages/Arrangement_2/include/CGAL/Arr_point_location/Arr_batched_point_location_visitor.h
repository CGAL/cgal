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
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef ARR_BATCHED_POINT_LOCATION_VISITOR_H
#define ARR_BATCHED_POINT_LOCATION_VISITOR_H

CGAL_BEGIN_NAMESPACE

//#include <CGAL/Arr_point_location/Arr_batched_point_location_event.h>
#include <CGAL/Sweep_line_2/Sweep_line_subcurve.h>
#include <CGAL/Sweep_line_2/Sweep_line_event.h>
#include <CGAL/Object.h>
#include <utility>


template<class Traits, class OutputIerator, class _Arrangement>
class Arr_batched_point_location_visitor
{
  typedef Arr_batched_point_location_visitor<Traits,
                                             OutputIerator,
                                             _Arrangement>        Self;   
  typedef _Arrangement                                            Arrangement;
  typedef Sweep_line_subcurve<Traits>                             Subcurve;
  typedef Sweep_line_event<Traits, Subcurve>                      Event;
  typedef typename Event::SubCurveIter                            SubCurveIter;

  typedef Basic_sweep_line_2<Traits,
                             Event,
                             Subcurve,
                             Self,
                             CGAL_ALLOCATOR(int)>                Sweep_line;

  typedef typename Sweep_line::StatusLineIter                   StatusLineIter;

   
  typedef typename Traits::X_monotone_curve_2            X_monotone_curve_2;
  typedef typename Traits::Point_2                       Point_2;
  typedef typename Arrangement::Halfedge_const_handle    Halfedge_const_handle;

  typedef typename Traits::Base_Point_2                  Base_Point_2;         
  typedef std::pair<Base_Point_2,Object>                      PL_Pair;


  public:


  Arr_batched_point_location_visitor(OutputIerator out,
                                     const Arrangement& arr):
                                     m_out(out),
                                     m_arr(arr)
  {}

  void attach(Sweep_line *sl)
  {
    m_sweep_line = sl;
  }


  void before_handle_event(Event* event){}

  //after_handle_event
  //(above_on_event is true iff 'above' subcurve is on the event
  bool after_handle_event(Event* event, StatusLineIter above, bool above_on_event)
  {
    if(! event->get_point().is_query())
      return true;

    // VERTEX
    if(event->has_right_curves() || event->has_left_curves())
    {
      Subcurve* sc = *(event->right_curves_begin());
      Halfedge_const_handle he = sc->get_last_curve().get_halfedge_handle();

      if(event->has_right_curves())
      {
        *m_out++ = std::make_pair(event->get_point(), make_object(he.target()));
        return true;
      }
      *m_out++ = std::make_pair(event->get_point(),make_object(he.source()));
      return true;
    }

     //UNBOUNDED_FACE
    if(above == m_sweep_line->status_line_end())
    {
      *m_out++ = std::make_pair(event->get_point(),
                                make_object(m_arr.unbounded_face()));
      return true;
    }

    // EDGE
    if(above_on_event)
    {
      Halfedge_const_handle he = (*above)->get_last_curve().get_halfedge_handle();
      *m_out++ = std::make_pair(event->get_point(),make_object(he));
      return true;
    }

    //FACE or UNBOUNDED_FACE
    Halfedge_const_handle he = (*above)->get_last_curve().get_halfedge_handle();
    *m_out++ = std::make_pair(event->get_point(),
                                make_object(he.face()));
    return true;
  }


 void add_subcurve(X_monotone_curve_2 cv,Subcurve* sc){}

  OutputIerator get_output_iterator()
  {
    return m_out;
  }


  void init_event(Event* event)
  {
   // event->set_query();
  }

  void after_sweep(){}
  void after_init(){}
  void update_event(Event* e,
                    const Point_2& end_point,
                    const X_monotone_curve_2& cv,
                    bool is_left_end)
  {}

  void update_event(Event* e,
                    Subcurve* sc1,
                    Subcurve* sc2,
                    bool created = false)
  {}

  void update_event(Event* e,
                    Subcurve* sc1)
  {}



protected:

  OutputIerator m_out;

  Sweep_line* m_sweep_line;

  const Arrangement& m_arr;
};



CGAL_END_NAMESPACE

#endif
