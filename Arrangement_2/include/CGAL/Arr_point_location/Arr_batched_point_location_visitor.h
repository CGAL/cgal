// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// $URL$
// $Id$
// 
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_ARR_BATCHED_POINT_LOCATION_VISITOR_H
#define CGAL_ARR_BATCHED_POINT_LOCATION_VISITOR_H

CGAL_BEGIN_NAMESPACE


#include <CGAL/Sweep_line_empty_visitor.h>
#include <CGAL/Object.h>
#include <CGAL/Arr_accessor.h>
#include <utility>


template< class Traits_, class OutputIerator, class Arrangement_ >
class Arr_batched_point_location_visitor : public Sweep_line_empty_visitor< Traits_ >
{
  typedef Arr_batched_point_location_visitor<Traits_,
                                             OutputIerator,
                                             Arrangement_>        Self;   
  typedef Arrangement_                                            Arrangement;
  typedef Traits_                                                 Traits;
  

  typedef Sweep_line_empty_visitor<Traits>             Base;
  typedef typename Base::Event                         Event;
  typedef typename Base::Subcurve                      Subcurve;
  typedef typename Base::SL_iterator                   SL_iterator;

   
  typedef typename Traits::X_monotone_curve_2            X_monotone_curve_2;
  typedef typename Traits::Point_2                       Point_2;
  typedef typename Traits::Base_Point_2                  Base_Point_2;
  typedef typename Arrangement::Halfedge_const_handle    Halfedge_const_handle;
  typedef typename Arrangement::Vertex_const_handle      Vertex_const_handle;
  typedef std::pair<Base_Point_2,Object>                 PL_Pair;
  typedef CGAL::Arr_accessor<Arrangement>                Arr_accessor;


  public:


  Arr_batched_point_location_visitor(OutputIerator out,
                                     const Arrangement& arr):
                                     m_out(out)
  {
    Arr_accessor             arr_access(const_cast<Arrangement&>(arr));
    //initialize m_top_fict
    Vertex_const_handle top_left_v = arr_access.top_left_fictitious_vertex();
    m_top_fict = top_left_v->incident_halfedges();
    if(m_top_fict->direction() == SMALLER)
      m_top_fict = m_top_fict->next()->twin();

    CGAL_assertion((m_top_fict->source() == 
                    arr_access.top_right_fictitious_vertex()) ||
                   ((m_top_fict->source()->boundary_in_x() == NO_BOUNDARY) &&
                    (m_top_fict->source()->boundary_in_y() == PLUS_INFINITY)));
    CGAL_assertion(m_top_fict->target() == 
                   arr_access.top_left_fictitious_vertex());
  }

  

  //after_handle_event
  //(above_on_event is true iff 'above' subcurve is on the event
  bool after_handle_event(Event* event, SL_iterator above, bool above_on_event)
  {
    if(!event->is_finite())
    {
      //its an event at infinity, we need to update m_top_fict in case its
      //in Y=+oo (vertical curve or curve with vertical asymptote)
      if(event->infinity_at_x() != NO_BOUNDARY)
        return true;

      Boundary_type y_inf = event->infinity_at_y();
      if(y_inf == PLUS_INFINITY)
        m_top_fict = m_top_fict->twin()->next()->twin();

      return true;
    }

    if(! event->is_query())
      return true;


    // ISOLATED VERTEX
    if(event->is_action())
    {
      Vertex_const_handle vh = event->get_point().get_vertex_handle();
      *m_out = std::make_pair(event->get_point().base_point(), 
                              CGAL::make_object(vh));
      ++m_out;
      return true;
    }

    // VERTEX
    if(event->has_right_curves() || event->has_left_curves())
    {
      if(event->has_right_curves())
      {
        Subcurve* sc = *(event->right_curves_begin());
        Halfedge_const_handle he = sc->get_last_curve().get_halfedge_handle();
        *m_out = std::make_pair(event->get_point().base_point(), 
                                CGAL::make_object(he->target()));
        ++m_out;
        return true;
      }
      Subcurve* sc = *(event->left_curves_begin());
      Halfedge_const_handle he = sc->get_last_curve().get_halfedge_handle();
      *m_out = std::make_pair (event->get_point().base_point(),
                               make_object(he->source()));
      ++m_out;
      return true;
    }

     //UNBOUNDED_FACE
    if(above == this ->status_line_end())
    {
      *m_out = std::make_pair (event->get_point().base_point(),
                               CGAL::make_object(m_top_fict->face()));
      ++m_out;
      return true;
    }

    // EDGE
    if(above_on_event)
    {
      Halfedge_const_handle he = 
        (*above)->get_last_curve().get_halfedge_handle();
      *m_out = std::make_pair (event->get_point().base_point(),
                               make_object(he));
      ++m_out;
      return true;
    }

    //FACE or UNBOUNDED_FACE
    Halfedge_const_handle he = 
      (*above)->get_last_curve().get_halfedge_handle();
    *m_out = std::make_pair (event->get_point().base_point(),
                             CGAL::make_object(he->face()));
    ++m_out;
    return true;
  }

  OutputIerator get_output_iterator()
  {
    return (m_out);
  }

protected:

  OutputIerator            m_out;
  Halfedge_const_handle    m_top_fict;
};



CGAL_END_NAMESPACE

#endif
