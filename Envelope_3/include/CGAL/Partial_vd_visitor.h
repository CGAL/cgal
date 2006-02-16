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
// Author(s)     : Michal Meyerovitch     <gorgymic@post.tau.ac.il>

#ifndef CGAL_PARTIAL_VD_VISITOR_H
#define CGAL_PARTIAL_VD_VISITOR_H

#include <CGAL/Sweep_line_2_empty_visitor.h>
#include <CGAL/Object.h>
#include <utility>

#include <iostream>


CGAL_BEGIN_NAMESPACE

template< class Traits_, class Arrangement_, class OutputIterator>
class Partial_vd_visitor : public Empty_visitor< Traits_ >
{
public:

  typedef Partial_vd_visitor<Traits_,
                             Arrangement_,
                             OutputIterator>      Self;
  typedef Arrangement_                            Arrangement;
  typedef Traits_                                 Traits;
  

  typedef Empty_visitor<Traits>                   Base;
  typedef typename Base::Event                    Event;
  typedef typename Base::Subcurve                 Subcurve;
  typedef typename Base::SL_iterator              SL_iterator;
  typedef typename Base::SubCurveIter             SubCurveIter;
  typedef typename Base::SubCurveRevIter          SubCurveRevIter;

   
  typedef typename Traits::X_monotone_curve_2            X_monotone_curve_2;
  typedef typename Traits::Point_2                       Point_2;
  typedef typename Traits::Base_Point_2                  Base_Point_2;
  typedef typename Arrangement::Halfedge_const_handle    Halfedge_const_handle;
  typedef typename Arrangement::Vertex_const_handle      Vertex_const_handle;
  typedef std::pair<Object, Object>                      Vd_Pair;

  Partial_vd_visitor(const Arrangement& arr, OutputIterator o) :
    m_arr(arr), m_out(o), m_last_event(NULL), m_last_above(NULL),
    m_is_last_pair(false), m_last_should_shoot_down(false),
    m_last_should_shoot_up(false)
  {}

  //after_handle_event
  //(above_on_event is true iff 'above' subcurve is on the event
  bool after_handle_event(Event* event, SL_iterator above, bool above_on_event)
  {
    #ifdef CGAL_DEBUG_PARTIAL_VD
      std::cout << "in after_handle_event " << event->get_point() << std::endl;
    #endif
    // we should never have above_on_event = true since we sweep an existing
    // arrangement
    CGAL_assertion(!above_on_event);
    
    Vertex_const_handle vh = event->get_point().get_vertex_handle();
    CGAL_assertion(vh != Vertex_const_handle(NULL));
    CGAL_assertion(m_arr.get_traits()->equal_2_object()(event->get_point(),
                                                        vh->point()));

    // these booleans are needed for "LESANEN" vertical edges in a partial
    // vertical decomposition
    bool event_should_shoot_down = should_shoot_down(event);
    bool event_should_shoot_up = should_shoot_up(event);
    
    // this bool is used for solving 2 problems:
    // 1. the redundant problem, where 2 vertices see each other, and we might
    //    discover this pair twice
    // 2. a problem where when shooting down what we see is not exactly right
    //    below the current event, since there were events in the middle that
    //    were removed from the status line (with same x coordinate)
    bool previous_event_see_current = false;
    // if last event is set, we should check if we interrupt its above sight
    if (m_is_last_pair)
    {
      CGAL_assertion(m_last_event != NULL);
      // special cases can happen if we stay in the same x-coordinate sweepline
      if(m_arr.get_traits()->compare_x_2_object()
                    (m_last_event->get_point(), event->get_point())
         == EQUAL &&
         ((above == this->status_line_end() && m_last_above == NULL) ||
         (above != this->status_line_end() && m_last_above == *above) ||
         (event == m_last_above->get_right_event()) // last above ends in the current event
         ))
      {
        // should change the previous pair - is sees the current event point
        // instead of the previous above curve
        #ifdef CGAL_DEBUG_PARTIAL_VD
          std::cout << "found a point pair: " << m_last_event->get_point()
                    << " and " << event->get_point() << std::endl;
        #endif
        Vertex_const_handle pvh = m_last_event->get_point().get_vertex_handle();
        if (event_should_shoot_down || m_last_should_shoot_up)
          *m_out++ = std::make_pair(CGAL::make_object(pvh), CGAL::make_object(vh));
        else
        {
          #ifdef CGAL_DEBUG_PARTIAL_VD
            std::cout << "point pair doesn't follow partial vd rules" 
                      << std::endl;
          #endif
        }
        previous_event_see_current = true;
      }
      else
      {
        // keep the previous pair, if existed
        if (m_last_above != NULL)
        {
          #ifdef CGAL_DEBUG_PARTIAL_VD
            std::cout << "found a pair: " << m_last_event->get_point()
                      << " and " << m_last_above->get_last_curve() << std::endl;
          #endif
          Halfedge_const_handle phe = m_last_above->get_last_curve().get_halfedge_handle();
          Vertex_const_handle pvh = m_last_event->get_point().get_vertex_handle();
          CGAL_assertion(m_arr.get_traits()->equal_2_object()(m_last_event->get_point(),
                                                            pvh->point()));
          if (m_last_should_shoot_up)
            *m_out++ = std::make_pair(CGAL::make_object(pvh), CGAL::make_object(phe));
          else
          {
            #ifdef CGAL_DEBUG_PARTIAL_VD
              std::cout << "pair doesn't follow partial vd rules"
                        << std::endl;
            #endif
          }
        }
      }
    }

    // todo: check if we should shoot from this point

    // if the event is an end of a vertical segment, we don't shoot down
    bool vertical_end = false;
    SubCurveIter lci = event->left_curves_begin();
    for(; lci != event->left_curves_end(); ++lci)
    {
      if (m_arr.get_traits()->is_vertical_2_object()((*lci)->get_last_curve()))
      {
        vertical_end = true;
        break;
      }

    }
    if (!vertical_end)          
    {
      // shoot down
      // if previous event sees the current event, then we don't have to shoot
      // down, since this pair was already discovered
      if (!previous_event_see_current)
      {
        // find who is right below the event point in the sweep line
        int num_of_right_curves = event->get_num_right_curves();
        SL_iterator below = above;
        int i;
        for(i=0; i<=num_of_right_curves && below != this->status_line_begin();
            ++i, --below);

        if (i == (num_of_right_curves+1))
        {
          // we see something below us
          // is it not an isolated point, since these were already removed from
          // the status line
          // it is also not an endpoint, since in this case, we would have
          // previous_event_see_current = true
          #ifdef CGAL_DEBUG_PARTIAL_VD
            std::cout << "found a pair: " << event->get_point()
                      << " and " << (*below)->get_last_curve() << std::endl;
          #endif
          CGAL_assertion(m_arr.get_traits()->equal_2_object()(event->get_point(),
                                                              vh->point()));

          // assert it is not an endpoint
          CGAL_assertion_code(
            Event* left_event = (*below)->get_left_event();
          );
          CGAL_assertion(m_arr.get_traits()->compare_x_2_object()
                         (left_event->get_point(), event->get_point()) != EQUAL);
                                  
          Halfedge_const_handle he = (*below)->get_last_curve().get_halfedge_handle();
          if (event_should_shoot_down)
            *m_out++ = std::make_pair(CGAL::make_object(vh), CGAL::make_object(he));
          else
          {
            #ifdef CGAL_DEBUG_PARTIAL_VD
              std::cout << "pair doesn't follow partial vd rules"
                        << std::endl;
            #endif
          }
        }
      }
    }

    // if the event is a start of a vertical segment, we don't shoot up
    bool vertical_start = false;
    lci = event->right_curves_begin();
    for(; lci != event->right_curves_end(); ++lci)
    {
      if (m_arr.get_traits()->is_vertical_2_object()((*lci)->get_last_curve()))
      {
        vertical_start = true;
        break;
      }

    }
    m_last_event = event;
    m_last_above = NULL;
    m_last_should_shoot_down = event_should_shoot_down;
    m_last_should_shoot_up = event_should_shoot_up;
    
    m_is_last_pair = false;
    if (!vertical_start)
    {
      m_is_last_pair = true;
      // shoot up
      if (above != this->status_line_end())
      {
        m_last_above = *above;
        // we see something above us, and keep it to be dealt with in the next event
        CGAL_assertion(m_arr.get_traits()->equal_2_object()(event->get_point(),
                                                            vh->point()));
        
      }
    }

    m_events.push_back(event);
    return false;
  }


  void after_sweep()
  {
    typename std::list<Event*>::iterator it = m_events.begin();
    for(;it != m_events.end(); ++it)
      this->deallocate_event(*it);
  }

protected:
  // check if should shoot down/up from the evnet's point
  // in partial vertical decomposition, we should check if the angle
  // formed by the edges which we insert the vertical segment between
  // is less than 180
  bool should_shoot_down(Event* event)
  {
    // if no curve ends in this event we always get angel > 180 in the
    // down direction
    if (!event->has_left_curves())
      return true;

    // we have at least one curve that ends in this event, we take the lowest
    // one
    // if we have a vertical curve downwards, then it should be it (and we
    // don't shoot down)
    // (we can use get_last_curve() since we work on existing arrangement
    // and don't split the curves)
    SubCurveIter lci = event->left_curves_begin();
    if (m_arr.get_traits()->is_vertical_2_object()((*lci)->get_last_curve()))
      return false;

    // if we don't have any curves that begin in this event, the angel is
    // always > 180 in the down direction
    if (!event->has_right_curves())
      return true;

    // get the lowest right curve
    SubCurveIter rci = event->right_curves_begin();

    // TODO: be careful - this check is only true for segments -
    // we should check the angle between curves and not between points
    // should define it more carefully within the traits
    
    Point_2 left = (*lci)->get_last_curve().source();
    Point_2 right = (*rci)->get_last_curve().target();
    // if (left, event->point(), right) is a left turn can shoot down
//    return m_arr.get_traits()->leftturn_2_object()(left, event->get_point(), right);
    return CGAL::left_turn(left, event->get_point(), right);
  }

  bool should_shoot_up(Event* event)
  {
    // if no curve starts in this event we always get angel > 180 in the
    // up direction
    if (!event->has_right_curves())
      return true;

    // we have at least one curve that starts in this event, we take the
    // highest one
    // if we have a vertical curve upwards, then it should be it (and we
    // don't shoot up)
    // (we can use get_last_curve() since we work on existing arrangement
    // and don't split the curves)
    SubCurveRevIter rci = event->right_curves_rbegin();
    if (m_arr.get_traits()->is_vertical_2_object()((*rci)->get_last_curve()))
      return false;

    // if we don't have any curves that end in this event, the angel is
    // always > 180 in the up direction
    if (!event->has_left_curves())
      return true;

    // get the highest left curve
    SubCurveRevIter lci = event->left_curves_rbegin();

    // TODO: be careful - this check is only true for segments -
    // we should check the angle between curves and not between points
    // should define it more carefully within the traits

    Point_2 left = (*lci)->get_last_curve().source();
    Point_2 right = (*rci)->get_last_curve().target();
    // if (left, event->point(), right) is a right turn can shoot up
//    return m_arr.get_traits()->right_turn_2_object()(left, event->get_point(), right);
    return CGAL::right_turn(left, event->get_point(), right);

  }
  
  const Arrangement& m_arr;
  OutputIterator     m_out;
  
  // the last event that was completely handled
  Event*             m_last_event;
  // the subcurve that was above the last event in the status line of that event
  Subcurve*          m_last_above;
  // indicate if in the last event, when we shoot up, we might have a
  // vertical pair
  bool               m_is_last_pair;
  // last values of should_shoot_down/up
  bool               m_last_should_shoot_down, m_last_should_shoot_up;

  // save all events here, and deallocate them in the end
  std::list<Event*>  m_events;
};



CGAL_END_NAMESPACE

#endif
