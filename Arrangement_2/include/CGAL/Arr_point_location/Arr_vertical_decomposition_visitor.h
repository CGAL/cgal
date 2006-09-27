// Copyright (c) 2006  Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein <wein@post.tau.ac.il>

#ifndef CGAL_ARR_VERTICAL_DECOMPOSITION_VISITOR_H
#define CGAL_ARR_VERTICAL_DECOMPOSITION_VISITOR_H

CGAL_BEGIN_NAMESPACE

#include <CGAL/Sweep_line_2_empty_visitor.h>
#include <CGAL/Arr_accessor.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Object.h>

/*!
 * A sweep-line visitor for the vertical decomposition operation.
 */
template< class Traits_, class Arrangement_>
class Arr_vertical_decomposition_visitor : public Empty_visitor<Traits_>
{
  typedef Traits_                                            Traits_2;
  typedef Arrangement_                                       Arrangement_2;

  typedef Arr_vertical_decomposition_visitor<Traits_2,
                                             Arrangement_2>  Self;
  
  typedef Empty_visitor<Traits_2>                            Base;
  typedef typename Base::Event                               Event;
  typedef typename Base::Subcurve                            Subcurve;
  typedef typename Base::SL_iterator                         SL_iterator;

   
  typedef typename Traits_2::X_monotone_curve_2          X_monotone_curve_2;
  typedef typename Traits_2::Point_2                     Point_2;
  typedef typename Traits_2::Base_Point_2                Base_Point_2;
  typedef typename Arrangement_2::Halfedge_const_handle  Halfedge_const_handle;
  typedef typename Arrangement_2::Vertex_const_handle    Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_around_vertex_const_circulator
                                       Halfedge_around_vertex_const_circulator;

  typedef Arr_accessor<Arrangement_2>                    Arr_accessor;

  typedef std::pair<CGAL::Object, CGAL::Object>            Vert_pair;
  typedef Unique_hash_map<Vertex_const_handle, Vert_pair>  Vert_map;
     
private:

  Halfedge_const_handle    m_top_he;
  Halfedge_const_handle    m_bottom_he;
  Vert_map                *m_vert_map;
  Vertex_const_handle      m_prev_vh;
  const typename Arrangement_2::Traits_2    *traits;

public:

  Arr_vertical_decomposition_visitor (const Arrangement_2& arr,
                                      Vert_map& map) :
    m_vert_map (&map),
    traits (arr.get_traits())
  {
    // Get a fictitious halfedge on the top edge of the bounding rectangle.
    Arr_accessor          arr_access (const_cast<Arrangement_2&>(arr));
    Vertex_const_handle   v_tl = arr_access.top_left_fictitious_vertex();
    
    m_top_he = v_tl->incident_halfedges();
    if (m_top_he->source()->infinite_in_y() != PLUS_INFINITY)
      m_top_he = m_top_he->next()->twin();

    // Get a fictitious halfedge on the bottom edge of the bounding rectangle.
    Vertex_const_handle   v_bl = arr_access.bottom_left_fictitious_vertex();
    
    m_bottom_he = v_bl->incident_halfedges();
    if (m_bottom_he->source()->infinite_in_y() != MINUS_INFINITY)
      m_bottom_he = m_bottom_he->next()->twin();
  }

  /*!
   * Notification after an event has been processed.
   * \param event The event point.
   * \param above A status-line iterator pointing on the subcurve above the
   *              event point.
   * \param above_on_event Is the subcurve above actually on the event.
   * \return Whether the event can be de-allocated.
   */
  bool after_handle_event(Event* event, SL_iterator above, bool above_on_event)
  {
    if (! event->is_finite())
    {
      // We only need to consider events with a finite x-cooridinates.
      if (event->infinity_at_x() != FINITE)
        return (true);

      Infinity_type y_inf = event->infinity_at_y();
      
      if (y_inf == PLUS_INFINITY)
      {
        // Update the fictitious top halfedge.
        m_top_he = m_top_he->twin()->next()->twin();
      }
      else
      {
        CGAL_assertion (y_inf == PLUS_INFINITY);

        // Update the fictitious bottom halfedge.
        m_bottom_he = m_bottom_he->prev();
      }
      return (true);
    }

    // Get the vertex handle associated with the current event.
    Vertex_const_handle vh = event->get_point().get_vertex_handle();

    // Check the feature from above.
    if (above == this->status_line_end())
    {
      // There is no concrete subcurve above the current event point, so the
      // ficitious top halfegde is above it.
      ((*m_vert_map)[vh]).second = CGAL::make_object(m_top_he);
    }
    else
    {
      // We have a valid subcurve above the event: get its halfedge handle
      // and associate it with the vertex.
      ((*m_vert_map)[vh]).second =
        CGAL::make_object((*above)->get_last_curve().get_halfedge_handle());
    }

    // Check if the previous vertex we handled has the same x-coordinate
    // as the current one (it lies vertically below the current vertex).
    const bool         prev_same_x =
      m_prev_vh != Vertex_const_handle() &&
      (traits->compare_x_2_object() (vh->point(),
                                     m_prev_vh->point()) == EQUAL);

    // Decrement the status-line iterator to reach the subcurve below the
    // event point. If the number of right subcurves associated with the
    // event is n_right, we decrement the iterator n_right times, then
    // check if it is possible to further decrement it.
    SL_iterator        below = above;
    int                n_right = event->get_num_right_curves();
    int                k;

    for (k = 0; k < n_right; k++)
      --below;

    if (below == this->status_line_begin())
    {
      if (prev_same_x)
      {
        // The previous vertex is vertically below the current vertex,
        // so we update their respective entries in the output map.
        // We first check if the two vertices are connected (by a vertical
        // segment) - if so, they "see" empty features.
        bool          vert_connected = false;
        
        if (! vh->is_isolated())
        {
          Halfedge_around_vertex_const_circulator   circ, first;

          first = circ = vh->incident_halfedges();
          do
          {
            if (circ->source() == m_prev_vh)
              vert_connected = true;
            ++circ;
          } while (! vert_connected && circ != first);
        }

        if (! vert_connected)
        {
          ((*m_vert_map)[vh]).first = CGAL::make_object(m_prev_vh);
          ((*m_vert_map)[m_prev_vh]).second = CGAL::make_object(vh);
        }
        else
        {
          ((*m_vert_map)[vh]).first = CGAL::Object();
          ((*m_vert_map)[m_prev_vh]).second = CGAL::Object();
        }
      }
      else
      {
        // There is no concrete subcurve below the current event point, so the
        // ficitious bottom halfegde is below it.
        ((*m_vert_map)[vh]).first = CGAL::make_object(m_bottom_he);
      }
    }
    else
    {
      // Decrement the iterator one more time, to reach the subcurve below the
      // event.
      --below;

      if (prev_same_x &&
          (traits->compare_y_at_x_2_object() (m_prev_vh->point(),
                                              (*below)->get_last_curve().
                                              base_curve()) != SMALLER))
      {
        // The previous vertex is vertically below the current vertex,
        // and above the subcurve that lies below vh. We therefore update
        // the respective entries in the output map accordingly.
        // We first check if the two vertices are connected (by a vertical
        // segment) - if so, they "see" empty features.
        bool          vert_connected = false;
        
        if (! vh->is_isolated())
        {
          Halfedge_around_vertex_const_circulator   circ, first;

          first = circ = vh->incident_halfedges();
          do
          {
            if (circ->source() == m_prev_vh)
              vert_connected = true;
            ++circ;
          } while (! vert_connected && circ != first);
        }

        if (! vert_connected)
        {
          ((*m_vert_map)[vh]).first = CGAL::make_object(m_prev_vh);
          ((*m_vert_map)[m_prev_vh]).second = CGAL::make_object(vh);
        }
        else
        {
          ((*m_vert_map)[vh]).first = CGAL::Object();
          ((*m_vert_map)[m_prev_vh]).second = CGAL::Object();
        }
      }
      else
      {
        // Get the halfedge handle of the subcurve below the current event and
        // associate it with its vertex.
        ((*m_vert_map)[vh]).first =
          CGAL::make_object((*below)->get_last_curve().get_halfedge_handle());
      }
    }

    // We are done with the current vertex, but we store it for future use.
    m_prev_vh = vh;

    // It is safe to deallocate the event.
    return (true);
  }
};

CGAL_END_NAMESPACE

#endif
