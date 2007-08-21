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

#include <CGAL/Sweep_line_empty_visitor.h>
#include <CGAL/Arr_accessor.h>
#include <CGAL/Object.h>

/*!
 * A sweep-line visitor for the vertical decomposition operation.
 */
template< class Traits_, class Arrangement_, class OutputIterator_>
class Arr_vertical_decomposition_visitor :
  public Sweep_line_empty_visitor<Traits_>
{
  typedef Traits_                                       Traits_2;
  typedef Arrangement_                                  Arrangement_2;
  typedef OutputIterator_                               OutputIterator;

  typedef Arr_vertical_decomposition_visitor<Traits_2,
                                             Arrangement_2,
                                             OutputIterator> Self;
  
  typedef Sweep_line_empty_visitor<Traits_2>            Base;
  typedef typename Base::Event                          Event;
  typedef typename Base::Subcurve                       Subcurve;
  typedef typename Base::Status_line_iterator           SL_iterator;
 
  typedef typename Traits_2::X_monotone_curve_2         X_monotone_curve_2;
  typedef typename Traits_2::Point_2                    Point_2;
  typedef typename Traits_2::Base_point_2               Base_point_2;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_around_vertex_const_circulator
                                       Halfedge_around_vertex_const_circulator;

  typedef CGAL::Arr_accessor<Arrangement_2>             Arr_accessor;

  typedef std::pair<CGAL::Object, CGAL::Object>         Vert_pair;
  typedef std::pair<Vertex_const_handle, Vert_pair>     Vert_entry;
  
private:

  OutputIterator           m_oi;
  Halfedge_const_handle    m_top_he;
  Halfedge_const_handle    m_bottom_he;
  Vertex_const_handle      m_prev_vh;
  CGAL::Object             m_prev_obj_below;
  CGAL::Object             m_prev_obj_above;
  const typename Arrangement_2::Traits_2    *traits;

public:

  Arr_vertical_decomposition_visitor (const Arrangement_2& arr,
                                      OutputIterator oi) :
    m_oi (oi),
    traits (arr.get_traits())
  {
    // Get a fictitious halfedge on the top edge of the bounding rectangle.
    Arr_accessor          arr_access (const_cast<Arrangement_2&>(arr));

    Vertex_const_handle   v_tl = arr_access.arrangement().top_left_fictitious_vertex();
    
    m_top_he = v_tl->incident_halfedges();
    if (m_top_he->source()->boundary_in_y() != PLUS_INFINITY)
      m_top_he = m_top_he->next()->twin();

    // Get a fictitious halfedge on the bottom edge of the bounding rectangle.
    Vertex_const_handle   v_bl = arr_access.arrangement().bottom_left_fictitious_vertex();
    
    m_bottom_he = v_bl->incident_halfedges();
    if (m_bottom_he->source()->boundary_in_y() != MINUS_INFINITY)
      m_bottom_he = m_bottom_he->next()->twin();
  }

  /*!
   * Notification before the sweep starts.
   */
  void before_sweep ()
  {
    m_prev_vh = Vertex_const_handle();     // Invalid previous vertex.
    return;
  }

  /*!
   * Notification after an event has been processed.
   * \param event The event point.
   * \param above A status-line iterator pointing on the subcurve above the
   *              event point.
   * \param above_on_event Is the subcurve above actually on the event.
   * \return Whether the event can be de-allocated.
   */
  bool after_handle_event(Event* event, SL_iterator above,
                          bool /* above_on_event */)
  {
    if (! event->is_finite())
    {
      // We only need to consider events with a finite x-cooridinates.
      if (event->boundary_in_x() != NO_BOUNDARY)
        return (true);
      Boundary_type y_inf = event->boundary_in_y();
      
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

    // Get the vertex handle associated with the current event, and insert
    // it into the output iterator.
    Vertex_const_handle vh = event->point().vertex_handle();
    CGAL::Object        obj_above, obj_below;

    // Check the feature from above.
    if (above == this->status_line_end())
    {
      // There is no concrete subcurve above the current event point, so the
      // ficitious top halfegde is above it.
      obj_above = CGAL::make_object(m_top_he);
    }
    else
    {
      // We have a valid subcurve above the event: get its halfedge handle
      // and associate it with the vertex.
      obj_above =
        CGAL::make_object((*above)->last_curve().halfedge_handle());
    }

    // Check if the previous vertex we handled has the same x-coordinate
    // as the current one (it lies vertically below the current vertex).
    const bool         prev_same_x =
      m_prev_vh != Vertex_const_handle() &&
      (traits->compare_x_2_object() (vh->point(), m_prev_vh->point()) == EQUAL);

    // Decrement the status-line iterator to reach the subcurve below the
    // event point. If the number of right subcurves associated with the
    // event is n_right, we decrement the iterator n_right times, then
    // check if it is possible to further decrement it.
    SL_iterator        below = above;
    int                n_right = event->number_of_right_curves();
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
          obj_below = CGAL::make_object(m_prev_vh);
          m_prev_obj_above = CGAL::make_object(vh);
        }
        else
        {
          obj_below = CGAL::Object();
          m_prev_obj_above = CGAL::Object();
        }
      }
      else
      {
        // There is no concrete subcurve below the current event point, so the
        // ficitious bottom halfegde is below it.
        obj_below = CGAL::make_object(m_bottom_he);
      }
    }
    else
    {
      // Decrement the iterator one more time, to reach the subcurve below the
      // event.
      --below;

      if (prev_same_x &&
          (traits->compare_y_at_x_2_object() (m_prev_vh->point(),
                                              (*below)->last_curve().base()) !=
           SMALLER))
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
          obj_below = CGAL::make_object(m_prev_vh);
          m_prev_obj_above = CGAL::make_object(vh);
        }
        else
        {
          obj_below = CGAL::Object();
          m_prev_obj_above = CGAL::Object();
        }
      }
      else
      {
        // Get the halfedge handle of the subcurve below the current event and
        // associate it with its vertex.
        obj_below =
          CGAL::make_object((*below)->last_curve().halfedge_handle());
      }
    }

    // We can now create the entry for the previous vertex, as we are not
    // going to change the identity of the features below or above it.
    if (m_prev_vh != Vertex_const_handle())
    {
      *m_oi = Vert_entry (m_prev_vh, Vert_pair (m_prev_obj_below,
                                                m_prev_obj_above));
      ++m_oi;
    }

    // We are done with the current vertex, but we cannot create the entry
    // yet - so we store it for the next event.
    m_prev_vh = vh;
    m_prev_obj_below = obj_below;
    m_prev_obj_above = obj_above;

    // It is safe to deallocate the event.
    return (true);
  }

  /*!
   * A notification issued when the sweep process is over.
   */
  void after_sweep ()
  {
    // Create an entry for the last vertex (the xy-largest one). 
    if (m_prev_vh != Vertex_const_handle())
    {
      *m_oi = Vert_entry (m_prev_vh, Vert_pair (m_prev_obj_below,
                                                m_prev_obj_above));
      ++m_oi;
    }

    return;
  }

  /*! Get the output iterator of vertices sorted xy-lexicographically. */
  OutputIterator output_iterator()
  {
    return (m_oi);
  }
};

CGAL_END_NAMESPACE

#endif
