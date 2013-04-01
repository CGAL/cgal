// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$ $Date$
// 
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ron Wein        <wein@post.tau.ac.il>

#ifndef CGAL_GSP_AGG_OP_SWEEP_H
#define CGAL_GSP_AGG_OP_SWEEP_H

#include <CGAL/Sweep_line_2.h>
#include <CGAL/Unique_hash_map.h>

namespace CGAL {

template <class Arrangement_,
          class MetaTraits_,
          class SweepVisitor,
          class CurveWrap,
          class SweepEvent,
          typename Allocator = CGAL_ALLOCATOR(int) >
class Gps_agg_op_sweep_line_2 :
  public Sweep_line_2<MetaTraits_,
                      SweepVisitor,
                      CurveWrap,
                      SweepEvent,
                      Allocator>
{
public:

  typedef Arrangement_                            Arrangement_2;
  typedef MetaTraits_                             Traits_2;
  typedef typename Traits_2::Point_2              Point_2;
  typedef typename Traits_2::X_monotone_curve_2   X_monotone_curve_2;

  typedef typename Arrangement_2::Vertex_handle       Vertex_handle;
  typedef typename Arrangement_2::Halfedge_handle     Halfedge_handle;

  typedef std::pair<Arrangement_2 *,
                    std::vector<Vertex_handle> *>     Arr_entry;

  typedef Sweep_line_2<Traits_2,
                       SweepVisitor,
                       CurveWrap,
                       SweepEvent,
                       Allocator>                 Base;

  typedef SweepEvent                              Event;


  typedef typename Base::Event_queue_iterator     EventQueueIter;
  typedef typename Event::Subcurve_iterator       EventCurveIter;

  typedef typename Base::Base_event               Base_event;
  typedef typename Base_event::Attribute          Attribute;

  typedef typename Base::Base_subcurve            Base_subcurve;
  
  typedef CurveWrap                               Subcurve;

  typedef std::list<Subcurve*>                    SubCurveList;
  typedef typename SubCurveList::iterator         SubCurveListIter; 


  typedef typename Base::Status_line_iterator     StatusLineIter;

public:

  /*!
   * Constructor.
   * \param visitor A pointer to a sweep-line visitor object.
   */
  Gps_agg_op_sweep_line_2 (SweepVisitor* visitor) : 
    Base (visitor)
  {}

  /*!
   * Constructor.
   * \param traits A pointer to a sweep-line traits object.
   * \param visitor A pointer to a sweep-line visitor object.
   */
  Gps_agg_op_sweep_line_2 (Traits_2 *traits, SweepVisitor* visitor) :
    Base(traits, visitor)
  {}

  /*! Perform the sweep. */
  template<class CurveInputIterator>
  void sweep (CurveInputIterator curves_begin,
              CurveInputIterator curves_end,
              unsigned int lower,
              unsigned int upper,
              unsigned int jump,
              std::vector<Arr_entry>& arr_vec)
  {
    CGAL_assertion (this->m_queue->empty() && 
                    this->m_statusLine.size() == 0);

    typedef Unique_hash_map<Vertex_handle, Event *>    Vertices_map;
    typedef typename Traits_2::Compare_xy_2            Compare_xy_2;

    this->m_visitor->before_sweep();
    // Allocate all of the Subcurve objects as one block.
    this->m_num_of_subCurves = static_cast<unsigned int>(std::distance (curves_begin, curves_end));
    this->m_subCurves = 
      this->m_subCurveAlloc.allocate (this->m_num_of_subCurves);

    this->m_curves_pair_set.resize (2 * this->m_num_of_subCurves);

    // Initialize the event queue using the vertices vectors. Note that these
    // vertices are already sorted, we simply have to merge them
    Vertices_map       vert_map;
    Vertex_handle      vh;
    Vertex_handle      invalid_v;
    unsigned int       i = lower;
    unsigned int       n = static_cast<unsigned int>((arr_vec[i].second)->size());
    unsigned int       j;
    EventQueueIter     q_iter;
    bool               first = true;
    Attribute          event_type;
    Event             *event;

    for (j = 0;
         j < n && (vh = (*(arr_vec[i].second))[j]) != invalid_v;
         j++)
    {
      // Insert the vertices of the first vector one after the other.
      event_type = _type_of_vertex (vh);
      if (event_type == Base_event::DEFAULT)
        continue;

      event = this->_allocate_event (vh->point(), event_type,
                                     ARR_INTERIOR, ARR_INTERIOR);
      // \todo When the boolean set operations are exteneded to support
      //       unbounded curves, we will need here a special treatment.

      #ifndef CGAL_ARRANGEMENT_ON_SURFACE_2_H
        event->set_finite();
      #endif
      
      if (! first)
      {
        q_iter = this->m_queue->insert_after (q_iter, event);
      }
      else
      {
        q_iter = this->m_queue->insert (event);
        first = false;
      }

      vert_map[vh] = event;
    }

    Comparison_result  res = LARGER;
    Compare_xy_2       comp_xy = this->m_traits->compare_xy_2_object();
    EventQueueIter     q_end = this->m_queue->end();

    for (i += jump; i <= upper; i += jump)
    {
      // Merge the vertices of the other vectors into the existing queue.
      q_iter = this->m_queue->begin();
      n = static_cast<unsigned int>((arr_vec[i].second)->size());

      for (j = 0;
           j < n && (vh = (*(arr_vec[i].second))[j]) != invalid_v;
           j++)
      {
        event_type = _type_of_vertex (vh);
        if (event_type == Base_event::DEFAULT)
          continue;

        while (q_iter != q_end &&
               (res = comp_xy (vh->point(), (*q_iter)->point())) == LARGER)
        {
          ++q_iter;
        }

        if (res == SMALLER || q_iter == q_end)
        {
          event = this->_allocate_event (vh->point(), event_type,
                                         ARR_INTERIOR, ARR_INTERIOR);
          // \todo When the boolean set operations are exteneded to support
          //       unbounded curves, we will need here a special treatment.
          
          #ifndef CGAL_ARRANGEMENT_ON_SURFACE_2_H
             event->set_finite();
          #endif

          this->m_queue->insert_before (q_iter, event);
          vert_map[vh] = event;
        }
        else if (res == EQUAL)
        {
          // In this case q_iter points to an event already associated with
          // the vertex, so we just update the map:
          vert_map[vh] = *q_iter;
        }
      }
    }

    // Go over all curves (which are associated with halfedges) and associate
    // them with the events we have just created.
    unsigned int           index = 0;
    CurveInputIterator     iter;
    Halfedge_handle        he;
    Event                 *e_left;
    Event                 *e_right;

    for (iter = curves_begin; iter != curves_end; ++iter, index++)
    {
      // Get the events associated with the end-vertices of the current
      // halfedge.
      he = iter->data().halfedge();

      CGAL_assertion (vert_map.is_defined (he->source()));
      CGAL_assertion (vert_map.is_defined (he->target()));

      if ((Arr_halfedge_direction)he->direction() == ARR_LEFT_TO_RIGHT)
      {
        e_left = vert_map[he->source()];
        e_right = vert_map[he->target()];
      }
      else
      {
        e_left = vert_map[he->target()];
        e_right = vert_map[he->source()];
      }

      // Create the subcurve object.
      this->m_subCurveAlloc.construct (this->m_subCurves + index,
                                       this->m_masterSubcurve);

      (this->m_subCurves + index)->init (*iter);
      (this->m_subCurves + index)->set_left_event(e_left);
      (this->m_subCurves + index)->set_right_event(e_right);
    
      e_right->add_curve_to_left (this->m_subCurves + index);  
      this->_add_curve_to_right (e_left, this->m_subCurves + index);
    }

    // Perform the sweep:
    this->_sweep();
    this->_complete_sweep();
    this->m_visitor->after_sweep();

    return;
  }
    
private:
   
  /*!
   * Check if the given vertex is an endpoint of an edge we are going
   * to use in the sweep.
   */
  Attribute _type_of_vertex (Vertex_handle v)
  {
    typename Arrangement_2::Halfedge_around_vertex_circulator  first, circ;

    circ = first = v->incident_halfedges();
    do
    {
      // Check if the current edge separates two faces with unequal
      // containment flags (otherwise we will simply not keep it).
      if (circ->face()->contained() != circ->twin()->face()->contained())
      {
        if ((Arr_halfedge_direction)circ->direction() == ARR_LEFT_TO_RIGHT)
          return (Base_event::RIGHT_END);
        else
          return (Base_event::LEFT_END);
      }
      ++circ;

    } while (circ != first);

    // If we reached here, we should not keep this vertex.
    return (Base_event::DEFAULT);
  }
};


} //namespace CGAL

#endif
