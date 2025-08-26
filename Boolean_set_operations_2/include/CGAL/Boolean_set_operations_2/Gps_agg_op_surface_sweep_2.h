// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Baruch Zukerman <baruchzu@post.tau.ac.il>
//             Ron Wein        <wein@post.tau.ac.il>

#ifndef CGAL_BSO_2_GSP_AGG_OP_SURFACE_SWEEP_2_H
#define CGAL_BSO_2_GSP_AGG_OP_SURFACE_SWEEP_2_H

#include <vector>

#include <CGAL/license/Boolean_set_operations_2.h>

#include <CGAL/Surface_sweep_2.h>
#include <CGAL/Unique_hash_map.h>

namespace CGAL {

namespace Ss2 = Surface_sweep_2;

template <typename Arrangement_, typename Visitor_>
class Gps_agg_op_surface_sweep_2 : public Ss2::Surface_sweep_2<Visitor_> {
public:
  using Arrangement_2 = Arrangement_;
  using Visitor = Visitor_;

  using Geometry_traits_2 = typename Visitor::Geometry_traits_2;

  using Arr = Arrangement_2;
  using Gt2 = Geometry_traits_2;

  using Point_2 = typename Gt2::Point_2;
  using X_monotone_curve_2 = typename Gt2::X_monotone_curve_2;

  using Vertex_handle = typename Arr::Vertex_handle;
  using Halfedge_handle = typename Arr::Halfedge_handle;

  using Arr_entry = std::pair<Arr*, std::vector<Vertex_handle> *>;

  using Base = Ss2::Surface_sweep_2<Visitor>;

  using Event = typename Visitor::Event;
  using Subcurve = typename Visitor::Subcurve;

  using EventQueueIter = typename Base::Event_queue_iterator;
  using EventCurveIter = typename Event::Subcurve_iterator;

  using Attribute = typename Event::Attribute;

  using SubCurveList = std::list<Subcurve*>;
  using SubCurveListIter = typename SubCurveList::iterator;

public:
  /*! Constructor.
   * \param visitor A pointer to a sweep-line visitor object.
   */
  Gps_agg_op_surface_sweep_2(Visitor* visitor) : Base(visitor) {}

  /*! Constructor.
   * \param traits A pointer to a sweep-line traits object.
   * \param visitor A pointer to a sweep-line visitor object.
   */
  Gps_agg_op_surface_sweep_2(Gt2* traits, Visitor* visitor) :
    Base(traits, visitor)
  {}

  /*! Perform the sweep. */
  template <typename CurveInputIterator>
  void sweep(CurveInputIterator curves_begin, CurveInputIterator curves_end,
             std::size_t lower, std::size_t upper, std::size_t jump,
             std::vector<Arr_entry>& arr_vec) {
    CGAL_assertion(this->m_queue->empty() && this->m_statusLine.size() == 0);

    using Vertices_map = Unique_hash_map<Vertex_handle, Event*>;
    using Compare_xy_2 = typename Gt2::Compare_xy_2;

    this->m_visitor->before_sweep();
    // Allocate all of the Subcurve objects as one block.
    this->m_num_of_subCurves = std::distance(curves_begin, curves_end);
    if (this->m_num_of_subCurves > 0)
      this->m_subCurves =
        this->m_subCurveAlloc.allocate(this->m_num_of_subCurves);


    // Initialize the event queue using the vertices vectors. Note that these
    // vertices are already sorted, we simply have to merge them
    Vertices_map vert_map;
    Vertex_handle vh;
    Vertex_handle invalid_v;
    std::size_t i = lower;
    auto n = (arr_vec[i].second)->size();
    std::size_t j;
    EventQueueIter q_iter;
    bool first = true;
    Attribute event_type;
    Event* event;

    for (j = 0; j < n && (vh = (*(arr_vec[i].second))[j]) != invalid_v; j++) {
      // Insert the vertices of the first vector one after the other.
      event_type = _type_of_vertex(vh);
      if (event_type == Event::DEFAULT) continue;

      event = this->_allocate_event(vh->point(), event_type,
                                    ARR_INTERIOR, ARR_INTERIOR);
      // \todo When the boolean set operations are extended to support
      //       unbounded curves, we will need here a special treatment.

      #ifndef CGAL_ARRANGEMENT_ON_SURFACE_2_H
        event->set_finite();
      #endif

      if (! first) {
        q_iter = this->m_queue->insert_after(q_iter, event);
      }
      else {
        q_iter = this->m_queue->insert(event);
        first = false;
      }

      vert_map[vh] = event;
    }

    Comparison_result res = LARGER;
    Compare_xy_2 comp_xy = this->m_traits->compare_xy_2_object();
    EventQueueIter q_end = this->m_queue->end();

    for (i += jump; i <= upper; i += jump) {
      // Merge the vertices of the other vectors into the existing queue.
      q_iter = this->m_queue->begin();
      n = (arr_vec[i].second)->size();

      for (j = 0; j < n && (vh = (*(arr_vec[i].second))[j]) != invalid_v; j++) {
        event_type = _type_of_vertex(vh);
        if (event_type == Event::DEFAULT) continue;

        while ((q_iter != q_end) &&
               (res = comp_xy(vh->point(), (*q_iter)->point())) == LARGER)
        {
          ++q_iter;
        }

        if (res == SMALLER || q_iter == q_end) {
          event = this->_allocate_event(vh->point(), event_type,
                                        ARR_INTERIOR, ARR_INTERIOR);
          // \todo When the boolean set operations are extended to support
          //       unbounded curves, we will need here a special treatment.

          #ifndef CGAL_ARRANGEMENT_ON_SURFACE_2_H
             event->set_finite();
          #endif

          this->m_queue->insert_before(q_iter, event);
          vert_map[vh] = event;
        }
        else if (res == EQUAL) {
          // In this case q_iter points to an event already associated with
          // the vertex, so we just update the map:
          vert_map[vh] = *q_iter;
        }
      }
    }

    // Go over all curves (which are associated with halfedges) and associate
    // them with the events we have just created.
    std::size_t index = 0;
    CurveInputIterator iter;
    Halfedge_handle he;
    Event* e_left;
    Event* e_right;

    for (iter = curves_begin; iter != curves_end; ++iter, index++) {
      // Get the events associated with the end-vertices of the current
      // halfedge.
      he = iter->data().halfedge();

      CGAL_assertion(vert_map.is_defined(he->source()));
      CGAL_assertion(vert_map.is_defined(he->target()));

      if ((Arr_halfedge_direction)he->direction() == ARR_LEFT_TO_RIGHT) {
        e_left = vert_map[he->source()];
        e_right = vert_map[he->target()];
      }
      else {
        e_left = vert_map[he->target()];
        e_right = vert_map[he->source()];
      }

      // Create the subcurve object.
      using Subcurve_alloc = decltype(this->m_subCurveAlloc);
      std::allocator_traits<Subcurve_alloc>::construct(this->m_subCurveAlloc,
                                                       this->m_subCurves + index,
                                                       this->m_masterSubcurve);
      (this->m_subCurves + index)->init(*iter);
      (this->m_subCurves + index)->set_left_event(e_left);
      (this->m_subCurves + index)->set_right_event(e_right);

      e_right->add_curve_to_left(this->m_subCurves + index);
      this->_add_curve_to_right(e_left, this->m_subCurves + index);
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
  Attribute _type_of_vertex(Vertex_handle v) {
    typename Arr::Halfedge_around_vertex_circulator first, circ;

    circ = first = v->incident_halfedges();
    do {
      // Check if the current edge separates two faces with unequal
      // containment flags (otherwise we will simply not keep it).
      if (circ->face()->contained() != circ->twin()->face()->contained()) {
        if ((Arr_halfedge_direction)circ->direction() == ARR_LEFT_TO_RIGHT)
          return (Event::RIGHT_END);
        else return (Event::LEFT_END);
      }
      ++circ;

    } while (circ != first);

    // If we reached here, we should not keep this vertex.
    return (Event::DEFAULT);
  }
};

} // namespace CGAL

#endif
