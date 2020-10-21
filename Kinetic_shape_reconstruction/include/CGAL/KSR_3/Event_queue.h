// Copyright (c) 2019 GeometryFactory SARL (France).
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
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_KSR_3_EVENT_QUEUE_H
#define CGAL_KSR_3_EVENT_QUEUE_H

// #include <CGAL/license/Kinetic_shape_reconstruction.h>

// Boost includes.
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/member.hpp>

// Internal includes.
#include <CGAL/KSR/utils.h>
#include <CGAL/KSR_3/Event.h>

namespace CGAL {
namespace KSR_3 {

// TODO: DOES NOT WORK WITH AN EXACT KERNEL! WHY?
// m_time for some reason evaluates to null ptr with no memory!
template<typename Data_structure>
class Event_queue {

public:
  // Kernel types.
  using Kernel = typename Data_structure::Kernel;
  using FT     = typename Kernel::FT;

  // Data structure types.
  using PVertex = typename Data_structure::PVertex;
  using PEdge   = typename Data_structure::PEdge;
  using PFace   = typename Data_structure::PFace;
  using IVertex = typename Data_structure::IVertex;
  using IEdge   = typename Data_structure::IEdge;

  // Event types.
  using Event = KSR_3::Event<Data_structure>;

  // Boost queue.
  using Queue = boost::multi_index_container<
  Event, boost::multi_index::indexed_by<
    boost::multi_index::ordered_non_unique
    <boost::multi_index::member<Event, FT, &Event::m_time> >,
    boost::multi_index::ordered_non_unique
    <boost::multi_index::member<Event, PVertex, &Event::m_pvertex> >,
    boost::multi_index::ordered_non_unique
    <boost::multi_index::member<Event, PVertex, &Event::m_pother> >,
    boost::multi_index::ordered_non_unique
    <boost::multi_index::member<Event, IEdge, &Event::m_iedge> >
    > >;

  using Queue_by_time        = typename Queue::template nth_index<0>::type;
  using Queue_by_pvertex_idx = typename Queue::template nth_index<1>::type;
  using Queue_by_pother_idx  = typename Queue::template nth_index<2>::type;
  using Queue_by_iedge_idx   = typename Queue::template nth_index<3>::type;

  // Size.
  const bool empty() const { return m_queue.empty(); }
  const KSR::size_t size() const { return m_queue.size(); }

  // Access.
  void push(const Event& event) {
    std::cout << "**** Pushing " << event << std::endl;
    m_queue.insert(event);
  }

  // Pop the event by the shortest time: short -> long
  const Event pop() {

    const auto event_iterator = queue_by_time().begin();
    const Event event = *event_iterator;
    m_queue.erase(event_iterator);

    if (queue_by_time().begin()->m_time == event.m_time)
      std::cerr << "WARNING: next Event is happening at the same time!" << std::endl;
    else if (CGAL::abs(queue_by_time().begin()->m_time - event.m_time) < 1e-15)
      std::cerr << "WARNING: next Event is happening at almost the same time!" << std::endl;
    return event;
  }

  // Get next event with the closest time.
  const Event next() {
    return *queue_by_time().begin();
  }

  // Erase all events of the iedge.
  void erase_vertex_events(const IEdge iedge) {

    // Erase by iedge.
    const auto pe = queue_by_iedge_idx().equal_range(iedge);
    const auto pe_range = CGAL::make_range(pe);

    for (const auto& event : pe_range)
      std::cout << "**** Erasing (by iedge) " << event << std::endl;
    queue_by_iedge_idx().erase(pe.first, pe.second);
  }

  // Erase all events of the pvertex.
  void erase_vertex_events(const PVertex pvertex) {

    // Erase by pvertex.
    const auto pv = queue_by_pvertex_idx().equal_range(pvertex);
    const auto pv_range = CGAL::make_range(pv);

    for (const auto& event : pv_range)
      std::cout << "**** Erasing (by pvertex) " << event << std::endl;
    queue_by_pvertex_idx().erase(pv.first, pv.second);

    // Erase by pother. TODO: Why is pother here?
    const auto po = queue_by_pother_idx().equal_range(pvertex);
    const auto po_range = CGAL::make_range(po);

    for (const auto& event : po_range)
      std::cout << "**** Erasing (by pother) " << event << std::endl;
    queue_by_pother_idx().erase(po.first, po.second);
  }

  // Sorting.
  const Queue_by_time& queue_by_time() const { return m_queue.template get<0>(); }
  const Queue_by_pvertex_idx& queue_by_pvertex_idx() const { return m_queue.template get<1>(); }
  const Queue_by_pother_idx& queue_by_pother_idx() const { return m_queue.template get<2>(); }
  const Queue_by_iedge_idx& queue_by_iedge_idx() const { return m_queue.template get<3>(); }

  Queue_by_time& queue_by_time() { return m_queue.template get<0>(); }
  Queue_by_pvertex_idx& queue_by_pvertex_idx() { return m_queue.template get<1>(); }
  Queue_by_pother_idx& queue_by_pother_idx() { return m_queue.template get<2>(); }
  Queue_by_iedge_idx& queue_by_iedge_idx() { return m_queue.template get<3>(); }

  // Helpers.
  void print() const {
    for (const auto& event : m_queue)
      std::cout << event << std::endl;
  }

private:
  Queue m_queue;
};

} // namespace KSR_3
} // namespace CGAL

#endif // CGAL_KSR_3_EVENT_QUEUE_H
