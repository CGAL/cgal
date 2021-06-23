// Copyright (c) 2019 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_KSR_3_EVENT_QUEUE_H
#define CGAL_KSR_3_EVENT_QUEUE_H

// #include <CGAL/license/Kinetic_shape_reconstruction.h>

// Boost includes.
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/composite_key.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/member.hpp>

// Internal includes.
#include <CGAL/KSR/utils.h>
#include <CGAL/KSR_3/Event.h>

namespace CGAL {
namespace KSR_3 {

template<typename Data_structure>
class Event_queue {

public:
  // Data structure types.
  using FT      = typename Data_structure::Kernel::FT;
  using PVertex = typename Data_structure::PVertex;
  using PEdge   = typename Data_structure::PEdge;
  using PFace   = typename Data_structure::PFace;
  using IVertex = typename Data_structure::IVertex;
  using IEdge   = typename Data_structure::IEdge;

  // Event types.
  using Event = KSR_3::Event<Data_structure>;
  using ETime = typename Event::ETime;

  // Boost queue.
  using Queue = boost::multi_index_container<
  Event, boost::multi_index::indexed_by<
    boost::multi_index::ordered_non_unique<
    boost::multi_index::member<Event, ETime, &Event::m_time> >,
    boost::multi_index::ordered_non_unique<
    boost::multi_index::member<Event, PVertex, &Event::m_pvertex> >,
    boost::multi_index::ordered_non_unique<
    boost::multi_index::member<Event, PVertex, &Event::m_pother> >,
    boost::multi_index::ordered_non_unique<
    boost::multi_index::composite_key<Event,
      boost::multi_index::member<Event, IEdge, &Event::m_iedge>,
      boost::multi_index::member<Event, std::size_t, &Event::m_support_plane_idx> > >
    > >;

  using Queue_by_time        = typename Queue::template nth_index<0>::type;
  using Queue_by_pvertex_idx = typename Queue::template nth_index<1>::type;
  using Queue_by_pother_idx  = typename Queue::template nth_index<2>::type;
  using Queue_by_iedge_idx   = typename Queue::template nth_index<3>::type;

  Event_queue(const bool verbose) :
  m_verbose(verbose)
  { }

  // Size.
  bool empty() const { return m_queue.empty(); }
  std::size_t size() const { return m_queue.size(); }
  void clear() { m_queue.clear(); }

  // Access.
  void push(const Event& event) {
    if (m_verbose) std::cout << "** pushing " << event << std::endl;
    m_queue.insert(event);
  }

  // Pop the event by the shortest time: short -> long
  Event pop() {

    // std::cout << "POPPING EVENTS: " << std::endl;
    // print();

    const auto event_iterator = queue_by_time().begin();
    const Event event = *event_iterator;
    m_queue.erase(event_iterator);

    const FT tol = KSR::tolerance<FT>();
    const FT time_diff = CGAL::abs(next().time() - event.time());
    if (time_diff < tol) {
      if (m_verbose) {
        std::cout << "WARNING: NEXT EVENT IS HAPPENNING AT THE SAME TIME!" << std::endl;
      }
    }
    return event;
  }

  // Get next event with the closest time.
  Event next() {
    return *queue_by_time().begin();
  }

  // Get next time within the range [min_time, max_time] that is greater
  // than curr_time by at least a tolerance.
  FT get_next_time(
    const FT min_time, const FT max_time, const FT curr_time) {

    const auto pother = Data_structure::null_pvertex();
    const auto ivertex = Data_structure::null_ivertex();
    const ETime e_min_time(min_time, pother, ivertex, true);
    const ETime e_max_time(max_time, pother, ivertex, true);

    const auto it_min = queue_by_time().lower_bound(e_min_time);
    const auto it_max = queue_by_time().upper_bound(e_max_time);
    const auto time_range = CGAL::make_range(it_min, it_max);

    for (const auto& event : time_range) {
      if (event.time() > (curr_time + KSR::tolerance<FT>())) {
        return event.time();
      }
    }
    CGAL_assertion(max_time > curr_time);
    return max_time;
  }

  // Erase all events of the iedge.
  void erase_vertex_events(
    const IEdge iedge,
    const std::size_t support_plane_idx) {

    // Erase by iedge.
    const auto pe = queue_by_iedge_idx().equal_range(
      boost::make_tuple(iedge, support_plane_idx));
    const auto pe_range = CGAL::make_range(pe);

    if (m_verbose) {
      for (const auto& event : pe_range)
        std::cout << "** erasing (by iedge) " << event << std::endl;
    }
    queue_by_iedge_idx().erase(pe.first, pe.second);
  }

  // Erase all events of the pvertex.
  void erase_vertex_events(const PVertex pvertex) {

    // Erase by pvertex.
    const auto pv = queue_by_pvertex_idx().equal_range(pvertex);
    const auto pv_range = CGAL::make_range(pv);

    if (m_verbose) {
      for (const auto& event : pv_range)
        std::cout << "** erasing (by pvertex) " << event << std::endl;
    }
    queue_by_pvertex_idx().erase(pv.first, pv.second);

    // Erase by pother.
    const auto po = queue_by_pother_idx().equal_range(pvertex);
    const auto po_range = CGAL::make_range(po);

    if (m_verbose) {
      for (const auto& event : po_range)
        std::cout << "** erasing (by pother) " << event << std::endl;
    }
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
  const bool m_verbose;
};

} // namespace KSR_3
} // namespace CGAL

#endif // CGAL_KSR_3_EVENT_QUEUE_H
