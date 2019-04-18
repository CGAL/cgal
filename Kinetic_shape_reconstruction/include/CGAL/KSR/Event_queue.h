// Copyright (c) 2019 GeometryFactory Sarl (France).
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

#ifndef CGAL_KSR_EVENT_QUEUE_H
#define CGAL_KSR_EVENT_QUEUE_H

//#include <CGAL/license/Kinetic_shape_reconstruction.h>

#include <CGAL/KSR/utils.h>
#include <CGAL/KSR/Event.h>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/member.hpp>

namespace CGAL
{

namespace KSR
{

template <typename GeomTraits>
class Event_queue
{
public:
  typedef GeomTraits Kernel;
  typedef typename Kernel::FT FT;

  typedef KSR::Event<GeomTraits> Event;

private:

  typedef boost::multi_index_container
  <Event,
   boost::multi_index::indexed_by<
     boost::multi_index::ordered_non_unique
     <boost::multi_index::member<Event, FT, &Event::m_time> >,
     boost::multi_index::ordered_non_unique
     <boost::multi_index::member<Event, KSR::size_t, &Event::m_vertex_idx> >
     > 
   > Queue;

  typedef typename Queue::iterator Queue_iterator;
  typedef typename Queue::template nth_index<0>::type Queue_by_time;
  typedef typename Queue_by_time::iterator Queue_by_time_iterator;
  typedef typename Queue::template nth_index<1>::type Queue_by_event_idx;
  typedef typename Queue_by_event_idx::iterator Queue_by_event_idx_iterator;
  
  Queue m_queue;

public:

  Event_queue() { }

  bool empty() const { return m_queue.empty(); }
  std::size_t size() const { return m_queue.size(); }

  void push (const Event& ev)
  {
    m_queue.insert (ev);
  }

  const Queue_by_time& queue_by_time() const { return m_queue.template get<0>(); }
  const Queue_by_event_idx& queue_by_event_idx() const { return m_queue.template get<1>(); }
  Queue_by_time& queue_by_time() { return m_queue.template get<0>(); }
  Queue_by_event_idx& queue_by_event_idx() { return m_queue.template get<1>(); }

  Event pop ()
  {
    Queue_iterator iter = queue_by_time().begin();
    Event out = *iter;
    m_queue.erase(iter);
    return out;
  }

  void print() const
  {
    for (const Event& e : m_queue)
      std::cerr << e << std::endl;
  }

  void erase_vertex_events (KSR::size_t vertex_idx, std::vector<Event>& events)
  {
    std::pair<Queue_by_event_idx_iterator, Queue_by_event_idx_iterator>
      range = queue_by_event_idx().equal_range(vertex_idx);

    std::copy (range.first, range.second, std::back_inserter (events));
    queue_by_event_idx().erase (range.first, range.second);
  }
     
};


}} // namespace CGAL::KSR


#endif // CGAL_KSR_EVENT_QUEUE_H
