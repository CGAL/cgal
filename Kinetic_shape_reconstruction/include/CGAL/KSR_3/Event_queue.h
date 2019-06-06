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

#ifndef CGAL_KSR_3_EVENT_QUEUE_H
#define CGAL_KSR_3_EVENT_QUEUE_H

//#include <CGAL/license/Kinetic_shape_reconstruction.h>

#include <CGAL/KSR/utils.h>
#include <CGAL/KSR_3/Event.h>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/member.hpp>

namespace CGAL
{

namespace KSR_3
{

template <typename Data>
class Event_queue
{
public:
  typedef typename Data::Kernel Kernel;
  typedef typename Kernel::FT FT;
  typedef typename Data::PVertex PVertex;
  typedef typename Data::PEdge PEdge;
  typedef typename Data::PFace PFace;
  typedef typename Data::IVertex IVertex;
  typedef typename Data::IEdge IEdge;

  typedef KSR_3::Event<Data> Event;

private:

  typedef boost::multi_index_container
  <Event,
   boost::multi_index::indexed_by<
     boost::multi_index::ordered_non_unique
     <boost::multi_index::member<Event, FT, &Event::m_time> >,
     boost::multi_index::ordered_non_unique
     <boost::multi_index::member<Event, PVertex, &Event::m_pvertex> >,
     boost::multi_index::ordered_non_unique
     <boost::multi_index::member<Event, PVertex, &Event::m_pother> >
     > 
   > Queue;

  typedef typename Queue::iterator Queue_iterator;
  typedef typename Queue::template nth_index<0>::type Queue_by_time;
  typedef typename Queue_by_time::iterator Queue_by_time_iterator;
  typedef typename Queue::template nth_index<1>::type Queue_by_pvertex_idx;
  typedef typename Queue_by_pvertex_idx::iterator Queue_by_pvertex_idx_iterator;
  typedef typename Queue::template nth_index<2>::type Queue_by_pother_idx;
  typedef typename Queue_by_pother_idx::iterator Queue_by_pother_idx_iterator;
  
  Queue m_queue;

public:

  Event_queue() { }

  bool empty() const { return m_queue.empty(); }
  std::size_t size() const { return m_queue.size(); }

  void push (const Event& ev)
  {
    CGAL_KSR_CERR(4) << "**** Pushing " << ev << std::endl;
    m_queue.insert (ev);
  }

  const Queue_by_time& queue_by_time() const { return m_queue.template get<0>(); }
  const Queue_by_pvertex_idx& queue_by_pvertex_idx() const { return m_queue.template get<1>(); }
  const Queue_by_pother_idx& queue_by_pother_idx() const { return m_queue.template get<2>(); }
  Queue_by_time& queue_by_time() { return m_queue.template get<0>(); }
  Queue_by_pvertex_idx& queue_by_pvertex_idx() { return m_queue.template get<1>(); }
  Queue_by_pother_idx& queue_by_pother_idx() { return m_queue.template get<2>(); }

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

  void erase_vertex_events (PVertex pvertex)
  {
    {
      std::pair<Queue_by_pvertex_idx_iterator, Queue_by_pvertex_idx_iterator>
        range = queue_by_pvertex_idx().equal_range(pvertex);

      for (const auto& ev : CGAL::make_range(range))
        CGAL_KSR_CERR(4) << "**** Erasing " << ev << std::endl;
    
      queue_by_pvertex_idx().erase (range.first, range.second);
    }
    {
      std::pair<Queue_by_pother_idx_iterator, Queue_by_pother_idx_iterator>
        range = queue_by_pother_idx().equal_range(pvertex);

      for (const auto& ev : CGAL::make_range(range))
        CGAL_KSR_CERR(4) << "**** Erasing " << ev << std::endl;
    
      queue_by_pother_idx().erase (range.first, range.second);
    }
  }

};


}} // namespace CGAL::KSR_3


#endif // CGAL_KSR_3_EVENT_QUEUE_H
