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

  typedef std::set<Event> Queue;
  typedef typename Queue::iterator iterator;
  typedef typename Queue::const_iterator const_iterator;
  
  Queue m_queue;

public:

  Event_queue() { }

  bool empty() const { return m_queue.empty(); }
  std::size_t size() const { return m_queue.size(); }
  iterator begin() { return m_queue.begin(); }
  iterator end() { return m_queue.end(); }
  const_iterator begin() const { return m_queue.begin(); }
  const_iterator end() const { return m_queue.end(); }

  void push (const Event& ev)
  {
    m_queue.insert (ev);
  }

  Event pop ()
  {
    Event out = *(m_queue.begin());
    m_queue.erase (m_queue.begin());
    return out;
  }

  void print() const
  {
    for (const Event& e : m_queue)
      std::cerr << e << std::endl;
  }

  void erase_vertex_events (KSR::size_t vertex_idx, std::vector<Event>& events)
  {
    iterator it = m_queue.begin();
    while (it != m_queue.end())
    {
      iterator current = it ++;
      if (current->vertex_idx() == vertex_idx)
      {
        events.push_back (*current);
        m_queue.erase(current);
      }
    }
  }
     
};


}} // namespace CGAL::KSR


#endif // CGAL_KSR_EVENT_QUEUE_H
