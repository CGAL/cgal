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
  
  typedef std::map<KSR::size_t, std::vector<iterator> > Map;
  typedef typename Map::iterator Map_iterator;

  Queue m_queue;
  Map m_map_vertices;

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
    iterator it = m_queue.insert (ev).first;
    save_vertex_event(it);
  }

  Event pop ()
  {
    Event out = *(m_queue.begin());
    remove_vertex_event (out.vertex_idx(), m_queue.begin());
    m_queue.erase (m_queue.begin());
//    remove_vertex_events(out.vertex());
    return out;
  }

  void print()
  {
    for (const Event& e : m_queue)
      std::cerr << e << std::endl;
  }

  void remove_vertex_events (KSR::size_t vertex)
  {
    Map_iterator mit = m_map_vertices.find (vertex);
    CGAL_assertion (mit != m_map_vertices.end());

    for (const iterator& it : mit->second)
      m_queue.erase(it);

    m_map_vertices.erase(mit);
  }

  void transfer_vertex_events (KSR::size_t old_vertex, KSR::size_t new_vertex)
  {
    std::vector<iterator>& vec = m_map_vertices[old_vertex];
    for (iterator iter : vec)
      push (Event (new_vertex, iter->intersection_line_idx(), iter->time()));

    remove_vertex_events (old_vertex);
  }
  
private:

  void remove_vertex_event (KSR::size_t vertex, iterator it)
  {
    std::vector<iterator>& vec = m_map_vertices[vertex];
    vec.erase (std::find(vec.begin(), vec.end(), it));
  }

  void save_vertex_event (iterator it)
  {
    KSR::size_t vertex = it->vertex_idx();
    Map_iterator mit = m_map_vertices.insert (std::make_pair (vertex, std::vector<iterator>())).first;
    mit->second.push_back(it);
  }


};


}} // namespace CGAL::KSR


#endif // CGAL_KSR_EVENT_QUEUE_H
