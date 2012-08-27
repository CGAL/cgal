// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KINETIC_DELAUNAY_SORT_WATCHER_EVENT_LOG_H
#define CGAL_KINETIC_DELAUNAY_SORT_WATCHER_EVENT_LOG_H
#include <CGAL/Kinetic/basic.h>
#include <string>
#include <sstream>
#include <vector>

namespace CGAL { namespace Kinetic {

struct Sort_event_log_visitor
{
    Sort_event_log_visitor(){}
   
    template <class Vertex_handle>
    void post_remove_vertex(Vertex_handle a) {
      typedef typename std::iterator_traits<Vertex_handle>::value_type Key;
      std::ostringstream out;
      out << "Removing vertex " << *a;
      log_.push_back(out.str());
      CGAL_LOG(Log::LOTS, "Logging: " << out.str());
    }

    template <class Vertex_handle>
    void post_insert_vertex(Vertex_handle a) {
      typedef typename std::iterator_traits<Vertex_handle>::value_type Key;
      std::ostringstream out;
      out << "Creating vertex " << *a;
      log_.push_back(out.str());
      CGAL_LOG(Log::LOTS, "Logging: " << out.str());
    }

    template <class Vertex_handle>
    void change_vertex(Vertex_handle a) {
      typedef typename std::iterator_traits<Vertex_handle>::value_type Key;
      std::ostringstream out;
      out << "Changing vertex " << *a;
      log_.push_back(out.str());
      CGAL_LOG(Log::LOTS, "Logging: " << out.str());
    }

    template <class Vertex_handle>
    void pre_swap(Vertex_handle a,Vertex_handle b) {
      typedef typename std::iterator_traits<Vertex_handle>::value_type Key;
      std::ostringstream out;
      out << "Before swap of " << *a << " and " << *b;
      log_.push_back(out.str());
      CGAL_LOG(Log::LOTS, "Logging: " << out.str());

    }
    template <class Vertex_handle>
    void post_swap(Vertex_handle a, Vertex_handle b) {
      typedef typename std::iterator_traits<Vertex_handle>::value_type Key;
      std::ostringstream out; 
      out << "After swap of " << *a << " and " << *b;
      log_.push_back(out.str());
      CGAL_LOG(Log::LOTS, "Logging: " << out.str());

    }

  typedef std::vector<std::string>::const_iterator Event_iterator;
  Event_iterator events_begin()  const
  {
    return log_.begin();
  }
  Event_iterator events_end()  const
  {
    return log_.end();
  }

  size_t size() const
  {
    return log_.size();
  }
  
  std::vector<std::string> log_;
};

} } //namespace CGAL::Kinetic
#endif
