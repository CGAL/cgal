// Copyright (c) 2005  Stanford University (USA).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KDS_DELAUNAY_SORT_WATCHER_EVENT_LOG_H
#define CGAL_KDS_DELAUNAY_SORT_WATCHER_EVENT_LOG_H
#include <CGAL/KDS/basic.h>
#include <string>
#include <sstream>
#include <vector>

CGAL_KDS_BEGIN_NAMESPACE

struct Sort_event_log_visitor
{
    Sort_event_log_visitor(){}
   
    template <class Vertex_handle>
    void remove_vertex(Vertex_handle a) {
      typedef typename std::iterator_traits<Vertex_handle>::value_type Key;
      std::ostringstream out;
      out << "Removing vertex " << *a;
      log_.push_back(out.str());
      CGAL_KDS_LOG(LOG_LOTS, "Logging: " << out.str());
    }

    template <class Vertex_handle>
    void create_vertex(Vertex_handle a) {
      typedef typename std::iterator_traits<Vertex_handle>::value_type Key;
      std::ostringstream out;
      out << "Creating vertex " << *a;
      log_.push_back(out.str());
      CGAL_KDS_LOG(LOG_LOTS, "Logging: " << out.str());
    }

    template <class Vertex_handle>
    void modify_vertex(Vertex_handle a) {
      typedef typename std::iterator_traits<Vertex_handle>::value_type Key;
      std::ostringstream out;
      out << "Changing vertex " << *a;
      log_.push_back(out.str());
      CGAL_KDS_LOG(LOG_LOTS, "Logging: " << out.str());
    }

    template <class Vertex_handle>
    void before_swap(Vertex_handle a,Vertex_handle b) {
      typedef typename std::iterator_traits<Vertex_handle>::value_type Key;
      std::ostringstream out;
      out << "Before swap of " << *a << " and " << *b;
      log_.push_back(out.str());
      CGAL_KDS_LOG(LOG_LOTS, "Logging: " << out.str());

    }
    template <class Vertex_handle>
    void after_swap(Vertex_handle a, Vertex_handle b) {
      typedef typename std::iterator_traits<Vertex_handle>::value_type Key;
      std::ostringstream out; 
      out << "After swap of " << *a << " and " << *b;
      log_.push_back(out.str());
      CGAL_KDS_LOG(LOG_LOTS, "Logging: " << out.str());

    }

  typedef std::vector<std::string>::const_iterator iterator;
  iterator begin()  const
  {
    return log_.begin();
  }
  iterator end()  const
  {
    return log_.end();
  }

  size_t size() const
  {
    return log_.size();
  }
  
  std::vector<std::string> log_;
};

CGAL_KDS_END_NAMESPACE
#endif
