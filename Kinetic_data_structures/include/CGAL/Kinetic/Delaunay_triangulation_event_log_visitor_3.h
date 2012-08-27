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

#ifndef CGAL_KINETIC_DELAUNAY_TRIANGULATION_3_WATCHER_ELV_H
#define CGAL_KINETIC_DELAUNAY_TRIANGULATION_3_WATCHER_ELV_H
#include <CGAL/Kinetic/basic.h>
#include <vector>
#include <CGAL/Kinetic/Delaunay_triangulation_visitor_base_3.h>
#include <CGAL/Kinetic/internal/triangulation_helpers_3.h>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>

namespace CGAL { namespace Kinetic {

struct Delaunay_triangulation_event_log_visitor_3: public Delaunay_triangulation_visitor_base_3
{
    Delaunay_triangulation_event_log_visitor_3(){}

    template <class Edge>
    void pre_edge_flip(Edge e) {
        typedef typename Edge::first_type::value_type::Vertex_handle::value_type::Point Point;
        std::ostringstream out;
        Point a= internal::vertex_of_edge(e, 0)->point();
        Point b= internal::vertex_of_edge(e, 1)->point();
        out << "Flipping away edge {" << (std::min)(a,b) << ", " << (std::max)(a,b) << "}";
        log_.push_back(out.str());
        CGAL_LOG(Log::LOTS, "Logging: " << out.str() << std::endl);
    }
    template <class Edge>
    void post_facet_flip(Edge e) {
        typedef typename Edge::first_type::value_type::Vertex_handle::value_type::Point Point;
        std::ostringstream out;
        Point a= internal::vertex_of_edge(e, 0)->point();
        Point b= internal::vertex_of_edge(e, 1)->point();
        out << "Flipping in edge {" << (std::min)(a,b) << ", " << (std::max)(a,b) << "}";
        log_.push_back(out.str());
        CGAL_LOG(Log::LOTS, "Logging: " << out.str() << std::endl);
    }

    template <class Facet>
    void pre_facet_flip(Facet e) {
        typedef typename Facet::first_type::value_type::Vertex_handle::value_type::Point Point;
        std::ostringstream out;
        Point pts[3];
        pts[0]= internal::vertex_of_facet(e, 0)->point();
        pts[1]= internal::vertex_of_facet(e, 1)->point();
        pts[2]= internal::vertex_of_facet(e, 2)->point();
        std::sort(&pts[0], &pts[0]+3);
        out << "Flipping away facet {" << pts[0] << ", " 
	    << pts[1] << ", " << pts[2] <<"}";
        log_.push_back(out.str());
        CGAL_LOG(Log::LOTS, "Logging: " << out.str() << std::endl);
    }

    template <class Facet>
    void post_edge_flip(Facet e) {
        typedef typename Facet::first_type::value_type::Vertex_handle::value_type::Point Point;
        std::ostringstream out;
        Point pts[3];
        pts[0]= internal::vertex_of_facet(e, 0)->point();
        pts[1]= internal::vertex_of_facet(e, 1)->point();
        pts[2]= internal::vertex_of_facet(e, 2)->point();
        std::sort(&pts[0], &pts[0]+3);
        out << "Flipping in facet {" << pts[0] << ", " 
	    << pts[1] << ", " << pts[2] <<"}";
        log_.push_back(out.str());
        CGAL_LOG(Log::LOTS, "Logging: " << out.str() << std::endl);
    }

  template <class Vertex_handle>
    void post_insert_vertex(Vertex_handle e) {
     std::ostringstream out;
    out << "Inserted vertex " << e->point();
    log_.push_back(out.str());
    CGAL_LOG(Log::LOTS, "Logging: " << out.str() << std::endl);
  }
  template <class Vertex_handle>
    void pre_remove_vertex(Vertex_handle e) {
     std::ostringstream out;
    out << "Removing vertex " << e->point();
    log_.push_back(out.str());
    CGAL_LOG(Log::LOTS, "Logging: " << out.str() << std::endl);
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
