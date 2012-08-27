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

#ifndef CGAL_KINETIC_DELAUNAY_TRIANGULATION_2_LOG_WATCHER_BASE_H
#define CGAL_KINETIC_DELAUNAY_TRIANGULATION_2_LOG_WATCHER_BASE_H
#include <CGAL/Kinetic/basic.h>
#include <vector>
#include <sstream>
#include <CGAL/Kinetic/Delaunay_triangulation_visitor_base_2.h>

namespace CGAL { namespace Kinetic {

struct Delaunay_triangulation_event_log_visitor_2: public Delaunay_triangulation_visitor_base_2
{
    Delaunay_triangulation_event_log_visitor_2(){}

    template <class Edge>
    void pre_flip(Edge e) {
        typedef typename Edge::first_type::value_type::Vertex_handle::value_type::Point Point;
        std::ostringstream out;
        Point a= e.first->vertex((e.second+1)%3)->point();
        Point b= e.first->vertex((e.second+2)%3)->point();
        out << "Flipping away edge {" << (std::min)(a,b) << ", " << (std::max)(a,b) << "}";
        log_.push_back(out.str());
        CGAL_LOG(Log::LOTS, "Logging: " << out.str() << std::endl);
    }

    template <class Edge>
    void post_flip(Edge e) {
        typedef typename Edge::first_type::value_type::Vertex_handle::value_type::Point Point;
        std::ostringstream out;
        Point a= e.first->vertex((e.second+1)%3)->point();
        Point b= e.first->vertex((e.second+2)%3)->point();
        out << "Flipping in edge {" << (std::min)(a,b) << ", " << (std::max)(a,b) << "}";
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
