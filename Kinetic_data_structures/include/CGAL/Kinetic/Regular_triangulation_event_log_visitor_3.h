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

#ifndef CGAL_KINETIC_REGULAR_TRIANGULATION_3_WATCHER_ELV_H
#define CGAL_KINETIC_REGULAR_TRIANGULATION_3_WATCHER_ELV_H
#include <CGAL/Kinetic/basic.h>
#include <CGAL/Kinetic/Delaunay_triangulation_event_log_visitor_3.h>

namespace CGAL { namespace Kinetic {

struct Regular_triangulation_event_log_visitor_3:
public Delaunay_triangulation_event_log_visitor_3
{
    typedef Delaunay_triangulation_event_log_visitor_3 P;
    Regular_triangulation_event_log_visitor_3(){}

    template <class Cell, class Str>
    void log_cell(Cell c, Str &out) const {
        typedef typename Cell::value_type::Vertex_handle::value_type::Point Point;
        Point pts[4];
        for (unsigned int i=0; i< 4; ++i){
            pts[i]= c->vertex(i)->point();
        }
        std::sort(&pts[0], &pts[0]+4);
        out << "{" << pts[0] << ", " << pts[1] << ", " 
	    << pts[2] << ", " << pts[3] << "}";
    }

    template <class Key, class Cell>
    void pre_move(Key k, Cell c) {
        std::ostringstream out;
        out << "Moving " << k << " from ";
        log_cell(c, out);
        log_.push_back(out.str());
        CGAL_LOG(Log::LOTS, "Logging: " << out.str() << std::endl);
    }

  template <class Key, class Cell>
    void post_move(Key k, Cell c) {
        std::ostringstream out;
        out << "Moved " << k << " from ";
        log_cell(c, out);
        log_.push_back(out.str());
        CGAL_LOG(Log::LOTS, "Logging: " << out.str() << std::endl);
    }

  /*template <class Key, class Cell>
    void pre_push(Key k, Cell c) {
        std::ostringstream out;
        out << "Pushing " << k << " into ";
        log_cell(c, out);
        log_.push_back(out.str());
        CGAL_LOG(Log::LOTS, "Logging: " << out.str() << std::endl);
    }

    template <class Vertex_handle>
    void post_push(Vertex_handle vh) {
        std::ostringstream out;
        out << "Pushed " << vh->point();
        log_.push_back(out.str());
        CGAL_LOG(Log::LOTS, "Logging: " << out.str() << std::endl);
    }

    template <class Vertex_handle>
    void pre_pop(Vertex_handle vh) {
        std::ostringstream out;
        out << "Popping " << vh->point();
        log_.push_back(out.str());
        CGAL_LOG(Log::LOTS, "Logging: " << out.str() << std::endl);
    }

    template <class Key, class Cell>
    void post_pop(Key k, Cell c) {
        std::ostringstream out;
        out << "Popped " << k << " from ";
        log_cell(c, out);
        log_.push_back(out.str());
        CGAL_LOG(Log::LOTS, "Logging: " << out.str() << std::endl);
	}*/
};

} } //namespace CGAL::Kinetic
#endif
