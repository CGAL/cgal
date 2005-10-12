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

#ifndef CGAL_KDS_DELAUNAY_TRIANGULATION_3_WATCHER_ELV_H
#define CGAL_KDS_DELAUNAY_TRIANGULATION_3_WATCHER_ELV_H
#include <CGAL/KDS/basic.h>
#include <vector>
#include <CGAL/KDS/Delaunay_triangulation_visitor_base_3.h>
#include <CGAL/KDS/internal/triangulation_helpers_3.h>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>

CGAL_KDS_BEGIN_NAMESPACE

struct Delaunay_triangulation_event_log_visitor_3: public Delaunay_triangulation_visitor_base_3  {
  Delaunay_triangulation_event_log_visitor_3(){}

  template <class Edge>
  void pre_edge_flip(Edge e){
    typedef typename Edge::first_type::value_type::Vertex_handle::value_type::Point Point;
    std::ostringstream out;
    Point a= internal::vertex_of_edge(e, 0)->point();
    Point b= internal::vertex_of_edge(e, 1)->point();
    out << "Flipping away edge {" << std::min(a,b) << ", " << std::max(a,b) << "}";
    log_.push_back(out.str());
    CGAL_KDS_LOG(LOG_LOTS, "Logging: " << out.str() << std::endl);
  }
  template <class Edge>
  void post_facet_flip(Edge e){
    typedef typename Edge::first_type::value_type::Vertex_handle::value_type::Point Point;
    std::ostringstream out;
    Point a= internal::vertex_of_edge(e, 0)->point();
    Point b= internal::vertex_of_edge(e, 1)->point();
    out << "Flipping in edge {" << std::min(a,b) << ", " << std::max(a,b) << "}";
    log_.push_back(out.str());
    CGAL_KDS_LOG(LOG_LOTS, "Logging: " << out.str() << std::endl);
  }

  template <class Facet>
  void pre_facet_flip(Facet e){
    typedef typename Facet::first_type::value_type::Vertex_handle::value_type::Point Point;
    std::ostringstream out;
    Point pts[3];
    pts[0]= internal::vertex_of_facet(e, 0)->point();
    pts[1]= internal::vertex_of_facet(e, 1)->point();
    pts[2]= internal::vertex_of_facet(e, 2)->point();
    std::sort(&pts[0], &pts[0]+3);
    out << "Flipping away facet {" << pts[0] << ", " << pts[1] << ", " << pts[2] <<"}";
    log_.push_back(out.str());
    CGAL_KDS_LOG(LOG_LOTS, "Logging: " << out.str() << std::endl);
  }
  
  template <class Facet>
  void post_edge_flip(Facet e){
    typedef typename Facet::first_type::value_type::Vertex_handle::value_type::Point Point;
    std::ostringstream out;
    Point pts[3];
    pts[0]= internal::vertex_of_facet(e, 0)->point();
    pts[1]= internal::vertex_of_facet(e, 1)->point();
    pts[2]= internal::vertex_of_facet(e, 2)->point();
    std::sort(&pts[0], &pts[0]+3);
    out << "Flipping in facet {" << pts[0] << ", " << pts[1] << ", " << pts[2] <<"}";
    log_.push_back(out.str());
    CGAL_KDS_LOG(LOG_LOTS, "Logging: " << out.str() << std::endl);
  }

  typedef std::vector<std::string>::const_iterator iterator;
  iterator begin()  const {
    return log_.begin();
  }

  iterator end()  const {
    return log_.end();
  }

  size_t size() const {
    return log_.size();
  }

  std::vector<std::string> log_;
};

CGAL_KDS_END_NAMESPACE

#endif
