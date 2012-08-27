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

#ifndef CGAL_KINETIC_DELAUNAY_TRIANGULATION_2_RE_WATCHER_BASE_H
#define CGAL_KINETIC_DELAUNAY_TRIANGULATION_2_RE_WATCHER_BASE_H
#include <CGAL/Kinetic/basic.h>
#include <CGAL/Kinetic/Delaunay_triangulation_visitor_base_2.h>
#include <set>

namespace CGAL { namespace Kinetic {

template <class Triangulation>
struct Delaunay_triangulation_recent_edges_visitor_2: public Delaunay_triangulation_visitor_base_2
{
  typedef typename Triangulation::Edge Edge;
  typedef typename Triangulation::Vertex_handle VH;
  Delaunay_triangulation_recent_edges_visitor_2(){}

  void pre_remove_vertex(VH) {
    recent_.clear();
  }
  void post_insert_vertex(VH) {
    recent_.clear();
  }

  void change_vertex(VH vh) {
    recent_.clear();
    typename Triangulation::Edge_circulator ec(vh), ef=ec;
    if (ec != NULL) {
      do {
	recent_.insert(*ec);
	++ec;
      } while (ec != ef);
    }
  }

  void pre_flip(Edge) {
    recent_.clear();
  }
  void post_flip(Edge e) {
    recent_.insert(e);
  }

  typedef typename std::set<Edge>::const_iterator iterator;
  iterator begin()  const
  {
    return recent_.begin();
  }
  iterator end()  const
  {
    return recent_.end();
  }

  bool contains(Edge e) const
  {
    return recent_.find(e) != recent_.end();
  }

  std::set<Edge> recent_;
};

} } //namespace CGAL::Kinetic
#endif
