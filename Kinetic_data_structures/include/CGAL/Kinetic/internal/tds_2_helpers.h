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

#ifndef CGAL_TRIANGULATION_DATA_STRUCTURE_HELPER_2_H
#define CGAL_TRIANGULATION_DATA_STRUCTURE_HELPER_2_H
#include <CGAL/Kinetic/basic.h>
#include <utility>
namespace CGAL { namespace Kinetic { namespace internal {

template <class TDS>
struct Triangulation_data_structure_helper_2
{
  typedef typename TDS::Edge Edge;
  typedef typename TDS::Vertex_handle Vertex_handle;
  typedef typename TDS::Face Face;
  //typedef typename Face::Edge_label Edge_label;

  static int low_degree(Vertex_handle vh, const TDS &tds) {
    unsigned int r=3;
    typename TDS::Face_circulator f= tds.incident_faces(vh), s=f;
    ++f; ++f; ++f;
    while (f != s) {
      ++r;
      if (r==5) break;
      ++f;
    };
    CGAL_postcondition(vh->degree() < 5 && vh->degree() == r 
		       || vh->degree() >=5 && r==5);
    return r;
  }

  static Vertex_handle origin(const Edge &e) {
    int o= e.first->ccw(e.second);
    return e.first->vertex(o);
  }

  static Vertex_handle destination(const Edge &e) {
    int o= e.first->cw(e.second);
    return e.first->vertex(o);
  }  

  static Vertex_handle third_vertex(const Edge &e) {
    return e.first->vertex(e.second);
  }
};



} } } //namespace CGAL::Kinetic::internal
#endif
