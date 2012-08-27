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

#ifndef CGAL_KINETIC_DELAUNAY_TRIANGULATION_VERTEX_BASE_H
#define CGAL_KINETIC_DELAUNAY_TRIANGULATION_VERTEX_BASE_H

#include <CGAL/Kinetic/basic.h>
#include <CGAL/Triangulation_vertex_base_2.h>

namespace CGAL { namespace Kinetic {

template < typename GT,
           typename Vb = Triangulation_vertex_base_2<GT> >
class Delaunay_triangulation_vertex_base_2
  : public Vb
{
  int degree_;
public:
  typedef typename Vb::Face_handle                   Face_handle;
  typedef typename Vb::Point                         Point;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other          Vb2;
    typedef Delaunay_triangulation_vertex_base_2<GT, Vb2>   Other;
  };

  Delaunay_triangulation_vertex_base_2()
    : Vb(), degree_(-1){}

  Delaunay_triangulation_vertex_base_2(const Point & p)
    : Vb(p), degree_(-1){}

  Delaunay_triangulation_vertex_base_2(const Point & p, Face_handle c)
    : Vb(p, c), degree_(-1) {}

  Delaunay_triangulation_vertex_base_2(Face_handle c)
    : Vb(c), degree_(-1) {}

  unsigned int neighbors() const { return degree_; } // was abs
  //bool neighbors_is_changed() const {return degree_<0;}
  /*void set_neighbors_is_changed(bool tf) {
    if (tf) degree_= -std::abs(degree_);
    else degree_= std::abs(degree_);
    CGAL_postcondition(neighbors_is_changed()==tf);
    } */
  void set_neighbors(int d) {degree_=d;}
};


} } //namespace CGAL::Kinetic

#endif
