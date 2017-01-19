// Copyright (c) 2003,2004  INRIA Sophia-Antipolis (France).
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
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_APOLLONIUS_GRAPH_HIERARCHY_VERTEX_BASE_2_H
#define CGAL_APOLLONIUS_GRAPH_HIERARCHY_VERTEX_BASE_2_H

#include <CGAL/license/Apollonius_graph_2.h>


#include <CGAL/basic.h>

namespace CGAL {

template < class Vbb >
class Apollonius_graph_hierarchy_vertex_base_2
 : public Vbb
{
  typedef Vbb                                                Base;
  typedef typename Base::Apollonius_graph_data_structure_2   Agds;

public:
  typedef typename Base::Site_2             Site_2;
  typedef Agds                              Apollonius_graph_data_structure_2;
  typedef typename Agds::Vertex_handle      Vertex_handle;
  typedef typename Agds::Face_handle        Face_handle;

  template < typename AGDS2 >
  struct Rebind_TDS {
    typedef typename Vbb::template Rebind_TDS<AGDS2>::Other   Vb2;
    typedef Apollonius_graph_hierarchy_vertex_base_2<Vb2>     Other;
  };

  Apollonius_graph_hierarchy_vertex_base_2()
    : Base(), _up(), _down()
    {}
  Apollonius_graph_hierarchy_vertex_base_2(const Site_2& p,
					   Face_handle f)
    : Base(p,f), _up(), _down()
    {}
  Apollonius_graph_hierarchy_vertex_base_2(const Site_2& p)
    : Base(p), _up(), _down()
    {}

  Vertex_handle up() {return _up;}
  Vertex_handle down() {return _down;}
  void set_up(Vertex_handle u) {_up=u;}
  void set_down(Vertex_handle d) {_down=d;}


 private:
  Vertex_handle  _up;    // same vertex one level above
  Vertex_handle  _down;  // same vertex one level below
};

} //namespace CGAL

#endif // CGAL_APOLLONIUS_GRAPH_HIERARCHY_VERTEX_BASE_2_H
