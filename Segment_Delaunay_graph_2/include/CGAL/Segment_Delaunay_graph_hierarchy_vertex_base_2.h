// Copyright (c) 2003,2004,2005  INRIA Sophia-Antipolis (France).
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


#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_HIERARCHY_VERTEX_BASE_2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_HIERARCHY_VERTEX_BASE_2_H

#include <CGAL/license/Segment_Delaunay_graph_2.h>


#include <CGAL/Segment_Delaunay_graph_2/basic.h>

namespace CGAL {


template <class Vbb>
class Segment_Delaunay_graph_hierarchy_vertex_base_2
 : public Vbb
{
public:
  typedef Vbb V_Base;
  typedef typename V_Base::Data_structure     D_S;

  typedef typename V_Base::Site_2             Site_2;
  typedef typename V_Base::Storage_site_2     Storage_site_2;

  typedef D_S                                  Data_structure;

  typedef typename D_S::Vertex_handle         Vertex_handle;
  typedef typename D_S::Face_handle           Face_handle;

  template < typename DS2 >
  struct Rebind_TDS {
    typedef typename Vbb::template Rebind_TDS<DS2>::Other         Vb2;
    typedef Segment_Delaunay_graph_hierarchy_vertex_base_2<Vb2>   Other;
  };

  Segment_Delaunay_graph_hierarchy_vertex_base_2()
    : V_Base(), up_( Vertex_handle() ), down_( Vertex_handle() ) {}

  Segment_Delaunay_graph_hierarchy_vertex_base_2(const Storage_site_2& ss,
						 Face_handle f)
    : V_Base(ss,f), up_( Vertex_handle() ), down_( Vertex_handle() ) {}

public:  // for use in hierarchy only
  Vertex_handle up()   { return up_; }
  Vertex_handle down() { return down_; }
  void set_up(Vertex_handle u)   { up_ = u; }
  void set_down(Vertex_handle d) { down_ = d; }

private:
  Vertex_handle up_;    // same vertex one level above
  Vertex_handle down_;  // same vertex one level below
};


} //namespace CGAL


#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_HIERARCHY_VERTEX_BASE_2_H
