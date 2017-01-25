// Copyright (c) 2004  INRIA Sophia-Antipolis (France).
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



#ifndef CGAL_TRIANGULATION_FACE_BASE_WITH_EDGES_2_H
#define CGAL_TRIANGULATION_FACE_BASE_WITH_EDGES_2_H

#include <CGAL/license/Segment_Delaunay_graph_2.h>



#include <CGAL/Triangulation_ds_face_base_2.h>
#include <CGAL/triangulation_assertions.h>



namespace CGAL {


template < class Gt,
	   class Fb = Triangulation_ds_face_base_2<> >
class Triangulation_face_base_with_edges_2
  : public Fb
{
protected:
  // local types
  typedef typename Fb::Triangulation_data_structure    TDS;

public:
  // TYPES
  //------
  typedef Gt                           Geom_traits;
  typedef Fb                           Base;
  typedef TDS                          Triangulation_data_structure_2;
  typedef typename TDS::Vertex_handle  Vertex_handle;
  typedef typename TDS::Face_handle    Face_handle;
  typedef typename TDS::Edge           Edge;


  template <typename TDS2>
  struct Rebind_TDS {
    typedef typename Fb::template Rebind_TDS<TDS2>::Other  Vb2;
    typedef Triangulation_face_base_with_edges_2<Gt,Vb2>   Other;
  }; 


public:
  // CREATION
  //---------
  Triangulation_face_base_with_edges_2() : Base()
  { init(); }

  Triangulation_face_base_with_edges_2(Vertex_handle v0,
				       Vertex_handle v1,
				       Vertex_handle v2)
    : Base(v0,v1,v2)
  { init(); }

  Triangulation_face_base_with_edges_2(Vertex_handle v0,
				       Vertex_handle v1,
				       Vertex_handle v2,
				       Face_handle n0,
				       Face_handle n1,
				       Face_handle n2)
    : Base(v0,v1,v2,n0,n1,n2)
  { init(); }

public:
  // OPERATIONS
  //-----------
  bool is_in_list(int i) const
  {
    CGAL_precondition( i >= 0 && i <= 2 );
    return ( next_edge_in_list[i].first != Face_handle() ||
	     prev_edge_in_list[i].first != Face_handle() );
  }

  void set_next(int i, const Edge& next)
  {
    CGAL_precondition( i >= 0 && i <= 2 );
    CGAL_precondition( next.first == Face_handle() ||
		       (next.second >= 0 && next.second <= 2) );
    next_edge_in_list[i] = next;
  }

  void set_previous(int i, const Edge&  prev)
  {
    CGAL_precondition( i >= 0 && i <= 2 );
    CGAL_precondition( prev.first == Face_handle() ||
		       (prev.second >= 0 && prev.second <= 2) );
    prev_edge_in_list[i] = prev;
  }

  Edge next(int i) const
  {
    CGAL_precondition( i >= 0 && i <= 2 );
    return next_edge_in_list[i];
  }

  Edge previous(int i) const
  {
    CGAL_precondition( i >= 0 && i <= 2 );
    return prev_edge_in_list[i];
  }


protected:
  // class variables
  Edge next_edge_in_list[3];
  Edge prev_edge_in_list[3];

protected:

  static int sentinel_index() { return -1; }

  static Edge sentinel_edge() {
    return Edge(Face_handle(), sentinel_index());
  }

  // initialization of in-place list pointers
  void init() {
    for (int i = 0; i < 3; i++) {
      next_edge_in_list[i] = sentinel_edge();
      prev_edge_in_list[i] = sentinel_edge();
    }
  }

};




} //namespace CGAL 

#endif // CGAL_TRIANGULATION_FACE_BASE_WITH_EDGES_2_H
