// Copyright (c) 2003  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>



#ifndef CGAL_APOLLONIUS_GRAPH_FACE_BASE_2_H
#define CGAL_APOLLONIUS_GRAPH_FACE_BASE_2_H


#include <CGAL/Apollonius_graph_short_names_2.h>


#include <CGAL/Triangulation_ds_face_base_2.h>
#include <CGAL/triangulation_assertions.h>



CGAL_BEGIN_NAMESPACE


template < class Gt,
	   class Fb = Triangulation_ds_face_base_2<> >
class Apollonius_graph_face_base_2
  : public Fb
{
protected:
  // local types
  typedef typename Fb::Triangulation_data_structure    AGDS;

public:
  // TYPES
  //------
  typedef Gt                            Geom_traits;
  typedef Fb                            Base;
  typedef AGDS                          Apollonius_graph_data_structure_2;
  typedef typename AGDS::Vertex_handle  Vertex_handle;
  typedef typename AGDS::Face_handle    Face_handle;
  typedef typename AGDS::Edge           Edge;


  template <typename AGDS2>
  struct Rebind_TDS {
    typedef typename Fb::template Rebind_TDS<AGDS2>::Other  Vb2;
    typedef Apollonius_graph_face_base_2<Gt,Vb2>            Other;
  }; 


public:
  // CREATION
  //---------
  Apollonius_graph_face_base_2() : Base()
  { init(); }

  Apollonius_graph_face_base_2(Vertex_handle v0,
			       Vertex_handle v1,
			       Vertex_handle v2)
    : Base(v0,v1,v2)
  { init(); }

  Apollonius_graph_face_base_2(Vertex_handle v0,
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
  inline bool is_in_list(int i) const
  {
    CGAL_precondition( i >= 0 && i <= 2 );
    return (next_edge_in_list[i].first != NULL ||
	    prev_edge_in_list[i].first != NULL);
  }

  inline void set_next(int i, const Edge& next)
  {
    CGAL_precondition( i >= 0 && i <= 2 );
    CGAL_precondition( next.first == NULL ||
				     (next.second >= 0 &&
				      next.second <= 2) );
    next_edge_in_list[i] = next;
  }

  inline void set_previous(int i, const Edge&  prev)
  {
    CGAL_precondition( i >= 0 && i <= 2 );
    CGAL_precondition( prev.first == NULL ||
		       (prev.second >= 0 && prev.second <= 2) );
    prev_edge_in_list[i] = prev;
  }

  inline Edge next(int i) const
  {
    CGAL_precondition( i >= 0 && i <= 2 );
    return next_edge_in_list[i];
  }

  inline Edge previous(int i) const
  {
    CGAL_precondition( i >= 0 && i <= 2 );
    return prev_edge_in_list[i];
  }


protected:
  // class variables
  Edge next_edge_in_list[3];
  Edge prev_edge_in_list[3];

protected:
  // initialization of in-place list pointers
  inline void init() {
    for (int i = 0; i < 3; i++) {
      next_edge_in_list[i] = Edge(Face_handle(NULL),-1);
      prev_edge_in_list[i] = Edge(Face_handle(NULL),-1);
    }
  }

};




CGAL_END_NAMESPACE 

#endif // CGAL_APOLLONIUS_GRAPH_FACE_BASE_2_H
