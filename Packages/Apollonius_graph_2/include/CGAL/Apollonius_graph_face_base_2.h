// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/Apollonius_graph_face_base_2.h
// package       : Apollonius_graph_2
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ======================================================================



#ifndef CGAL_APOLLONIUS_GRAPH_FACE_BASE_2_H
#define CGAL_APOLLONIUS_GRAPH_FACE_BASE_2_H

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
  typedef AGDS                          Apollonius_graph_data_structure;
  typedef typename AGDS::Vertex_handle  Vertex_handle;
  typedef typename AGDS::Face_handle    Face_handle;
  typedef typename AGDS::Edge           Edge;

  // Borland requires this.
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
    CGAL_triangulation_assertion( i >= 0 && i <= 2 );
    return (next_face_in_list[i] != NULL ||
	    prev_face_in_list[i] != NULL);
  }

  inline void set_next(int i, const Edge& next)
  {
    CGAL_triangulation_precondition( i >= 0 && i <= 2 );
    CGAL_triangulation_precondition( next.first == NULL ||
				     (next.second >= 0 &&
				      next.second <= 2) );
    next_face_in_list[i] = next.first;
    next_indx_in_list[i] = next.second;
  }

  inline void set_previous(int i, const Edge&  prev)
  {
    CGAL_triangulation_precondition( i >= 0 && i <= 2 );
    CGAL_triangulation_precondition( prev.first == NULL ||
				     (prev.second >= 0 &&
				      prev.second <= 2) );
    prev_face_in_list[i] = prev.first;
    prev_indx_in_list[i] = prev.second;
  }

  inline Edge next(int i) const
  {
    CGAL_triangulation_precondition( i >= 0 && i <= 2 );
    return Edge(next_face_in_list[i], next_indx_in_list[i]);
  }

  inline Edge previous(int i) const
  {
    CGAL_triangulation_precondition( i >= 0 && i <= 2 );
    return Edge(prev_face_in_list[i], prev_indx_in_list[i]);
  }


protected:
  // class variables
  Face_handle next_face_in_list[3];
  int         next_indx_in_list[3];
  Face_handle prev_face_in_list[3];
  int         prev_indx_in_list[3];

protected:
  // initialization of in-place list pointers
  inline void init() {
    for (int i = 0; i < 3; i++) {
      next_face_in_list[i] = Face_handle(NULL);
      prev_face_in_list[i] = Face_handle(NULL);
    }
  }

};




CGAL_END_NAMESPACE 

#endif // CGAL_APOLLONIUS_GRAPH_FACE_BASE_2_H
