// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
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
// $URL$
// $Id$
// 
//
// Author(s)     : Mariette Yvinec

// Copyright notice to be modified

#ifndef CGAL_TRIANGULATION_FACE_BASE_WITH_IN_PLACE_EDGE_LIST_2_H
#define CGAL_TRIANGULATION_FACE_BASE_WITH_IN_PLACE_EDGE_LIST_2_H

#include <CGAL/basic.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_ds_face_base_2.h>

CGAL_BEGIN_NAMESPACE 

template < typename Gt, typename Fb = Triangulation_ds_face_base_2<> >
class Triangulation_face_base_with_in_place_edge_list_2 
  : public Fb
{
public:
  typedef Gt                                           Geom_traits;
  typedef typename Fb::Vertex_handle                   Vertex_handle;
  typedef typename Fb::Face_handle                     Face_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Fb::template Rebind_TDS<TDS2>::Other  Fb2;
    typedef Triangulation_face_base_with_in_place_edge_list_2<Gt, Fb2> Other;
  };

public:
  Triangulation_face_base_with_in_place_edge_list_2()
    : Fb()
  {
    initialize_in_place_edge_list();
  }

  Triangulation_face_base_with_in_place_edge_list_2(Vertex_handle v0, 
						    Vertex_handle v1, 
						    Vertex_handle v2)
    : Fb(v0,v1,v2)
  {
    initialize_in_place_edge_list();
  }

  Triangulation_face_base_with_in_place_edge_list_2(Vertex_handle v0, 
						    Vertex_handle v1, 
						    Vertex_handle v2,
						    Face_handle n0, 
						    Face_handle n1, 
						    Face_handle n2)
    : Fb(v0,v1,v2,n0,n1,n2)
  {
    initialize_in_place_edge_list();
  }

  static int ccw(int i) {return Triangulation_cw_ccw_2::ccw(i);}
  static int  cw(int i) {return Triangulation_cw_ccw_2::cw(i);}

  // in-place edge list stuff -- start
private:
  typedef std::pair<Face_handle,int>  Edge;

  bool selected[3];
  Edge next_[3], prev_[3];

  void initialize_in_place_edge_list()
  {
    static Edge SENTINEL_QUEUE_EDGE = Edge(Face_handle(), -1);

    for (int i = 0; i < 3; ++i) {
      selected[i] = false;
      next_[i] = SENTINEL_QUEUE_EDGE;
      prev_[i] = SENTINEL_QUEUE_EDGE;
    }
  }
public:
  Edge next(unsigned int i)     const { return next_[i]; }
  Edge previous(unsigned int i) const { return prev_[i]; }

  void set_next(unsigned int i, const Edge& next) { next_[i] = next; }
  void set_previous(unsigned int i, const Edge& prev) { prev_[i] = prev; }
  
  bool is_in_list(unsigned int i) const {
    return next_[i].second != -1 || prev_[i].second != -1;
  }

  void mark_selected(unsigned int i) { selected[i] = true; }
  void mark_unselected(unsigned int i) { selected[i] = false; }
  bool is_selected(unsigned int i) const { return selected[i]; }
  // in-place edge list stuff -- end

#ifndef CGAL_NO_DEPRECATED_CODE
  Vertex_handle mirror_vertex(int i) const;
  int mirror_index(int i) const;
#endif

};

#ifndef CGAL_NO_DEPRECATED_CODE
template < class Gt, class Fb >
inline
typename Triangulation_face_base_with_in_place_edge_list_2<Gt,Fb>::Vertex_handle
Triangulation_face_base_with_in_place_edge_list_2<Gt,Fb>::
mirror_vertex(int i) const
{
  CGAL_triangulation_precondition ( this->neighbor(i) != Face_handle()
				    && this->dimension() >= 1);
  //return neighbor(i)->vertex(neighbor(i)->index(this->handle()));
  return this->neighbor(i)->vertex(mirror_index(i));
}

template < class Gt, class Fb >
inline int
Triangulation_face_base_with_in_place_edge_list_2<Gt,Fb>::
mirror_index(int i) const
{
  // return the index of opposite vertex in neighbor(i);
  CGAL_triangulation_precondition (this->neighbor(i) != Face_handle() &&
	                           this->dimension() >= 1);
  if (this->dimension() == 1) {
    return 1 - (this->neighbor(i)->index(this->vertex(1-i)));
  }
  return this->ccw( this->neighbor(i)->index(this->vertex(this->ccw(i))));
}
#endif

CGAL_END_NAMESPACE 

#endif //CGAL_TRIANGULATION_FACE_BASE_WITH_IN_PLACE_EDGE_LIST_2_H
