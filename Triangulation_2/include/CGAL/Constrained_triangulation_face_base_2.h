// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Mariette Yvinec

#ifndef CGAL_CONSTRAINED_TRIANGULATION_FACE_BASE_2_H
#define CGAL_CONSTRAINED_TRIANGULATION_FACE_BASE_2_H

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_face_base_2.h>

namespace CGAL { 

template <class Gt, class Fb = Triangulation_face_base_2<Gt> >
class Constrained_triangulation_face_base_2
  :  public Fb
{
  typedef Fb                                           Base;
  typedef typename Fb::Triangulation_data_structure    TDS;
public:
  typedef Gt                                   Geom_traits;
  typedef TDS                                  Triangulation_data_structure;
  typedef typename TDS::Vertex_handle          Vertex_handle;
  typedef typename TDS::Face_handle            Face_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Fb::template Rebind_TDS<TDS2>::Other    Fb2;
    typedef Constrained_triangulation_face_base_2<Gt,Fb2>    Other;
  };


protected:
  bool C[3];
 
public:
  Constrained_triangulation_face_base_2()
    : Base()
  {
    set_constraints(false,false,false);
  }

  Constrained_triangulation_face_base_2(Vertex_handle v0, 
					Vertex_handle v1, 
					Vertex_handle v2)
    : Base(v0,v1,v2)
  {
    set_constraints(false,false,false);
  }

  Constrained_triangulation_face_base_2(Vertex_handle v0, 
					Vertex_handle v1, 
					Vertex_handle v2,
					Face_handle n0, 
					Face_handle n1, 
					Face_handle n2)
    : Base(v0,v1,v2,n0,n1,n2)
  {
    set_constraints(false,false,false);
  }


  Constrained_triangulation_face_base_2(Vertex_handle v0, 
					Vertex_handle v1, 
					Vertex_handle v2,
					Face_handle n0, 
					Face_handle n1, 
					Face_handle n2,
					bool c0, 
					bool c1, 
					bool c2 )
    : Base(v0,v1,v2,n0,n1,n2)
  {
    set_constraints(c0,c1,c2);
  }


  bool is_constrained(int i) const ;
  void set_constraints(bool c0, bool c1, bool c2) ;
  void set_constraint(int i, bool b);
  void reorient();
  void ccw_permute();
  void cw_permute();
  
};

template <class Gt, class Fb>
inline void
Constrained_triangulation_face_base_2<Gt,Fb>::
set_constraints(bool c0, bool c1, bool c2)
{
  C[0]=c0;
  C[1]=c1;
  C[2]=c2;
}

template <class Gt, class Fb>
inline void
Constrained_triangulation_face_base_2<Gt,Fb>::
set_constraint(int i, bool b)
{
  CGAL_triangulation_precondition( i == 0 || i == 1 || i == 2);
  C[i] = b;
}
    
template <class Gt, class Fb>
inline bool
Constrained_triangulation_face_base_2<Gt,Fb>::
is_constrained(int i) const
{
  return(C[i]);
}

template <class Gt, class Fb>
inline void
Constrained_triangulation_face_base_2<Gt,Fb>::
reorient()
{
  Base::reorient();
  set_constraints(C[1],C[0],C[2]);
}

template <class Gt, class Fb>
inline void
Constrained_triangulation_face_base_2<Gt,Fb>::
ccw_permute()
{
  Base::ccw_permute();
  set_constraints(C[2],C[0],C[1]);
}

template <class Gt, class Fb>
inline void
Constrained_triangulation_face_base_2<Gt,Fb>::
cw_permute()
{
  Base::cw_permute();
  set_constraints(C[1],C[2],C[0]);
}
  
} //namespace CGAL 
  
#endif //CGAL_CONSTRAINED_TRIANGULATION_FACE_BASE_2_H
