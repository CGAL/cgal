// Copyright (c) 2016  INRIA Nancy (France).
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
// $URL: 
// $Id: 
// 
//
// Author(s)     : Monique Teillaud <Monique.Teillaud@inria.fr>
//                 Mikhail Bogdanov

#ifndef CGAL_HYPERBOLIC_TRIANGULATION_FACE_BASE_2_H
#define CGAL_HYPERBOLIC_TRIANGULATION_FACE_BASE_2_H

#include <CGAL/basic.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_ds_face_base_2.h>

namespace CGAL { 

template < typename Gt, typename Fb = Triangulation_ds_face_base_2<> >
class Hyperbolic_triangulation_face_base_2 
  : public Fb
{
public:
  typedef Gt                                           Geom_traits;
  typedef typename Fb::Vertex_handle                   Vertex_handle;
  typedef typename Fb::Face_handle                     Face_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Fb::template Rebind_TDS<TDS2>::Other  Fb2;
    typedef Hyperbolic_triangulation_face_base_2<Gt, Fb2>             Other;
  };

public:
  Hyperbolic_triangulation_face_base_2()
    : Fb(), _is_finite_non_hyperbolic(false), _non_hyperbolic_edge(UCHAR_MAX)
  {}

  Hyperbolic_triangulation_face_base_2(Vertex_handle v0, 
			               Vertex_handle v1, 
			               Vertex_handle v2)
    : Fb(v0,v1,v2), _is_finite_non_hyperbolic(false), _non_hyperbolic_edge(UCHAR_MAX)
  {}

  Hyperbolic_triangulation_face_base_2(Vertex_handle v0, 
			               Vertex_handle v1, 
			               Vertex_handle v2,
			               Face_handle n0, 
			               Face_handle n1, 
			               Face_handle n2)
    : Fb(v0,v1,v2,n0,n1,n2), _is_finite_non_hyperbolic(false), _non_hyperbolic_edge(UCHAR_MAX)
  {}

  static int ccw(int i) {return Triangulation_cw_ccw_2::ccw(i);}
  static int  cw(int i) {return Triangulation_cw_ccw_2::cw(i);}

  bool get_flag() const
  {
    return _is_finite_non_hyperbolic;
  }
  
  void set_flag(bool flag)
  {
    _is_finite_non_hyperbolic = flag;
  }
  
  // Supposed to be called before "get_non_hyperbolic_edge"
  // iiordanov: leaving here for documentation reasons
  // bool has_non_hyperbolic_edge() const
  // {
  //   return _non_hyperbolic_edge <= 2;
  // }
  
  // Higly recommended to call "has_non_hyperbolic_edge" before 
  // iiordanov: why? the result of "has_non_hyperbolic_edge" is in the second precondition
  unsigned char get_char() const
  {
    CGAL_triangulation_precondition(_is_finite_non_hyperbolic);
    CGAL_triangulation_precondition(_non_hyperbolic_edge <= 2);
    
    return _non_hyperbolic_edge;
  }
  
  void set_char(unsigned char uschar)
  {
    CGAL_triangulation_precondition(_is_finite_non_hyperbolic);
    CGAL_triangulation_precondition(uschar <= 2); 
    
    _non_hyperbolic_edge = uschar;
  }

#ifndef CGAL_NO_DEPRECATED_CODE
  Vertex_handle mirror_vertex(int i) const;
  int mirror_index(int i) const;
#endif

private:
  // a finite face is non_hyperbolic if its circumscribing circle intersects the circle at infinity
  bool _is_finite_non_hyperbolic;
  
  // defined only if the face is finite and non_hyperbolic
  unsigned char _non_hyperbolic_edge;
};

#ifndef CGAL_NO_DEPRECATED_CODE
template < class Gt, class Fb >
inline
typename Hyperbolic_triangulation_face_base_2<Gt,Fb>::Vertex_handle
Hyperbolic_triangulation_face_base_2<Gt,Fb>::
mirror_vertex(int i) const
{
  CGAL_triangulation_precondition ( this->neighbor(i) != Face_handle()
				    && this->dimension() >= 1);
  //return neighbor(i)->vertex(neighbor(i)->index(this->handle()));
  return this->neighbor(i)->vertex(mirror_index(i));
}

template < class Gt, class Fb >
inline int
Hyperbolic_triangulation_face_base_2<Gt,Fb>::
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

} //namespace CGAL 

#endif //CGAL_HYPERBOLIC_TRIANGULATION_FACE_BASE_2_H
