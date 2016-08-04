// Copyright (c) 2016  INRIA Nancy (France).
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
// $URL: 
// $Id: 
// 
//
// Author(s)     : Monique Teillaud <Monique.Teillaud@inria.fr>

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

  bool is_finite_non_hyperbolic() const
  {
    return _is_finite_non_hyperbolic;
  }
  
  void set_finite_non_hyperbolic(bool is_finite_non_hyperbolic)
  {
    _is_finite_non_hyperbolic = is_finite_non_hyperbolic;
  }
  
  // Supposed to be called before "get_non_hyperbolic_edge"
  bool has_non_hyperbolic_edge() const
  {
    return _non_hyperbolic_edge <= 2;
  }
  
  // Higly recommended to call "has_non_hyperbolic_edge" before 
  unsigned char get_non_hyperbolic_edge() const
  {
    assert(_is_finite_non_hyperbolic);
    assert(_non_hyperbolic_edge <= 2);
    
    return _non_hyperbolic_edge;
  }
  
  void set_non_hyperbolic_edge(unsigned char non_hyperbolic_edge)
  {
    assert(_is_finite_non_hyperbolic);
    assert(non_hyperbolic_edge <= 2); 
    
    _non_hyperbolic_edge = non_hyperbolic_edge;
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
