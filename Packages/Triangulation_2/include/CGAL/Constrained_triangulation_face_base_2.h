// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Constrained_triangulation_face_base_2.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_CONSTRAINED_TRIANGULATION_FACE_BASE_2_H
#define CGAL_CONSTRAINED_TRIANGULATION_FACE_BASE_2_H

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_ds_face_base_2.h>

CGAL_BEGIN_NAMESPACE 

template <class Fb = Triangulation_ds_face_base_2<> >
class Constrained_triangulation_face_base_2
  :  public Fb
{
  typedef Fb                                           Base;
  typedef typename Fb::Triangulation_data_structure    TDS;
public:
  typedef TDS                                  Triangulation_data_structure;
  typedef typename TDS::Vertex_handle          Vertex_handle;
  typedef typename TDS::Face_handle            Face_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Fb::template Rebind_TDS<TDS2>::Other    Fb2;
    typedef Constrained_triangulation_face_base_2<Fb2>         Other;
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

template <class Gt>
inline void
Constrained_triangulation_face_base_2<Gt>::
set_constraints(bool c0, bool c1, bool c2)
{
  C[0]=c0;
  C[1]=c1;
  C[2]=c2;
}

template <class Gt>
inline void
Constrained_triangulation_face_base_2<Gt>::
set_constraint(int i, bool b)
{
  CGAL_triangulation_precondition( i == 0 || i == 1 || i == 2);
  C[i] = b;
}
    
template <class Gt>
inline bool
Constrained_triangulation_face_base_2<Gt>::
is_constrained(int i) const
{
  return(C[i]);
}

template <class Gt>
inline void
Constrained_triangulation_face_base_2<Gt>::
reorient()
{
  Base::reorient();
  set_constraints(C[1],C[0],C[2]);
}

template <class Gt>
inline void
Constrained_triangulation_face_base_2<Gt>::
ccw_permute()
{
  Base::ccw_permute();
  set_constraints(C[2],C[0],C[1]);
}

template <class Gt>
inline void
Constrained_triangulation_face_base_2<Gt>::
cw_permute()
{
  Base::cw_permute();
  set_constraints(C[1],C[2],C[0]);
}
  
CGAL_END_NAMESPACE 
  
#endif //CGAL_CONSTRAINED_TRIANGULATION_FACE_BASE_2_H







