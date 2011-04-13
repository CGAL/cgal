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
#include <CGAL/Triangulation_face_base_2.h>

CGAL_BEGIN_NAMESPACE 

template <class Gt>
class Constrained_triangulation_face_base_2
  :  public Triangulation_face_base_2<Gt>
{
public:
  typedef Gt Geom_traits;
  typedef Triangulation_face_base_2<Gt> Fab;
  typedef Constrained_triangulation_face_base_2<Gt> Constrained_face_base;
  typedef typename Gt::Point_2  Point;

protected:
  bool C[3];
 
public:
  Constrained_triangulation_face_base_2()
    : Fab()
  {
    set_constraints(false,false,false);
  }

  Constrained_triangulation_face_base_2(void* v0, void* v1, void* v2)
    : Fab(v0,v1,v2)
  {
    set_constraints(false,false,false);
  }

  Constrained_triangulation_face_base_2(void* v0, void* v1, void* v2,
					void* n0, void* n1, void* n2)
    : Fab(v0,v1,v2,n0,n1,n2)
  {
    set_constraints(false,false,false);
  }


  Constrained_triangulation_face_base_2(void* v0, void* v1, void* v2,
					void* n0, void* n1, void* n2,
					bool c0, bool c1, bool c2 )
    : Fab(v0,v1,v2,n0,n1,n2)
  {
    set_constraints(c0,c1,c2);
  }


  bool is_constrained(int i) const ;
  void set_constraints(bool c0, bool c1, bool c2) ;
  void set_constraint(int i, bool b);
  void reorient();
  void ccw_permute();
  void cw_permute();
  bool is_valid(bool verbose = false, int level = 0) const;
  
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
  Fab::reorient();
  set_constraints(C[1],C[0],C[2]);
}

template <class Gt>
inline void
Constrained_triangulation_face_base_2<Gt>::
ccw_permute()
{
  Fab::ccw_permute();
  set_constraints(C[2],C[0],C[1]);
}

template <class Gt>
inline void
Constrained_triangulation_face_base_2<Gt>::
cw_permute()
{
  Fab::cw_permute();
  set_constraints(C[1],C[2],C[0]);
}
  
template <class Gt>
inline bool
Constrained_triangulation_face_base_2<Gt>::
is_valid(bool verbose, int level) const
{
  bool result = Fab::is_valid(verbose, level);
  CGAL_triangulation_assertion(result);
  if (dimension() == 2) {
    for(int i = 0; i < 3; i++) {
      Constrained_face_base*  n = 
	static_cast<Constrained_face_base*>(neighbor(i));
      if(n != NULL){
	int ni = n->face_index(this);
	result = result && ( is_constrained(i) == n->is_constrained(ni));
      }
    }
  }
  return (result);
}

CGAL_END_NAMESPACE 
  
#endif //CGAL_CONSTRAINED_TRIANGULATION_FACE_BASE_2_H







