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
// file          : Triangulation_face_base_2.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_FACE_BASE_2_H
#define CGAL_TRIANGULATION_FACE_BASE_2_H


#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_utils_2.h>

CGAL_BEGIN_NAMESPACE 

template < class Gt >
class Triangulation_face_base_2 
  : public Triangulation_cw_ccw_2
{
public:
  // typedef typename Gt::Triangle Triangle;
  typedef Triangulation_face_base_2  Face_base;

  inline
  Triangulation_face_base_2()
  {
    set_vertices();
    set_neighbors();
  }

  inline
  Triangulation_face_base_2( void* v0, void* v1, void* v2)
  {
    set_vertices(v0, v1, v2);
    set_neighbors();
  }

  inline
  Triangulation_face_base_2(void* v0, void* v1, void* v2,
				 void* n0, void* n1, void* n2)
  {
    set_vertices(v0, v1, v2);
    set_neighbors(n0, n1, n2);
  }


  inline 
  void* vertex(int i) const
  {
    CGAL_triangulation_precondition( i == 0 || i == 1 || i == 2);
    return V[i];
  } 

 inline 
 bool has_vertex(const void* v) const
  {
    return (V[0] == v) || (V[1] == v) || (V[2]== v);
  }
    
    
  inline 
  bool has_vertex(const void* v, int& i) const
  {
    if (v == V[0]) {
      i = 0;
      return true;
    }
    if (v == V[1]) {
      i = 1;
      return true;
    }
    if (v == V[2]) {
      i = 2;
      return true;
    }
    return false;
  }
    
    
  inline 
  int vertex_index(const void* v) const
  {
    if (v == V[0]) return 0;
    if (v == V[1]) return 1;
    CGAL_triangulation_assertion( v == V[2] );
    return 2;
  }

  inline 
  void* neighbor(int i) const
  {
    CGAL_triangulation_precondition( i == 0 || i == 1 || i == 2);
    return N[i];
  }
    
    
  inline 
  bool has_neighbor(const void* n) const
  {
    return (N[0] == n) || (N[1] == n) || (N[2] == n);
  }
    
    
  inline 
  bool has_neighbor(const void* n, int& i) const
  {
    if(n == N[0]){
      i = 0;
      return true;
    }
    if(n == N[1]){
      i = 1;
      return true;
    }
    if(n == N[2]){
      i = 2;
      return true;
    }
    return false;
  }
    
    
  inline 
  int face_index(const void* n) const
  {
    if (n == N[0]) return 0;
    if (n == N[1]) return 1;
    CGAL_triangulation_assertion( n == N[2] );
    return 2;
  }
    
 

  inline 
  void set_vertex(int i, void* v)
  {
    CGAL_triangulation_precondition( i == 0 || i == 1 || i == 2);
    V[i] = v;
  }
    
    
  inline 
  void set_neighbor(int i, void* n)
  {
    CGAL_triangulation_precondition( i == 0 || i == 1 || i == 2);
    N[i] = n;
  }

  inline 
  void set_vertices()
  {
    V[0] = V[1] = V[2] = NULL;
  }
    
    
  inline 
  void set_vertices(void* v0,
		    void* v1,
		    void* v2)
  {
    V[0] = v0;
    V[1] = v1;
    V[2] = v2;
  }
    
  inline 
  void set_neighbors()
  {
    N[0] = N[1] = N[2] = NULL;
  }
    
  inline
  void set_neighbors(void* n0,
		     void* n1,
		     void* n2)
  {
    N[0] = n0;
    N[1] = n1;
    N[2] = n2;
  }
 
   
 
  //the following trivial is_valid to allow
  // the user of derived face base classes 
  // to add their own purpose checking
  bool is_valid(bool /* verbose */ = false, int /* level */ = 0) const
  {return true;}


private:
  void* V[3];
  void* N[3];
};

CGAL_END_NAMESPACE 

#endif //CGAL_TRIANGULATION_FACE_BASE_2_H
