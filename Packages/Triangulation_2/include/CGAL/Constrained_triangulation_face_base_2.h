#ifndef CGAL_CONSTRAINED_TRIANGULATION_FACE_BASE_2_H
#define CGAL_CONSTRAINED_TRIANGULATION_FACE_BASE_2_H

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h


template <class Gt>
class CGAL_Constrained_triangulation_face_base_2
  :  public CGAL_triangulation_face_base_2<Gt>
{
protected:
    bool C[3];

public:
  typedef Gt Geom_traits;
  typedef CGAL_Triangulation_face_base_2<Gt> Fb;
  typedef CGAL_Constrained_triangulation_face_base_2<Gt> Constrained_face_base;

  CGAL_Constrained_triangulation_face_base_2()
    : Fb()
  {
    set_constrained(false,false,false);
  }

  CGAL_Constrained_triangulation_face_base_2(void* v0, void* v1, void* v2)
    : Fb(v0,v1,v2)
  {
    set_constrained(false,false,false);
  }

  CGAL_Constrained_triangulation_face_base_2(void* v0, void* v1, void* v2
					     void* n0, void* n1, void* n2)
    : Fb(v0,v1,v2,n0,n1,n2)
  {
    set_constrained(false,false,false);
  }


  CGAL_Constrained_triangulation_face_base_2(void* v0, void* v1, void* v2
					     void* n0, void* n1, void* n2
					     bool c0, bool c1, bool c2 )
    : Fb(v0,v1,v2,n0,n1,n2)
  {
    set_constrained(c0,c1,c2);
  }

  void set_constrained(bool c0, bool c1, bool c2)
  {
    C[0]=c0;
    C[1]=c1;
    C[2]=c2;
  }

  void set_constrained(int i, bool b)
  {
    CGAL_triangulation_precondition( i == 0 || i == 1 || i == 2);
    C[i] = b;
  }
    
  bool is_constrained(int i)
  {
    return(C[i]);
  }
  
  bool is_valid(bool verbose = false, int level = 0)
  {
    bool result = Fb::is_valid();
    CGAL_triangulation_assertion(result);
    for(int i = 0; i < 3; i++) {
      Face_handle n = neighbor(i);
      if(n != NULL){
	int ni = n->index(handle());
	result = result && ( is_constrained(i) == n->is_constrained(ni));
      }
    }
    return (result);
  }
  
#endif CGAL_CONSTRAINED_TRIANGULATION_FACE_BASE_2_H







