#ifndef CGAL_CONSTRAINED_TRIANGULATION_FACE_BASE_2_H
#define CGAL_CONSTRAINED_TRIANGULATION_FACE_BASE_2_H

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <list.h>

template <class Gt>
class CGAL_Constrained_triangulation_face_base_2
  :  public CGAL_Triangulation_face_base_2<Gt>
{
public:
  typedef Gt Geom_traits;
  typedef CGAL_Triangulation_face_base_2<Gt> Fb;
  typedef CGAL_Regular_triangulation_face_base_2<Gt> Regular_face_base;
  typedef typename Gt::Point  Point;
  typedef list<Point> Point_list;

protected;
 Point_list* point_list;

  CGAL_Constrained_triangulation_face_base_2()
    : Fb()
  {
    set_constraints(false,false,false);
  }

  CGAL_Constrained_triangulation_face_base_2(void* v0, void* v1, void* v2)
    : Fb(v0,v1,v2)
  {
    set_constraints(false,false,false);
  }

  CGAL_Constrained_triangulation_face_base_2(void* v0, void* v1, void* v2,
					     void* n0, void* n1, void* n2)
    : Fb(v0,v1,v2,n0,n1,n2)
  {
    set_constraints(false,false,false);
  }


  CGAL_Constrained_triangulation_face_base_2(void* v0, void* v1, void* v2,
					     void* n0, void* n1, void* n2,
					     bool c0, bool c1, bool c2 )
    : Fb(v0,v1,v2,n0,n1,n2)
  {
    set_constraints(c0,c1,c2);
  }

  void set_constraints(bool c0, bool c1, bool c2)
  {
    C[0]=c0;
    C[1]=c1;
    C[2]=c2;
  }

  void set_constraint(int i, bool b)
  {
    CGAL_triangulation_precondition( i == 0 || i == 1 || i == 2);
    C[i] = b;
  }
    
  bool is_constrained(int i) const
  {
    return(C[i]);
  }
  
  bool is_valid(bool verbose = false, int level = 0) const
  {
    bool result = Fb::is_valid();
    CGAL_triangulation_assertion(result);
    for(int i = 0; i < 3; i++) {
      Constrained_face_base*  n = (Constrained_face_base*)neighbor(i);
      if(n != NULL){
	 // The following seems natural, but it may fail if the faces
	// this and n are neighbors on two edges (1-dim triangulation,
	// with infinite faces
	// int ni = n->index(this);//
	int ni = cw(n->vertex_index(vertex(cw(i))));
	result = result && ( is_constrained(i) == n->is_constrained(ni));
      }
    }
    return (result);
  }
};

  
#endif CGAL_CONSTRAINED_TRIANGULATION_FACE_BASE_2_H







