#ifndef CGAL_REGULAR_TRIANGULATION_FACE_BASE_2_H
#define CGAL_REGULAR_TRIANGULATION_FACE_BASE_2_H

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h>


template <class Gt>
class CGAL_Regular_triangulation_face_base_2
  :  public CGAL_Triangulation_face_base_2<Gt>
{
public:
  typedef Gt Geom_traits;
  typedef CGAL_Triangulation_face_base_2<Gt> Fb;
  typedef CGAL_Regular_triangulation_face_base_2<Gt> Regular_face_base;
  typedef typename Gt::Point  Point;
  typedef list<Point> Point_list;

protected:
 Point_list  plist;

public:
 CGAL_Regular_triangulation_face_base_2()
    : Fb(), plist()
  {}

  CGAL_Regular_triangulation_face_base_2(void* v0, void* v1, void* v2)
    : Fb(v0,v1,v2),plist() 
  { }

  CGAL_Regular_triangulation_face_base_2(void* v0, void* v1, void* v2,
					     void* n0, void* n1, void* n2)
    : Fb(v0,v1,v2,n0,n1,n2), plist()
  { }

  ~CGAL_Regular_triangulation_face_base_2()
  { 
    plist.clear();
  }

  Point_list& point_list()
  {
    return plist;
  }

};



#endif // CGAL_REGULAR_TRIANGULATION_FACE_BASE_2_H
