#include <CGAL/Cartesian.h>

#include <iostream>
#include <fstream>
#include <list>
#include <vector>

#include <CGAL/Bbox_3.h>

#include <CGAL/Halfedge_data_structure_polyhedron_default_3.h>
#include <CGAL/Polyhedron_default_traits_3.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/IO/Postscript_stream_3.h>
#include "PS_Stream_3.C"
#include "PS_Stream.C"

typedef CGAL::Bbox_3            PS_BBox3;
typedef CGAL::Cartesian<double> R;
typedef R::Direction_3          Direction;
typedef R::Point_3              Point3;
typedef R::Vector_3             Vector3;
typedef CGAL::Halfedge_data_structure_polyhedron_default_3<R> HDS;
typedef CGAL::Polyhedron_default_traits_3<R> Traits;
typedef CGAL::Polyhedron_3<Traits,HDS> Polyhedron;

int main()
{
  PS_BBox3 bb3(-1,-1,-1,1,1,1);
  Direction dir(-2,-4,-2);
  // Point3 a2(2,0,0);
//   Point3 b2(1,1,0);
//   Point3 c2(1,0,1);
//   Point3 d2(2,1,1);

  Point3 a2(1,0,0);
  Point3 b2(0,1,0);
  Point3 c2(0,0,1);
  Point3 d2(0,0,0);

// Point3 a(-2,0,0);
// Point3 b(-1,0,0);
// Point3 c(-1,1,0);
// Point3 d(-2,1,0);
// Point3 e(-2,1,-1);
// Point3 f(-1,1,-1);
// Point3 g(-1,0,-1);
// Point3 h(-2,0,-1);

  CGAL::PS_Stream_3 ps(bb3,dir,300,"tetra1.ps",CGAL::PS_Stream::READABLE_EPS);

  Polyhedron P2;

// P.make_triangle(a,b,c);
// P.make_triangle(c,a,d);
// P.make_triangle(d,a,h);
// P.make_triangle(h,d,e);
// P.make_triangle(e,d,c);
// P.make_triangle(c,e,f);
// P.make_triangle(f,e,h);
// P.make_triangle(h,f,g);
// P.make_triangle(g,f,c);
// P.make_triangle(c,g,b);
// P.make_triangle(b,g,h);
// P.make_triangle(h,b,a);
 
  P2.make_tetrahedron(a2,b2,c2,d2);

  ps <<P2;
 
  return 0;
}
