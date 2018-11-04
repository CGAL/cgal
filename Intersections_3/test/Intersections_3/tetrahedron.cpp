#include <cassert>

#include <CGAL/Cartesian.h>
#include <CGAL/intersections.h>

typedef CGAL::Cartesian<double> K;

typedef K::Point_3         Point;
typedef K::Tetrahedron_3   Tetrahedron;

typedef K::Segment_3       Segment;
typedef K::Triangle_3      Triangle;
typedef K::Iso_cuboid_3    Iso_cuboid;
typedef K::Sphere_3        Sphere;

typedef K::Plane_3         Plane;
typedef K::Line_3          Line;
typedef K::Ray_3           Ray;

typedef CGAL::Bbox_3       Bbox;

int main()
{
  Point p(0,0,0), q(10,0,0), r(10,10,0), s(0, 10,10);
  Point p2(1,1,1), q2(20,20,20), r2(0,0,20);
  Tetrahedron tet(p,q,r,s);
  Sphere sp(p2,1.0);
  
  CGAL::do_intersect(tet,Triangle(p2,q2,r2));
  CGAL::do_intersect(tet,Segment(p2,q2));
  CGAL::do_intersect(tet,Iso_cuboid(p2,q2));
  CGAL::do_intersect(tet,sp);
  CGAL::do_intersect(tet,Plane(p2,q2,r2));
  CGAL::do_intersect(tet,Line(p2,q2));
  CGAL::do_intersect(tet,Ray(p2,q2));  
  CGAL::do_intersect(tet,tet);
  CGAL::do_intersect(tet,sp.bbox());
  CGAL::do_intersect(sp, Line(p2,q2));
  CGAL::do_intersect(sp, Ray(p2,q2));
  CGAL::do_intersect(sp, Segment(p2,q2));
  
  
  CGAL::do_intersect(Triangle(p2,q2,r2), tet);
  CGAL::do_intersect(Segment(p2,q2), tet);
  CGAL::do_intersect(Iso_cuboid(p2,q2), tet);
  CGAL::do_intersect(sp, tet);
  CGAL::do_intersect(Plane(p2,q2,r2), tet);
  CGAL::do_intersect(Line(p2,q2), tet);
  CGAL::do_intersect(Ray(p2,q2), tet);  
  CGAL::do_intersect(sp.bbox(), tet);
  
  CGAL::do_intersect(Line(p2,q2), sp);
  CGAL::do_intersect(Ray(p2,q2), sp);
  CGAL::do_intersect(Segment(p2,q2), sp);

  return 0;
}
