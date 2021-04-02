#include <cassert>

#include <CGAL/Cartesian.h>

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
  Tetrahedron tet(p,q,r,s),
      tet2(p, Point(9,0,0), Point(15, 15, 0), Point(0, 15, 10)),
      v_v(p,Point(-10,0,0),Point(-10,-10,0), Point(0,-10,-10)),
      v_e(Point(-10,0,0), Point(0,0,10), Point(0,0,-10), Point(-10,10,0)),
      v_f(Point(-10,0,0), Point(0,-10,-10), Point(0,-10,10), Point(0,10,0)),
      e_e(Point(-10,0,0), Point(-15,0,0), Point(0,5,10), Point(0,10,5)),
      e_f(Point(0,15,15), Point(-15,0,0), Point(0,-11,10), Point(0,10,-11)),
      f_f(Point(10,10,10), q, r, s),
      tet3(Point(-1,0,0), Point(-10,0,0), Point(-10,-10,0), Point(0,-10,-10));

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

  CGAL_assertion(CGAL::do_intersect(tet, v_e));
  CGAL_assertion(CGAL::do_intersect(tet, v_f));
  CGAL_assertion(CGAL::do_intersect(tet, v_v));
  CGAL_assertion(CGAL::do_intersect(tet, tet2));
  CGAL_assertion(CGAL::do_intersect(tet, e_e));
  CGAL_assertion(CGAL::do_intersect(tet, e_f));
  CGAL_assertion(CGAL::do_intersect(tet, f_f));

  CGAL_assertion(!CGAL::do_intersect(tet, tet3));


  return 0;
}
