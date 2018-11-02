#include <cassert>

#include <CGAL/Cartesian.h>
#include <CGAL/intersections.h>

typedef CGAL::Cartesian<double> K;

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
  Tetrahedron tet;
  Sphere sp;

  CGAL::do_intersect(tet,Triangle());
  CGAL::do_intersect(tet,Segment());
  CGAL::do_intersect(tet,Iso_cuboid());
  CGAL::do_intersect(tet,Sphere());
  CGAL::do_intersect(tet,Plane());
  CGAL::do_intersect(tet,Line());
  CGAL::do_intersect(tet,Ray());
  CGAL::do_intersect(tet,tet);
  CGAL::do_intersect(tet,Bbox());
  CGAL::do_intersect(sp, Line());
  CGAL::do_intersect(sp, Ray());
  CGAL::do_intersect(sp, Segment());

  CGAL::do_intersect(Triangle(), tet);
  CGAL::do_intersect(Segment(), tet);
  CGAL::do_intersect(Iso_cuboid(), tet);
  CGAL::do_intersect(Sphere(), tet);
  CGAL::do_intersect(Plane(), tet);
  CGAL::do_intersect(Line(), tet);
  CGAL::do_intersect(Ray(), tet);
  CGAL::do_intersect(Bbox(), tet);

  CGAL::do_intersect(Line(), sp);
  CGAL::do_intersect(Ray(), sp);
  CGAL::do_intersect(Segment(), sp);


  return 0;
}
