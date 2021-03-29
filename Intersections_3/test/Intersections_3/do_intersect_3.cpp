#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Triangle_3 Triangle_3;
typedef K::Segment_3 Segment_3;
typedef K::Ray_3 Ray_3;
typedef K::Line_3 Line_3;
typedef K::Tetrahedron_3 Tetrahedron_3;
typedef K::Point_3 Point_3;
typedef K::Plane_3 Plane_3;
typedef K::Iso_cuboid_3 Iso_cuboid_3;
typedef K::Sphere_3 Sphere_3;

int main()
{
  Point_3 p(0,0,0), q(1,1,1), r(1,0,0), s(1,1,0);
  Segment_3 s0(p,q),s1(p,q);
  Ray_3 ray(p,q);
  Line_3 lin(p,q);
  Triangle_3 t0(p,q,r), t1(p,q,r);
  Tetrahedron_3 tet(p,q,r,s);
  Plane_3 pla(p,q,r);
  Iso_cuboid_3 ic(p,q);
  Sphere_3 sp(p,1);
  CGAL::Bbox_3 bb;

  do_intersect(p,sp);
  do_intersect(bb,t0);
  do_intersect(t0,t1);
  do_intersect(t0,s0);
  do_intersect(t0,ray);
  do_intersect(t0,lin);

  do_intersect(s0,s1);
  do_intersect(pla,s0);
  do_intersect(pla,lin);
  do_intersect(pla,ray);
  do_intersect(pla,pla);
  do_intersect(pla,pla,pla);
  do_intersect(tet,p);
  do_intersect(ic,s0);
  do_intersect(ic,ray);
  do_intersect(ic,ic);
  do_intersect(ic,lin);
  do_intersect(ic,pla);
  do_intersect(ic,ic);
  do_intersect(sp,sp);


  {
    Plane_3 p0(Point_3(0,0,0), Point_3(1,0,0), Point_3(1,1,0));
    Plane_3 p1(Point_3(0,0,1), Point_3(1,0,1), Point_3(1,1,1));
    Plane_3 p2(Point_3(0,0,2), Point_3(1,0,2), Point_3(1,1,2));

    assert(do_intersect(p0,p1,p2) == false); // three parallel
  }

  {
    Plane_3 p0(Point_3(0,0,0), Point_3(1,0,0), Point_3(1,1,0));
    Plane_3 p2(Point_3(2,2,0), Point_3(4,2,0), Point_3(4,4,0));
    Plane_3 p1(Point_3(10,10,0), Point_3(10,20,0), Point_3(20,20,0));

    assert(do_intersect(p0,p1,p2) == true); // three equal
  }

  {
    Plane_3 p0(Point_3(0,0,0), Point_3(1,0,0), Point_3(1,1,0));
    Plane_3 p2(Point_3(2,2,0), Point_3(4,2,0), Point_3(4,4,0));
    Plane_3 p1(Point_3(10,10,1), Point_3(10,20,1), Point_3(20,20,1));

    assert(do_intersect(p0,p1,p2) == false); // two equal, one parallel
  }

  {
    Plane_3 p0(Point_3(0,0,0), Point_3(0,1,0), Point_3(0,0,1));

    Plane_3 p2(Point_3(2,2,0), Point_3(4,2,0), Point_3(4,4,0));
    Plane_3 p1(Point_3(10,10,1), Point_3(10,20,1), Point_3(20,20,1));

    assert(do_intersect(p0,p1,p2) == false); // two parallel, one intersecting both
  }

  {
    Plane_3 p0(Point_3(0,0,0), Point_3(0,1,0), Point_3(0,0,1));
    Plane_3 p1(Point_3(0,0,0), Point_3(1,0,0), Point_3(1,1,0));
    Plane_3 p2(Point_3(2,2,0), Point_3(4,2,0), Point_3(4,4,0));
    assert(do_intersect(p0,p1,p2) == true); // two equal, one intersecting
  }
  {
    Plane_3 p0(Point_3(0,0,0), Point_3(1,1,1), Point_3(0,0,1));
    Plane_3 p1(Point_3(2,2,2), Point_3(6,6,6), Point_3(0,1,0));
    Plane_3 p2(Point_3(3,3,3), Point_3(5,5,5), Point_3(1,0,0));

    assert(do_intersect(p0,p1,p2) == true); // three in a line
  }


  {
    Plane_3 p0(Point_3(0,0,0), Point_3(1,0,0), Point_3(0,1,0)); // xy
    Plane_3 p1(Point_3(0,0,0), Point_3(0,0,1), Point_3(0,1,1)); // yz
    Plane_3 p2(Point_3(0,0,10), Point_3(0,1,10), Point_3(1,0,0)); //

    assert(do_intersect(p0,p1,p2) == false); // pairwise intersections, i.e. 3 lines
  }

  {
    Plane_3 p0(Point_3(0,0,0), Point_3(1,0,0), Point_3(0,1,0));
    Plane_3 p1(Point_3(0,0,0), Point_3(1,1,1), Point_3(1,1,0));
    Plane_3 p2(Point_3(1,0,1), Point_3(1,1,5), Point_3(0,1,1));

    assert(do_intersect(p0,p1,p2) == true); // three in a point
  }
  return 0;
}
