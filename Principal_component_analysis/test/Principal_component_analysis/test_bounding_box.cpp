// Test program for the bounding_box() function.
// Sylvain Pion.

#include <vector>
#include <cassert>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/bounding_box.h>

typedef double                      FT;
typedef CGAL::Simple_cartesian<FT>  K;
typedef K::Point_2                  Point_2;
typedef K::Point_3                  Point_3;
typedef K::Segment_2                Segment_2;
typedef K::Triangle_3               Triangle_3;

void test_2()
{
  // Test bounding_box() for Point_2
  Point_2 p0 (1, 2);
  Point_2 p1 (2, 1);

  std::vector<Point_2> pts;
  pts.push_back(p0);
  pts.push_back(p0);
  pts.push_back(p1);

  assert( CGAL::bounding_box(pts.begin(), pts.end())
          == K::Iso_rectangle_2(p0, p1));

#if 0 // Does not work yet
  // Test bounding_box() for Segment_2
  Segment_2 s0 (p0, p1);
  Segment_2 s1 (p0, p0);

  std::vector<Segment_2> segs;
  segs.push_back(s0);
  segs.push_back(s0);
  segs.push_back(s1);

  assert( CGAL::bounding_box(segs.begin(), segs.end())
          == K::Iso_rectangle_2(p0, p1));
#endif
}

void test_3()
{
  // Test bounding_box() for Point_3
  Point_3 p0 (1, 2, 5);
  Point_3 p1 (2, 1, 6);

  std::vector<Point_3> pts;
  pts.push_back(p0);
  pts.push_back(p1);
  pts.push_back(p1);

  assert( CGAL::bounding_box(pts.begin(), pts.end())
          == K::Iso_cuboid_3(p0, p1));

#if 0 // Does not work yet
  // Test bounding_box() for Triangle_3
  Triangle_3 t0 (p0, p1);
  Triangle_3 t1 (p0, p0);

  std::vector<Triangle_3> trs;
  trs.push_back(t0);
  trs.push_back(t0);
  trs.push_back(t1);

  assert( CGAL::bounding_box(trs.begin(), trs.end())
          == K::Iso_cuboid_3(p0, p1));
#endif
}


int main()
{
  test_2();
  test_3();
  return 0;
}
