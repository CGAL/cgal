// Test program for Geomview_stream with kernel objects.  It's mostly a
// compilation test, I do not verify the output automatically.
//
//  Sylvain Pion, 2000.

#include <CGAL/basic.h>

#include <CGAL/Cartesian.h>

#include <CGAL/Point_2.h>
#include <CGAL/Point_3.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Segment_3.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Sphere_3.h>

#include <CGAL/IO/Geomview_stream.h>

#include <unistd.h> // for sleep()

typedef CGAL::Cartesian<double> K;

void test_parse_point()
{
  const char *test_point="( 123 456 789 1 )";
  K::Point_3 p;
  CGAL::parse_point(test_point, p);
  CGAL_assertion(p == K::Point_3(123, 456, 789));
}

int main()
{
  test_parse_point();

  CGAL::Geomview_stream gv(CGAL::Bbox_3(0, 0, 0, 350, 350, 350));
  gv.set_trace(true);

  gv.clear(); // remove pickplane.
  gv << K::Point_2 (200, 100);
  gv << K::Point_3 (200, 100, 100);
  gv << K::Segment_2 (K::Point_2(200, 100),
                      K::Point_2(300, 100));
  gv << K::Segment_3 (K::Point_3(200, 100, 100),
                      K::Point_3(300, 100, 200));
  gv << K::Sphere_3 (K::Point_3(100, 100, 100), 1000);
  gv << K::Triangle_2 (K::Point_2(200, 200),
                       K::Point_2(220, 220),
                       K::Point_2(180, 220));
  gv << K::Triangle_3 (K::Point_3(200, 200, 50),
                       K::Point_3(220, 220, 80),
                       K::Point_3(180, 220, 100));
  gv << K::Tetrahedron_3 (K::Point_3(100, 100, 180),
                          K::Point_3(120,  70, 220),
                          K::Point_3(100, 100, 220),
                          K::Point_3(120, 150, 250));
  gv << CGAL::Bbox_2(10, 10, 30, 30);
  gv << CGAL::Bbox_3(10, 10, 10, 30, 30, 30);

  gv.look_recenter();
  sleep(10);

  return 0;
}
