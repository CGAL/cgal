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

typedef CGAL::Cartesian<double> K;

int main()
{
  CGAL::Geomview_stream gv(CGAL::Bbox_3(0, 0, 0, 350, 350, 350));

  K::Point_3 P3(200, 100, 100);
  K::Point_2 P2(200, 100);
  K::Segment_2 S2(K::Point_2(200, 100),
                  K::Point_2(300, 100));
  K::Segment_3 S3(K::Point_3(200, 100, 100),
                  K::Point_3(300, 100, 200));
  K::Sphere_3 Sph3(K::Point_3(100, 100, 100), 1000);
  K::Triangle_2 Tr2(K::Point_2(200, 200),
                    K::Point_2(220, 220),
                    K::Point_2(180, 220));
  K::Triangle_3 Tr3(K::Point_3(200, 200, 50),
                    K::Point_3(220, 220, 80),
                    K::Point_3(180, 220, 100));
  K::Tetrahedron_3 T3(K::Point_3(100, 100, 180),
                      K::Point_3(120,  70, 220),
                      K::Point_3(100, 100, 220),
                      K::Point_3(120, 150, 250));

  gv << Sph3;
  gv << P3;
  gv << P2;
  gv << S2;
  gv << S3;
  gv << Tr2;
  gv << Tr3;
  gv << T3;

  sleep(30);

  return 0;
}
