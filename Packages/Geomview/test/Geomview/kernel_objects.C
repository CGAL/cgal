// Test program for Geomview_stream with kernel objects.  It's mostly a
// compilation test, I do not verify the output.
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

typedef CGAL::Cartesian<double> Kernel;

int main()
{
  CGAL::Geomview_stream gv(CGAL::Bbox_3(0, 0, 0, 350, 350, 350));

  Kernel::Point_3 P3(200, 100, 100);
  Kernel::Sphere_3 S3(Kernel::Point_3(100, 100, 100), 1000);

  gv << S3;
  gv << P3;

  sleep(10);

  return 0;
}
