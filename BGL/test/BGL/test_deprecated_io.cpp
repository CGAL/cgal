#include <CGAL/Installation/internal/disable_deprecation_warnings_and_errors.h>

#include <fstream>
#include <iostream>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/generators.h>

#include <CGAL/boost/graph/IO/OFF.h>
#include <CGAL/boost/graph/IO/VTK.h>
#include <CGAL/boost/graph/IO/WRL.h>


typedef CGAL::Simple_cartesian<double>         Kernel;
typedef Kernel::Point_3                        Point_3;
typedef CGAL::Surface_mesh<Point_3>            SM;

int main()
{
  // OFF
  SM sm_in, sm_out;
  Point_3 p0(0,0,0), p1(1,0,0), p2(0,1,0);
  CGAL::make_triangle(p0, p1, p2, sm_out);
  bool ok = CGAL::write_off("tmp_deprecated.off", sm_out);
  assert(ok);
  ok = CGAL::read_off("tmp_deprecated.off", sm_in);
  assert(ok);
  assert(num_vertices(sm_in) == 3 && num_faces(sm_in) == 1);
  sm_in.clear();

  std::ofstream os("tmp_deprecated.off");
  ok = CGAL::write_off(os, sm_out);
  assert(ok);
  os.close();
  std::ifstream is("tmp_deprecated.off");
  ok = CGAL::read_off(is, sm_in);
  assert(ok);
  assert(num_vertices(sm_in) == 3 && num_faces(sm_in) == 1);
  is.close();
  sm_in.clear();
#ifdef CGAL_USE_VTK
  //vtk
  os.open("tmp_deprecated.vtp");
  ok = CGAL::write_vtp(os, sm_out);
  assert(ok);
  os.close();

  ok = CGAL::IO::read_VTP("tmp_deprecated.vtp", sm_in);
  assert(ok);
  assert(num_vertices(sm_in) == 3 && num_faces(sm_in) == 1);
  sm_in.clear();
#endif
  //wrl
  os.open("tmp_deprecated.wrl");
  ok = CGAL::write_wrl(os, sm_out);
  assert(ok);
  os.close();
  return EXIT_SUCCESS;
}
