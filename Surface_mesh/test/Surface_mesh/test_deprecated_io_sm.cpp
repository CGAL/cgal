#include <CGAL/internal/disable_deprecation_warnings_and_errors.h>
#include <fstream>
#include <iostream>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/generators.h>
#include <CGAL/Surface_mesh/IO/3MF.h>
#include <CGAL/Surface_mesh/IO/OFF.h>
#include <CGAL/Surface_mesh/IO/PLY.h>


typedef CGAL::Simple_cartesian<double>         Kernel;
typedef Kernel::Point_3                        Point_3;
typedef CGAL::Surface_mesh<Point_3>            SM;

int main()
{
  // OFF
  SM sm_in, sm_out;
  Point_3 p0(0,0,0), p1(1,0,0), p2(0,1,0);
  CGAL::make_triangle(p0, p1, p2, sm_out);
  bool ok = CGAL::write_off(sm_out, "tmp.off");
  assert(ok);
  ok = CGAL::read_off(sm_in, "tmp.off");
  assert(ok);
  assert(num_vertices(sm_in) == 3 && num_faces(sm_in) == 1);
  sm_in.clear();

  std::ofstream os("tmp.off");
  ok = CGAL::write_off(os, sm_out);
  assert(ok);
  os.close();
  std::ifstream is("tmp.off");
  ok = CGAL::read_off(is, sm_in);
  assert(ok);
  assert(num_vertices(sm_in) == 3 && num_faces(sm_in) == 1);
  is.close();
  sm_in.clear();

  //PLY
  os.open("tmp.ply");
  std::string comments;
  ok = CGAL::write_ply(os, sm_out, comments);
  assert(ok);
  os.close();
  is.open("tmp.ply");
  ok = CGAL::read_ply(is, sm_in, comments);
  assert(ok);
  assert(num_vertices(sm_in) == 3 && num_faces(sm_in) == 1);
  is.close();
  sm_in.clear();

#ifdef CGAL_LINKED_WITH_3MF
  // 3mf
  std::vector<SM> output_3mf;
  ok = CGAL::read_3mf("test.3mf", output_3mf);
  assert(ok);
  assert(output_3mf.size() == 2);
  sm_in.clear();
#endif

  //others
  ok = CGAL::write_mesh(sm_out, "tmp.off");
  assert(ok);
  ok = CGAL::read_mesh(sm_in, "tmp.ply");
  assert(ok);
  assert(num_vertices(sm_in) == 3 && num_faces(sm_in) == 1);
  sm_in.clear();
  return EXIT_SUCCESS;
}
