#include <CGAL/internal/disable_deprecation_warnings_and_errors.h>

#include <CGAL/Simple_cartesian.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO/LAS.h>
#include <CGAL/Point_set_3/IO/OFF.h>
#include <CGAL/Point_set_3/IO/PLY.h>

#include <iostream>
#include <fstream>


typedef CGAL::Simple_cartesian<double>         Kernel;
typedef Kernel::Point_3                        Point_3;
typedef Kernel::Vector_3 Vector_3;

int main()
{
  CGAL::Point_set_3<Point_3, Vector_3> ps, ps2;
  std::ifstream is("data/oni.pwn");
  std::ofstream os;

  if(!CGAL::read_xyz_point_set(is, ps))
  {
    std::cerr<<"Error while reading input."<<std::endl;
    return EXIT_FAILURE;
  }
  is.close();
  bool ok = false;
#ifdef CGAL_LINKED_WITH_LASLIB
  //LAS
  os.open("tmp.las", std::ios::binary);
  ok = CGAL::write_las_point_set(os, ps);
  os.close();
  assert (ok);

  is.open("tmp.las", std::ios::binary);
  ok = CGAL::read_las_point_set(is, ps2);
  is.close();
  assert(ok);
  ps2.clear();
#endif

  //OFF
  os.open("tmp.off", std::ios::binary);
  ok = CGAL::write_off_point_set(os, ps);
  os.close();
  assert (ok);

  is.open("tmp.off", std::ios::binary);
  ok = CGAL::read_off_point_set(is, ps2);
  is.close();
  assert(ok);
  ps2.clear();

  //PLY
  os.open("tmp.ply", std::ios::binary);
  ok = CGAL::write_ply_point_set(os, ps);
  os.close();
  assert (ok);

  is.open("tmp.ply", std::ios::binary);
  ok = CGAL::read_ply_point_set(is, ps2);
  is.close();
  assert(ok);
  ps2.clear();

  //XYZ
  os.open("tmp.xyz", std::ios::binary);
  ok = CGAL::write_xyz_point_set(os, ps);
  os.close();
  assert (ok);

  is.open("tmp.xyz", std::ios::binary);
  ok = CGAL::read_xyz_point_set(is, ps2);
  is.close();
  assert(ok);
  ps2.clear();

  return EXIT_SUCCESS;
}
