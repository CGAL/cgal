#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/OBJ.h>
#include <CGAL/IO/OFF.h>
#include <CGAL/IO/PLY.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/IO/STL.h>
#include <CGAL/IO/XYZ/read_xyz_points.h>
#include <CGAL/Surface_mesh.h>

#include <cstddef>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 P;
typedef std::vector<std::vector<std::size_t> > Faces;

int return_code = 0;

template <class F>
void run(const char* label, const std::string& content, F f)
{
  std::istringstream is(content);
  std::vector<P> pts; Faces faces;
  bool ok = f(is, pts, faces);
  std::cout << label << ": ok=" << ok << " npts=" << pts.size() << " nfaces=" << faces.size() << "\n";
  if(!ok) {
    std::cerr << "ERROR: Failed to read " << label << "\n";
    return_code = 1;
  }
}

int main()
{
  std::cout << std::unitbuf;

  auto xyz = [](std::istream& is, std::vector<P>& pts, Faces&) {
    return CGAL::IO::read_XYZ(is, std::back_inserter(pts)); };
  auto obj = [](std::istream& is, std::vector<P>& pts, Faces& f) {
    return CGAL::IO::read_OBJ(is, pts, f); };
  auto stl = [](std::istream& is, std::vector<P>& pts, Faces& f) {
    return CGAL::IO::read_STL(is, pts, f); };
  auto ply = [](std::istream& is, std::vector<P>& pts, Faces& f) {
    return CGAL::IO::read_PLY(is, pts, f); };
  auto off = [](std::istream& is, std::vector<P>& pts, Faces& f) {
    return CGAL::IO::read_OFF(is, pts, f); };

  run("XYZ plain           ", "0 0 0\n1 0 0\n0 1 0", xyz);
  run("XYZ leading +       ", "+0 0 0\n1 0 0\n0 1 0\n", xyz);
  run("OBJ plain           ", "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n", obj);
  run("OBJ leading +       ", "v +0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n", obj);
  run("STL ascii plain     ", "solid s\nfacet normal 0 0 1\nouter loop\nvertex 0 0 0\nvertex 1 0 0\nvertex 0 1 0\nendloop\nendfacet\nendsolid s\n", stl);
  run("STL ascii leading + ", "solid s\nfacet normal +0 +0 +1\nouter loop\nvertex +0 0 0\nvertex 1 0 0\nvertex 0 1 0\nendloop\nendfacet\nendsolid s\n", stl);
  run("PLY ascii plain     ", "ply\nformat ascii 1.0\nelement vertex 3\nproperty double x\nproperty double y\nproperty double z\nelement face 1\nproperty list uchar int vertex_indices\nend_header\n0 0 0\n1 0 0\n0 1 0\n3 0 1 2\n", ply);
  run("PLY ascii leading + ", "ply\nformat ascii 1.0\nelement vertex 3\nproperty double x\nproperty double y\nproperty double z\nelement face 1\nproperty list uchar int vertex_indices\nend_header\n+0 0 0\n1 0 0\n0 1 0\n3 0 1 2\n", ply);
  run("OFF plain           ", "OFF\n3 1 0\n0 0 0\n1 0 0\n0 1 0\n3 0 1 2\n", off);
  run("OFF 1e+00 exponent  ", "OFF\n3 1 0\n0e+00 0 0\n1 0 0\n0 1 0\n3 0 1 2\n", off);
  run("OFF inline comment  ", "OFF\n3 1 0\n0 0 0 # c\n1 0 0\n0 1 0\n3 0 1 2\n", off);
  return return_code;
}
