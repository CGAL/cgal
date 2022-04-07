#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h>

#include <iostream>
#include <sstream>
#include <fstream>

typedef CGAL::Exact_predicates_exact_constructions_kernel Exact_kernel;
typedef CGAL::Polyhedron_3<Exact_kernel> Polyhedron;
typedef CGAL::Surface_mesh<Exact_kernel::Point_3> Surface_mesh;
typedef CGAL::Nef_polyhedron_3<Exact_kernel> Nef_polyhedron;


void fill_cube_1(Polyhedron& poly)
{
  std::string input =
"OFF\n\
8 12 0\n\
-1 -1 -1\n\
-1 1 -1\n\
1 1 -1\n\
1 -1 -1\n\
-1 -1 1\n\
-1 1 1\n\
1 1 1\n\
1 -1 1\n\
3  0 1 3\n\
3  3 1 2\n\
3  0 4 1\n\
3  1 4 5\n\
3  3 2 7\n\
3  7 2 6\n\
3  4 0 3\n\
3  7 4 3\n\
3  6 4 7\n\
3  6 5 4\n\
3  1 5 6\n\
3  2 1 6";

  std::stringstream ss;
  ss << input;
  ss >> poly;
}

void fill_cube_2(Polyhedron& poly)
{
  std::string input =
"OFF\n\
8 12 0\n\
-0.5 -0.5 -0.5\n\
-0.5 0.5 -0.5\n\
0.5 0.5 -0.5\n\
0.5 -0.5 -0.5\n\
-0.5 -0.5 0.5\n\
-0.5 0.5 0.5\n\
0.5 0.5 0.5\n\
0.5 -0.5 0.5\n\
3  0 1 3\n\
3  3 1 2\n\
3  0 4 1\n\
3  1 4 5\n\
3  3 2 7\n\
3  7 2 6\n\
3  4 0 3\n\
3  7 4 3\n\
3  6 4 7\n\
3  6 5 4\n\
3  1 5 6\n\
3  2 1 6";

  std::stringstream ss;
  ss << input;
  ss >> poly;
}

int main()
{
  // create two nested cubes
  Polyhedron cube1, cube2;
  fill_cube_1(cube1);
  fill_cube_2(cube2);

  Nef_polyhedron nef1(cube1);
  Nef_polyhedron nef2(cube2);

  // compute the difference of the nested cubes
  Nef_polyhedron nef=nef1-nef2;

  // output the result into a Surface_mesh
  Surface_mesh output;
  CGAL::convert_nef_polyhedron_to_polygon_mesh(nef, output);
  std::ofstream out;
  out.open("out.off");
  out << output;
  out.close();
}
