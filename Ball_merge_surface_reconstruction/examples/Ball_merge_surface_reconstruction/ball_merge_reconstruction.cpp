#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Ball_merge_surface_reconstruction.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/IO/read_points.h>

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;



int main(int argc, char **argv)
{
  double tlen;
  const char *inFilename = argv[1];//Filename
  std::ifstream inStream(inFilename);//Read the file
  std::vector<std::array<unsigned char, 3>> meshVertexColors;//Variable to hold if the input has textures
  Point p;
  double par = atof(argv[2]);//Parameter to check IR
  int option = atoi(argv[3]);//Option to select global and local variants - 1 for global and 0 for local
  if (argc >= 5)//If local, there is an option to give extra optional parameter to filter long triangles
    tlen = atof(argv[4]);
  else
    tlen = 200.0;//If not provided, default 200 will be taken

  std::vector<Point> points;
  CGAL::IO::read_points(argv[1], std::back_inserter(points));

  CGAL::Ball_merge_surface_reconstruction<K, CGAL::Parallel_tag> bmsr;
  bmsr.option=option;
  bmsr(points, par, tlen);

  // AF: In case of duplicated points the colors of vertices will be wrong

  std::vector<std::vector<int>> meshFaceIndices;
  bmsr.set_triangle_indices_hull1(meshFaceIndices);
  std::string st = "BMOut1.ply";

  CGAL::IO::write_polygon_soup(st, points, meshFaceIndices);
  meshFaceIndices.clear();

  if (option == 1){//Sometimes, in the gloabl case, the largest group would be a mould created by the convex hull, just to avoid it, we will write the second largest group as well
    bmsr.set_triangle_indices_hull2(meshFaceIndices);
    st="BMOut2.ply";
    CGAL::IO::write_polygon_soup(st, points, meshFaceIndices);
  }

  return 0;
}
