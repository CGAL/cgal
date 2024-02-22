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
  const char *inFilename = argv[1];//Filename
  std::ifstream inStream(inFilename);//Read the file
  std::vector<std::array<unsigned char, 3>> meshVertexColors;//Variable to hold if the input has textures
  Point p;
  double par = atof(argv[2]);//Parameter to check IR
  int option = atoi(argv[3]);//Option to select global and local variants - 1 for global and 0 for local


  std::vector<Point> points;
  CGAL::IO::read_points(argv[1], std::back_inserter(points));

  if (option == 0)
  {
    double tlen = (argc >= 5) ? atof(argv[4]) : 200.;
    std::vector<std::array<int,3>> meshFaceIndices;
    CGAL::ball_merge_surface_reconstruction_local<CGAL::Parallel_if_available_tag>(points, meshFaceIndices, par, tlen);
    CGAL::IO::write_polygon_soup("BMOut1.ply", points, meshFaceIndices);
  }
  else{
    std::vector<std::array<int,3>> meshFaceIndices1, meshFaceIndices2;
    CGAL::ball_merge_surface_reconstruction_global<CGAL::Parallel_if_available_tag>(points, meshFaceIndices1, meshFaceIndices2, par);
    CGAL::IO::write_polygon_soup("BMOut1.ply", points, meshFaceIndices1);
    CGAL::IO::write_polygon_soup("BMOut2.ply", points, meshFaceIndices2);
  }

  return 0;
}
