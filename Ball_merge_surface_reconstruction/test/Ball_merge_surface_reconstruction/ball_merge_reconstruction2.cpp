#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Ball_merge_surface_reconstruction.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/IO/read_points.h>

// #include <fstream>
// #include <sstream>
// #include <iostream>
// #include <string>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;

int main(int argc, char **argv)
{
  const char *inFilename = CGAL::data_file_path("points_3/kitten.xyz"));//Filename
  std::ifstream inStream(inFilename);//Read the file
  // std::vector<std::array<unsigned char, 3>> meshVertexColors;//Variable to hold if the input has textures
  // double par = atof(argv[2]);//Parameter to check IR
  // int option = atoi(argv[3]);//Option to select global and local variants - 1 for global and 0 for local

  std::vector<Point> points;
  CGAL::IO::read_points(argv[1], std::back_inserter(points));

  std::vector<std::array<int,3>> meshFaceIndices1, meshFaceIndices2;
    CGAL::ball_merge_surface_reconstruction_global<CGAL::Parallel_if_available_tag>(points, meshFaceIndices1, meshFaceIndices2, 1.7);
    CGAL::IO::write_polygon_soup("BMOut1.ply", points, meshFaceIndices1);
    CGAL::IO::write_polygon_soup("BMOut2.ply", points, meshFaceIndices2);

  // if (option == 0)
  // {
  //   double eta = (argc >= 5) ? atof(argv[4]) : 200.;
  //   std::vector<std::array<int,3>> meshFaceIndices;
  //   CGAL::ball_merge_surface_reconstruction_local<CGAL::Parallel_if_available_tag>(points, meshFaceIndices, par, eta);
  //   CGAL::IO::write_polygon_soup("BMOut1.ply", points, meshFaceIndices);
  // }
  // else{
  //   std::vector<std::array<int,3>> meshFaceIndices1, meshFaceIndices2;
  //   CGAL::ball_merge_surface_reconstruction_global<CGAL::Parallel_if_available_tag>(points, meshFaceIndices1, meshFaceIndices2, par);
  //   CGAL::IO::write_polygon_soup("BMOut1.ply", points, meshFaceIndices1);
  //   CGAL::IO::write_polygon_soup("BMOut2.ply", points, meshFaceIndices2);
  // }

  return 0;
}
