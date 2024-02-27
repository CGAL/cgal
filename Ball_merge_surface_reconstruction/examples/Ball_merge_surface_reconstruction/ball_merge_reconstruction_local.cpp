#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Ball_merge_surface_reconstruction.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/IO/read_points.h>

#include <string>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;

int main(int argc, char **argv)
{
  std::string inFilename=argc>1?argv[1]:CGAL::data_file_path("points_3/hippo2.ply"); //Test input file
  std::vector<Point> points;
  CGAL::IO::read_points(inFilename, std::back_inserter(points)); //Reading the input points

  std::vector<std::array<int,3>> meshFaceIndices;
  CGAL::ball_merge_surface_reconstruction_local<CGAL::Parallel_if_available_tag>(points, meshFaceIndices, 1.8, 40);
  CGAL::IO::write_polygon_soup("BMOut1.ply", points, meshFaceIndices);

  return 0;
}
