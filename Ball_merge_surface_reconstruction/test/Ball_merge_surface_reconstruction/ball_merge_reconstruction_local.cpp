#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Ball_merge_surface_reconstruction.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/IO/read_points.h>

#include <string>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;

namespace params=CGAL::parameters;

int main(int argc, char **argv)
{
  // reading input points and parameters
  std::string inFilename=argc>1?argv[1]:CGAL::data_file_path("points_3/hippo2.ply"); //Test input file
  double delta = argc>2?atof(argv[2]):1.8;
  double eta= argc>3?atof(argv[3]):40.;
  std::vector<Point> points;
  CGAL::IO::read_points(inFilename, std::back_inserter(points)); //Reading the input points

  if (points.empty())
  {
    std::cerr << inFilename << " cannot be read correctly or is empty.\n";
    return 1;
  }

  //vector for storing output triangle indices
  std::vector<std::array<int,3>> meshFaceIndices;

  // run the reconstruction with the given parameters
  CGAL::ball_merge_surface_reconstruction_local(
    points, meshFaceIndices, params::delta(delta)
                                    .eta(eta)
                                    .concurrency_tag(CGAL::Parallel_if_available_tag()));

  // write output triangle soups
  CGAL::IO::write_polygon_soup("BMOut1.ply", points, meshFaceIndices);

  return 0;
}
