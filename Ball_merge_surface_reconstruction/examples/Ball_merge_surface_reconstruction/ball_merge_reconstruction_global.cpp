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
  const std::string inFilename=argc>1?argv[1]:CGAL::data_file_path("points_3/kitten.xyz");
  double delta = argc>2?atof(argv[2]):1.7;
  std::vector<Point> points;
  CGAL::IO::read_points(inFilename, std::back_inserter(points));

  if (points.empty())
  {
    std::cerr << inFilename << " cannot be read correctly or is empty.\n";
    return 1;
  }

  //vectors for storing output triangle indices
  std::vector<std::array<int,3>> meshFaceIndices1, meshFaceIndices2;

  // run the reconstruction with the given parameter
  CGAL::ball_merge_surface_reconstruction_global(
    points,
    meshFaceIndices1, meshFaceIndices2,
    params::delta(delta).concurrency_tag(CGAL::Parallel_if_available_tag()));

  // write output triangle soups
  CGAL::IO::write_polygon_soup("BMOut1.ply", points, meshFaceIndices1); //The first resulting mesh
  CGAL::IO::write_polygon_soup("BMOut2.ply", points, meshFaceIndices2); //The second resulting mesh

  return 0;
}
