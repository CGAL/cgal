#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Ball_merge_surface_reconstruction.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/IO/read_points.h>

#include <string>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;

int main(int argc, char **argv)
{
  // reading input points and parameters
  std::string inFilename = argc>1?argv[1]:CGAL::data_file_path("points_3/kitten.xyz");//Filename
  double delta = argc>2?atof(argv[2]):1.7;

  //Select to select between global and local variants 0 for local and any other value for global
  int option = argc>3?atoi(argv[3]):1;

  std::vector<Point> points;
  CGAL::IO::read_points(inFilename, std::back_inserter(points));

  if (points.empty())
  {
    std::cerr << inFilename << " cannot be read correctly or is empty.\n";
    return 1;
  }

  if (option == 0)
  {
    // reading an extra parameter
    double eta = (argc >= 5) ? atof(argv[4]) : 200.;

    //vector for storing output triangle indices
    std::vector<std::array<int,3>> meshFaceIndices;

    // run the reconstruction with the given parameters
    CGAL::ball_merge_surface_reconstruction_local<CGAL::Parallel_if_available_tag>(points, meshFaceIndices, delta, eta);

    // write output triangle soup
    CGAL::IO::write_polygon_soup("BMOut1.ply", points, meshFaceIndices);
  }
  else{
    //vectors for storing output triangle indices
    std::vector<std::array<int,3>> meshFaceIndices1, meshFaceIndices2;

    // run the reconstruction with the given parameter
    CGAL::ball_merge_surface_reconstruction_global<CGAL::Parallel_if_available_tag>(points, meshFaceIndices1, meshFaceIndices2, delta);

    // write output triangle soups
    CGAL::IO::write_polygon_soup("BMOut1.ply", points, meshFaceIndices1);
    CGAL::IO::write_polygon_soup("BMOut2.ply", points, meshFaceIndices2);
  }

  return 0;
}
