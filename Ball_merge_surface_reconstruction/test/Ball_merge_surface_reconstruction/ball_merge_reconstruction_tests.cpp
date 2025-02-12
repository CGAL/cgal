#include <CGAL/Ball_merge_surface_reconstruction.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/IO/read_points.h>
#include <boost/container/small_vector.hpp>
#if defined(CGAL_LINKED_WITH_TBB)
#include <tbb/concurrent_vector.h>
#endif

#include <string>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;

namespace params = CGAL::parameters;
CGAL::Parallel_if_available_tag ctag;

void test1()
{
  std::string inFilename=CGAL::data_file_path("points_3/half.xyz"); //Test input file
  std::vector<Point> points;
  CGAL::IO::read_points(inFilename, std::back_inserter(points)); //Reading the input points

  std::vector<std::array<int,3>> meshFaceIndices1, meshFaceIndices2; //For saving the results
  CGAL::ball_merge_surface_reconstruction_global(points, meshFaceIndices1, meshFaceIndices2, params::concurrency_tag(ctag)); //Calling global BallMerge with parameter=1.7
  CGAL::IO::write_polygon_soup("BMOut1.ply", points, meshFaceIndices1); //The first resulting mesh
  CGAL::IO::write_polygon_soup("BMOut2.ply", points, meshFaceIndices2); //The second resulting mesh
}

void test2()
{
  const std::string inFilename = CGAL::data_file_path("points_3/kitten.xyz");//Filename
  std::ifstream inStream(inFilename);//Read the file
  // std::vector<std::array<unsigned char, 3>> meshVertexColors;//Variable to hold if the input has textures
  // double par = atof(argv[2]);//Parameter to check IR
  // int option = atoi(argv[3]);//Option to select global and local variants - 1 for global and 0 for local

  std::vector<Point> points;
  CGAL::IO::read_points(inFilename, std::back_inserter(points));

  std::vector<std::array<int,3>> meshFaceIndices1, meshFaceIndices2;
    CGAL::ball_merge_surface_reconstruction_global(points, meshFaceIndices1, meshFaceIndices2, params::concurrency_tag(ctag));
    CGAL::IO::write_polygon_soup("BMOut1.ply", points, meshFaceIndices1);
    CGAL::IO::write_polygon_soup("BMOut2.ply", points, meshFaceIndices2);
}

void test_containers()
{
  const std::string inFilename = CGAL::data_file_path("points_3/kitten.xyz");//Filename
  std::ifstream inStream(inFilename);//Read the file

  std::vector<Point> points;
  CGAL::IO::read_points(inFilename, std::back_inserter(points));

  {
  std::vector<std::vector<std::size_t>> meshFaceIndices, meshFaceIndices1, meshFaceIndices2;
  CGAL::ball_merge_surface_reconstruction_local(points, meshFaceIndices, params::concurrency_tag(ctag));
  CGAL::ball_merge_surface_reconstruction_global(points, meshFaceIndices1, meshFaceIndices2, params::concurrency_tag(ctag));
  }
  {
  std::list<boost::container::small_vector<std::size_t,3>> meshFaceIndices, meshFaceIndices1, meshFaceIndices2;
  CGAL::ball_merge_surface_reconstruction_local(points, meshFaceIndices, params::concurrency_tag(ctag));
  CGAL::ball_merge_surface_reconstruction_global(points, meshFaceIndices1, meshFaceIndices2, params::concurrency_tag(ctag));
  }

#if defined(CGAL_LINKED_WITH_TBB)
  {
  tbb::concurrent_vector<Point> points;
  CGAL::IO::read_points(inFilename, std::back_inserter(points));
  tbb::concurrent_vector<tbb::concurrent_vector<std::size_t>> meshFaceIndices, meshFaceIndices1, meshFaceIndices2;
  CGAL::ball_merge_surface_reconstruction_local(points, meshFaceIndices, params::concurrency_tag(ctag));
  CGAL::ball_merge_surface_reconstruction_global(points, meshFaceIndices1, meshFaceIndices2, params::concurrency_tag(ctag));
  }
#endif
}

int main()
{
  test1();
  test2();
  test_containers();

  return 0;
}
