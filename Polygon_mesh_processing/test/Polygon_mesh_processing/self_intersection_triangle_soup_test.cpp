#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/tags.h>

#include <CGAL/Timer.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     EPICK;
typedef CGAL::Exact_predicates_exact_constructions_kernel       EPECK;

namespace PMP = ::CGAL::Polygon_mesh_processing;
namespace CP = ::CGAL::parameters;

template <typename K>
int test_self_intersections(const std::string filename,
                            const bool expected)
{
  std::vector<typename K::Point_3> points;
  std::vector< std::array<std::size_t, 3> > triangles;

  bool read_ok = CGAL::IO::read_polygon_soup(filename, points, triangles);

  if ( !read_ok ) {
    std::cerr << "Error: cannot read file: " << filename << std::endl;
    return 1;
  }

  std::cout << "Reading file: " << filename << std::endl;

  CGAL::Timer timer;
  timer.start();

  std::vector<std::pair<std::size_t, std::size_t> > intersected_tris;

  PMP::triangle_soup_self_intersections<CGAL::Parallel_if_available_tag>(points, triangles, std::back_inserter(intersected_tris));
  bool intersecting_1 = !intersected_tris.empty();

  std::cout << "self_intersections test took " << timer.time() << " sec." << std::endl;
  std::cout << intersected_tris.size() << " pairs of triangles are intersecting." << std::endl;

  timer.reset();

  bool intersecting_2 =
    PMP::does_triangle_soup_self_intersect<CGAL::Parallel_if_available_tag>(points, triangles);

  std::cout << "does_self_intersect test took " << timer.time() << " sec." << std::endl;
  std::cout << (intersecting_2 ? "There is a self-intersection." :
                                 "There are no self-intersections.") << std::endl;

  assert(intersecting_1 == intersecting_2);
  assert(intersecting_1 == expected);

  std::cout << filename << " passed the tests." << std::endl << std::endl;

  return 0;
}

template <typename K>
int test_limited_self_intersections(const std::string& filename)
{
  std::vector<typename K::Point_3> points;
  std::vector< std::array<std::size_t, 3> > triangles;

  bool read_ok = CGAL::IO::read_polygon_soup(filename, points, triangles);

  if ( !read_ok ) {
    std::cerr << "Error: cannot read file: " << filename << std::endl;
    return 1;
  }

  CGAL::Timer timer;
  timer.start();

  std::vector<std::pair<std::size_t, std::size_t> > intersected_tris;

#ifdef CGAL_LINKED_WITH_TBB
  PMP::triangle_soup_self_intersections()<CGAL::Parallel_if_available_tag>(
        points, triangles,
        std::back_inserter(intersected_tris), CGAL::parameters::maximum_number(40));
  std::cout << "self_intersections test for 40 SI took " << timer.time() << " sec." << std::endl;
  std::cout << "Found " << intersected_tris.size() << " SIs." << std::endl;
  if(intersected_tris.size() < 40)
  {
    std::cerr<<"Not enough intersections found in parallel."<<std::endl;
    return 1;
  }
  intersected_tris.clear();
  timer.reset();
#endif
  PMP::triangle_soup_self_intersections<CGAL::Sequential_tag>(
    points, triangles,
    std::back_inserter(intersected_tris), CGAL::parameters::maximum_number(40));
  std::cout << "self_intersections test for 40 SI took " << timer.time() << " sec." << std::endl;
  timer.reset();
  if(intersected_tris.size() != 40)
  {
    std::cerr<<"Too many intersections found in sequential"<<std::endl;
    return 1;
  }
  return 0;
}

int main(int argc, char** argv)
{
  // If file(s) are provided, the associated expected result must also be provided.
  // Note that this expected value is a Boolean that is passed in command line
  // with either 'true' or 'false' (and not integers), that is for example:
  // > self_intersection_surface_mesh_test data/U.off false

  // First test ----------------------------------------------------------------
  bool expected = false;
  std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/elephant.off");
  if(argc > 1) {
    assert(argc > 2);
    std::stringstream ss(argv[2]);
    ss >> std::boolalpha >> expected;
    assert(!ss.fail()); // make sure that argv[2] is either 'true' or 'false'
  }

  std::cout << "First test (EPICK):" << std::endl;
  int r = test_self_intersections<EPICK>(filename, expected);

  std::cout << "First test (EPECK):" << std::endl;
  r += test_self_intersections<EPECK>(filename, expected);

  // Second test ---------------------------------------------------------------
  expected = true;
  filename = (argc > 3) ? argv[3] : CGAL::data_file_path("meshes/mannequin-devil.off");
  if(argc > 3) {
    assert(argc > 4);
    std::stringstream ss(argv[4]);
    ss >> std::boolalpha >> expected;
    assert(!ss.fail());
  }

  std::cout << "Second test (EPICK):" << std::endl;
  r += test_self_intersections<EPICK>(filename, expected);

  std::cout << "Second test (EPECK):" << std::endl;
  r += test_self_intersections<EPECK>(filename, expected);

  // Third test ----------------------------------------------------------------
  expected = true;
  filename = (argc > 5) ? argv[5] : "data/overlapping_triangles.off";
  if(argc > 5) {
    assert(argc > 6);
    std::stringstream ss(argv[6]);
    ss >> std::boolalpha >> expected;
    assert(!ss.fail());
  }

  std::cout << "Third test (EPICK):" << std::endl;
  r += test_self_intersections<EPICK>(filename, expected);

  std::cout << "Third test (EPECK):" << std::endl;
  r += test_self_intersections<EPECK>(filename, expected);

  // Fourth test ----------------------------------------------------------------
  expected = true;
  filename = (argc > 7) ? argv[7] : "data_degeneracies/degtri_single.off";
  if(argc > 7) {
    assert(argc > 8);
    std::stringstream ss(argv[8]);
    ss >> std::boolalpha >> expected;
    assert(!ss.fail());
  }

  std::cout << "Fourth test (EPICK):" << std::endl;
  r += test_self_intersections<EPICK>(filename, expected);

  std::cout << "Fourth test (EPECK):" << std::endl;
  r += test_self_intersections<EPECK>(filename, expected);

  filename = (argc > 9) ? argv[9] : CGAL::data_file_path("meshes/mannequin-devil.off");

  std::cout << "Test with maximum_number (EPICK):" << std::endl;
  r += test_limited_self_intersections<EPICK>(filename);

  std::cout << "Test with maximum_number (EPECK):" << std::endl;
  r += test_limited_self_intersections<EPECK>(filename);

  // extra hand written tests
  {
    // shared edge with same point: no self-intersection
    typedef EPICK::Point_3 Point_3;
    std::vector<EPICK::Point_3> points = {Point_3(0,0,0), Point_3(1,0,0), Point_3(0,1,0), Point_3(0,-1,0)};
    std::vector< std::array<std::size_t, 3> > triangles={{0,1,2}, {0,1,3}};
    assert(!PMP::does_triangle_soup_self_intersect(points, triangles));
  }
  {
    // shared edge with duplicated points: self-intersection
    typedef EPICK::Point_3 Point_3;
    std::vector<EPICK::Point_3> points = {Point_3(0,0,0), Point_3(1,0,0), Point_3(0,1,0),Point_3(0,0,0), Point_3(1,0,0), Point_3(0,-1,0)};
    std::vector< std::array<std::size_t, 3> > triangles={{0,1,2}, {3,4,5}};
    assert(PMP::does_triangle_soup_self_intersect(points, triangles));
  }
  {
    // shared vertex with same point: no self-intersection
    typedef EPICK::Point_3 Point_3;
    std::vector<EPICK::Point_3> points = {Point_3(0,0,0), Point_3(1,0,0), Point_3(0,1,0), Point_3(0,2,0), Point_3(1,2,0)};
    std::vector< std::array<std::size_t, 3> > triangles={{0,1,2}, {3,4,2}};
    assert(!PMP::does_triangle_soup_self_intersect(points, triangles));
  }
  {
    // shared vertex with duplicated points: self-intersection
    typedef EPICK::Point_3 Point_3;
    std::vector<EPICK::Point_3> points = {Point_3(0,0,0), Point_3(1,0,0), Point_3(0,1,0), Point_3(0,2,0), Point_3(1,2,0), Point_3(0,1,0)};
    std::vector< std::array<std::size_t, 3> > triangles={{0,1,2}, {3,4,5}};
    assert(PMP::does_triangle_soup_self_intersect(points, triangles));
  }
  {
    // 4 triangles around a shared edge: no self-intersection
    typedef EPICK::Point_3 Point_3;
    std::vector<EPICK::Point_3> points = {Point_3(0,0,0), Point_3(1,0,0), Point_3(0,1,0), Point_3(0,-1,0), Point_3(0,0,1), Point_3(0,0,-1)};
    std::vector< std::array<std::size_t, 3> > triangles={{0,1,2},{0,1,3},{0,1,4},{0,1,5}};
    assert(!PMP::does_triangle_soup_self_intersect(points, triangles));
  }
  {
    // 4 triangles around a shared edge but two triangles intersecting: self-intersection
    typedef EPICK::Point_3 Point_3;
    std::vector<EPICK::Point_3> points = {Point_3(0,0,0), Point_3(1,0,0), Point_3(0,1,0), Point_3(0,0.5,0), Point_3(0,0,1), Point_3(0,0,-1)};
    std::vector< std::array<std::size_t, 3> > triangles={{0,1,2},{0,1,3},{0,1,4},{0,1,5}};
    assert(PMP::does_triangle_soup_self_intersect(points, triangles));
  }

  return r;
}
