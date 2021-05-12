#define CGAL_OPTIMAL_BOUNDING_BOX_DEBUG

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/boost/graph/generators.h>
#include <CGAL/optimal_bounding_box.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

#include <array>
#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef K::FT                                                       FT;
typedef K::Point_3                                                  Point;

typedef CGAL::Polyhedron_3<K>                                       Mesh;

typedef CGAL::Oriented_bounding_box_traits_3<K>                     Traits;
typedef Traits::Matrix                                              Matrix;

bool is_equal(const FT d1, const FT d2)
{
  const FT epsilon = 1e-3;

  bool ok;
  if(std::is_floating_point<FT>::value)
    ok = CGAL::abs(d1 - d2) < (std::max)(epsilon * d1, epsilon);
  else
    ok = (d1 == d2);

  if(!ok)
  {
    std::cout << "Got " << d1 << " but expected: " << d2 << std::endl;
    return false;
  }

  return true;
}

template <typename PointRange>
void test_OBB_data(const PointRange& points,
                   const double expected_vol,
                   const bool with_convex_hull = true)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  // the algorithm is allowed to fail, but not too often
  int failure_count = 0;
  for(int i=0; i<10; ++i)
  {
    std::cout << "iter #" << i << std::endl;

    CGAL::Surface_mesh<Point> obb_mesh;
    CGAL::oriented_bounding_box(points, obb_mesh, CGAL::parameters::use_convex_hull(with_convex_hull));
    PMP::triangulate_faces(obb_mesh);

    // the triangulate algorithm might fail if the algorithm manages
    // to fit perfectly the box to have a true 0 volume
    if(CGAL::is_triangle_mesh(obb_mesh))
    {
      double vol = PMP::volume(obb_mesh);
      std::cout << "  volume is: " << vol << ", expected: " << expected_vol << std::endl;
      if(!is_equal(vol, expected_vol))
      {
        std::cout << "Failure!" << std::endl;
        ++failure_count;
      }
    }
  }

  std::cout << "failures: " << failure_count << std::endl;
  assert(failure_count < 2); // 10% failure
}

void test_OBB_of_mesh(const std::string fname,
                      const double expected_vol)
{
  std::cout << "Test: " << fname << std::endl;

  std::ifstream input(fname);
  Mesh mesh;
  if(!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << fname << " is not a valid input file." << std::endl;
    std::exit(1);
  }

  std::vector<Point> points;
  for(const auto& v : vertices(mesh))
    points.push_back(v->point());

  test_OBB_data(points, expected_vol);
}

void test_OBB_of_point_set(const std::string fname,
                           const double expected_vol)
{
  std::cout << "Test: " << fname << std::endl;

  std::ifstream input(fname);
  if(!input)
  {
    std::cerr << fname << " is not a valid input file." << std::endl;
    std::exit(1);
  }

  std::deque<Point> points;
  double x, y, z;
  while(input >> x >> y >> z)
    points.emplace_back(x, y, z);

  test_OBB_data(points, expected_vol, false /*no convex hull due to degenerate data*/);
}

int main()
{
  std::cout.precision(17);

  test_OBB_of_mesh("data/elephant.off", 0.294296);
  test_OBB_of_mesh("data/long_tetrahedron.off", 0.04);
  test_OBB_of_mesh("data/reference_tetrahedron.off", 1);

  // degenerate cases, disabled because
  // - some testsuite platforms are too slow in debug, and they timeout
  // - the algorithm does not fully support degenerate data: it returns a box which is optimal
  //   in terms of volume (0), but is not necessarily optimal in the lower dimensions (i.e., the base
  //   of the OBB is not optimal).
//  test_OBB_of_mesh("data/triangles.off", 0); // 2D data set
//  test_OBB_of_mesh("data/flat_mesh.off", 0); // 2D data set
//  test_OBB_of_point_set("data/points_2D.xyz", 0); // 2D data set
//  test_OBB_of_point_set("data/points_1D.xyz", 0); // 1D data set
//  test_OBB_of_point_set("data/points_0D.xyz", 0); // 0D data set

  std::cout << "Done!" << std::endl;

  return 0;
}
