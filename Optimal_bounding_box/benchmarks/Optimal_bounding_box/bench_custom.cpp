#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/Eigen_linear_algebra_traits.h>
#include <CGAL/Optimal_bounding_box/optimal_bounding_box.h>
#include <CGAL/subdivision_method_3.h>
#include <CGAL/Timer.h>

#include <iostream>
#include <fstream>
#include <vector>

//#define CGAL_OPTIMAL_BOUNDING_BOX_DEBUG_BENCHMARK

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

void bench(const char* fname)
{
  std::vector<K::Point_3> sm_points, obb_points;
  std::ifstream in(fname);

  K::Point_3 p;
  int i = 0;
  while(in >> p)
  {
    if(i % 2 == 0) // avoid normals
      sm_points.push_back(p);

    ++i;
  }

  std::cout << "input data (points + normals)= " << i << std::endl;
  std::cout << "number of points= " << sm_points.size() << std::endl;

  CGAL::Eigen_linear_algebra_traits la_traits;

  // use convex hull - > true
  // no convex hull - > false
  CGAL::Optimal_bounding_box::compute_optimal_bounding_box(sm_points, obb_points, la_traits, false);

  std::cout << "done" << '\n';

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG_BENCHMARK
  std::cout.precision(17);
  for(int i =0; i < obb_points.size(); i++)
    std::cout << obb_points[i] << std::endl;

  CGAL::Surface_mesh<K::Point_3> mesh;
  CGAL::make_hexahedron(obb_points[0], obb_points[1], obb_points[2], obb_points[3],
                        obb_points[4], obb_points[5], obb_points[6], obb_points[7], mesh);

  std::ofstream out("/tmp/result_obb.off");
  out << mesh;
  out.close();
#endif
}

int main(int argc, char* argv[])
{
  bench(argv[1]);

  return 0;
}
