#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Optimal_bounding_box/optimization_algorithms.h>
#include <CGAL/Optimal_bounding_box/population.h>
#include <CGAL/Optimal_bounding_box/obb.h>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>

#include <CGAL/subdivision_method_3.h>
#include <CGAL/Timer.h>

//#define OBB_DEBUG_BENCHMARK

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

void
bench(const char* fname)
{
  std::vector<K::Point_3> sm_points, obb_points;
  std::ifstream in(fname);

  K::Point_3 p;
  int i = 0;
  while(in >> p){

    if(i % 2 == 0) // avoid normals
    {
      sm_points.push_back(p);
    }
    ++i;
  }

  std::cout << "input data (points + normals)= " << i << std::endl;
  std::cout << "number of points= " << sm_points.size() << std::endl;

  CGAL::Optimal_bounding_box::find_obb(sm_points, obb_points, false);

  std::cout << "done" << '\n';

#ifdef OBB_DEBUG_BENCHMARK
  std::cout.precision(17);
  for(int i =0; i < obb_points.size(); i++)
  {
    std::cout << obb_points[i] << std::endl;
  }

  CGAL::Surface_mesh<K::Point_3> mesh;
  CGAL::make_hexahedron(obb_points[0], obb_points[1], obb_points[2], obb_points[3], obb_points[4], obb_points[5],
      obb_points[6], obb_points[7], mesh);

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
