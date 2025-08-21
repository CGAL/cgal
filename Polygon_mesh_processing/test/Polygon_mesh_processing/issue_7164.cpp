#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Random.h>
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = Kernel::Point_3;
using Mesh = CGAL::Surface_mesh<Point_3>;
namespace PMP = CGAL::Polygon_mesh_processing;

int main(/*int argc, char** argv*/)
{
  // A simple triangle
  std::vector<Point_3> pts_A;
  std::vector<std::vector<size_t>> trs_A;
  pts_A.emplace_back( 0.26641936088212415, 0.2664193608821242, 0.73358063911787585);
  pts_A.emplace_back(-0.14011519816541251, 0.6017979969632727, 1.1810107045967466);
  pts_A.emplace_back(-0.14011519816541279,-0.1810107045967464, 0.39820200303672726);
  trs_A.emplace_back(std::vector<size_t>{0,1,2});
  Mesh A;
  PMP::polygon_soup_to_polygon_mesh(pts_A, trs_A, A);

  // An open tetrahedron
  std::vector<Point_3> pts_B;
  std::vector<std::vector<size_t>> trs_B;
  pts_B.emplace_back(0,0,0);
  pts_B.emplace_back(1,1,0);
  pts_B.emplace_back(1,0,1);
  pts_B.emplace_back(0,1,1);
  trs_B.emplace_back(std::vector<size_t>{0,1,2});
  trs_B.emplace_back(std::vector<size_t>{3,1,0});
  trs_B.emplace_back(std::vector<size_t>{3,2,1});
  Mesh B;
  PMP::polygon_soup_to_polygon_mesh(pts_B, trs_B, B);

  double bound = 0.01 * 0.42149467833714593;
  PMP::bounded_error_Hausdorff_distance<CGAL::Sequential_tag>(A, B, bound);
  PMP::bounded_error_Hausdorff_distance<CGAL::Sequential_tag>(B, A, bound);

  return EXIT_SUCCESS;
}
