#define PHI 1.6180339887498948482

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/generators.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#include <CGAL/Barycentric_coordinates_3/Discrete_harmonic_coordinates_3.h>

using Kernel =  CGAL::Exact_predicates_inexact_constructions_kernel;
using FT = Kernel::FT;
using Point_3 =  Kernel::Point_3;
using Surface_mesh = CGAL::Surface_mesh<Point_3>;
namespace PMP = CGAL::Polygon_mesh_processing;
using CP3 = CGAL::Barycentric_coordinates::Computation_policy_3;

int main() {
  Surface_mesh icosahedron;
  CGAL::make_icosahedron(icosahedron, Point_3(0.0, 0.0, 0.0), 2.0);
  PMP::triangulate_faces(faces(icosahedron), icosahedron);

  std::vector<FT> coords;
  std::vector<Point_3> queries{
    Point_3(-1, 1 + PHI, PHI), Point_3(0.5, (1+3*PHI)/2, PHI/2), Point_3(1, 1+PHI, -PHI), //Boundary
    Point_3(-1, 1, 1), Point_3(0, 0, 1), Point_3(0, 2, 1),                                //Interior
    Point_3(0, 2*PHI, 4), Point_3(0, 3, 2*PHI), Point_3(4, 0, 0)};                        //Exterior

  std::cout << std::endl << "Discrete harmonic coordinates : " << std::endl << std::endl;

  for (const auto& query : queries){

    coords.clear();
    CGAL::Barycentric_coordinates::discrete_harmonic_coordinates_3(
      icosahedron, query, std::back_inserter(coords), CP3::FAST);

    // Output discrete harmonics coordinates.
    for (std::size_t i = 0; i < coords.size() -1; ++i) {
      std::cout << coords[i] << ", ";
    }
    std::cout << coords[coords.size() -1] << std::endl;
  }
  std::cout << std::endl;

  return EXIT_SUCCESS;
}
