#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3>             Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/blobby.off";
  std::ifstream input(filename);

  Mesh mesh;
  if (!input || !(input >> mesh))
  {
    std::cerr << "Input mesh is not a valid off file." << std::endl;
    return 1;
  }
  input.close();

  std::cout << "Test surface_self_intersection\n";
  std::vector< std::vector<K::Point_3> >polylines;

  PMP::experimental::surface_self_intersection(mesh, std::back_inserter(polylines));

  //dump polylines
  std::ofstream output("intersection_polylines.cgal");
  for(const std::vector<K::Point_3>& polyline : polylines)
  {
    output << polyline.size() << " ";
    std::copy(polyline.begin(), polyline.end(),std::ostream_iterator<K::Point_3>(output," "));
    output << "\n";
  }
  output.close();

  std::cout << "Number of vertices before autorefinement " << mesh.number_of_vertices() << "\n";
  PMP::experimental::autorefine(mesh);
  std::cout << "Number of vertices after autorefinement " << mesh.number_of_vertices() << "\n";

  output.open("mesh_autorefined.off");
  output << mesh;
  output.close();

  input.open(filename);
  mesh.clear();
  input >> mesh;
  std::cout << "Number of vertices before self-intersection removal " << mesh.number_of_vertices() << "\n";
  if (!PMP::experimental::autorefine_and_remove_self_intersections(mesh))
    std::cout << "WARNING: Cannot remove all self-intersections\n";
  std::cout << "Number of vertices after self-intersection removal " << mesh.number_of_vertices() << "\n";

  output.open("mesh_fixed.off");
  output << std::setprecision(17) << mesh;
  output.close();

  return 0;
}
