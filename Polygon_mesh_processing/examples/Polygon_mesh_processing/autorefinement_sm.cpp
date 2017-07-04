#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/intersection.h>

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

  
  std::cout << "Test surface_autointersection\n";
  std::vector< std::vector<K::Point_3> >polylines;

  PMP::surface_autointersection(mesh, std::back_inserter(polylines));
  
  //dump polylines
  std::ofstream output("intersection_polylines.cgal");
  BOOST_FOREACH(const std::vector<K::Point_3>& polyline, polylines)
  {
    output << polyline.size() << " ";
    std::copy(polyline.begin(), polyline.end(),std::ostream_iterator<K::Point_3>(output," "));
    output << "\n";
  }
  output.close();

  return 0;
}
