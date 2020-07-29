#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/iterator.h>

#include <CGAL/Polygon_mesh_processing/intersection.h>

namespace PMP=CGAL::Polygon_mesh_processing;

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Mesh;

typedef boost::graph_traits<Mesh>::halfedge_descriptor    halfedge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor       face_descriptor;
typedef boost::graph_traits<Mesh>::vertex_descriptor      vertex_descriptor;

int main(int argc, char* argv[])
{
  const char* filename1 = (argc > 1) ? argv[1] : "data/blobby.off";
  std::ifstream input(filename1);

  Mesh mesh1;
  if ( !input || !(input >> mesh1) ) {
    std::cerr << filename1 << " is not a valid off file." << std::endl;
    return 1;
  }
  input.close();

  const char* filename2 = (argc > 2) ? argv[2] : "data/eight.off";
  input.open(filename2);

  Mesh mesh2;
  if ( !input || !(input >> mesh2) ) {
    std::cerr << filename2 << " is not a valid off file." << std::endl;
    return 1;
  }

  std::vector< std::vector<Point> > polylines;
  PMP::surface_intersection(mesh1, mesh2, std::back_inserter(polylines));

  //dump polylines
  std::ofstream output("intersection_polylines.cgal");
  output.precision(17);
  for(const std::vector<Point>& polyline : polylines)
  {
    output << polyline.size() << " ";
    std::copy(polyline.begin(), polyline.end(),std::ostream_iterator<Point>(output," "));
    output << "\n";
  }

  return 0;
}
