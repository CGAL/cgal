#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

#include <CGAL/boost/graph/dijkstra_shortest_path.h>
#include <CGAL/IO/polygon_mesh_io.h>


#include <string>
#include <vector>
#include <fstream>
#include <exception>
#include <algorithm>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = K::Point_3;
using Mesh = CGAL::Surface_mesh<Point>;

using vertex_descriptor = boost::graph_traits<Mesh>::vertex_descriptor;
using edge_descriptor = boost::graph_traits<Mesh>::edge_descriptor;
using halfedge_descriptor = boost::graph_traits<Mesh>::halfedge_descriptor;

namespace params = CGAL::parameters;


// Example main
int main(int argc, char** argv)
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/elephant.off");

  // Try building a surface_mesh
  Mesh sm;
  bool ok = CGAL::IO::read_polygon_mesh(filename, sm);
  if (!ok || !sm.is_valid() || sm.is_empty())
  {
    std::cerr << "Error: Invalid facegraph" << std::endl;
    std::cerr << "Filename = " << filename << std::endl;
    return EXIT_FAILURE;
  }

  const std::size_t i0 = 0;
  const std::size_t i1 = num_vertices(sm) / 2;


  // Get the vertex descriptors of the source and target vertices
  const vertex_descriptor vs = *vertices(sm).first;
  vertex_descriptor vt;
  std::size_t vid = 0;
  for (const vertex_descriptor v : vertices(sm))
  {
    if (vid++ == i1)
    {
      vt = v;
      break;
    }
  }

  std::list<halfedge_descriptor> halfedge_sequence;
  CGAL::dijkstra_shortest_path(vs, vt, sm,
    std::back_inserter(halfedge_sequence));

  assert(source(halfedge_sequence.front(), sm)==vs);
  assert(target(halfedge_sequence.back(), sm)==vt);

  // dump
  std::cout << "Shortest path between vertices " << i0 << " and " << i1
            << " is made of " << halfedge_sequence.size() << " halfedges." << std::endl;

  // Get the property map of the points of the mesh
  auto vpmap = get(CGAL::vertex_point, sm);

  std::ofstream out("shortest_path.polylines.txt");
  out << halfedge_sequence.size()+1 << " " << get(vpmap, source(halfedge_sequence.front(),sm));
  for (const halfedge_descriptor he : halfedge_sequence)
  {
    const vertex_descriptor v = target(he, sm);
    out << " " << get(vpmap, v);
  }
  out << std::endl;
  out.close();

  return EXIT_SUCCESS;
}
