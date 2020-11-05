#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/boost/graph/iterator.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/OFF_reader.h>


#include <fstream>
#include <iostream>
#include <iterator>
#include <vector>

namespace PMP = CGAL::Polygon_mesh_processing;
namespace NP = CGAL::parameters;

typedef CGAL::Exact_predicates_inexact_constructions_kernel          K;
typedef CGAL::Surface_mesh<K::Point_3>                               Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor                 vertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor               halfedge_descriptor;

void merge_vertices(vertex_descriptor v_keep, vertex_descriptor v_rm, Mesh& mesh)
{
  std::cout << "merging vertices " << v_keep << " and " << v_rm << std::endl;

  for(halfedge_descriptor h : CGAL::halfedges_around_target(v_rm, mesh)){
    set_target(h, v_keep, mesh); // to ensure that no halfedge points at the deleted vertex
  }

  remove_vertex(v_rm, mesh);
}

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/blobby.off";
  std::ifstream input(filename);

  Mesh mesh;
  if(!input || !(input >> mesh) || num_vertices(mesh) == 0)
  {
    std::cerr << filename << " is not a valid off file." << std::endl;
    return EXIT_FAILURE;
  }

  // Artificially create non-manifoldness for the sake of the example by merging some vertices
  vertex_descriptor v0 = *(vertices(mesh).begin());
  vertex_descriptor v1 = *(--(vertices(mesh).end()));
  merge_vertices(v0, v1, mesh);

  // Count non manifold vertices
  int counter = 0;
  for(vertex_descriptor v : vertices(mesh))
  {
    if(PMP::is_non_manifold_vertex(v, mesh))
    {
      std::cout << "vertex " << v << " is non-manifold" << std::endl;
      ++counter;
    }
  }

  std::cout << counter << " non-manifold occurrence(s)" << std::endl;

  // Fix manifoldness by splitting non-manifold vertices
  std::vector<std::vector<vertex_descriptor> > duplicated_vertices;
  std::size_t new_vertices_nb = PMP::duplicate_non_manifold_vertices(mesh,
                                                                     NP::output_iterator(
                                                                       std::back_inserter(duplicated_vertices)));

  std::cout << new_vertices_nb << " vertices have been added to fix mesh manifoldness" << std::endl;

  for(std::size_t i=0; i<duplicated_vertices.size(); ++i)
  {
    std::cout << "Non-manifold vertex " << duplicated_vertices[i].front() << " was fixed by creating";
    for(std::size_t j=1; j<duplicated_vertices[i].size(); ++j)
      std::cout << " " << duplicated_vertices[i][j];
    std::cout << std::endl;
  }

  return EXIT_SUCCESS;
}
