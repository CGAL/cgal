#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>

#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using Mesh = CGAL::Polyhedron_3<Kernel>;

using Vertex_handle = Mesh::Vertex_handle;
using Halfedge_handle = Mesh::Halfedge_handle;
using Facet_handle = Mesh::Facet_handle;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/mech-holes-shark.off");

  Mesh poly;
  if(!PMP::IO::read_polygon_mesh(filename, poly))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  // Incrementally fill the holes
  unsigned int nb_holes = 0;
  for(Halfedge_handle h : halfedges(poly))
  {
    if(h->is_border())
    {
      std::vector<Facet_handle>  patch_facets;
      std::vector<Vertex_handle> patch_vertices;
      bool success = std::get<0>(PMP::triangulate_refine_and_fair_hole(poly,
                                                                       h,
                                                                       CGAL::parameters::face_output_iterator(std::back_inserter(patch_facets))
                                                                                        .vertex_output_iterator(std::back_inserter(patch_vertices))
                                                                                        .vertex_point_map(get(CGAL::vertex_point, poly))
                                                                                        .geom_traits(Kernel())));

      std::cout << " Number of facets in constructed patch: " << patch_facets.size() << std::endl;
      std::cout << " Number of vertices in constructed patch: " << patch_vertices.size() << std::endl;
      std::cout << " Fairing : " << (success ? "succeeded" : "failed") << std::endl;
      ++nb_holes;
    }
  }

  std::cout << std::endl;
  std::cout << nb_holes << " holes have been filled" << std::endl;

  std::ofstream out("filled.off");
  out.precision(17);
  out << poly << std::endl;

  return 0;
}
