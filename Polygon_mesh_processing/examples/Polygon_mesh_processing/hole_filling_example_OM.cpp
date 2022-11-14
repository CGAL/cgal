#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>
#include <CGAL/boost/graph/helpers.h>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

#include <cassert>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   Kernel;

typedef OpenMesh::PolyMesh_ArrayKernelT< >                    Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor          vertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor        halfedge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor            face_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/mech-holes-shark.off");

  Mesh mesh;
  OpenMesh::IO::read_mesh(mesh, filename);

  // Incrementally fill the holes
  unsigned int nb_holes = 0;
  for(halfedge_descriptor h : halfedges(mesh))
  {
    if(CGAL::is_border(h,mesh))
    {
      std::vector<face_descriptor>  patch_facets;
      std::vector<vertex_descriptor> patch_vertices;
      bool success = std::get<0>(PMP::triangulate_refine_and_fair_hole(mesh, h,
                                                                       std::back_inserter(patch_facets),
                                                                       std::back_inserter(patch_vertices),
                                                                       CGAL::parameters::vertex_point_map(get(CGAL::vertex_point, mesh))
                                                                                        .geom_traits(Kernel())));

      assert(CGAL::is_valid_polygon_mesh(mesh));

      std::cout << "* FILL HOLE NUMBER " << ++nb_holes << std::endl;
      std::cout << "  Number of facets in constructed patch: " << patch_facets.size() << std::endl;
      std::cout << "  Number of vertices in constructed patch: " << patch_vertices.size() << std::endl;
      std::cout << "  Is fairing successful: " << success << std::endl;
    }
  }

  assert(CGAL::is_valid_polygon_mesh(mesh));
  std::cout << std::endl;
  std::cout << nb_holes << " holes have been filled" << std::endl;

  OpenMesh::IO::write_mesh(mesh, "filled_OM.off");

  return 0;
}
