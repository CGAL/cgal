#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>

#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <iostream>
#include <fstream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel                              Kernel;
typedef Kernel::Point_3                                                                  Point;
typedef CGAL::Linear_cell_complex_traits<3, Kernel>                                      MyTraits;
typedef CGAL::Linear_cell_complex_for_bgl_combinatorial_map_helper<2, 3, MyTraits>::type LCC;

typedef boost::graph_traits<LCC>::vertex_descriptor                                      vertex_descriptor;
typedef boost::graph_traits<LCC>::halfedge_descriptor                                    halfedge_descriptor;
typedef boost::graph_traits<LCC>::face_descriptor                                        face_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/mech-holes-shark.off");

  LCC mesh;
  if(!PMP::IO::read_polygon_mesh(filename, mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  // Incrementally fill the holes
  unsigned int nb_holes = 0;
  for ( halfedge_descriptor h : halfedges(mesh))
  {
    if(is_border(h,mesh))
    {
      std::vector<face_descriptor>  patch_facets;
      std::vector<vertex_descriptor> patch_vertices;
      bool success = std::get<0>(PMP::triangulate_refine_and_fair_hole(mesh,
                                                                       h,
                                                                       std::back_inserter(patch_facets),
                                                                       std::back_inserter(patch_vertices),
                                                                       PMP::parameters::vertex_point_map(get(CGAL::vertex_point, mesh))
                                                                                       .geom_traits(Kernel())));

      std::cout << "* Number of facets in constructed patch: " << patch_facets.size() << std::endl;
      std::cout << "  Number of vertices in constructed patch: " << patch_vertices.size() << std::endl;
      std::cout << "  Is fairing successful: " << success << std::endl;
      ++nb_holes;
    }
  }

  std::cout << std::endl;
  std::cout << nb_holes << " holes have been filled" << std::endl;

  std::ofstream out("filled_LCC.off");
  out.precision(17);
  CGAL::IO::write_OFF(out, mesh);

  return 0;
}
