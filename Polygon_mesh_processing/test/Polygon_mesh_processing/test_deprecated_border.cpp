#include <CGAL/Installation/internal/disable_deprecation_warnings_and_errors.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/border.h>

#include <CGAL/Iterator_range.h>

#include <iostream>
#include <iterator>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel                K;
typedef CGAL::Surface_mesh<K::Point_3>                                     Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

typedef boost::graph_traits<Mesh>::halfedge_descriptor                     halfedge_descriptor;
typedef boost::graph_traits<Mesh>::face_iterator                           face_iterator;

int main(int, char**)
{
  Mesh mesh;
  if(!CGAL::IO::read_polygon_mesh(CGAL::data_file_path("meshes/quads_to_stitch.off"), mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return EXIT_FAILURE;
  }

  std::vector<halfedge_descriptor> border_hedges;
  PMP::extract_boundary_cycles(mesh, std::back_inserter(border_hedges));
  std::cout << border_hedges.size() << " boundary cycle(s)" << std::endl;
  assert(border_hedges.size() == 2);

  assert(PMP::number_of_borders(mesh) == 2);

  std::vector<halfedge_descriptor> bhs;
  PMP::border_halfedges(mesh, std::back_inserter(bhs));
  std::cout << bhs.size() << " border halfedge(s)" << std::endl;
  assert(bhs.size() == 20);

  // now with only the first four faces
  bhs.clear();
  CGAL::Iterator_range<face_iterator> fr(faces(mesh).begin(), std::next(faces(mesh).begin(), 4));
  PMP::border_halfedges(fr, mesh, std::back_inserter(bhs));
  std::cout << bhs.size() << " border halfedge(s) [4 faces]" << std::endl;
  assert(bhs.size() == 10);

  std::cout << "Done" << std::endl;
  return EXIT_SUCCESS;
}
