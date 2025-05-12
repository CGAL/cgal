#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/boost/graph/copy_face_graph.h>

#include <CGAL/make_conforming_constrained_Delaunay_triangulation_3.h>

#include <CGAL/draw_constrained_triangulation_3.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = K::Point_3;
using Surface_mesh = CGAL::Surface_mesh<Point>;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const auto filename = (argc > 1) ? argv[1]
                       : CGAL::data_file_path("meshes/mpi_and_sphere.off");

  CGAL::Surface_mesh<K::Point_3> mesh;
  std::ifstream in(filename);
  if(!in || !(in >> mesh)) {
    std::cerr << "Error: cannot read file " << filename << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Number of facets in " << filename << ": "
            << mesh.number_of_faces() << "\n";

  // Check if the mesh is a triangle mesh
  bool triangle_mesh = CGAL::is_triangle_mesh(mesh);
  if(!triangle_mesh)
  {
    std::cout << "Mesh is not a triangle mesh, triangulate faces"
              << " to check self-intersections...\n";

    CGAL::Surface_mesh<K::Point_3> trimesh;
    CGAL::copy_face_graph(mesh, trimesh);
    PMP::triangulate_faces(trimesh);

    if(PMP::does_self_intersect(trimesh))
    {
      std::cout << "Mesh self-intersects, let's keep the triangulated version"
                << " for future autorefinement\n";
      CGAL::copy_face_graph(trimesh, mesh);
      mesh = std::move(trimesh);
      triangle_mesh = true;
    }
  }

  if(triangle_mesh && PMP::does_self_intersect(mesh))
  {
    std::cout << "Mesh is a self-intersecting triangle mesh, perform autorefinement...\n";
    PMP::autorefine(mesh);
  }

  std::cout << "Number of facets after preprocessing: "
            << mesh.number_of_faces() << "\n";

  auto ccdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3(mesh);

  std::cout << "Number of constrained facets in the CDT: "
            << ccdt.number_of_constrained_facets() << '\n';

  CGAL::draw(ccdt);

  return EXIT_SUCCESS;
}
