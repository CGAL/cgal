#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/VSA/approximate_mesh.h>

#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

int main()
{
  // read input surface triangle mesh
  Polyhedron input;
  std::ifstream file("data/mask.off");
  file >> input;

  // output indexed triangles
  std::vector<Kernel::Point_3> anchors;
  std::vector<CGAL::cpp11::array<std::size_t, 3> > triangles; // triplets of indices

  // free function interface with named parameters
  bool is_manifold = CGAL::VSA::approximate_mesh(input,
    CGAL::parameters::seeding_method(CGAL::VSA::HIERARCHICAL). // hierarchical seeding
    max_number_of_proxies(200). // seeding with maximum number of proxies
    number_of_iterations(30). // number of clustering iterations after seeding
    anchors(std::back_inserter(anchors)). // anchor vertices
    triangles(std::back_inserter(triangles))); // indexed triangles

  std::cout << "#anchor vertices: " << anchors.size() << std::endl;
  std::cout << "#triangles: " << triangles.size() << std::endl;

  if (is_manifold) {
    std::cout << "oriented, 2-manifold output." << std::endl;

    // convert from soup to polyhedron mesh
    CGAL::Polygon_mesh_processing::orient_polygon_soup(anchors, triangles);
    Polyhedron mesh;
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(anchors, triangles, mesh);
    if (CGAL::is_closed(mesh) && (!CGAL::Polygon_mesh_processing::is_outward_oriented(mesh)))
      CGAL::Polygon_mesh_processing::reverse_face_orientations(mesh);

    std::ofstream out("mask-dump.off");
    out << mesh;
    out.close();
  }

  return EXIT_SUCCESS;
}
