#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>

#include <CGAL/make_conforming_constrained_Delaunay_triangulation_3.h>

#include <CGAL/draw_constrained_triangulation_3.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/IO/write_MEDIT.h>

#include <vector>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = K::Point_3;
using Surface_mesh = CGAL::Surface_mesh<Point>;

namespace PMP = CGAL::Polygon_mesh_processing;

// Function to verify preconditions for a mesh input
bool verify_preconditions_mesh(const Surface_mesh& mesh)
{
  if(CGAL::is_triangle_mesh(mesh))
    return !CGAL::Polygon_mesh_processing::does_self_intersect(mesh);

  Surface_mesh triangle_mesh;
  CGAL::copy_face_graph(mesh, triangle_mesh);
  bool tri_ok = CGAL::Polygon_mesh_processing::triangulate_faces(triangle_mesh);
  if(!tri_ok)
    return false;

  return !CGAL::Polygon_mesh_processing::does_self_intersect(triangle_mesh);
}

// Function to verify preconditions for a polygon soup input
template <typename PointRange, typename PolygonRange>
bool verify_preconditions_soup(const PointRange& points, const PolygonRange& polygons)
{
  return !CGAL::Polygon_mesh_processing::does_polygon_soup_self_intersect(points, polygons);
}

int main(int argc, char* argv[])
{
  const auto filename = (argc > 1) ? argv[1]
                       : CGAL::data_file_path("meshes/mpi_and_sphere.off");

  CGAL::Surface_mesh<K::Point_3> mesh;
  if(!CGAL::IO::read_polygon_mesh(filename, mesh))
  {
    std::cerr << "Error: cannot read file " << filename << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Number of facets in " << filename << ": "
            << mesh.number_of_faces() << "\n";

  // Verify preconditions for the mesh input
  if(!verify_preconditions_mesh(mesh))
  {
    // If the mesh is not a valid triangle mesh or has self-intersections,
    // convert it to a polygon soup and verify the preconditions for the soup
    std::vector<K::Point_3> points;
    std::vector<std::vector<std::size_t>> polygons;
    CGAL::Polygon_mesh_processing::polygon_mesh_to_polygon_soup(mesh, points, polygons);
    if(!verify_preconditions_soup(points, polygons))
    {
      std::cerr << "Error: input polygon soup is not a valid input for CCDT_3\n";
      return EXIT_FAILURE;
    }

    std::cerr << "Error: input mesh is not a valid input for CCDT_3\n";
    return EXIT_FAILURE;
  }

  auto ccdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3(mesh);

  if(ccdt.number_of_constrained_facets() == 0)
  {
    std::cerr << "Error: no constrained facets in the CDT.\n";
    std::cerr << "Invalid input.\n";
    return EXIT_SUCCESS;
  }

  std::cout << "Number of constrained facets in the CDT: "
            << ccdt.number_of_constrained_facets() << '\n';

  std::ofstream ofs(argc > 2 ? argv[2] : "out.mesh");
  ofs.precision(17);
  CGAL::IO::write_MEDIT(ofs, ccdt);

  CGAL::draw(ccdt);

  return EXIT_SUCCESS;
}
