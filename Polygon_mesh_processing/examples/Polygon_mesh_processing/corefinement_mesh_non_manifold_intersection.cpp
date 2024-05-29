#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <iostream>
#include <string>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Surface_mesh<K::Point_3>                        Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const std::string filename1 = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/blobby.off");
  const std::string filename2 = (argc > 2) ? argv[2] : CGAL::data_file_path("meshes/eight.off");

  Mesh mesh1, mesh2;
  if(!PMP::IO::read_polygon_mesh(filename1, mesh1) || !PMP::IO::read_polygon_mesh(filename2, mesh2))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }


  PMP::Corefinement::Non_manifold_output_visitor<Mesh> visitor(mesh1, mesh2);

  Mesh out;
  bool valid_inter = PMP::corefine_and_compute_intersection(mesh1, mesh2, out, CGAL::parameters::visitor(visitor));

  if(valid_inter)
  {
    std::cout << "Intersection was successfully computed as a manifold mesh\n";
    CGAL::IO::write_polygon_mesh("inter.off", out, CGAL::parameters::stream_precision(17));
  }
  else
  {
    std::cout << "Intersection was successfully computed but is non-manifold, exporting a triangle soup\n";

    std::vector<K::Point_3> points;
    std::vector< std::array<std::size_t, 3> > polygons;

    visitor.extract_intersection(points, polygons);
    CGAL::IO::write_polygon_soup("inter_soup.off", points, polygons, CGAL::parameters::stream_precision(17));
    // make the soup topologically manifold (but geometrically self-intersecting)
    PMP::orient_polygon_soup(points, polygons);
    // fill a mesh with the intersection
    PMP::polygon_soup_to_polygon_mesh(points, polygons, out);
    CGAL::IO::write_polygon_mesh("inter.off", out, CGAL::parameters::stream_precision(17));
  }

  return 0;
}
