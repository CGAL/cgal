#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/generators.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/Barycentric_coordinates_3/Mean_value_coordinates_3.h>
#include <fstream>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using Point_3 =  Kernel::Point_3;
using FT =  Kernel::FT;

using Surface_mesh =  CGAL::Surface_mesh<Point_3>;
namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char** argv) {

  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/sphere.off");
  Surface_mesh sm;

  if(!CGAL::IO::read_polygon_mesh(filename, sm))
  {
    std::cerr << "Invalid input: " << filename  << std::endl;
    return 1;
  }

  Surface_mesh deformed;
  deformed = sm;

  Surface_mesh quad_cage;

  const Point_3 p0(2, -2, -2), p0_new(5, -5, -5);
  const Point_3 p1(2, 2, -2), p1_new(3, 3, -3);
  const Point_3 p2(-2, 2, -2), p2_new(-2, 2, -2);
  const Point_3 p3(-2, -2, -2), p3_new(-3, -3, -3);

  const Point_3 p4(-2, -2, 2), p4_new(-3, -3, 3);
  const Point_3 p5(2, -2, 2), p5_new(4, -4, 4);
  const Point_3 p6(2, 2, 2), p6_new(2, 2, 3);
  const Point_3 p7(-2, 2, 2), p7_new(-3, 3, 3);

  CGAL::make_hexahedron(p0, p1, p2, p3, p4, p5, p6, p7, quad_cage,
                        CGAL::parameters::do_not_triangulate_faces(false));

  CGAL::Barycentric_coordinates::Mean_value_coordinates_3<Surface_mesh, Kernel> mv(quad_cage);
  auto vertex_point_map = get_property_map(CGAL::vertex_point, deformed);

  std::vector<FT> coords;
  std::vector<Point_3> target_cube{p0_new, p1_new, p2_new, p3_new,
                                   p4_new, p5_new, p6_new, p7_new};

  for(Surface_mesh::Vertex_index v : vertices(deformed)){

    const Point_3 vertex_val = get(vertex_point_map, v);
    coords.clear();
    mv(vertex_val, std::back_inserter(coords));

    FT x = FT(0), y = FT(0), z = FT(0);
    for(std::size_t i = 0; i < 8; i++){

      x += target_cube[i].x() * coords[i];
      y += target_cube[i].y() * coords[i];
      z += target_cube[i].z() * coords[i];
    }

    put(vertex_point_map, v, Point_3(x, y, z));
  }

  std::ofstream out_original("sphere.off");
  out_original << sm << std::endl;

  std::ofstream out_deformed("deformed_sphere.off");
  out_deformed << deformed << std::endl;

  return EXIT_SUCCESS;
}
