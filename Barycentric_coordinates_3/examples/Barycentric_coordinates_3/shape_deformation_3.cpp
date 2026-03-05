#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/generators.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/Barycentric_coordinates_3.h>
#include <fstream>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 =  Kernel::Point_3;
using Surface_mesh =  CGAL::Surface_mesh<Point_3>;

int main(int argc, char** argv) {
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/armadillo.off");
  Surface_mesh sm;

  if(!CGAL::IO::read_polygon_mesh(filename, sm)) {
    std::cerr << "Invalid input: " << filename  << std::endl;
    return 1;
  }

  Surface_mesh deformed;
  deformed = sm;

  Surface_mesh quad_cage;

  const Point_3 p0( 200, -200, -200), p0_new( 500, -500, -500);
  const Point_3 p1( 200,  200, -200), p1_new( 300,  300, -300);
  const Point_3 p2(-200,  200, -200), p2_new(-200,  200, -200);
  const Point_3 p3(-200, -200, -200), p3_new(-300, -300, -300);

  const Point_3 p4(-200, -200, 200), p4_new(-300, -300, 300);
  const Point_3 p5( 200, -200, 200), p5_new( 400, -400, 400);
  const Point_3 p6( 200,  200, 200), p6_new( 200,  200, 300);
  const Point_3 p7(-200,  200, 200), p7_new(-300,  300, 300);

  CGAL::make_hexahedron(p0, p1, p2, p3, p4, p5, p6, p7, quad_cage,
                        CGAL::parameters::do_not_triangulate_faces(false));

  CGAL::Barycentric_coordinates::Mean_value_coordinates_3<Surface_mesh> mv(quad_cage);
  auto vertex_point_map = get_property_map(CGAL::vertex_point, deformed);

  std::vector<double> coords;
  std::vector<Point_3> target_cube{p0_new, p1_new, p2_new, p3_new,
                                   p4_new, p5_new, p6_new, p7_new};

  for(Surface_mesh::Vertex_index v : vertices(deformed)) {
    const Point_3 vertex_val = get(vertex_point_map, v);
    coords.clear();
    mv(vertex_val, std::back_inserter(coords));

    Point_3 p = CGAL::Barycentric_coordinates::apply_barycentric_coordinates(target_cube, coords);

    put(vertex_point_map, v, p);
  }

  std::ofstream out_original("armadillo.off");
  out_original << sm << std::endl;

  std::ofstream out_deformed("deformed_armadillo_mv.off");
  out_deformed << deformed << std::endl;

  return EXIT_SUCCESS;
}
