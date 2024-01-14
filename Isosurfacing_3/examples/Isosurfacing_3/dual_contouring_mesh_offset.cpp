#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Isosurfacing_3/dual_contouring_3.h>
#include <CGAL/Isosurfacing_3/Implicit_Cartesian_grid_domain_3.h>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Side_of_triangle_mesh.h>

#include <CGAL/boost/graph/IO/OFF.h>

#include <iostream>
#include <string>
#include <vector>

using Kernel = CGAL::Simple_cartesian<double>;
using FT = typename Kernel::FT;
using Point = typename Kernel::Point_3;
using Vector = typename Kernel::Vector_3;

using Mesh = CGAL::Surface_mesh<Point>;

using Primitive = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
using Traits = CGAL::AABB_traits<Kernel, Primitive>;
using Tree = CGAL::AABB_tree<Traits>;

using Point_range = std::vector<Point>;
using Polygon_range = std::vector<std::vector<std::size_t> >;

int main(int, char**)
{
  const std::string input_name = CGAL::data_file_path("meshes/cross.off");
  const Vector grid_spacing(0.1, 0.1, 0.1);
  const FT offset_value = 0.2;

  Mesh mesh_input;
  if(!CGAL::IO::read_OFF(input_name, mesh_input))
  {
    std::cerr << "Could not read input mesh" << std::endl;
    return EXIT_FAILURE;
  }

  // compute bounding box
  CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh_input);
  Vector aabb_increase_vec = Vector(offset_value + 0.01, offset_value + 0.01, offset_value + 0.01);
  bbox += (Point(bbox.xmax(), bbox.ymax(), bbox.zmax()) + aabb_increase_vec).bbox();
  bbox += (Point(bbox.xmin(), bbox.ymin(), bbox.zmin()) - aabb_increase_vec).bbox();

  // construct AABB tree
  Tree tree(mesh_input.faces_begin(), mesh_input.faces_end(), mesh_input);

  CGAL::Side_of_triangle_mesh<Mesh, CGAL::GetGeomTraits<Mesh>::type> sotm(mesh_input);

  // functors for addressing distance and normal queries
  auto mesh_distance = [&tree](const Point& p)
  {
    const Point x = tree.closest_point(p);
    return sqrt((p - x).squared_length());
  };

  auto mesh_normal = [&tree](const Point& p)
  {
    const Point x = tree.closest_point(p);
    const Vector n = p - x; // @todo address case where norm = zero
    return n / sqrt(n.squared_length()); // normalize output vector
  };

  // create a domain with given bounding box and grid spacing
  auto domain = CGAL::Isosurfacing::create_implicit_Cartesian_grid_domain(bbox, grid_spacing,
                                                                          mesh_distance, mesh_normal);
  // containers for output indexed surface mesh
  Point_range points;
  Polygon_range polygons;

  // run dual contouring
  CGAL::Isosurfacing::dual_contouring(domain, offset_value, points, polygons);

  // save indexed face set to file, in the OFF format
  CGAL::IO::write_OFF("result.off", points, polygons);

  return EXIT_SUCCESS;
}
