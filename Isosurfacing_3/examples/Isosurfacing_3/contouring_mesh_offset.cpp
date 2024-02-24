#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_3/dual_contouring_3.h>
#include <CGAL/Isosurfacing_3/Dual_contouring_domain_3.h>
#include <CGAL/Isosurfacing_3/Value_function_3.h>
#include <CGAL/Isosurfacing_3/Finite_difference_gradient_3.h>
#include <CGAL/Isosurfacing_3/marching_cubes_3.h>
#include <CGAL/Isosurfacing_3/Marching_cubes_domain_3.h>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Side_of_triangle_mesh.h>

#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <iostream>
#include <string>
#include <vector>

using Kernel = CGAL::Simple_cartesian<double>;
using FT = typename Kernel::FT;
using Point = typename Kernel::Point_3;
using Vector = typename Kernel::Vector_3;

using Mesh = CGAL::Surface_mesh<Point>;

using Grid = CGAL::Isosurfacing::Cartesian_grid_3<Kernel>;
using Values = CGAL::Isosurfacing::Value_function_3<Grid>;
using Gradients = CGAL::Isosurfacing::Finite_difference_gradient_3<Kernel>;

using Point_range = std::vector<Point>;
using Polygon_range = std::vector<std::vector<std::size_t> >;

struct Offset_oracle
{
  using Primitive = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
  using Traits = CGAL::AABB_traits<Kernel, Primitive>;
  using Tree = CGAL::AABB_tree<Traits>;

private:
  const bool is_closed;
  const Tree tree;
  CGAL::Side_of_triangle_mesh<Mesh, Kernel> sotm;

public:
  Offset_oracle(const Mesh& mesh)
    : is_closed(CGAL::is_closed(mesh)), tree(mesh.faces_begin(), mesh.faces_end(), mesh), sotm(mesh)
  { }

  FT distance(const Point& p) const
  {
    const Point cp = tree.closest_point(p);
    FT d = sqrt((p - cp).squared_length());

    if(is_closed && sotm(p) == (CGAL::ON_BOUNDED_SIDE))
      d *= -1;

    return d;
  }
};

void run_marching_cubes(const Grid& grid,
                         const FT offset_value,
                         const Offset_oracle& offset_oracle)
{
  using Domain = CGAL::Isosurfacing::Marching_cubes_domain_3<Grid, Values>;

  std::cout << "\n ---- " << std::endl;
  std::cout << "Running Marching Cubes with offset value = " << offset_value << std::endl;

  // fill up values
  auto mesh_distance = [&offset_oracle](const Point& p) { return offset_oracle.distance(p); };
  Values values { mesh_distance, grid };
  Domain domain { grid, values };

  Point_range points;
  Polygon_range triangles;

  // run marching cubes
  std::cout << "Running Marching Cubes with isovalue = " << offset_value << std::endl;
  CGAL::Isosurfacing::marching_cubes(domain, offset_value, points, triangles,
                                     CGAL::parameters::do_not_triangulate_faces(true));

  std::cout << "Output #vertices (MC): " << points.size() << std::endl;
  std::cout << "Output #triangles (MC): " << triangles.size() << std::endl;
  CGAL::IO::write_polygon_soup("marching_cubes_offsets.off", points, triangles);
}

void run_dual_contouring(const Grid& grid,
                         const FT offset_value,
                         const Offset_oracle& offset_oracle)
{
  using Domain = CGAL::Isosurfacing::Dual_contouring_domain_3<Grid, Values, Gradients>;

  std::cout << "\n ---- " << std::endl;
  std::cout << "Running Dual Contouring with offset value = " << offset_value << std::endl;

  // fill up values and gradients
  auto mesh_distance = [&offset_oracle](const Point& p) { return offset_oracle.distance(p); };

  Values values { mesh_distance, grid };
  Gradients gradients { values };
  Domain domain { grid, values, gradients };

  // output containers
  Point_range points;
  Polygon_range triangles;

  // run dual contouring
  std::cout << "Running Dual Contouring with isovalue = " << offset_value << std::endl;
  CGAL::Isosurfacing::dual_contouring(domain, offset_value, points, triangles);

  std::cout << "Output #vertices (DC): " << points.size() << std::endl;
  std::cout << "Output #triangles (DC): " << triangles.size() << std::endl;
  CGAL::IO::write_polygon_soup("dual_contouring_mesh_offset.off", points, triangles);
}

int main(int argc, char** argv)
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/cross.off");
  const FT offset_value = (argc > 2) ? std::stod(argv[2]) : 0.2;

  if(offset_value < 0)
  {
    std::cerr << "Offset value must be positive" << std::endl;
    return EXIT_FAILURE;
  }

  Mesh mesh;
  if(!CGAL::IO::read_polygon_mesh(filename, mesh) || is_empty(mesh))
  {
    std::cerr << "Could not read input mesh" << std::endl;
    return EXIT_FAILURE;
  }

  if(CGAL::is_closed(mesh))
    std::cout << "Input mesh is closed - using signed distance offset" << std::endl;
  else
    std::cout << "Input mesh is not closed - using unsigned distance offset" << std::endl;

  // construct loose bounding box from input mesh
  CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh);

  const FT diag_length = sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                              CGAL::square(bbox.ymax() - bbox.ymin()) +
                              CGAL::square(bbox.zmax() - bbox.zmin()));
  const FT loose_offset = offset_value + 0.1 * diag_length;

  Vector aabb_increase_vec = Vector(loose_offset, loose_offset, loose_offset);
  bbox += (Point(bbox.xmax(), bbox.ymax(), bbox.zmax()) + aabb_increase_vec).bbox();
  bbox += (Point(bbox.xmin(), bbox.ymin(), bbox.zmin()) - aabb_increase_vec).bbox();

  const int nv = 15;

  bbox = CGAL::Bbox_3{ -5, -5, -5, 15, 15, 15 };
  Grid grid { bbox, CGAL::make_array<std::size_t>(nv, nv, nv) };

  std::cout << "Bbox: " << grid.bbox() << std::endl;
  std::cout << "Cell dimensions: " << grid.spacing()[0] << " " << grid.spacing()[1] << " " << grid.spacing()[2] << std::endl;
  std::cout << "Cell #: " << grid.xdim() << ", " << grid.ydim() << ", " << grid.zdim() << std::endl;

  Offset_oracle offset_oracle(mesh);

  run_marching_cubes(grid, offset_value, offset_oracle);

  run_dual_contouring(grid, offset_value, offset_oracle);

  return EXIT_SUCCESS;
}
