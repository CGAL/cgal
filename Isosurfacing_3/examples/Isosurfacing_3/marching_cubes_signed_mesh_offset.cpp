#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_3/Explicit_Cartesian_grid_domain_3.h>
#include <CGAL/Isosurfacing_3/marching_cubes_3.h>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Side_of_triangle_mesh.h>

#include <CGAL/IO/polygon_soup_io.h>

#include <iostream>
#include <vector>

using Kernel = CGAL::Simple_cartesian<double>;
using FT = typename Kernel::FT;
using Point = typename Kernel::Point_3;
using Vector = typename Kernel::Vector_3;

using Grid = CGAL::Isosurfacing::Cartesian_grid_3<Kernel>;

using Mesh = CGAL::Surface_mesh<Point>;

using Primitive = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
using Traits = CGAL::AABB_traits<Kernel, Primitive>;
using Tree = CGAL::AABB_tree<Traits>;

using Point_range = std::vector<Point>;
using Polygon_range = std::vector<std::vector<std::size_t> >;

// computes the Euclidean distance from query point p to the mesh
// via the AABB tree data structure
inline Kernel::FT distance_to_mesh(const Tree& tree,
                                   const Point& p)
{
  const Point cp = tree.closest_point(p);
  return sqrt((p - cp).squared_length());
}

int main(int argc, char* argv[])
{
  const std::string input_name = (argc > 1) ? argv[1] :CGAL::data_file_path("meshes/cross.off");
  const int n_voxels = 20;
  const FT offset_value = 0.2;

  // load input mesh
  Mesh mesh_input;
  if(!CGAL::IO::read_OFF(input_name, mesh_input))
  {
    std::cerr << "Could not read input mesh" << std::endl;
    return EXIT_FAILURE;
  }

  // compute loose bounding box of the mesh
  CGAL::Bbox_3 aabb_grid = CGAL::Polygon_mesh_processing::bbox(mesh_input);
  const FT loose_offset = offset_value + 0.01;
  Vector aabb_increase_vec = Vector(loose_offset, loose_offset, loose_offset);
  aabb_grid += (Point(aabb_grid.xmax(), aabb_grid.ymax(), aabb_grid.zmax()) + aabb_increase_vec).bbox();
  aabb_grid += (Point(aabb_grid.xmin(), aabb_grid.ymin(), aabb_grid.zmin()) - aabb_increase_vec).bbox();

  // construct AABB tree and functor to address inside/outside point queries
  Tree tree(mesh_input.faces_begin(), mesh_input.faces_end(), mesh_input);
  CGAL::Side_of_triangle_mesh<Mesh, CGAL::GetGeomTraits<Mesh>::type> sotm(mesh_input);

  // create grid
  Grid grid { n_voxels, n_voxels, n_voxels, aabb_grid };

  for(std::size_t k=0; k<grid.zdim(); ++k) {
    for(std::size_t j=0; j<grid.ydim(); ++j) {
      for(std::size_t i=0; i<grid.xdim(); ++i)
      {
        const FT pos_x = i * grid.spacing()[0] + grid.bbox().xmin();
        const FT pos_y = j * grid.spacing()[1] + grid.bbox().ymin();
        const FT pos_z = k * grid.spacing()[2] + grid.bbox().zmin();
        const Point p(pos_x, pos_y, pos_z);

        // compute unsigned distance to input mesh
        grid.value(i, j, k) = distance_to_mesh(tree, p);

        // sign distance so that it is negative inside the mesh
        const bool is_inside = (sotm(p) == CGAL::ON_BOUNDED_SIDE);
        if(is_inside)
          grid.value(i, j, k) *= -1.0;
      }
    }
  }

  // create domain from the grid
  auto domain = CGAL::Isosurfacing::create_explicit_Cartesian_grid_domain(grid);

  // containers for the triangle soup output
  Point_range points;
  Polygon_range polygons;

  // execute marching cubes with an isovalue equating offset
  CGAL::Isosurfacing::marching_cubes(domain, offset_value, points, polygons);

  // save the output
  CGAL::IO::write_polygon_soup("marching_cubes_signed_mesh_offset.off", points, polygons);

  return EXIT_SUCCESS;
}
