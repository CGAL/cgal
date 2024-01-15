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
  const Point x = tree.closest_point(p);
  return std::sqrt((p - x).squared_length());
}

// Usage : marching_cubes_multiple_mesh_offsets input.off
int main(int argc, char *argv[])
{
  const std::string input_name = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/cross.off");

  std::cout << "Input file: " << input_name << std::endl;

  // load input mesh
  Mesh mesh_input;
  if(!CGAL::IO::read_OFF(input_name, mesh_input))
  {
    std::cerr << "Could not read input mesh" << std::endl;
    return EXIT_FAILURE;
  }

  // construct loose bounding box from input mesh
  CGAL::Bbox_3 aabb_grid = CGAL::Polygon_mesh_processing::bbox(mesh_input);
  const FT loose_offset = 3.0;
  Vector aabb_increase_vec = Vector(loose_offset, loose_offset, loose_offset);
  aabb_grid += (Point(aabb_grid.xmax(), aabb_grid.ymax(), aabb_grid.zmax()) + aabb_increase_vec).bbox();
  aabb_grid += (Point(aabb_grid.xmin(), aabb_grid.ymin(), aabb_grid.zmin()) - aabb_increase_vec).bbox();

  // construct AABB tree and functor to address inside/outside point queries
  Tree tree(mesh_input.faces_begin(), mesh_input.faces_end(), mesh_input);
  CGAL::Side_of_triangle_mesh<Mesh, CGAL::GetGeomTraits<Mesh>::type> sotm(mesh_input);

  // create grid
  const int n_voxels = 250;
  Grid grid { n_voxels, n_voxels, n_voxels, aabb_grid };

  for(std::size_t z = 0; z < grid.zdim(); ++z) {
    for(std::size_t y = 0; y < grid.ydim(); ++y) {
      for(std::size_t x = 0; x < grid.xdim(); ++x)
      {
        const FT pos_x = x * grid.spacing()[0] + grid.bbox().xmin();
        const FT pos_y = y * grid.spacing()[1] + grid.bbox().ymin();
        const FT pos_z = z * grid.spacing()[2] + grid.bbox().zmin();
        const Point p(pos_x, pos_y, pos_z);

        // compute unsigned distance to input mesh
        grid.value(x, y, z) = distance_to_mesh(tree, p);

        // sign distance so that it is negative inside the mesh
        const bool is_inside = (sotm(p) == CGAL::ON_BOUNDED_SIDE);
        if(is_inside)
          grid.value(x, y, z) *= -1.0;
      }
    }
  }

  // create domain from the grid
  auto domain = CGAL::Isosurfacing::create_explicit_Cartesian_grid_domain(grid);

  // run Marching cubes with a range of offsets,
  // and save all output meshes to files "output-index.off"
  int index = 0;
  for(FT offset = 0.0; offset < 0.3; offset += 0.01, index++)
  {
          // containers for the triangle soup output
          Point_range points;
          Polygon_range triangles;

          // execute marching cubes with an isovalue equating offset
      std::cout << "Marching cubes with offset " << offset << "...";
      CGAL::Isosurfacing::marching_cubes(domain, offset, points, triangles);
      std::cout << "done" << std::endl;

          // save the output
      std::string filename("output-");
      filename.append(std::to_string(index));
      filename.append(std::string(".off"));
      std::cout << "Save to file " << filename << "...";
          CGAL::IO::write_polygon_soup(filename, points, triangles);
      std::cout << "done" << std::endl;
  }

  return EXIT_SUCCESS;
}
