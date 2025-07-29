#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Side_of_triangle_mesh.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef K::Point_3 Point_3;
typedef CGAL::Bbox_3 Bbox_3;
typedef CGAL::Surface_mesh<Point_3> Mesh;

typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
typedef CGAL::AABB_traits_3<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const std::string in_filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/tetrahedron.off");

  Mesh mesh;
  if(!CGAL::IO::read_polygon_mesh(in_filename, mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return EXIT_FAILURE;
  }

  const bool is_closed = CGAL::is_closed(mesh);
  if(!is_closed) {
    std::cerr << "Warning: the mesh is not closed." << std::endl;
  }

  auto vol_id_map = mesh.add_property_map<face_descriptor, std::size_t>().first;
  std::size_t vccn = PMP::volume_connected_components(mesh, vol_id_map,
    CGAL::parameters::do_orientation_tests(false));
  std::size_t ccn = PMP::internal::number_of_connected_components(mesh);
  if (vccn != ccn) {
    std::cerr << "Error: input has nested connected components" << std::endl;
    return EXIT_FAILURE;
  }

  if (!PMP::is_outward_oriented(mesh)) {
    std::cout << "Mesh was not outward oriented, reorienting it." << std::endl;
    PMP::reverse_face_orientations(mesh);
  }

  // save input to build dir for convenience
  CGAL::IO::write_polygon_mesh("input_mesh.off", mesh, CGAL::parameters::stream_precision(17));

  CGAL::Bbox_3 bbox = PMP::bbox(mesh);
  std::cout << "Bounding box: " << bbox << std::endl;

  // Scale the bounding box to have significant outside values
  const FT scale_factor = 1.5;
  bbox.scale(scale_factor);

  // Create an AABB tree from the mesh
  Tree tree(faces(mesh).first, faces(mesh).second, mesh);
  std::cout << "AABB tree created (" << tree.size() << " faces)" << std::endl;

  // Side of triangle mesh
  CGAL::Side_of_triangle_mesh<Mesh, K> sotm(mesh);

  // sample the bbox with a grid of n points, and compute the distance to the mesh
  // If the mesh is closed, that distance is signed (negative inside, positive outside)
  // write the result to a CSV file with format: x,y,z,distance
  const int n = (argc > 2) ? std::atoi(argv[2]) : 20;
  const char* out_filename = (argc > 3) ? argv[3] : "distances.csv";

  // let's also dump it as a colored PLY point set for visualization
  // the color will be based on the distance and whether the point is inside or outside the mesh
  // Outside: shades of red (positive distance)
  // Inside: shades of blue (negative distance)
  // Color intensity is based on the absolute distance value

  // The max distance is achieved at the corners of the bounding box
  // so query the first eight corners to find the maximum distance
  FT max_distance = 0;
  for(const Point_3& corner : {Point_3(bbox.xmin(), bbox.ymin(), bbox.zmin()),
                               Point_3(bbox.xmax(), bbox.ymin(), bbox.zmin()),
                               Point_3(bbox.xmin(), bbox.ymax(), bbox.zmin()),
                               Point_3(bbox.xmax(), bbox.ymax(), bbox.zmin()),
                               Point_3(bbox.xmin(), bbox.ymin(), bbox.zmax()),
                               Point_3(bbox.xmax(), bbox.ymin(), bbox.zmax()),
                               Point_3(bbox.xmin(), bbox.ymax(), bbox.zmax()),
                               Point_3(bbox.xmax(), bbox.ymax(), bbox.zmax())}) {
    FT d = CGAL::approximate_sqrt(tree.squared_distance(corner));
    max_distance = std::max(max_distance, d);
  }

  const int num_vertices_per_dir = n + 1;
  const int total_vertices = std::pow(num_vertices_per_dir, 3);

  std::ofstream ply_out("distances.ply");
  ply_out << "ply\n"
          << "format ascii 1.0\n"
          << "element vertex " << total_vertices << "\n"
          << "property double x\n"
          << "property double y\n"
          << "property double z\n"
          << "property uchar red\n"
          << "property uchar green\n"
          << "property uchar blue\n"
          << "end_header\n";
  ply_out.precision(17);

  std::ofstream out(out_filename);
  out.precision(17);

  if (!out) {
    std::cerr << "Could not open output file: " << out_filename << std::endl;
    return EXIT_FAILURE;
  }

  // Loop from 0 to n (inclusive) for n+1 vertices per direction
  for(int i=0; i<=n; ++i)
  {
    for(int j=0; j<=n; ++j)
    {
      for(int k=0; k<=n; ++k)
      {
        Point_3 p(bbox.xmin() + (bbox.xmax() - bbox.xmin()) * i / n,
                  bbox.ymin() + (bbox.ymax() - bbox.ymin()) * j / n,
                  bbox.zmin() + (bbox.zmax() - bbox.zmin()) * k / n);

        FT d = CGAL::approximate_sqrt(tree.squared_distance(p));
        // std::cout << "Point: " << p << " distance to mesh: " << d << std::endl;

        if (is_closed) {
          d *= (sotm(p) == CGAL::ON_BOUNDED_SIDE) ? -1 : 1;
        }

        out << p.x() << "," << p.y() << "," << p.z() << "," << d << "\n";

        // Color based on distance
        // The intensity is maximal for points close to the surface (small distance), and minimal for points far away (max_distance)
        unsigned char red = 0, green = 0, blue = 0;
        if (d > 0) {
          // Outside the mesh
          red = static_cast<unsigned char>(255 * (1 - d / max_distance));
        } else {
          // Inside the mesh
          blue = static_cast<unsigned char>(255 * (1 + d / max_distance));
        }

        ply_out << p.x() << " " << p.y() << " " << p.z() << " "
                << static_cast<int>(red) << " "
                << static_cast<int>(green) << " "
                << static_cast<int>(blue) << "\n";
      }
    }
  }

  std::cout << "Done" << std::endl;
  return EXIT_SUCCESS;
}
