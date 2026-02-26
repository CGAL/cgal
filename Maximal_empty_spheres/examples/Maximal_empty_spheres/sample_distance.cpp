#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Side_of_triangle_mesh.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef K::Point_3 Point_3;
typedef K::Vector_3 Vector_3;
typedef CGAL::Bbox_3 Bbox_3;
typedef CGAL::Surface_mesh<Point_3> Mesh;

typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
typedef CGAL::AABB_traits_3<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

void translate(Mesh& mesh, const Vector_3& translation_vector)
{
  for(auto v : vertices(mesh)) {
    mesh.point(v) = mesh.point(v) + translation_vector;
  }
}

void scale(Mesh& mesh, const FT scale_factor)
{
  for(auto v : vertices(mesh)) {
    const Point_3& p = mesh.point(v);
    mesh.point(v) = Point_3(p.x() * scale_factor,
                            p.y() * scale_factor,
                            p.z() * scale_factor);
  }
}

void resize_to_bbox(Mesh& mesh,
                    const Bbox_3& bbox)
{
  // Compute current bounding box
  Bbox_3 curr_bbox = PMP::bbox(mesh);
  Point_3 curr_min{curr_bbox.xmin(), curr_bbox.ymin(), curr_bbox.zmin()};
  Point_3 curr_max{curr_bbox.xmax(), curr_bbox.ymax(), curr_bbox.zmax()};
  Point_3 curr_center{(curr_bbox.xmin() + curr_bbox.xmax()) / 2.0,
                      (curr_bbox.ymin() + curr_bbox.ymax()) / 2.0,
                      (curr_bbox.zmin() + curr_bbox.zmax()) / 2.0};

  // Target bbox min, max, center
  Point_3 target_min{bbox.xmin(), bbox.ymin(), bbox.zmin()};
  Point_3 target_max{bbox.xmax(), bbox.ymax(), bbox.zmax()};
  Point_3 target_center{(bbox.xmin() + bbox.xmax()) / 2.0,
                        (bbox.ymin() + bbox.ymax()) / 2.0,
                        (bbox.zmin() + bbox.zmax()) / 2.0};

  // Compute sizes
  FT curr_size_x = curr_bbox.xmax() - curr_bbox.xmin();
  FT curr_size_y = curr_bbox.ymax() - curr_bbox.ymin();
  FT curr_size_z = curr_bbox.zmax() - curr_bbox.zmin();
  FT target_size_x = bbox.xmax() - bbox.xmin();
  FT target_size_y = bbox.ymax() - bbox.ymin();
  FT target_size_z = bbox.zmax() - bbox.zmin();

  // Compute scale factor (uniform, fit largest axis)
  FT scale_factor = std::numeric_limits<double>::max();
  if (curr_size_x > 0) scale_factor = std::min(scale_factor, target_size_x / curr_size_x);
  if (curr_size_y > 0) scale_factor = std::min(scale_factor, target_size_y / curr_size_y);
  if (curr_size_z > 0) scale_factor = std::min(scale_factor, target_size_z / curr_size_z);

  // Translate mesh to origin
  translate(mesh, Vector_3(-curr_center.x(), -curr_center.y(), -curr_center.z()));

  // Scale mesh
  if (scale_factor != 1.0) {
    scale(mesh, scale_factor);
  }

  // Translate mesh to target center
  translate(mesh, Vector_3(CGAL::ORIGIN, target_center));
}

int main(int argc, char* argv[])
{
  const std::string in_filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/tetrahedron.off");

  Mesh mesh;
  if(!PMP::IO::read_polygon_mesh(in_filename, mesh))
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

#if 0
  CGAL::Bbox_3 bbox = PMP::bbox(mesh);
#else
  // scale the input such that it fits in a cube [-1,1]^3
  CGAL::Bbox_3 bbox = Bbox_3(-1, -1, -1, 1, 1, 1);
  resize_to_bbox(mesh, bbox);
#endif
  std::cout << "Bounding box: " << bbox << std::endl;

  // save input to build dir for convenience
  CGAL::IO::write_polygon_mesh("input_mesh.off", mesh, CGAL::parameters::stream_precision(17));

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
