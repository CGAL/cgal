// Minimal test for user Fast_winding_number implementation (orders 1,2,3)
// Steps only: load mesh -> build AABB tree -> build face property maps -> create FWN(1,2,3)
// sample points in bbox -> compute approximations and exact winding number -> print comparison.

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <boost/property_map/property_map.hpp>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <vector>

#include <CGAL/Fast_winding_number.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = Kernel::Point_3;
using Vector = Kernel::Vector_3;
using FT = Kernel::FT;
using Mesh = CGAL::Surface_mesh<Point>;
using Primitive = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
using Traits = CGAL::AABB_traits_3<Kernel, Primitive>;
using Tree = CGAL::AABB_tree<Traits>;

int main() {
  Mesh mesh;
  const std::string filename = CGAL::data_file_path("meshes/elephant.off");
  if(!CGAL::IO::read_polygon_mesh(filename, mesh)) {
    std::cerr << "Failed to read mesh: " << filename << '\n';
    return EXIT_FAILURE;
  }
  if(!CGAL::is_triangle_mesh(mesh)) {
    std::cerr << "Input mesh not triangulated" << std::endl;
    return EXIT_FAILURE;
  }
  Tree tree(faces(mesh).begin(), faces(mesh).end(), mesh);
  tree.build();

  // Property maps (face -> normal, area, centroid)
  using face_descriptor = boost::graph_traits<Mesh>::face_descriptor;

  using Face_normal_tag = CGAL::dynamic_face_property_t<Vector>;
  using Face_normal_map = typename boost::property_map<Mesh, Face_normal_tag>::type;
  using Face_area_tag = CGAL::dynamic_face_property_t<FT>;
  using Face_area_map = typename boost::property_map<Mesh, Face_area_tag>::type;
  using Face_centroid_tag = CGAL::dynamic_face_property_t<Point>;
  using Face_centroid_map = typename boost::property_map<Mesh, Face_centroid_tag>::type;

  auto vpm = get_const_property_map(CGAL::vertex_point, mesh);
  // Instantiate Fast_winding_number for orders 1,2,3 (beta chosen as 2.0)
  const FT beta = FT(2.0);
  Face_normal_map face_normal_map_ = get(Face_normal_tag(), mesh, Vector(0., 0., 0.));
  Face_area_map face_area_map_ = get(Face_area_tag(), mesh, 0.);
  Face_centroid_map face_centroid_map_ = get(Face_centroid_tag(), mesh, Point(0., 0., 0.));
  // Compute face normals, areas, centroids
  namespace PMP = CGAL::Polygon_mesh_processing;
  PMP::compute_face_normals(mesh, face_normal_map_);
  for(face_descriptor f : faces(mesh)) {
    double area = PMP::face_area(f, mesh);

    const auto hedge = halfedge(f, mesh);
    const auto vertices = vertices_around_face(hedge, mesh);
    CGAL_precondition(vertices.size() > 0);

    FT sum = FT(0), x = FT(0), y = FT(0), z = FT(0);
    for(const auto vertex : vertices) {
      const Point& point = get(vpm, vertex);
      x += point.x();
      y += point.y();
      z += point.z();
      sum += FT(1);
    }
    CGAL_precondition(sum > FT(0));
    x /= sum;
    y /= sum;
    z /= sum;
    put(face_centroid_map_, f, Point(x, y, z));
    put(face_area_map_, f, area);
  }
  CGAL::Fast_winding_number<Mesh, Face_normal_map, Face_area_map, Face_centroid_map, Tree, int(1)> fwn1(
      mesh, face_normal_map_, face_area_map_, face_centroid_map_, tree, CGAL::parameters::fast_winding_number_beta(beta));
  CGAL::Fast_winding_number<Mesh, Face_normal_map, Face_area_map, Face_centroid_map, Tree, int(2)> fwn2(
      mesh, face_normal_map_, face_area_map_, face_centroid_map_, tree, CGAL::parameters::fast_winding_number_beta(beta));
  CGAL::Fast_winding_number<Mesh, Face_normal_map, Face_area_map, Face_centroid_map, Tree, int(3)> fwn3(
      mesh, face_normal_map_, face_area_map_, face_centroid_map_, tree, CGAL::parameters::fast_winding_number_beta(beta));

  // Sample grid points in bbox
  CGAL::Bbox_3 bb = CGAL::Polygon_mesh_processing::bbox(mesh);
  std::vector<Point> samples;
  const int N = 10; // 10^3 = 1000 points
  for(int ix = 0; ix < N; ++ix)
    for(int iy = 0; iy < N; ++iy)
      for(int iz = 0; iz < N; ++iz) {
        double fx = double(ix) / (N - 1);
        double fy = double(iy) / (N - 1);
        double fz = double(iz) / (N - 1);
        double x = bb.xmin() + fx * (bb.xmax() - bb.xmin());
        double y = bb.ymin() + fy * (bb.ymax() - bb.ymin());
        double z = bb.zmin() + fz * (bb.zmax() - bb.zmin());
        samples.emplace_back(x, y, z);
      }

  double sum_err1 = 0, sum_err2 = 0, sum_err3 = 0;
  double max_err1 = 0, max_err2 = 0, max_err3 = 0;
  std::size_t count = 0;
  for(const Point& p : samples) {
    double exact = CGAL::to_double(fwn3.exact_winding_number(p));
    double w1 = CGAL::to_double(fwn1.fast_winding_number(p));
    double w2 = CGAL::to_double(fwn2.fast_winding_number(p));
    double w3 = CGAL::to_double(fwn3.fast_winding_number(p));
    double e1 = std::abs(w1 - exact);
    double e2 = std::abs(w2 - exact);
    double e3 = std::abs(w3 - exact);
    sum_err1 += e1;
    sum_err2 += e2;
    sum_err3 += e3;
    if(e1 > max_err1)
      max_err1 = e1;
    if(e2 > max_err2)
      max_err2 = e2;
    if(e3 > max_err3)
      max_err3 = e3;
    ++count;
    std::cout << p << " : exact=" << exact << "  w1=" << w1 << "  w2=" << w2 << "  w3=" << w3 << "\n";
  }
  std::cout << "Samples: " << count << "\n";
  std::cout << "Order  MeanError  MaxError\n";
  std::cout << "1      " << (sum_err1 / count) << "  " << max_err1 << "\n";
  std::cout << "2      " << (sum_err2 / count) << "  " << max_err2 << "\n";
  std::cout << "3      " << (sum_err3 / count) << "  " << max_err3 << "\n";

  return EXIT_SUCCESS;
}
