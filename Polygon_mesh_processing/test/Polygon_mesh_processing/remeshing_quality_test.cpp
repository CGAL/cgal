#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/Adaptive_sizing_field.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <fstream>
#include <limits>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/hand.off");

  Mesh mesh;
  if (!PMP::IO::read_polygon_mesh(filename, mesh) || !CGAL::is_triangle_mesh(mesh))
  {
    std::cerr << "Not a valid input file." << std::endl;
    return 1;
  }

  std::cout << "Start remeshing of " << filename
    << " (" << num_faces(mesh) << " faces)..." << std::endl;

  const double tol = 0.001;
  const std::pair edge_min_max{0.001, 0.5};
  PMP::Adaptive_sizing_field sizing_field(tol, edge_min_max, faces(mesh), mesh);
  const unsigned int nb_iter = 5;

  PMP::isotropic_remeshing(
      faces(mesh),
      sizing_field,
      mesh,
      PMP::parameters::number_of_iterations(nb_iter).number_of_relaxation_steps(3)
      );

  /*
   * More information on quality metrics can be found here: https://ieeexplore.ieee.org/document/9167456
   */
  std::cout << "Remeshing done, checking triangle quality...\n" << std::endl;
  double qmin = std::numeric_limits<double>::max();      // minimum triangle quality
  double qavg = 0.;           // average quality
  double min_angle = std::numeric_limits<double>::max(); // minimum angle
  double avg_min_angle = 0.; // average minimum angle
  for (auto face : mesh.faces())
  {
    if (PMP::is_degenerate_triangle_face(face, mesh))
    {
      std::cout << "Found degenerate triangle!" << std::endl;
      continue;
    }

    // Calculate Q(t) triangle quality indicator
    std::vector<Mesh::Point> pts; pts.reserve(3);
    for (auto vx : mesh.vertices_around_face(mesh.halfedge(face)))
      pts.push_back(mesh.point(vx));

    double half_perim = 0.; // half-perimeter
    double lmax = 0.;       // longest edge
    for (int i = 0; i < 3; ++i)
    {
      const double length = CGAL::sqrt(CGAL::squared_distance(pts[i], pts[(i + 1) % 2]));

      half_perim += length;
      if (length > lmax) lmax = length;
    }
    half_perim /= 2.;
    const double area = CGAL::sqrt(CGAL::squared_area(pts[0], pts[1], pts[2]));

    const double face_quality = 6. / CGAL::sqrt(3.) * area / half_perim / lmax;

    qavg += face_quality;
    if (face_quality < qmin) qmin = face_quality;

    // Calculate minimum angle
    const auto v0 = pts[1] - pts[0];
    const auto v1 = pts[2] - pts[0];
    const auto v2 = pts[2] - pts[1];

    const double dotp0 = CGAL::scalar_product(v0,v1);
    const double angle0 = acos(dotp0 / (sqrt(v0.squared_length()) * sqrt(v1.squared_length())));
    const double dotp1 = CGAL::scalar_product(-v0, v2);
    const double angle1 = acos(dotp1 / (sqrt(v0.squared_length()) * sqrt(v2.squared_length())));
    const double angle2 = CGAL_PI - (angle0 + angle1);

    double curr_min_angle = angle1;
    if (angle1 < curr_min_angle) curr_min_angle = angle1;
    if (angle2 < curr_min_angle) curr_min_angle = angle2;

    avg_min_angle += curr_min_angle;
    if (curr_min_angle < min_angle) min_angle = curr_min_angle;
  }
  qavg /= mesh.number_of_faces();
  avg_min_angle /= mesh.number_of_faces();

  std::cout << "Mesh size: " << mesh.number_of_faces() << std::endl;
  std::cout << "Average quality: " << qavg << std::endl;
  std::cout << "Minimum quality: " << qmin << std::endl;
  std::cout << "Average minimum angle: " << avg_min_angle * 180. / CGAL_PI
    << " deg" << std::endl;
  std::cout << "Minimum angle: " << min_angle * 180. / CGAL_PI
    << " deg" << std::endl;

  return 0;
}