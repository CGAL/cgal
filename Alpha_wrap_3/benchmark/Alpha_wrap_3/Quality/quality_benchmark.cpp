#include <distance_utils.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/alpha_wrap_3.h>

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <array>
#include <cmath>
#include <iostream>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = Kernel::Point_3;
using Vector_3 = Kernel::Vector_3;
using Triangle_3 = Kernel::Triangle_3;
using FT = Kernel::FT;

using Mesh = CGAL::Surface_mesh<Point_3>;
using face_descriptor = boost::graph_traits<Mesh>::face_descriptor;

using Oracle = CGAL::Alpha_wraps_3::internal::Triangle_mesh_oracle<Kernel>;
using Dt = CGAL::Alpha_wraps_3::internal::Alpha_wrapper_3<Oracle>::Triangulation;

namespace PMP = CGAL::Polygon_mesh_processing;

std::array<FT, 3> triangle_angles(const Triangle_3& tr)
{
  FT sq_a = CGAL::squared_distance(tr[0], tr[1]);
  FT sq_b = CGAL::squared_distance(tr[1], tr[2]);
  FT sq_c = CGAL::squared_distance(tr[2], tr[0]);

  FT two_ab = 2. * CGAL::sqrt(sq_a) * CGAL::sqrt(sq_b);
  FT two_bc = 2. * CGAL::sqrt(sq_b) * CGAL::sqrt(sq_c);
  FT two_ca = 2. * CGAL::sqrt(sq_c) * CGAL::sqrt(sq_a);

  FT angle_a = (sq_b + sq_c - sq_a) / two_bc;
  FT angle_b = (sq_c + sq_a - sq_b) / two_ca;
  FT angle_c = (sq_a + sq_b - sq_c) / two_ab;
  if(angle_a < -1.) angle_a = -1.;
  if(angle_b < -1.) angle_b = -1.;
  if(angle_c < -1.) angle_c = -1.;
  if(angle_a > 1.) angle_a = 1.;
  if(angle_b > 1.) angle_b = 1.;
  if(angle_c > 1.) angle_c = 1.;
  angle_a = std::acos(angle_a);
  angle_b = std::acos(angle_b);
  angle_c = std::acos(angle_c);

  return {angle_a, angle_b, angle_c};
}

bool is_almost_degenerate(const Triangle_3& tr,
                          double threshold)
{
  FT sq_area = tr.squared_area();
  return (CGAL::sqrt(CGAL::to_double(sq_area)) < threshold);
}

auto surface_mesh_face_to_triangle(const face_descriptor fd,
                                   const Mesh& sm)
{
  typename boost::graph_traits<Mesh>::halfedge_descriptor hd = halfedge(fd,sm);
   return Triangle_3(sm.point(target(hd,sm)),
                     sm.point(target(next(hd,sm),sm)),
                     sm.point(target(next(next(hd,sm),sm),sm)));
}

double mean_min_angle(const Mesh& mesh)
{
  double mean_min_angle = 0.;
  for(const face_descriptor f : faces(mesh))
  {
    const Triangle_3 tr = surface_mesh_face_to_triangle(f, mesh);
    std::array<FT, 3> angles = triangle_angles(tr);

    FT min_angle = std::min({angles[0], angles[1], angles[2]});

    min_angle = min_angle * (180.0 / CGAL_PI);
    mean_min_angle += min_angle;
  }

  mean_min_angle /= static_cast<double>(mesh.number_of_faces());
  return mean_min_angle;
}

double mean_max_angle(const Mesh& mesh)
{
  double mean_max_angle = 0.;
  for(const face_descriptor f : faces(mesh))
  {
    const Triangle_3 tr = surface_mesh_face_to_triangle(f, mesh);
    std::array<FT, 3> angles = triangle_angles(tr);

    FT max_angle = std::max({angles[0], angles[1], angles[2]});

    max_angle = max_angle * (180.0 / CGAL_PI);
    mean_max_angle += max_angle;
  }

  mean_max_angle /= static_cast<double>(mesh.number_of_faces());
  return mean_max_angle;
}

double mean_radius_ratio(const Mesh& mesh,
                         double degenerate_threshold)
{
  double mean_radius_ratio = 0.;
  size_t num_almost_degenerate_tri = 0;
  for(const face_descriptor f : faces(mesh))
  {
    const Triangle_3 tr = surface_mesh_face_to_triangle(f, mesh);
    if(is_almost_degenerate(tr, degenerate_threshold))
    {
      ++num_almost_degenerate_tri;
      continue;
    }

    FT circumsphere_radius = std::sqrt(CGAL::squared_radius(tr[0], tr[1], tr[2]));

    FT a = std::sqrt(CGAL::squared_distance(tr[0], tr[1]));
    FT b = std::sqrt(CGAL::squared_distance(tr[1], tr[2]));
    FT c = std::sqrt(CGAL::squared_distance(tr[2], tr[0]));
    FT s = 0.5 * (a + b + c);
    FT inscribed_radius = std::sqrt((s * (s - a) * (s - b) * (s - c)) / s);
    FT radius_ratio = circumsphere_radius / inscribed_radius;
    radius_ratio /= 2.;  // normalized
    mean_radius_ratio += radius_ratio;
  }

  mean_radius_ratio /= static_cast<double>(mesh.number_of_faces() - num_almost_degenerate_tri);
  return mean_radius_ratio;
}

double mean_edge_ratio(const Mesh& mesh,
                       double degenerate_threshold)
{
  double mean_edge_ratio = 0.;
  size_t num_almost_degenerate_tri = 0;

  for(const face_descriptor f : faces(mesh))
  {
    const Triangle_3 tr = surface_mesh_face_to_triangle(f, mesh);
    if(is_almost_degenerate(tr, degenerate_threshold))
    {
      ++num_almost_degenerate_tri;
      continue;
    }

    FT a = std::sqrt(CGAL::squared_distance(tr[0], tr[1]));
    FT b = std::sqrt(CGAL::squared_distance(tr[1], tr[2]));
    FT c = std::sqrt(CGAL::squared_distance(tr[2], tr[0]));
    FT min_edge = std::min({a, b, c});
    FT max_edge = std::max({a, b, c});
    FT edge_ratio = max_edge / min_edge;

    mean_edge_ratio += edge_ratio;
  }

  mean_edge_ratio /= static_cast<double>(mesh.number_of_faces() - num_almost_degenerate_tri);
  return mean_edge_ratio;
}

double mean_aspect_ratio(const Mesh& mesh,
                         double degenerate_threshold)
{
  double mean_aspect_ratio = 0.;
  size_t num_almost_degenerate_tri = 0;
  for(const face_descriptor f : faces(mesh))
  {
    const Triangle_3 tr = surface_mesh_face_to_triangle(f, mesh);
    if(is_almost_degenerate(tr, degenerate_threshold))
    {
      ++num_almost_degenerate_tri;
      continue;
    }

    FT a = std::sqrt(CGAL::squared_distance(tr[0], tr[1]));
    FT b = std::sqrt(CGAL::squared_distance(tr[1], tr[2]));
    FT c = std::sqrt(CGAL::squared_distance(tr[2], tr[0]));
    FT s = 0.5 * (a + b + c);
    FT inscribed_radius = std::sqrt((s * (s - a) * (s - b) * (s - c)) / s);
    FT max_edge = std::max({a, b, c});
    FT aspect_ratio = max_edge / inscribed_radius;
    aspect_ratio /= (2. * std::sqrt(3.));  // normalized
    mean_aspect_ratio += aspect_ratio;
  }

  mean_aspect_ratio /= static_cast<double>(mesh.number_of_faces() - num_almost_degenerate_tri);
  return mean_aspect_ratio;
}

size_t num_almost_degenerate_tri(const Mesh& mesh,
                                 double degenerate_threshold)
{
  size_t num_almost_degenerate_tri = 0;
  for(const face_descriptor f : faces(mesh))
  {
    const Triangle_3 tr = surface_mesh_face_to_triangle(f, mesh);
    if(is_almost_degenerate(tr, degenerate_threshold))
    {
      ++num_almost_degenerate_tri;
    }
  }
  return num_almost_degenerate_tri;
}

int main(int argc, char** argv)
{
  const int argc_check = argc - 1;
  char *entry_name_ptr = nullptr;
  double relative_alpha_ratio = 20.;
  double relative_offset_ratio = 600.;

  for(int i=1; i<argc; ++i)
  {
    if(!strcmp("-i", argv[i]) && i < argc_check) {
      entry_name_ptr = argv[++i];
    } else if(!strcmp("-a", argv[i]) && i < argc_check) {
      relative_alpha_ratio = std::stod(argv[++i]);
    } else if(!strcmp("-d", argv[i]) && i < argc_check) {
      relative_offset_ratio = std::stod(argv[++i]);
    }
  }

  if(argc < 3 || relative_alpha_ratio <= 0.)
  {
    std::cerr << "Error: bad input parameters." << std::endl;
    return EXIT_FAILURE;
  }

  Mesh input_mesh;
  if(!PMP::IO::read_polygon_mesh(entry_name_ptr, input_mesh) ||
     is_empty(input_mesh) ||
     !is_triangle_mesh(input_mesh))
  {
    return EXIT_FAILURE;
  }

  CGAL::Bbox_3 bbox = PMP::bbox(input_mesh);
  const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                       CGAL::square(bbox.ymax() - bbox.ymin()) +
                                       CGAL::square(bbox.zmax() - bbox.zmin()));
  const double alpha = diag_length / relative_alpha_ratio;
  const double offset = diag_length / relative_offset_ratio;

  Mesh wrap;
  CGAL::alpha_wrap_3(input_mesh, alpha, offset, wrap);

  double degenerate_threshold = 0.;
  for(const face_descriptor f : faces(wrap))
  {
    const Triangle_3 tr = surface_mesh_face_to_triangle(f, wrap);
    degenerate_threshold += CGAL::sqrt(CGAL::to_double(tr.squared_area()));
  }

  degenerate_threshold /= wrap.number_of_faces();
  degenerate_threshold /= 1000;

  std::cout << mean_min_angle(wrap) << "\n";
  std::cout << mean_max_angle(wrap) << "\n";
  std::cout << mean_radius_ratio(wrap, degenerate_threshold) << "\n";
  std::cout << mean_edge_ratio(wrap, degenerate_threshold) << "\n";
  std::cout << mean_aspect_ratio(wrap, degenerate_threshold) << "\n";
  std::cout << wrap.number_of_faces() << "\n";
  std::cout << num_almost_degenerate_tri(wrap, degenerate_threshold) << "\n";
  std::cout << 100. * approximate_distance_relative_to_bbox(wrap, input_mesh, Aw3i::HAUSDORFF) << "\n";

  return 0;
}
