#define CGAL_AW3_DEBUG

#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/alpha_wrap_3.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/Point_set_3.h>

#include <boost/property_map/function_property_map.hpp>

#include <cassert>
#include <deque>
#include <fstream>
#include <iostream>
#include <vector>

using namespace CGAL::Alpha_wraps_3::internal;

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = Kernel::Point_3;
using Segment_3 = Kernel::Segment_3;
using Vector_3 = Kernel::Vector_3;

using Mesh = CGAL::Surface_mesh<Point_3>;

using Points = std::vector<Point_3>;
using Segments = std::deque<Segment_3>;
using Face = std::array<std::size_t, 3>;
using Faces = std::vector<Face>;

typedef CGAL::Point_set_3<Point_3> Point_set;

void test_mesh_API(const std::string& tm_filename,
                   const double alpha_rel = 10,
                   const double offset_rel = 300)
{
  std::cout << "test_mesh_API()" << std::endl;

  Mesh input_mesh;
  bool res = CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(tm_filename, input_mesh);
  assert(res);
  assert(!is_empty(input_mesh) && is_triangle_mesh(input_mesh));

  std::cout << "Input: " << num_vertices(input_mesh) << " vertices, " << num_faces(input_mesh) << " faces" << std::endl;
  std::cout << "Alpha rel: " << alpha_rel << " Offset rel: " << offset_rel << std::endl;

  Mesh wrap;
  auto in_vpm = get(CGAL::vertex_point, input_mesh);
  auto out_vpm = get(CGAL::vertex_point, wrap);

  // No alpha
//  CGAL::alpha_wrap_3(input_mesh, wrap);

//  CGAL::alpha_wrap_3(input_mesh, wrap,
//                     CGAL::parameters::vertex_point_map(in_vpm));

//  CGAL::alpha_wrap_3(input_mesh, wrap,
//                     CGAL::parameters::vertex_point_map(in_vpm),
//                     CGAL::parameters::vertex_point_map(out_vpm));

  // No offset
  CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(input_mesh);
  const Vector_3 longest_diag = Point_3(bbox.xmax(), bbox.ymax(), bbox.zmax()) -
                                Point_3(bbox.xmin(), bbox.ymin(), bbox.zmin());
  double longest_diag_length = CGAL::to_double(CGAL::approximate_sqrt(longest_diag.squared_length()));

  const double alpha = longest_diag_length / alpha_rel;

  CGAL::alpha_wrap_3(input_mesh, alpha, wrap);

  CGAL::alpha_wrap_3(input_mesh, alpha, wrap,
                     CGAL::parameters::vertex_point_map(in_vpm));

  CGAL::alpha_wrap_3(input_mesh, alpha, wrap,
                     CGAL::parameters::vertex_point_map(in_vpm),
                     CGAL::parameters::vertex_point_map(out_vpm));

  // Alpha and offset
  const double offset = longest_diag_length / offset_rel;

  CGAL::alpha_wrap_3(input_mesh, alpha, offset, wrap);

  CGAL::alpha_wrap_3(input_mesh, alpha, offset, wrap,
                     CGAL::parameters::vertex_point_map(in_vpm));

  CGAL::alpha_wrap_3(input_mesh, alpha, offset, wrap,
                     CGAL::parameters::vertex_point_map(in_vpm),
                     CGAL::parameters::vertex_point_map(out_vpm));
}

void test_triangles_API(const std::string& ts_filename,
                        const double alpha_rel = 10,
                        const double offset_rel = 300)
{
  std::cout << "test_triangles_API()" << std::endl;

  Points points;
  Faces faces;
  bool res = CGAL::IO::read_polygon_soup(ts_filename, points, faces);
  assert(res && !faces.empty());

  std::cout << "Input: " << points.size() << " vertices, " << faces.size() << " faces" << std::endl;
  std::cout << "Alpha rel: " << alpha_rel << " Offset rel: " << offset_rel << std::endl;

  Mesh wrap;
  auto in_pm = CGAL::Identity_property_map<Point_3>();
  auto out_vpm = get(CGAL::vertex_point, wrap);

  // No alpha
//  CGAL::alpha_wrap_3(points, faces, wrap);

//  CGAL::alpha_wrap_3(points, faces, wrap,
//                     CGAL::parameters::vertex_point_map(in_pm));

//  CGAL::alpha_wrap_3(points, faces, wrap,
//                     CGAL::parameters::vertex_point_map(in_pm),
//                     CGAL::parameters::vertex_point_map(out_vpm));

  // No offset
  CGAL::Bbox_3 bbox;
  for(const auto& f : faces)
    for(int i=0; i<3; ++i)
      bbox += points[f[i]].bbox();

  const Vector_3 longest_diag = Point_3(bbox.xmax(), bbox.ymax(), bbox.zmax()) -
                                Point_3(bbox.xmin(), bbox.ymin(), bbox.zmin());
  double longest_diag_length = CGAL::to_double(CGAL::approximate_sqrt(longest_diag.squared_length()));

  const double alpha = longest_diag_length / alpha_rel;

  CGAL::alpha_wrap_3(points, faces, alpha, wrap);

  CGAL::alpha_wrap_3(points, faces, alpha, wrap,
                     CGAL::parameters::vertex_point_map(in_pm));

  CGAL::alpha_wrap_3(points, faces, alpha, wrap,
                     CGAL::parameters::vertex_point_map(in_pm),
                     CGAL::parameters::vertex_point_map(out_vpm));

  // Alpha and offset
  const double offset = longest_diag_length / offset_rel;

  CGAL::alpha_wrap_3(points, faces, alpha, offset, wrap);

  CGAL::alpha_wrap_3(points, faces, alpha, offset, wrap,
                     CGAL::parameters::vertex_point_map(in_pm));

  CGAL::alpha_wrap_3(points, faces, alpha, offset, wrap,
                     CGAL::parameters::vertex_point_map(in_pm),
                     CGAL::parameters::vertex_point_map(out_vpm));
}

#if 0
void test_segments_API(const std::string& ss_filename,
                       const double alpha_rel = 10,
                       const double offset_rel = 300)
{
  Segments segments;
  std::ifstream ifs(ss_filename);
  int len = 0;
  while(ifs >> len)
  {
    std::vector<Point_3> polyline;
    while(len--)
    {
      Point_3 point;
      ifs >> point;
      polyline.push_back(point);
    }

    if(polyline.size() >= 2)
    {
      for(std::size_t i=0; i<polyline.size() - 1; ++i)
        segments.emplace_back(polyline[i], polyline[i+1]);
    }
  }

  std::cout << segments.size() << " segments (segment soup)" << std::endl;

  Mesh wrap;
  auto in_pm = CGAL::Identity_property_map<Point_3>();
  auto out_vpm = get(CGAL::vertex_point, wrap);

  // No alpha
//  CGAL::alpha_wrap_3(segments, wrap);

//  CGAL::alpha_wrap_3(segments, wrap,
//                     CGAL::parameters::vertex_point_map(in_pm));

//  CGAL::alpha_wrap_3(segments, wrap,
//                     CGAL::parameters::vertex_point_map(in_pm),
//                     CGAL::parameters::vertex_point_map(out_vpm));

  // No offset
  CGAL::Bbox_3 bbox;
  for(const Segment_3& s : segments)
  {
    bbox += s.source().bbox();
    bbox += s.target().bbox();
  }

  const Vector_3 longest_diag = Point_3(bbox.xmax(), bbox.ymax(), bbox.zmax()) -
                                Point_3(bbox.xmin(), bbox.ymin(), bbox.zmin());
  double longest_diag_length = CGAL::to_double(CGAL::approximate_sqrt(longest_diag.squared_length()));

  const double alpha = longest_diag_length / alpha_rel;

  CGAL::alpha_wrap_3(segments, alpha, wrap);

  CGAL::alpha_wrap_3(segments, alpha, wrap,
                     CGAL::parameters::vertex_point_map(in_pm));

  CGAL::alpha_wrap_3(segments, alpha, wrap,
                     CGAL::parameters::vertex_point_map(in_pm),
                     CGAL::parameters::vertex_point_map(out_vpm));

  // Alpha and offset
  const double offset = longest_diag_length / offset_rel;

  CGAL::alpha_wrap_3(segments, alpha, offset, wrap);

  CGAL::alpha_wrap_3(segments, alpha, offset, wrap,
                     CGAL::parameters::vertex_point_map(in_pm));

  CGAL::alpha_wrap_3(segments, alpha, offset, wrap,
                     CGAL::parameters::vertex_point_map(in_pm),
                     CGAL::parameters::vertex_point_map(out_vpm));
}
#endif

void test_points_API(const std::string& ps_filename,
                     const double alpha_rel = 10,
                     const double offset_rel = 300)
{
  std::cout << "test_points_API()" << std::endl;

  Point_set point_set;
  point_set.add_normal_map();
  bool res = CGAL::IO::read_points(ps_filename, point_set.index_back_inserter(),
                                   CGAL::parameters::point_map(point_set.point_push_map())
                                                    .normal_map(point_set.normal_push_map()));
  assert(res && !point_set.empty());

  std::cout << point_set.size() << " points (Point Set)" << std::endl;

  struct To_p : public CGAL::cpp98::unary_function<int, Point_3>
  {
    To_p(const Point_set& ps) : ps(ps) { }
    Point_3 operator()(const int i) const { return ps.point(i); }
    Point_set ps;
  };

  auto p3r = point_set.points();
  auto in_pm = CGAL::Identity_property_map<Point_3>();
  Mesh wrap;
  auto out_vpm = get(CGAL::vertex_point, wrap);

  // No alpha
//  CGAL::alpha_wrap_3(p3r, wrap);

//  CGAL::alpha_wrap_3(p3r, wrap,
//                     CGAL::parameters::vertex_point_map(in_pm));

//  CGAL::alpha_wrap_3(p3r, wrap,
//                     CGAL::parameters::vertex_point_map(in_pm),
//                     CGAL::parameters::vertex_point_map(out_vpm));

  // No offset
  CGAL::Bbox_3 bbox;
  for(const auto& i : point_set)
      bbox += point_set.point(i).bbox();

  const Vector_3 longest_diag = Point_3(bbox.xmax(), bbox.ymax(), bbox.zmax()) -
                                Point_3(bbox.xmin(), bbox.ymin(), bbox.zmin());
  double longest_diag_length = CGAL::to_double(CGAL::approximate_sqrt(longest_diag.squared_length()));

  const double alpha = longest_diag_length / alpha_rel;

  CGAL::alpha_wrap_3(p3r, alpha, wrap);

  CGAL::alpha_wrap_3(p3r, alpha, wrap,
                     CGAL::parameters::vertex_point_map(in_pm));

  CGAL::alpha_wrap_3(p3r, alpha, wrap,
                     CGAL::parameters::vertex_point_map(in_pm),
                     CGAL::parameters::vertex_point_map(out_vpm));

  // Alpha and offset
  const double offset = longest_diag_length / offset_rel;

  CGAL::alpha_wrap_3(p3r, alpha, offset, wrap);

  CGAL::alpha_wrap_3(p3r, alpha, offset, wrap,
                     CGAL::parameters::vertex_point_map(in_pm));

  CGAL::alpha_wrap_3(p3r, alpha, offset, wrap,
                     CGAL::parameters::vertex_point_map(in_pm),
                     CGAL::parameters::vertex_point_map(out_vpm));
}

int main(int argc, char** argv)
{
  const std::string tm_filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/building.off"); // mesh
  const std::string ts_filename = (argc > 2) ? argv[2] : CGAL::data_file_path("meshes/oblong-shuffled.off"); // triangle soup
  const std::string ss_filename = (argc > 3) ? argv[3] : CGAL::data_file_path("images/420.polylines.txt"); // segment soup
  const std::string ps_filename = (argc > 4) ? argv[4] : CGAL::data_file_path("points_3/b9_training.ply"); // point set

  test_mesh_API(tm_filename);
  test_triangles_API(ts_filename);
//  test_segments_API(ss_filename); // segment soup is not supported in free functions yet
  test_points_API(ps_filename);

  std::cout << "Done!" << std::endl;
  return EXIT_SUCCESS;
}
