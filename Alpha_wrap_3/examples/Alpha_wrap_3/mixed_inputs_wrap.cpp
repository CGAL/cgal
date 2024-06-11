#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/alpha_wrap_3.h>
#include <CGAL/Real_timer.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <iostream>
#include <string>

namespace AW3 = CGAL::Alpha_wraps_3;
namespace PMP = CGAL::Polygon_mesh_processing;

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;
using Segment_3 = K::Segment_3;

using Segments = std::vector<Segment_3>;
using Points = std::vector<Point_3>;
using Face = std::array<std::size_t, 3>;
using Faces = std::vector<Face>;

using Mesh = CGAL::Surface_mesh<Point_3>;

int main(int argc, char** argv)
{
  // Read the inputs
  const std::string ts_filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/armadillo.off"); // triangle soup
  const std::string ss_filename = (argc > 2) ? argv[2] : CGAL::data_file_path("images/420.polylines.txt"); // segment soup
  const std::string ps_filename = (argc > 3) ? argv[3] : CGAL::data_file_path("points_3/ball.ply"); // point set

  std::cout << "Triangle soup input: " << ts_filename << std::endl;
  std::cout << "Segment soup input: " << ss_filename << std::endl;
  std::cout << "Point set input: " << ps_filename << std::endl;

  // = read the soup
  Points points;
  Faces faces;
  if(!CGAL::IO::read_polygon_soup(ts_filename, points, faces) || points.empty() || faces.empty())
  {
    std::cerr << "Invalid soup input: " << ts_filename << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << points.size() << " points (triangle soup)" << std::endl;

  // = read the polylines
  std::ifstream ifs(ss_filename);
  Segments segments;
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

  // = read the points
  Points ps_points;
  if(!CGAL::IO::read_points(ps_filename, std::back_inserter(ps_points)))
  {
    std::cerr << "Invalid point set input: " << ps_filename << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << ps_points.size() << " points (point set)" << std::endl;

  const double relative_alpha = (argc > 4) ? std::stod(argv[4]) : 15.;
  const double relative_offset = (argc > 5) ? std::stod(argv[5]) : 450.;

  CGAL::Bbox_3 bbox = bbox_3(std::cbegin(points), std::cend(points));
  CGAL::Bbox_3 ps_bbox = bbox_3(std::cbegin(ps_points), std::cend(ps_points));
  bbox += ps_bbox;

  const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                       CGAL::square(bbox.ymax() - bbox.ymin()) +
                                       CGAL::square(bbox.zmax() - bbox.zmin()));
  const double alpha = diag_length / relative_alpha;
  const double offset = diag_length / relative_offset;

  CGAL::Real_timer t;
  t.start();

  using TS_Oracle = CGAL::Alpha_wraps_3::internal::Triangle_soup_oracle<K>;
  using SS_Oracle = CGAL::Alpha_wraps_3::internal::Segment_soup_oracle<K, TS_Oracle>;
  using Oracle = CGAL::Alpha_wraps_3::internal::Point_set_oracle<K, SS_Oracle>;

  TS_Oracle ts_oracle(K{});
  SS_Oracle ss_oracle(ts_oracle);
  Oracle oracle(ss_oracle);

  oracle.add_triangle_soup(points, faces, CGAL::parameters::default_values());
  oracle.add_segment_soup(segments, CGAL::parameters::default_values());
  oracle.add_point_set(ps_points, CGAL::parameters::default_values());

  CGAL::Alpha_wraps_3::internal::Alpha_wrapper_3<Oracle> aw3(oracle);

  Mesh wrap;
  aw3(alpha, offset, wrap);

  t.stop();
  std::cout << "Took " << t.time() << std::endl;

  std::string ts_name = std::string(ts_filename);
  ts_name = ts_name.substr(ts_name.find_last_of("/") + 1, ts_name.length() - 1);
  ts_name = ts_name.substr(0, ts_name.find_last_of("."));
  std::string ss_name = std::string(ss_filename);
  ss_name = ss_name.substr(ss_name.find_last_of("/") + 1, ss_name.length() - 1);
  ss_name = ss_name.substr(0, ss_name.find_last_of("."));
  std::string ps_name = std::string(ps_filename);
  ps_name = ps_name.substr(ps_name.find_last_of("/") + 1, ps_name.length() - 1);
  ps_name = ps_name.substr(0, ps_name.find_last_of("."));
  std::string output_name = ts_name + "_" + ss_name + "_"  + ps_name + "_"
                                  + std::to_string(static_cast<int>(relative_alpha)) + "_"
                                  + std::to_string(static_cast<int>(relative_offset)) + ".off";
  std::cout << "Writing to " << output_name << std::endl;
  CGAL::IO::write_polygon_mesh(output_name, wrap, CGAL::parameters::stream_precision(17));

  return EXIT_SUCCESS;
}
