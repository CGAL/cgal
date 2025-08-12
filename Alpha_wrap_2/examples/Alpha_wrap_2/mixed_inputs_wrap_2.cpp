
#define CGAL_AW2_DEBUG
#define CGAL_AW2_DEBUG_DUMP_INTERMEDIATE_WRAPS
#define CGAL_AW2_DEBUG_QUEUE
#define CGAL_AW2_DEBUG_DUMP_EVERY_STEP
// #define CGAL_AW2_DEBUG_QUEUE_PP
#define CGAL_AW2_DEBUG_STEINER_COMPUTATION
// #define CGAL_AW2_DEBUG_SPHERE_MARCHING
#define CGAL_AW2_USE_SORTED_PRIORITY_QUEUE
#define CGAL_AW2_DEBUG_EDGE_STATUS
#define CGAL_AW2_DEBUG_TRAVERSABILITY

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/alpha_wrap_2.h>

#include <CGAL/point_generators_2.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Polygon_2.h>

#include <iostream>
#include <string>
#include <fstream>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = K::Point_2;
using Segment_2 = K::Segment_2;
using Triangle_2 = K::Triangle_2;

using Segments = std::vector<Segment_2>;
using Points = std::vector<Point_2>;
using Triangles = std::vector<Triangle_2>;

using Polygon_2 = CGAL::Polygon_2<K>;

int main(int argc, char** argv)
{
  // Generate random points in a unit square
  const int num_points = (argc > 1) ? std::stoi(argv[1]) : 20;
  const int num_segments = (argc > 2) ? std::stoi(argv[2]) : 5;
  const int num_triangles = (argc > 3) ? std::stoi(argv[3]) : 5;

  // Generate random points
  Points points;
  CGAL::Random_points_in_square_2<Point_2> point_gen(0.5); // radius 0.5 for unit square centered at origin
  std::copy_n(point_gen, num_points, std::back_inserter(points));
  std::cout << points.size() << " random points generated" << std::endl;

  // Generate random segments between random points in the square
  Segments segments;
  CGAL::Random_points_in_square_2<Point_2> segment_point_gen(0.5);
  for(int i = 0; i < num_segments; ++i)
  {
    Point_2 p = *segment_point_gen++;
    Point_2 q = *segment_point_gen++;
    segments.emplace_back(p, q);
  }
  std::cout << segments.size() << " random segments generated" << std::endl;

  // Generate random triangles in the square
  Triangles triangles;
  CGAL::Random_points_in_square_2<Point_2> triangle_point_gen(0.5);
  for(int i = 0; i < num_triangles; ++i)
  {
    Point_2 p = *triangle_point_gen++;
    Point_2 q = *triangle_point_gen++;
    Point_2 r = *triangle_point_gen++;
    triangles.emplace_back(p, q, r);
  }
  std::cout << triangles.size() << " random triangles generated" << std::endl;

  std::ofstream out_points("mixed_input_points.xyz");
  out_points.precision(17);
  for(const auto& p : points)
    out_points << p.x() << " " << p.y() << " 0\n";
  out_points.close();

  std::ofstream out_segments("mixed_input_segments.txt");
  out_segments.precision(17);
  for(const auto& s : segments)
    out_segments << "2 " << s.source().x() << " " << s.source().y() << " 0 "
                         << s.target().x() << " " << s.target().y() << " 0\n";
  out_segments.close();

  std::ofstream out_triangles("mixed_input_triangles.txt");
  out_triangles.precision(17);
  for(const auto& t : triangles)
    out_triangles << "4 " << t.vertex(0).x() << " " << t.vertex(0).y() << " 0 "
                          << t.vertex(1).x() << " " << t.vertex(1).y() << " 0 "
                          << t.vertex(2).x() << " " << t.vertex(2).y() << " 0 "
                          << t.vertex(0).x() << " " << t.vertex(0).y() << " 0\n";
  out_triangles.close();

  const double relative_alpha = (argc > 4) ? std::stod(argv[4]) : 15.;
  const double relative_offset = (argc > 5) ? std::stod(argv[5]) : 450.;

  // Since we're using a unit square, diagonal length is sqrt(2)
  const double diag_length = std::sqrt(2.0);
  const double alpha = diag_length / relative_alpha;
  const double offset = diag_length / relative_offset;

  CGAL::Real_timer t;
  t.start();

  using Triangle_Oracle = CGAL::Alpha_wraps_2::internal::Triangle_soup_oracle<K>;
  using Segment_Oracle = CGAL::Alpha_wraps_2::internal::Segment_soup_oracle<K, Triangle_Oracle>;
  using Oracle = CGAL::Alpha_wraps_2::internal::Point_set_oracle<K, Segment_Oracle>;

  Triangle_Oracle triangle_oracle(K{});
  Segment_Oracle segment_oracle(triangle_oracle);
  Oracle oracle(segment_oracle);

  oracle.add_point_set(points, CGAL::parameters::default_values());
  oracle.add_segment_soup(segments, CGAL::parameters::default_values());
  oracle.add_triangle_soup(triangles, CGAL::parameters::default_values());

  CGAL::Alpha_wraps_2::internal::Alpha_wrapper_2<Oracle> aw2(oracle);

  std::vector<Polygon_2> wrap;
  aw2(alpha, offset, wrap);

  t.stop();
  std::cout << "Took " << t.time() << " seconds" << std::endl;

  std::string output_name = "mixed_wrap_" +
                           std::to_string(num_points) + "_" +
                           std::to_string(num_segments) + "_" +
                           std::to_string(num_triangles) + "_" +
                           std::to_string(static_cast<int>(relative_alpha)) + "_" +
                           std::to_string(static_cast<int>(relative_offset)) + ".txt";
  std::cout << "Writing to " << output_name << std::endl;

  std::ofstream out(output_name);
  out.precision(17);
  for(const auto& p : wrap)
    out << p << std::endl;

  return EXIT_SUCCESS;
}
