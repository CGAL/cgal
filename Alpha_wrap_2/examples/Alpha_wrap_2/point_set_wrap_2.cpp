#define CGAL_AW2_DEBUG
#define CGAL_AW2_DEBUG_DUMP_INTERMEDIATE_WRAPS
#define CGAL_AW2_DEBUG_QUEUE
#define CGAL_AW2_DEBUG_DUMP_EVERY_STEP
// #define CGAL_AW2_DEBUG_QUEUE_PP
#define CGAL_AW2_DEBUG_STEINER_COMPUTATION
// #define CGAL_AW2_DEBUG_SPHERE_MARCHING
// #define CGAL_AW2_USE_SORTED_PRIORITY_QUEUE
#define CGAL_AW2_DEBUG_EDGE_STATUS
#define CGAL_AW2_DEBUG_TRAVERSABILITY

#include "output_helper.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/alpha_wrap_2.h>
#include <CGAL/Polygon_2.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/Real_timer.h>

#include <iostream>
#include <string>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = K::Point_2;
using Point_3 = K::Point_3;
using Vector_2 = K::Vector_2;
using Vector_3 = K::Vector_3;

using Point_set_2 = CGAL::Point_set_3<Point_2, Vector_2>;
using Point_set_3 = CGAL::Point_set_3<Point_3, Vector_3>;

using Polygon_2 = CGAL::Polygon_2<K>;

int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  const std::string filename = argc > 1 ? argv[1] : CGAL::data_file_path("points_3/circles.ply");

  // This code reads a _3D point file_
  Point_set_3 point_set_3;
  if(!CGAL::IO::read_points(filename,
                            point_set_3.index_back_inserter(),
                            CGAL::parameters::point_map(point_set_3.point_push_map())))
  {
  std::cerr << "Can't read input file " << filename << std::endl;
  return EXIT_FAILURE;
  }

  // Project onto the xy plane
  Point_set_2 point_set_2;
  for (const auto& p3 : point_set_3.points())
    point_set_2.insert(Point_2(p3.x(), p3.y()));

  std::cout << point_set_2.size() << " points" << std::endl;

  // Compute the alpha and offset values
  const double relative_alpha = (argc > 2) ? std::stod(argv[2]) : 10.;
  const double relative_offset = (argc > 3) ? std::stod(argv[3]) : 300.;

  CGAL::Bbox_2 bbox;
  for (const auto& p2 : point_set_2.points())
    bbox += p2.bbox();

  std::ofstream out("projected.xyz");
  out.precision(17);
  for (const auto& p2 : point_set_2.points())
    out << p2 << " 0\n";
  out.close();

  const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                       CGAL::square(bbox.ymax() - bbox.ymin()));
  const double alpha = diag_length / relative_alpha;
  const double offset = diag_length / relative_offset;
  std::cout << "absolute alpha = " << alpha << " absolute offset = " << offset << std::endl;

  // Construct the wrap
  CGAL::Real_timer t;
  t.start();

  std::vector<Polygon_2> wrap;
  CGAL::alpha_wrap_2(point_set_2.points(), alpha, offset, wrap);

  t.stop();
  std::cout << "Result: " << wrap.size() << " polygons" << std::endl;
  std::cout << "Took " << t.time() << " s." << std::endl;

  // Save the result
  const std::string output_name = generate_output_name(filename, relative_alpha, relative_offset);
  std::cout << "Writing to " << output_name << std::endl;
  write_polygons(output_name, wrap);

  std::cout << "Done." << std::endl;

  return EXIT_SUCCESS;
}
