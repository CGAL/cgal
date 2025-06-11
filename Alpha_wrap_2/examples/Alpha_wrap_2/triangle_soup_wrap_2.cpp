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

#include "output_helper.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/alpha_wrap_2.h>
#include <CGAL/Polygon_2.h>

#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Real_timer.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <array>
#include <iostream>
#include <string>
#include <vector>

namespace AW3 = CGAL::Alpha_wraps_2;

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = K::Point_2;
using Point_3 = K::Point_3;

using Polygon_2 = CGAL::Polygon_2<K>;

int main(int argc, char** argv)
{
  std::cout.precision(17);

  // Read the input
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/blobby-shuffled.off");
  std::cout << "Reading " << filename << "..." << std::endl;

  std::vector<Point_3> points;
  std::vector<std::array<std::size_t, 3> > faces;
  if(!CGAL::IO::read_polygon_soup(filename, points, faces) || faces.empty())
  {
    std::cerr << "Invalid input:" << filename << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Input: " << points.size() << " points, " << faces.size() << " faces" << std::endl;

  // Project onto the XY plane
  std::vector<Point_2> points_2(points.size());
  for(std::size_t i = 0; i < points.size(); ++i)
  points_2[i] = Point_2(points[i].x(), points[i].y());

  std::vector<Point_3> points_3(points_2.size());
  for(std::size_t i = 0; i < points_2.size(); ++i)
    points_3[i] = Point_3(points_2[i].x(), points_2[i].y(), 0.);
  CGAL::IO::write_polygon_soup("projected.off", points_3, faces, CGAL::parameters::stream_precision(17));

  // Compute the alpha and offset values
  const double relative_alpha = (argc > 2) ? std::stod(argv[2]) : 20.;
  const double relative_offset = (argc > 3) ? std::stod(argv[3]) : 600.;

  CGAL::Bbox_2 bbox;
  for(const Point_2& p : points_2)
    bbox += p.bbox();

  const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                       CGAL::square(bbox.ymax() - bbox.ymin()));

  const double alpha = diag_length / relative_alpha;
  const double offset = diag_length / relative_offset;

  // Construct the wrap
  CGAL::Real_timer t;
  t.start();

  std::vector<Polygon_2> wrap;
  CGAL::alpha_wrap_2(points_2, faces, alpha, offset, wrap);

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
