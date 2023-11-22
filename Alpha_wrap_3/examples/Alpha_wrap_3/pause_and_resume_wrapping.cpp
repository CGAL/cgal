// This example demonstrates how to interrupt the wrapping process before it has terminated,
// and how to resume afterwards.
//
// -------------------------------- !! Warning !! --------------------------------------------------
// By default, the wrapper uses an unsorted LIFO queue of faces to refine. This means that
// the intermediate result is not very useful because the algorithm carves deep and not wide
// (somewhat like a DFS vs a BFS).
//
// The sorted queue option is enabled with the macro below to make the refinement algorithm
// more uniform. The downside is that it is slower.
// -------------------------------------------------------------------------------------------------
#define CGAL_AW3_USE_SORTED_PRIORITY_QUEUE

#include "output_helper.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/alpha_wrap_3.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Random.h>
#include <CGAL/Timer.h>

#include <iostream>
#include <string>

namespace AW3 = CGAL::Alpha_wraps_3;
namespace PMP = CGAL::Polygon_mesh_processing;

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;

using Points = std::vector<Point_3>;
using Face = std::array<std::size_t, 3>;
using Faces = std::vector<Face>;

using Mesh = CGAL::Surface_mesh<Point_3>;
using face_descriptor = boost::graph_traits<Mesh>::face_descriptor;

struct Interrupter_visitor
  : public AW3::internal::Wrapping_default_visitor
{
  using Base = AW3::internal::Wrapping_default_visitor;

  CGAL::Real_timer timer;
  double max_time = -1; // in seconds

public:
  void set_max_time(double t) { max_time = t; }

public:
  template <typename AlphaWrapper>
  void on_flood_fill_begin(const AlphaWrapper&)
  {
    std::cout << "Starting timer..." << std::endl;
    timer.start();
  }

  template <typename Wrapper>
  bool go_further(const Wrapper&)
  {
    if(timer.time() > max_time)
    {
      timer.stop();
      std::cout << "Paused after " << timer.time() << " s." << std::endl;
      return false;
    }

    return true;
  }
};

int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  CGAL::Random rng;
  std::cout << "Random seed = " << rng.get_seed() << std::endl;

  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/armadillo.off");

  // = read the soup
  Points points;
  Faces faces;
  if(!CGAL::IO::read_polygon_soup(filename, points, faces) || faces.empty())
  {
    std::cerr << "Invalid soup input: " << filename << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Input: " << points.size() << " points, " << faces.size() << " faces" << std::endl;

  // Compute the alpha and offset values
  const double relative_alpha = (argc > 2) ? std::stod(argv[2]) : rng.get_double(150., 200.);
  const double relative_offset = (argc > 3) ? std::stod(argv[3]) : 600.;
  std::cout << "relative_alpha = " << relative_alpha << std::endl;

  CGAL::Bbox_3 bbox;
  for(const Point_3& p : points)
    bbox += p.bbox();

  const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                       CGAL::square(bbox.ymax() - bbox.ymin()) +
                                       CGAL::square(bbox.zmax() - bbox.zmin()));

  const double alpha = diag_length / relative_alpha;
  const double offset = diag_length / relative_offset;

  // Build the wrapper
  using Oracle = CGAL::Alpha_wraps_3::internal::Triangle_soup_oracle<K>;
  Oracle oracle(alpha);
  oracle.add_triangle_soup(points, faces, CGAL::parameters::default_values());
  CGAL::Alpha_wraps_3::internal::Alpha_wrapper_3<Oracle> aw3(oracle);

  // --- Launch the wrapping, and pause when the algorithm has spent 1s flooding
  Interrupter_visitor interrupter;
  interrupter.set_max_time(1.);

  Mesh wrap;
  aw3(alpha, offset, wrap, CGAL::parameters::visitor(interrupter));
  std::cout << ">>> The current wrap has " << num_vertices(wrap) << " vertices" << std::endl;
  CGAL::IO::write_polygon_mesh("stopped_1.off", wrap, CGAL::parameters::stream_precision(17));

  // --- Restart from the previous state, and pause a bit further
  interrupter.set_max_time(2.);
  aw3(alpha, offset, wrap, CGAL::parameters::visitor(interrupter)
                                            .refine_triangulation(true));
  std::cout << ">>> The current wrap has " << num_vertices(wrap) << " vertices" << std::endl;
  CGAL::IO::write_polygon_mesh("stopped_2.off", wrap, CGAL::parameters::stream_precision(17));

  // --- Restart from the previous state, and let it finish
  aw3(alpha, offset, wrap, CGAL::parameters::refine_triangulation(true));
  std::cout << ">>> The final (resumed) wrap has " << num_vertices(wrap) << " vertices" << std::endl;
  std::string output_name = generate_output_name(filename, relative_alpha, relative_offset);
  std::cout << "Writing to " << "resumed_" + output_name << std::endl;
  CGAL::IO::write_polygon_mesh("resumed_" + output_name, wrap, CGAL::parameters::stream_precision(17));

  // --- Get the final wrap, in one go:
  Mesh single_pass_wrap;
  CGAL::alpha_wrap_3(points, faces, alpha, offset, single_pass_wrap);
  std::cout << ">>> The final (from scratch) wrap has " << num_vertices(single_pass_wrap) << " vertices" << std::endl;

  output_name = generate_output_name(filename, relative_alpha, relative_offset);
  std::cout << "Writing to " << output_name << std::endl;
  CGAL::IO::write_polygon_mesh(output_name, single_pass_wrap, CGAL::parameters::stream_precision(17));

  // --- Compare the results to ensure both approaches yield identical meshes
  std::vector<std::pair<face_descriptor, face_descriptor> > common;
  std::vector<face_descriptor> m1_only;
  std::vector<face_descriptor> m2_only;
  PMP::match_faces(wrap, single_pass_wrap,
                   std::back_inserter(common),
                   std::back_inserter(m1_only),
                   std::back_inserter(m2_only));
  if(!m1_only.empty() || !m2_only.empty())
  {
    std::cerr << "Error: The two wraps should have been identical!" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
