#define CGAL_AW3_TIMER

#include "output_helper.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/alpha_wrap_3.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Real_timer.h>

#include <iostream>
#include <string>

namespace PMP = CGAL::Polygon_mesh_processing;

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT = K::FT;
using Point_3 = K::Point_3;

using Mesh = CGAL::Surface_mesh<Point_3>;

// We want decreasing alphas, and these are relative ratios, so they need to be increasing
const std::vector<FT> relative_alphas = { 1, 2/*50, 100, 150, 200, 250*/ };
const FT relative_offset = 600;

int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  // Read the input
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/cube.off");
  std::cout << "Reading " << filename << "..." << std::endl;

  Mesh mesh;
  if(!PMP::IO::read_polygon_mesh(filename, mesh) || is_empty(mesh) || !is_triangle_mesh(mesh))
  {
    std::cerr << "Invalid input:" << filename << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Input: " << num_vertices(mesh) << " vertices, " << num_faces(mesh) << " faces" << std::endl;

  const CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh);
  const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                       CGAL::square(bbox.ymax() - bbox.ymin()) +
                                       CGAL::square(bbox.zmax() - bbox.zmin()));

  // ===============================================================================================
  // Naive approach:

  CGAL::Real_timer t;
  double total_time = 0.;

  for(std::size_t i=0; i<relative_alphas.size(); ++i)
  {
    t.reset();
    t.start();

    const double alpha = diag_length / relative_alphas[i];
    const double offset = diag_length / relative_offset;
    std::cout << ">>> [" << i << "] alpha: " << alpha << " offset: " << offset << std::endl;

    Mesh wrap;
    CGAL::alpha_wrap_3(mesh, alpha, offset, wrap,
                       CGAL::parameters::do_enforce_manifoldness(false));

    t.stop();
    std::cout << "  Result: " << num_vertices(wrap) << " vertices, " << num_faces(wrap) << " faces" << std::endl;
    std::cout << "  Elapsed time: " << t.time() << " s." << std::endl;

    const std::string output_name = generate_output_name(filename, relative_alphas[i], relative_offset);
    std::cout << "Writing to " << output_name << std::endl;
    CGAL::IO::write_polygon_mesh(output_name, wrap, CGAL::parameters::stream_precision(17));

    total_time += t.time();
  }

  std::cout << "Total elapsed time (naive): " << total_time << " s.\n" << std::endl;

  // ===============================================================================================
  // Re-use approach
  //
  // Here, we restart from the triangulation of the previous state, and carve according
  // to a (smaller) alpha value. This enables considerable speed-up: the cumulated time taken
  // to run `n` successive instances of `{alpha_wrap(alpha_i)}_(i=1...n)` will be equal to the
  // time taken to run alpha_wrap(alpha_n) from scratch.
  //
  // For example:
  // naive:
  //   alpha_wrap(alpha_1, ...) ---> 2s
  //   alpha_wrap(alpha_2, ...) ---> 4s
  //   alpha_wrap(alpha_3, ...) ---> 8s
  // will become with reusability:
  //   alpha_wrap(alpha_1, ..., reuse) ---> 2s
  //   alpha_wrap(alpha_2, ..., reuse) ---> 2s // 2+2 = 4s = naive alpha_2
  //   alpha_wrap(alpha_3, ..., reuse) ---> 4s // 2+2+4 = 8s = naive alpha_3
  // Thus, if we care about the intermediate results, we save 6s (8s instead of 14s).
  // The speed-up increases with the number of intermediate results, and if the alpha values
  // are close.
  //
  // !! Warning !!
  // The result of alpha_wrap(alpha_1, ...) followed by alpha_wrap(alpha_2, ...) with alpha_2
  // smaller than alpha_1 is going to be close but NOT exactly equal to that produced by calling
  // alpha_wrap(alpha_2, ...) directly.

  total_time = 0.;
  t.reset();

  using Oracle = CGAL::Alpha_wraps_3::internal::Triangle_mesh_oracle<K>;
  using Wrapper = CGAL::Alpha_wraps_3::internal::Alpha_wrapper_3<Oracle>;
  Wrapper wrapper; // contains the triangulation that is being refined iteratively

  for(std::size_t i=0; i<relative_alphas.size(); ++i)
  {
    t.reset();
    t.start();

    const double alpha = diag_length / relative_alphas[i];
    const double offset = diag_length / relative_offset;
    std::cout << ">>> [" << i << "] alpha: " << alpha << " offset: " << offset << std::endl;

    // The triangle mesh oracle can be initialized with alpha to internally perform a split
    // of too-big facets while building the AABB Tree. This split in fact yields a significant
    // speed-up for meshes with elements that are large compared to alpha. This speed-up makes it
    // faster to re-build the AABB tree for every value of alpha than to use a non-optimized tree.
    Oracle oracle(alpha);
    oracle.add_triangle_mesh(mesh, CGAL::parameters::default_values());
    wrapper.oracle() = oracle;

    Mesh wrap;
    wrapper(alpha, offset, wrap,
            CGAL::parameters::do_enforce_manifoldness(false)
                             .refine_triangulation((i != 0)));

    t.stop();
    std::cout << "  Result: " << num_vertices(wrap) << " vertices, " << num_faces(wrap) << " faces" << std::endl;
    std::cout << "  Elapsed time: " << t.time() << " s." << std::endl;

    total_time += t.time();
  }

  std::cout << "Total elapsed time (successive): " << total_time << " s." << std::endl;

  return EXIT_SUCCESS;
}
