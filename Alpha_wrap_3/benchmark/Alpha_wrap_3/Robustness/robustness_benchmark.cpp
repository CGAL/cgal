#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/alpha_wrap_3.h>
#include <CGAL/Alpha_wrap_3/internal/validation.h>

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = Kernel::Point_3;
using Mesh = CGAL::Surface_mesh<Point_3>;

namespace CGAL {
namespace Alpha_wraps_3 {
namespace internal {
namespace {

enum Robustness_benchmark_exit_code
{
  // Success
  VALID_SOLID_OUTPUT = 0,

  // Failure
  INPUT_IS_INVALID = 1,
  OUTPUT_IS_NOT_TRIANGLE_MESH = 2,
  OUTPUT_IS_COMBINATORIAL_NON_MANIFOLD = 3,
  OUTPUT_HAS_BORDERS = 4,
  OUTPUT_HAS_DEGENERATED_FACES = 5,
  OUTPUT_HAS_GEOMETRIC_SELF_INTERSECTIONS = 6,
  OUTPUT_DOES_NOT_BOUND_VOLUME = 7,
  OUTPUT_DOES_NOT_CONTAIN_INPUT = 8,
  OUTPUT_DISTANCE_IS_TOO_LARGE = 9,
};

} // namespace
} // namespace internal
} // namespace Alpha_wraps_3
} // namespace CGAL

namespace PMP = CGAL::Polygon_mesh_processing;
namespace AW3i = CGAL::Alpha_wraps_3::internal;

int main(int argc, char** argv)
{
  const int argc_check = argc - 1;
  char* entry_name_ptr = nullptr;
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
    return AW3i::INPUT_IS_INVALID;

  Mesh input_mesh;
  if(!PMP::IO::read_polygon_mesh(entry_name_ptr, input_mesh) ||
     is_empty(input_mesh) ||
     !is_triangle_mesh(input_mesh)
#ifndef CGAL_ALPHA_WRAP_3_TOLERATE_DEGENERACIES
     || AW3i::has_degenerated_faces(input_mesh)
#endif
     )
  {
    return AW3i::INPUT_IS_INVALID;
  }

  const CGAL::Bbox_3 bbox = PMP::bbox(input_mesh);
  const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                       CGAL::square(bbox.ymax() - bbox.ymin()) +
                                       CGAL::square(bbox.zmax() - bbox.zmin()));
  const double alpha = diag_length / relative_alpha_ratio;
  const double offset = diag_length / relative_offset_ratio;

  Mesh wrap;
  alpha_wrap_3(input_mesh, alpha, offset, wrap);

  if(!is_triangle_mesh(wrap))
    return AW3i::OUTPUT_IS_NOT_TRIANGLE_MESH;

  if(!is_closed(wrap))
    return AW3i::OUTPUT_HAS_BORDERS;

  if(AW3i::has_degenerated_faces(wrap))
    return AW3i::OUTPUT_HAS_DEGENERATED_FACES;

  if(AW3i::is_combinatorially_non_manifold(wrap))
    return AW3i::OUTPUT_IS_COMBINATORIAL_NON_MANIFOLD;

  if(PMP::does_self_intersect(wrap))
    return AW3i::OUTPUT_HAS_GEOMETRIC_SELF_INTERSECTIONS;

  if(!PMP::does_bound_a_volume(wrap))
    return AW3i::OUTPUT_DOES_NOT_BOUND_VOLUME;

  if(!AW3i::is_outer_wrap_of_triangle_mesh(wrap, input_mesh))
    return AW3i::OUTPUT_DOES_NOT_CONTAIN_INPUT;

  if(!AW3i::has_expected_Hausdorff_distance(wrap, input_mesh, alpha, offset))
     return AW3i::OUTPUT_DISTANCE_IS_TOO_LARGE;

  return AW3i::VALID_SOLID_OUTPUT;
}
