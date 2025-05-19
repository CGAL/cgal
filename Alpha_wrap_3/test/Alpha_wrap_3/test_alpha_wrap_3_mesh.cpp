#define CGAL_AW3_TIMER
//#define CGAL_AW3_DEBUG
#define CGAL_AW3_DEBUG_MANIFOLDNESS
//#define CGAL_AW3_DEBUG_STEINER_COMPUTATION
//#define CGAL_AW3_DEBUG_INITIALIZATION
//#define CGAL_AW3_DEBUG_QUEUE

#include <CGAL/alpha_wrap_3.h>
#include <CGAL/Alpha_wrap_3/internal/validation.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/repair_degeneracies.h>
#include <CGAL/Random.h>

using namespace CGAL::Alpha_wraps_3::internal;

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT = Kernel::FT;
using Point_3 = Kernel::Point_3;
using Vector_3 = Kernel::Vector_3;

using Mesh = CGAL::Surface_mesh<Point_3>;

void alpha_wrap_triangle_mesh(Mesh& input_mesh,
                              const double alpha,
                              const double offset)
{
  namespace AW3 = CGAL::Alpha_wraps_3;
  namespace PMP = CGAL::Polygon_mesh_processing;

  std::cout << "Input: " << num_vertices(input_mesh) << " vertices, " << num_faces(input_mesh) << " faces" << std::endl;
  std::cout << "Alpha: " << alpha << " Offset: " << offset << std::endl;

  bool has_degeneracies = !PMP::remove_degenerate_faces(input_mesh);
  if(has_degeneracies)
    std::cerr << "Warning: Failed to remove some degenerate faces." << std::endl;

  std::cout << "Processed input: " << vertices(input_mesh).size() << " vertices, " << faces(input_mesh).size() << " faces" << std::endl;

//  CGAL::IO::write_polygon_mesh("input.off", input_mesh, CGAL::parameters::stream_precision(17));

  const bool enforce_manifoldness = true;

  Mesh wrap;
  CGAL::alpha_wrap_3(input_mesh, alpha, offset, wrap,
                     CGAL::parameters::do_enforce_manifoldness(enforce_manifoldness));

  std::cout << "Result: " << vertices(wrap).size() << " vertices, " << faces(wrap).size() << " faces" << std::endl;

//  CGAL::IO::write_polygon_mesh("last.off", wrap, CGAL::parameters::stream_precision(17));

  if(!has_degeneracies)
  {
    assert(AW3::internal::is_valid_wrap(wrap, enforce_manifoldness));
    assert(AW3::internal::is_outer_wrap_of_triangle_mesh(wrap, input_mesh));

    if(!enforce_manifoldness)
      assert(AW3::internal::has_expected_Hausdorff_distance(wrap, input_mesh, alpha, offset));
  }

  if(!enforce_manifoldness)
    assert(AW3::internal::check_edge_length(wrap, alpha));
}

void alpha_wrap_triangle_mesh(const std::string& filename,
                              const double alpha_rel,
                              const double offset_rel)
{
  Mesh input_mesh;
  bool res = CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(filename, input_mesh);
  assert(res);
  assert(!is_empty(input_mesh) && is_triangle_mesh(input_mesh));

  CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(input_mesh);
  const Vector_3 longest_diag = Point_3(bbox.xmax(), bbox.ymax(), bbox.zmax()) -
                                Point_3(bbox.xmin(), bbox.ymin(), bbox.zmin());
  double longest_diag_length = CGAL::to_double(CGAL::approximate_sqrt(longest_diag.squared_length()));

  const double alpha = longest_diag_length / alpha_rel;
  const double offset = longest_diag_length / offset_rel;

  alpha_wrap_triangle_mesh(input_mesh, alpha, offset);
}

void alpha_wrap_triangle_mesh(const std::string& filename)
{
  CGAL::Random r;

  Mesh input_mesh;
  bool res = CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(filename, input_mesh);
  assert(res);
  assert(!is_empty(input_mesh) && is_triangle_mesh(input_mesh));

  CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(input_mesh);
  const Vector_3 longest_diag = Point_3(bbox.xmax(), bbox.ymax(), bbox.zmax()) -
                                Point_3(bbox.xmin(), bbox.ymin(), bbox.zmin());
  double longest_diag_length = CGAL::to_double(CGAL::approximate_sqrt(longest_diag.squared_length()));

  for(int i=0; i<2; ++i)
  {
    const double alpha_expo = r.get_double(0., 6); // to have alpha_rel between 1 and 64
    const double offset_expo = r.get_double(0., 6);
    const double alpha_rel = std::pow(2, alpha_expo);
    const double offset_rel = std::pow(2, offset_expo);
    const double alpha = longest_diag_length / alpha_rel;
    const double offset = longest_diag_length / offset_rel;

    std::cout << "===================================================" << std::endl;
    std::cout << filename << " " << alpha << " (rel " << alpha_rel << ")"
                          << " " << offset << " (rel " << offset_rel << ")" << std::endl;
    std::cout << "Random seed = " << r.get_seed() << std::endl;

    alpha_wrap_triangle_mesh(input_mesh, alpha, offset);
  }
}

int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  // For convenience to do manual testing
  if(argc > 1)
  {
    const double relative_alpha = (argc > 2) ? std::stod(argv[2]) : 20.;
    const double relative_offset = (argc > 3) ? std::stod(argv[3]) : 600.;

    alpha_wrap_triangle_mesh(argv[1], relative_alpha, relative_offset);
    return EXIT_SUCCESS;
  }

  alpha_wrap_triangle_mesh("data/tetrahedron.off");
  alpha_wrap_triangle_mesh("data/plane.off");
  alpha_wrap_triangle_mesh("data/open_tetrahedron.off");
  alpha_wrap_triangle_mesh("data/tetrahedron_open_tip.off");
  alpha_wrap_triangle_mesh("data/sphere_one_hole.off");
  alpha_wrap_triangle_mesh("data/tetrahedron_duplicates.off");
  alpha_wrap_triangle_mesh("data/tetrahedron_degenerencies.off");
  alpha_wrap_triangle_mesh("data/non_manifold.off");
  alpha_wrap_triangle_mesh("data/combinatorial_manifold.off");
  alpha_wrap_triangle_mesh("data/combinatorial_manifold_multiple_components.off");
  alpha_wrap_triangle_mesh("data/tetrahedron_self_intersection_tip.off");
  alpha_wrap_triangle_mesh("data/tetrahedron_twisted_tip.off");
  alpha_wrap_triangle_mesh("data/tetrahedron_random_perturbation.off");
  alpha_wrap_triangle_mesh("data/overlay_triangulation.off");
  alpha_wrap_triangle_mesh("data/two_knives.off");
  alpha_wrap_triangle_mesh("data/three_knives.off");
  alpha_wrap_triangle_mesh("data/bunny_random_perturbation.off");

  std::cout << "Done!" << std::endl;
  return EXIT_SUCCESS;
}
