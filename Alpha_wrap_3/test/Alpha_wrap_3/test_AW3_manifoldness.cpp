#define CGAL_AW3_TIMER
#define CGAL_AW3_DEBUG
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
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Random.h>

using namespace CGAL::Alpha_wraps_3::internal;

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT = Kernel::FT;
using Point_3 = Kernel::Point_3;
using Vector_3 = Kernel::Vector_3;

using Mesh = CGAL::Surface_mesh<Point_3>;

void alpha_wrap_triangle_manifoldness(Mesh& input_mesh,
                                      const double alpha,
                                      const double offset)
{
  namespace AW3 = CGAL::Alpha_wraps_3;
  namespace PMP = CGAL::Polygon_mesh_processing;

  std::cout << "Input: " << num_vertices(input_mesh) << " vertices, " << num_faces(input_mesh) << " faces" << std::endl;
  std::cout << "Alpha: " << alpha << " Offset: " << offset << std::endl;

  const bool has_degeneracies = !PMP::remove_degenerate_faces(input_mesh);
  if(has_degeneracies)
    std::cerr << "Warning: Failed to remove some degenerate faces." << std::endl;

  std::cout << "Processed input: " << vertices(input_mesh).size() << " vertices, " << faces(input_mesh).size() << " faces" << std::endl;

//  CGAL::IO::write_polygon_mesh("input.off", input_mesh, CGAL::parameters::stream_precision(17));

  Mesh nm_wrap;
  CGAL::alpha_wrap_3(input_mesh, alpha, offset, nm_wrap,
                     CGAL::parameters::do_enforce_manifoldness(false));

  std::cout << "Result: " << vertices(nm_wrap).size() << " vertices, " << faces(nm_wrap).size() << " faces" << std::endl;

  if(!has_degeneracies)
  {
    assert(AW3::internal::is_valid_wrap(nm_wrap, false));
    assert(AW3::internal::is_outer_wrap_of_triangle_mesh(nm_wrap, input_mesh));
    assert(AW3::internal::has_expected_Hausdorff_distance(nm_wrap, input_mesh, alpha, offset));
  }

  assert(AW3::internal::check_edge_length(nm_wrap, alpha));

  FT base_vol = 0;
  if(!is_closed(nm_wrap))
    std::cerr << "W: non-manifold wrap is not closed" << std::endl;
  else
    base_vol = PMP::volume(nm_wrap);

  Mesh m_wrap;
  CGAL::alpha_wrap_3(input_mesh, alpha, offset, m_wrap,
                     CGAL::parameters::do_enforce_manifoldness(true));

//  CGAL::IO::write_polygon_mesh("last.off", wrap, CGAL::parameters::stream_precision(17));

  if(!has_degeneracies)
  {
    assert(AW3::internal::is_valid_wrap(m_wrap, true /*manifoldness*/));
    assert(AW3::internal::is_outer_wrap_of_triangle_mesh(m_wrap, input_mesh));
  }

  const FT final_vol = PMP::volume(m_wrap);

  if(base_vol != 0)
  {
    const FT ratio = final_vol / base_vol;

    std::cout << "Volumes post-manifoldness fix:\n"
              << "before: " << base_vol << "\n"
              << "after:  " << final_vol << "\n"
              << "ratio:  " << ratio << std::endl;
    if(ratio > 1.1) // more than 10% extra volume
      std::cerr << "W: large increase of volume after manifoldness resolution" << std::endl;
  }
}

void alpha_wrap_triangle_manifoldness(const std::string& filename,
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

  alpha_wrap_triangle_manifoldness(input_mesh, alpha, offset);
}

void alpha_wrap_triangle_manifoldness(const std::string& filename)
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

    alpha_wrap_triangle_manifoldness(input_mesh, alpha, offset);
  }
}

int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  if(argc > 1)
  {
    const double relative_alpha = (argc > 2) ? std::stod(argv[2]) : 20.;
    const double relative_offset = (argc > 3) ? std::stod(argv[3]) : 600.;

    alpha_wrap_triangle_manifoldness(argv[1], relative_alpha, relative_offset);
    return EXIT_SUCCESS;
  }

  alpha_wrap_triangle_manifoldness("data/tetrahedron.off");
  alpha_wrap_triangle_manifoldness("data/plane.off");
  alpha_wrap_triangle_manifoldness("data/open_tetrahedron.off");
  alpha_wrap_triangle_manifoldness("data/tetrahedron_open_tip.off");
  alpha_wrap_triangle_manifoldness("data/sphere_one_hole.off");
  alpha_wrap_triangle_manifoldness("data/tetrahedron_duplicates.off");
  alpha_wrap_triangle_manifoldness("data/tetrahedron_degenerencies.off");
  alpha_wrap_triangle_manifoldness("data/non_manifold.off");
  alpha_wrap_triangle_manifoldness("data/combinatorial_manifold.off");
  alpha_wrap_triangle_manifoldness("data/combinatorial_manifold_multiple_components.off");
  alpha_wrap_triangle_manifoldness("data/tetrahedron_self_intersection_tip.off");
  alpha_wrap_triangle_manifoldness("data/tetrahedron_twisted_tip.off");
  alpha_wrap_triangle_manifoldness("data/tetrahedron_random_perturbation.off");
  alpha_wrap_triangle_manifoldness("data/overlay_triangulation.off");
  alpha_wrap_triangle_manifoldness("data/two_knives.off");
  alpha_wrap_triangle_manifoldness("data/three_knives.off");
  alpha_wrap_triangle_manifoldness("data/bunny_random_perturbation.off");

  std::cout << "Done!" << std::endl;
  return EXIT_SUCCESS;
}
