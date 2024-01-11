#define CGAL_AW3_TIMER
#define CGAL_AW3_DEBUG
#define CGAL_AW3_DEBUG_MANIFOLDNESS
// #define CGAL_AW3_DEBUG_INITIALIZATION

#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/alpha_wrap_3.h>
#include <CGAL/Alpha_wrap_3/internal/validation.h>

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

using namespace CGAL::Alpha_wraps_3::internal;

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT = Kernel::FT;
using Point_3 = Kernel::Point_3;
using Vector_3 = Kernel::Vector_3;

using Mesh = CGAL::Surface_mesh<Point_3>;

using Seeds = std::vector<Point_3>;

template <typename Oracle>
void generate_random_seeds(const Oracle& oracle,
                           const double offset,
                           Seeds& seeds,
                           CGAL::Random& r)
{
  const auto bbox = CGAL::Alpha_wraps_3::internal::Alpha_wrapper_3<Oracle>(oracle).construct_bbox(offset);
  const double sq_offset = CGAL::square(offset);

  while(seeds.size() < 3)
  {
    const Point_3 seed (bbox.xmin() + r.get_double(0., 1.) * (bbox.xmax() - bbox.xmin()),
                        bbox.ymin() + r.get_double(0., 1.) * (bbox.ymax() - bbox.ymin()),
                        bbox.zmin() + r.get_double(0., 1.) * (bbox.zmax() - bbox.zmin()));

    const FT sqd = oracle.squared_distance(seed);
#ifdef CGAL_AW3_DEBUG
    std::cout << "Generate " << seed << " at squared distance " << sqd << " (sqo: " << sq_offset << ")" << std::endl;
#endif

    if(sqd > sq_offset)
      seeds.push_back(seed);
  }
}

void alpha_wrap_triangle_mesh(Mesh& input_mesh,
                              const double alpha,
                              const double offset,
                              Seeds& seeds,
                              CGAL::Random& r)
{
  namespace AW3 = CGAL::Alpha_wraps_3;
  namespace PMP = CGAL::Polygon_mesh_processing;

  using Geom_traits = typename CGAL::GetGeomTraits<Mesh>::type;
  using Oracle = AW3::internal::Triangle_mesh_oracle<Geom_traits>;

  std::cout << "Input: " << num_vertices(input_mesh) << " vertices, " << num_faces(input_mesh) << " faces" << std::endl;

  bool has_degeneracies = !PMP::remove_degenerate_faces(input_mesh);
  if(has_degeneracies)
    std::cerr << "Warning: Failed to remove some degenerate faces." << std::endl;

  std::cout << "Processed input: " << vertices(input_mesh).size() << " vertices, " << faces(input_mesh).size() << " faces" << std::endl;
//  CGAL::IO::write_polygon_mesh("input.off", input_mesh, CGAL::parameters::stream_precision(17));

  Oracle oracle;
  oracle.add_triangle_mesh(input_mesh);
  AW3::internal::Alpha_wrapper_3<Oracle> aw3(oracle);

  if(seeds.empty())
    generate_random_seeds(oracle, offset, seeds, r);

  std::cout << "Seeds:" << std::endl;
  for(const Point_3& s : seeds)
    std::cout << s << std::endl;

  assert(!seeds.empty());

  const bool enforce_manifoldness = true;

  Mesh wrap;
  aw3(alpha, offset, wrap,
      CGAL::parameters::seed_points(std::ref(seeds))
                       .do_enforce_manifoldness(enforce_manifoldness));

  std::cout << "Result: " << vertices(wrap).size() << " vertices, " << faces(wrap).size() << " faces" << std::endl;

  // Tolerate failed initialization since we use random seeds and it's difficult to guarantee it
  if(is_empty(wrap))
    return;

//  CGAL::IO::write_polygon_mesh("last.off", wrap, CGAL::parameters::stream_precision(17));

  if(!has_degeneracies)
  {
    assert(AW3::internal::is_valid_wrap(wrap, enforce_manifoldness));

    if(!enforce_manifoldness)
      assert(AW3::internal::has_expected_Hausdorff_distance(wrap, input_mesh, alpha, offset));
  }

  if(!enforce_manifoldness)
    assert(AW3::internal::check_edge_length(wrap, alpha));
}

void alpha_wrap_triangle_mesh(Mesh& input_mesh,
                              const double alpha,
                              const double offset,
                              CGAL::Random& r)
{
  Seeds seeds;
  return alpha_wrap_triangle_mesh(input_mesh, alpha, offset, seeds, r);
}

void alpha_wrap_triangle_mesh(const std::string& filename)
{
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
    CGAL::Random r;

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

    alpha_wrap_triangle_mesh(input_mesh, alpha, offset, r);
  }
}

void alpha_wrap_triangle_mesh(const std::string& filename,
                              Seeds& seeds)
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
    const double alpha_expo = r.get_double(-3., 3.);
    const double offset_expo = r.get_double(-3., 3.);
    const double alpha = longest_diag_length / std::pow(2, alpha_expo);
    const double offset = longest_diag_length / std::pow(2, offset_expo);

    std::cout << "===================================================" << std::endl;
    std::cout << filename << " " << alpha << " " << offset << std::endl;
    std::cout << "Random seed = " << r.get_seed() << std::endl;

    alpha_wrap_triangle_mesh(input_mesh, alpha, offset, seeds, r);
  }
}

int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  if(argc > 1)
  {
    alpha_wrap_triangle_mesh(argv[1]);
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
