#define CGAL_AW3_TIMER
#define CGAL_AW3_DEBUG
#define CGAL_AW3_DEBUG_MANIFOLDNESS

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/alpha_wrap_3.h>
#include <CGAL/Alpha_wrap_3/internal/validation.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>

#include <iostream>
#include <vector>

using namespace CGAL::Alpha_wraps_3::internal;

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT = Kernel::FT;
using Point_3 = Kernel::Point_3;
using Vector_3 = Kernel::Vector_3;

using Points = std::vector<Point_3>;
using Face = std::vector<std::size_t>;
using Faces = std::vector<Face>;

using Mesh = CGAL::Surface_mesh<Point_3>;

void alpha_wrap_triangle_soup(Points& pr,
                              Faces& fr,
                              double alpha,
                              double offset)
{
  namespace AW3 = CGAL::Alpha_wraps_3;
  namespace PMP = CGAL::Polygon_mesh_processing;

  using Oracle = AW3::internal::Triangle_soup_oracle<Kernel, int, false /*subdivide*/>;

  std::cout << "Input: " << pr.size() << " points, " << fr.size() << " faces" << std::endl;

  PMP::repair_polygon_soup(pr, fr);
  std::cout << "Processed input: " << pr.size() << " points, " << fr.size() << " faces" << std::endl;
//  CGAL::IO::write_polygon_soup("input.off", pr, fr, CGAL::parameters::stream_precision(17));

  Mesh input_mesh; // only required for Hausdorff
  PMP::orient_polygon_soup(pr, fr);
  assert(PMP::is_polygon_soup_a_polygon_mesh(fr));
  PMP::polygon_soup_to_polygon_mesh(pr, fr, input_mesh);

  // AW3
  Oracle oracle;
  oracle.add_triangle_soup(pr, fr);
  AW3::internal::Alpha_wrapper_3<Oracle> aw3(oracle);

  Mesh wrap;
  aw3(alpha, offset, wrap, CGAL::parameters::do_enforce_manifoldness(false));

  std::cout << "First call result: " << vertices(wrap).size() << " vertices, " << faces(wrap).size() << " faces" << std::endl;

//  CGAL::IO::write_polygon_mesh("last.off", wrap, CGAL::parameters::stream_precision(17));

  assert(AW3::internal::is_valid_wrap(wrap, false /*manifoldness*/));
  assert(AW3::internal::is_outer_wrap_of_triangle_soup(wrap, pr, fr));
  assert(AW3::internal::has_expected_Hausdorff_distance(wrap, input_mesh, alpha, offset));
  assert(AW3::internal::check_edge_length(wrap, alpha));

  alpha *= 2;
  offset *= 2;

  Mesh wrap_2;
  aw3(alpha, offset, wrap_2, CGAL::parameters::do_enforce_manifoldness(false));

  std::cout << "Second call result: " << vertices(wrap_2).size() << " vertices, " << faces(wrap_2).size() << " faces" << std::endl;

//  CGAL::IO::write_polygon_mesh("last.off", wrap_2, CGAL::parameters::stream_precision(17));

  // Might fail for very small meshes
//  assert(num_vertices(wrap_2) <= num_vertices(wrap) && num_faces(wrap_2) <= num_faces(wrap));

  assert(AW3::internal::is_valid_wrap(wrap_2, false /*manifoldness*/));
  assert(AW3::internal::is_outer_wrap_of_triangle_soup(wrap_2, pr, fr));
  assert(AW3::internal::has_expected_Hausdorff_distance(wrap_2, input_mesh, alpha, offset));
  assert(AW3::internal::check_edge_length(wrap_2, alpha));
}

void alpha_wrap_triangle_soup(const std::string& filename)
{
  Points points;
  Faces faces;
  bool res = CGAL::IO::read_polygon_soup(filename, points, faces);
  assert(res);
  assert(!faces.empty());

  CGAL::Bbox_3 bbox;
  for(const auto& f : faces)
    for(int i=0; i<3; ++i)
      bbox += points[f[i]].bbox();

  const Vector_3 longest_diag = Point_3(bbox.xmax(), bbox.ymax(), bbox.zmax()) -
                                Point_3(bbox.xmin(), bbox.ymin(), bbox.zmin());
  double longest_diag_length = CGAL::to_double(CGAL::approximate_sqrt(longest_diag.squared_length()));

  CGAL::Random r;

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

    alpha_wrap_triangle_soup(points, faces, alpha, offset);
  }
}

int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  if(argc > 1)
  {
    alpha_wrap_triangle_soup(argv[1]);
    return EXIT_SUCCESS;
  }

  alpha_wrap_triangle_soup("data/tetrahedron.off");
  alpha_wrap_triangle_soup("data/plane.off");
  alpha_wrap_triangle_soup("data/open_tetrahedron.off");
  alpha_wrap_triangle_soup("data/tetrahedron_open_tip.off");
  alpha_wrap_triangle_soup("data/sphere_one_hole.off");
  alpha_wrap_triangle_soup("data/tetrahedron_duplicates.off");
  alpha_wrap_triangle_soup("data/tetrahedron_degenerencies.off");
  alpha_wrap_triangle_soup("data/non_manifold.off");
  alpha_wrap_triangle_soup("data/combinatorial_manifold.off");
  alpha_wrap_triangle_soup("data/combinatorial_manifold_multiple_components.off");
  alpha_wrap_triangle_soup("data/tetrahedron_self_intersection_tip.off");
  alpha_wrap_triangle_soup("data/tetrahedron_twisted_tip.off");
  alpha_wrap_triangle_soup("data/tetrahedron_random_perturbation.off");
//  alpha_wrap_triangle_soup("data/overlay_triangulation.off"); // due to geometrically degenerate faces + using soups here
  alpha_wrap_triangle_soup("data/two_knives.off");
  alpha_wrap_triangle_soup("data/three_knives.off");
  alpha_wrap_triangle_soup("data/bunny_random_perturbation.off");

  std::cout << "Done!" << std::endl;
  return EXIT_SUCCESS;
}
