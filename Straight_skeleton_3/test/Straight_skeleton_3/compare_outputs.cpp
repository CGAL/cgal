#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/boost/graph/helpers.h>

#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/repair_degeneracies.h>
#include <CGAL/Polygon_mesh_processing/shape_predicates.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>

#include <iostream>
#include <filesystem>
#include <fstream>
#include <type_traits>
#include <vector>

namespace PMP = ::CGAL::Polygon_mesh_processing;

// Kernel used in reading operations
// -- Don't change this: we only want to read and write doubles
using IK = CGAL::Exact_predicates_inexact_constructions_kernel;

// Kernel used in testing operations
// using K = CGAL::Exact_predicates_exact_constructions_kernel;
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
//using K = CGAL::Single_precision_epick;

template <typename Mesh>
bool has_degenerate_faces(const Mesh& mesh)
{
  for (auto f : faces(mesh))
    if (PMP::is_degenerate_triangle_face(f, mesh))
      return true;

  return false;
}

template <typename Mesh>
bool compare(const Mesh& ours,
             const Mesh& theirs)
{
  using VPM = typename CGAL::GetVertexPointMap<Mesh>::const_type;
  using P = typename boost::property_traits<VPM>::value_type;
  using FT = typename CGAL::Kernel_traits<P>::type::FT;
  static_assert((std::is_same_v<FT, K::FT>));

  bool are_similar = true;

  // volume
  const FT vol_ours = PMP::volume(ours);
  const FT vol_theirs = PMP::volume(theirs);
  std::cout << "vol_ours = " << vol_ours << std::endl;
  std::cout << "vol_theirs = " << vol_theirs << std::endl;
  std::cout << "volume ratio = " << vol_ours / vol_theirs << std::endl;

  if(CGAL::abs(vol_ours - vol_theirs) > 0.001 * vol_ours) {
    std::cerr << "Error: difference in volume" << std::endl;
    if(CGAL::abs(vol_ours - vol_theirs) > 0.01 * vol_ours) {
      std::cerr << "Error: significant difference in volume" << std::endl;
      if(CGAL::abs(vol_ours - vol_theirs) > 0.1 * vol_ours) {
        std::cerr << "Error: REALLY significant difference in volume" << std::endl;
      }
    }
    are_similar = false;
  }

  // area
  const FT area_ours = PMP::area(ours);
  const FT area_theirs = PMP::area(theirs);
  std::cout << "area_ours = " << area_ours << std::endl;
  std::cout << "area_theirs = " << area_theirs << std::endl;
  std::cout << "area ratio = " << area_ours / area_theirs << std::endl;

  if(CGAL::abs(area_ours - area_theirs) > 0.01 * area_ours) {
    std::cerr << "Error: significant difference in area" << std::endl;
    are_similar = false;
  }

  // Hausdorff
  CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(ours);
  const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                       CGAL::square(bbox.ymax() - bbox.ymin()) +
                                       CGAL::square(bbox.zmax() - bbox.zmin()));

  FT hd = PMP::approximate_symmetric_Hausdorff_distance<CGAL::Parallel_if_available_tag>(
                ours, theirs,
                CGAL::parameters::number_of_points_on_edges(10)
                                 .number_of_points_on_faces(100)
                                 .random_seed(0));
  std::cout << "symmetric Hausdorff distance: " << hd << std::endl;
  std::cout << "symmetric Hausdorff ratio: " << hd / diag_length << std::endl;
  if(hd > 0.001 * diag_length) {
    std::cerr << "Error: large Hausdorff distance" << std::endl;
    are_similar = false;
  }

  return are_similar;
}

int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  if(argc < 2)
  {
    std::cerr << "Usage: "
              << argv[0] << "\n"
              << "\tours_filename\n"
              << "\ttheirs_filename\n"
              << "\tsave_path\n"
              << std::endl;
    std::cerr << "Input format: any format readable with CGAL::IO::read_polygon_mesh()" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Reading with Kernel: " << typeid(IK).name() << std::endl;
  std::cout << "Testing with Kernel: " << typeid(K).name() << std::endl;

  const char* ours_filename = argv[1];
  const char* theirs_filename = argv[2];

  CGAL::Surface_mesh<IK::Point_3> in_ours;
  if(!CGAL::IO::read_polygon_mesh(ours_filename, in_ours)) {
    std::cerr << "Error: failed to read input (ours)" << std::endl;
    return EXIT_FAILURE;
  }

  CGAL::Surface_mesh<IK::Point_3> in_theirs;
  if(!CGAL::IO::read_polygon_mesh(theirs_filename, in_theirs)) {
    std::cerr << "Error: failed to read input (theirs)" << std::endl;
    return EXIT_FAILURE;
  }

  // convert to testing kernels
  CGAL::Surface_mesh<K::Point_3> ours, theirs;
  CGAL::copy_face_graph(in_ours, ours);
  CGAL::copy_face_graph(in_theirs, theirs);

  std::cout << "Ours: " << num_vertices(ours) << " NV " << num_faces(theirs) << " NF" << std::endl;
  std::cout << "Theirs: " << num_vertices(theirs) << " NV " << num_faces(theirs) << " NF" << std::endl;

  std::filesystem::path save_path = (argc > 3) ? std::filesystem::path(argv[3]) : std::filesystem::current_path();

  CGAL::IO::write_polygon_mesh(save_path / "ours.off", ours,
                               CGAL::parameters::stream_precision(std::numeric_limits<K::FT>::max_digits10));

#define CGAL_SS3_ATTEMPT_REPAIR
#ifdef CGAL_SS3_ATTEMPT_REPAIR
  if (!CGAL::is_closed(ours) ||
      has_degenerate_faces(ours) ||
      PMP::does_self_intersect(ours))
  {
    std::cout << "attempting repair..." << std::endl;

    CGAL::Bbox_3 bbox = PMP::bbox(ours);
    const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                        CGAL::square(bbox.ymax() - bbox.ymin()) +
                                        CGAL::square(bbox.zmax() - bbox.zmin()));
    const double threshold = 1e-6 * diag_length;

    PMP::remove_almost_degenerate_faces(ours, CGAL::parameters::collapse_length_threshold(threshold)
                                                               .flip_triangle_height_threshold(threshold));

    std::vector<K::Point_3> ours_points;
    std::vector<std::vector<std::size_t> > ours_polygons;
    PMP::polygon_mesh_to_polygon_soup(ours, ours_points, ours_polygons);
    bool is_sane = PMP::autorefine_triangle_soup(ours_points, ours_polygons, CGAL::parameters::apply_iterative_snap_rounding(true));
    if (!is_sane) {
      std::cout << "Still has SI after autoref" << std::endl;
    }

    PMP::orient_polygon_soup(ours_points, ours_polygons);
    CGAL_assertion(PMP::is_polygon_soup_a_polygon_mesh(ours_polygons));

    clear(ours);
    PMP::polygon_soup_to_polygon_mesh(ours_points, ours_polygons, ours);

    std::cout << PMP::internal::number_of_connected_components(ours) << " CCs" << std::endl;
    PMP::keep_largest_connected_components(ours, 1);

    std::cout << "is closed? " << CGAL::is_closed(ours) << std::endl;
    PMP::stitch_borders(ours);

    CGAL::IO::write_polygon_mesh(save_path / "ours_processed.off", ours,
                                 CGAL::parameters::stream_precision(std::numeric_limits<K::FT>::max_digits10));
  }
#endif

  const bool closed_ours = CGAL::is_closed(ours);
  const bool closed_theirs = CGAL::is_closed(theirs);
  const bool has_degen_ours = has_degenerate_faces(ours);
  const bool has_degen_theirs = has_degenerate_faces(theirs);
  const bool SI_ours = PMP::does_self_intersect(ours);
  const bool SI_theirs = PMP::does_self_intersect(theirs);

  std::cout << "Health checks (ours) " << closed_ours << " " << !has_degen_ours << " " << !SI_ours << std::endl;
  std::cout << "Health checks (theirs) " << closed_theirs << " " << !has_degen_theirs << " " << !SI_theirs << std::endl;

  if(!closed_ours || has_degen_ours || SI_ours) {
    std::cerr << "Error: 'ours' is BAD" << std::endl;
    return EXIT_FAILURE;
  }

  if(!closed_theirs || has_degen_theirs || SI_theirs) {
    std::cerr << "Error: 'theirs' is BAD" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Valid inputs" << std::endl;

  if(!compare(ours, theirs)) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
