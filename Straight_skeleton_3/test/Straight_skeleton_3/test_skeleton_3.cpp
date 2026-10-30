#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include "CGAL/_test_pipeline.h"

#include <CGAL/property_map.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <list>
#include <vector>

namespace SS3 = CGAL::Straight_skeletons_3;

bool failed = false;

using K = CGAL::Exact_predicates_exact_constructions_kernel;
using FT = K::FT;
using Point_3 = K::Point_3;

using Mesh = CGAL::Surface_mesh<Point_3>;
using face_descriptor = typename boost::graph_traits<Mesh>::face_descriptor;

void test(const std::string& mesh_filename,
          const std::vector<FT>& save_times,
          const std::filesystem::path& save_path = std::filesystem::current_path(),
          const std::filesystem::path& config_path = {},
          bool redict_stdout = true)
{
  std::cout << "\n==== TEST: " << mesh_filename << " ====" << std::endl;
  std::cout << "save_times: ";
  for (const auto& t : save_times) {
    std::cout << t << " ";
  }
  std::cout << std::endl;

  // do nothing if filename is commented (used in scripts)
  if (std::string_view(mesh_filename) == "#" || std::string_view(mesh_filename) == "//") {
    std::cerr << mesh_filename << " is not a valid path" << std::endl;
    failed = true;
    return;
  }

  // triangulate if necessary, and check if this is a sane input
  Mesh sm;
  if (!SS3::utils::preprocess_input(mesh_filename, sm, save_times)) {
    std::cerr << mesh_filename << " is not a valid input" << std::endl;
    failed = true;
    return;
  } else {
    std::cout << mesh_filename << " is a valid input" << std::endl;
  }

#if 0
  if (num_faces(sm) > 500) {
    std::cerr << mesh_filename << " is too large" << std::endl;
    return;
  }
#endif

  std::ofstream log;

  if (redict_stdout) {
    std::string mesh_file_stem = std::filesystem::path(mesh_filename).stem().string();
    std::string log_path = "logs/" + mesh_file_stem + ".log";
    std::cout << "Writing to " << log_path << std::endl;
    log.open(log_path, std::ios::out | std::ios::trunc);
    log.precision(17);
    CGAL::Straight_skeletons_3::internal::set_log_stream(log);
  }

  // weights
  auto fwm = get(CGAL::dynamic_face_property_t<FT>{}, sm);

  auto inner_fwm = sm.template add_property_map<face_descriptor, double>("f:weight");
  if (!inner_fwm.second)
  {
    std::cout << "Input mesh had an inner f:weight property map" << std::endl;
    for (face_descriptor f : faces(sm)) {
      std::cout << "  " << f << ": " << get(inner_fwm.first, f) << std::endl;
      put(fwm, f, get(inner_fwm.first, f));
    }

    SS3::utils::clamp_weights(sm, fwm);
  } else {
    // constant property map
    for (face_descriptor f : faces(sm)) {
      put(fwm, f, 1);
    }
  }

  bool result = SS3::utils::run(sm, fwm, save_times, save_path, config_path);
  std::cout << "result: " << result << std::endl;
  if (!result)
    failed = true;
  assert(result);
}

void test_events(const std::vector<FT>& save_times = {})
{
  test("data/events/DblEdgeMergeEvent.obj", save_times);
  test("data/events/DblEdgeMergeEvent_proper.obj", save_times);
  test("data/events/DblTriangleEvent.obj", save_times);
  test("data/events/DblTriangleEvent_reflex.obj", save_times);
  test("data/events/EdgeEvent_noflip.obj", save_times);
  test("data/events/EdgeEvent.obj", save_times);
  test("data/events/EdgeMergeEvent_inverted.obj", save_times);
  test("data/events/EdgeMergeEvent.obj", save_times);
  test("data/events/EdgeMergeEvent_proper.obj", save_times);
  test("data/events/EdgeMergeEvent_reflex.obj", save_times);
  test("data/events/EdgeSplitEvent.obj", save_times);
  test("data/events/FlipVertexEvent.obj", save_times);
  test("data/events/misc_0.obj", save_times);
  test("data/events/misc_1.obj", save_times);
  test("data/events/PierceEvent2.obj", save_times);
  test("data/events/PierceEvent3bis.obj", save_times);
  test("data/events/PierceEvent3.obj", save_times);
  test("data/events/PierceEvent3quater.obj", save_times);
  test("data/events/PierceEvent3ter.obj", save_times);
  test("data/events/PierceEvent4.obj", save_times);
  test("data/events/PierceEvent.obj", save_times);
  test("data/events/PolyhedronSplitEvent-deg3.obj", save_times);
  test("data/events/PolyhedronSplitEvent.obj", save_times);
  test("data/events/SplitMergeEvent.obj", save_times);
  test("data/events/SplitMergeEvent_reflex.obj", save_times);
  test("data/events/SplitMergeEvent_reflex_proper.obj", save_times);
  test("data/events/SurfaceEvent2.obj", save_times);
  test("data/events/SurfaceEvent.obj", save_times);
  test("data/events/SurfaceEvent_reflex.obj", save_times);
  test("data/events/SurfaceEvent_topol.obj", save_times);
  test("data/events/TetrahedronEvent.obj", save_times);
  test("data/events/TriangleEvent2.obj", save_times);
  test("data/events/TriangleEvent.obj", save_times);
  test("data/events/VertexEvent_EdgeSplit.obj", save_times);
  test("data/events/VertexEvent.obj", save_times);
  test("data/events/VertexEvent_topol.obj", save_times);
}

void test_custom_models(const std::vector<FT>& save_times = {})
{
  test("data/4xtreme.obj", save_times);
  test("data/724_star-Butterfly.obj", save_times);
  test("data/726_star-PNSplit.obj", save_times);
  test("data/aoi_asteroid.obj", save_times);
  test("data/Armadillo_000096.obj", save_times);
  test("data/Armadillo_000194.obj", save_times);
  test("data/arrangement.obj", save_times);
  test("data/auren_8.obj", save_times);
  test("data/box_2.obj", save_times);
  test("data/box_3.obj", save_times);
  test("data/box.obj", save_times);
  test("data/bunny_small.obj", save_times);
  test("data/chess_bauer.obj", save_times);
  test("data/chinese_lion_174.obj", save_times);
  test("data/connected_tunnel.obj", save_times);
  test("data/convex_piece_2.obj", save_times);
  test("data/cubet.obj", save_times);
  test("data/cube_with_hole.off", save_times);
  test("data/culver_iron_maiden.obj", save_times);
  test("data/discontinuous.obj", save_times);
  test("data/disproof_spherical.obj", save_times);
  test("data/doublebox_n.obj", save_times);
  test("data/doublebox_n-tmp.obj", save_times);
  test("data/double_hole_n.obj", save_times);
  test("data/EdgeEvent_old.obj", save_times);
  test("data/golem_tunnel_n.obj", save_times);
  test("data/graph_test_2.obj", save_times);
  test("data/graph_test_auren.obj", save_times);
  test("data/hand_small.obj", save_times);
  test("data/held_convex.obj", save_times);
  test("data/idea_lower_bound.obj", save_times);
  test("data/iron_maiden.obj", save_times);
  test("data/journal_armadillo.obj", save_times);
  test("data/journal_asteroid.obj", save_times);
  test("data/journal_bunny.obj", save_times);
  test("data/journal_lion.obj", save_times);
  test("data/journal_venus.obj", save_times);
  test("data/journal_verworrtakelt.obj", save_times);
  test("data/kiev.obj", save_times);
  test("data/krampus_tunnel_n.obj", save_times);
  test("data/mini_maiden.obj", save_times);
  test("data/mpi.off", save_times);
  test("data/prisms5_n.obj", save_times);
  test("data/prisms5.obj", save_times);
  test("data/pyramid_3c2r.obj", save_times);
  test("data/pyramid_5c1r.obj", save_times);
  test("data/pyramid.off", save_times);
  test("data/pyramids_hole_n.obj", save_times);
  test("data/pyramids_hole.obj", save_times);
  test("data/saddle_2.obj", save_times);
  test("data/saddle_3.obj", save_times);
  test("data/saddle6.obj", save_times);
  test("data/saddle.obj", save_times);
  test("data/Schoenhardt.obj", save_times);
  test("data/seastar.obj", save_times);
  test("data/small_sphere_shake.obj", save_times);
  test("data/splitConvexVertex5.obj", save_times);
  test("data/splitConvexVertex.obj", save_times);
  test("data/splitEdgeEvent_noflip.obj", save_times);
  test("data/splitReflexVertex.obj", save_times);
  test("data/splitTermination.obj", save_times);
  test("data/split_wedge_tabletop.obj", save_times);
  test("data/Star.obj", save_times);
  test("data/StarVertex.obj", save_times);
  test("data/sticking.obj", save_times);
  test("data/ted_mit_loch.obj", save_times);
  test("data/testWeightVertexSplitter.obj", save_times);
  test("data/torus.obj", save_times);
  test("data/torus_void_n.obj", save_times);
  test("data/torus_void.obj", save_times);
  test("data/transformer_clean_n.obj", save_times);
  test("data/transformer_last_n.obj", save_times);
  test("data/transformer_n.obj", save_times);
  test("data/tunnel.obj", save_times);
  test("data/venus_115.obj", save_times);
  test("data/venus_267.obj", save_times);
  test("data/verworrtakelt_fixed.obj", save_times);
  test("data/verworrtakelt_n.obj", save_times);
  test("data/verworrtakelt.obj", save_times);
  test("data/wedge_tabletop.obj", save_times);
}

void test_defaults()
{
  test_events();
  test_events({1, 2, 3});
  test_custom_models();
  test_custom_models({1, 2, 3});
}

int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  if (argc == 1) {
    test_defaults();
  } else {
    // Argument parsing
    std::filesystem::path mesh_filename;
    std::filesystem::path weights_filename;
    std::filesystem::path save_path = std::filesystem::current_path();
    std::filesystem::path config_path;
    std::vector<FT> save_times;

    if (!SS3::utils::parse_args(argc, argv, mesh_filename, weights_filename, save_path, config_path, save_times)) {
      return EXIT_FAILURE;
    }

    test(mesh_filename, save_times, save_path, config_path, false /*do not redirect stdout*/);
  }

  std::cout << "Done" << std::endl;
  return failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
