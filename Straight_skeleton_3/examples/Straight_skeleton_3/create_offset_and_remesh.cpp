#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Straight_skeleton_3/create_offset_polyhedra.h>

#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/remesh_planar_patches.h>
#include <CGAL/Polygon_mesh_processing/region_growing.h>
#include <CGAL/Polygon_mesh_processing/repair_degeneracies.h>

#include <CGAL/Real_timer.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/Color.h>

#include <iostream>
#include <sstream>
#include <vector>

namespace PMP = CGAL::Polygon_mesh_processing;
namespace SS3 = CGAL::Straight_skeletons_3;

using K = CGAL::Exact_predicates_exact_constructions_kernel;
using FT = K::FT;
using Point_3 = K::Point_3;

using Mesh = CGAL::Surface_mesh<Point_3>;
using face_descriptor = boost::graph_traits<Mesh>::face_descriptor;
using edge_descriptor = boost::graph_traits<Mesh>::edge_descriptor;
using halfedge_descriptor = boost::graph_traits<Mesh>::halfedge_descriptor;

// Generate a unique color based on face index
CGAL::IO::Color generate_color(const std::size_t index)
{
  // Use a simple hash-like coloring scheme to generate distinct colors
  unsigned char r = (index * 73) % 256;
  unsigned char g = (index * 127) % 256;
  unsigned char b = (index * 181) % 256;
  return { r, g, b };
}

// visitor to preserve face colors during triangulation
template <typename ColorMap>
struct Color_preserving_visitor
  : public PMP::Triangulate_faces::Default_visitor<Mesh>
{
  ColorMap color_map;
  CGAL::IO::Color original_color;

  Color_preserving_visitor(ColorMap cm) : color_map(cm) {}

  void before_subface_creations(face_descriptor f_old)
  {
    original_color = get(color_map, f_old);
  }

  void after_subface_created(face_descriptor f_new)
  {
    put(color_map, f_new, original_color);
  }
};

int main(int argc, char** argv)
{
  if (argc == 1) {
    std::cout << "Usage: " << argv[0] << " <input filename>" << std::endl;
    std::cout << "  <input filename>: Path to input mesh file (required)" << std::endl;
    return EXIT_FAILURE;
  }

  const char* mesh_filename = argv[1];

  Mesh sm;
  if(!CGAL::IO::read_polygon_mesh(mesh_filename, sm) || CGAL::is_empty(sm) || !is_valid_face_graph(sm)) {
    std::cerr << "Error: failed to read input " << mesh_filename << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Input mesh: " << vertices(sm).size() << " NV " << faces(sm).size() << " NF" << std::endl;

  CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(sm);
  const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                       CGAL::square(bbox.ymax() - bbox.ymin()) +
                                       CGAL::square(bbox.zmax() - bbox.zmin()));
  const double threshold_lgth = 0.01 * diag_length;

  // Color the original mesh with unique colors *before* triangulation
  auto input_color_map = get(CGAL::dynamic_face_property_t<CGAL::IO::Color>(), sm);
  std::size_t face_index = 0;
  for (auto f : faces(sm))
    put(input_color_map, f, generate_color(face_index++));

  CGAL::IO::write_OFF("colored_input.off", sm, CGAL::parameters::face_color_map(input_color_map).stream_precision(17));

  // let's tolerate untriangulated inputs
  if (!CGAL::is_triangle_mesh(sm)) {
    std::cout << "Triangulating input mesh..." << std::endl;
    Color_preserving_visitor color_visitor(input_color_map);
    PMP::triangulate_faces(sm, CGAL::parameters::visitor(color_visitor));
  }

  std::cout << "Input mesh: " << vertices(sm).size() << " NV " << faces(sm).size() << " NF (triangulated)" << std::endl;

  std::vector<FT> save_times = { -25, -50, -75 };

  // every face has weight 3.0, meaning faces move at thrice the default speed
  CGAL::Constant_property_map<face_descriptor, double> face_weight_map(3.);

  std::vector<Mesh> results;
  results.reserve(save_times.size());

  CGAL::Real_timer timer;
  timer.start();

  // correspondence maps from output faces to input faces
  std::vector<CGAL::unordered_flat_map<face_descriptor, face_descriptor>> offset_faces_to_input_faces(save_times.size());

  // Main call
  bool success = CGAL::create_straight_skeleton_and_offset_polyhedra_3(
    sm, save_times, results,
    CGAL::parameters::face_weight_map(face_weight_map));

  timer.stop();
  std::cout << "Elapsed: " << timer.time() << std::endl;

  // save the results
  for (std::size_t i=0; i<results.size(); ++i) {
    const FT save_time = save_times[i];
    Mesh& result_sm = results[i];

    std::cout << "Result mesh @ " << save_time << " : " << vertices(result_sm).size() << " NV " << faces(result_sm).size() << " NF (triangulated)" << std::endl;

    std::stringstream out_ss;
    out_ss << "result_" << save_time << ".off";

    if (!CGAL::IO::write_OFF(out_ss.str(), result_sm, CGAL::parameters::stream_precision(17))) {
      std::cerr << "Error: failed to write result " << out_ss.str() << std::endl;
      return EXIT_FAILURE;
    }

    // Now, let's simplify the mesh to get polygonal faces...

    // Do some preprocessing to remove nasty elements that were created by perturbation of almost parallel faces
    PMP::remove_almost_degenerate_faces(result_sm,
                                        CGAL::parameters::cap_threshold(std::cos(175. / 180 * CGAL_PI))
                                                         .needle_threshold(10)
                                                         .collapse_length_threshold(threshold_lgth)
                                                         .flip_triangle_height_threshold(threshold_lgth));

    std::stringstream out_css;
    out_css << "cleaned_result_" << save_time << ".off";

    if (!CGAL::IO::write_OFF(out_css.str(), result_sm, CGAL::parameters::stream_precision(17))) {
      std::cerr << "Error: failed to write result " << out_css.str() << std::endl;
      return EXIT_FAILURE;
    }

    std::vector<std::size_t> region_ids(num_faces(result_sm));
    std::vector<std::size_t> corner_id_map(num_vertices(result_sm), -1); // corner status of vertices
    std::vector<bool> ecm(num_edges(result_sm), false); // mark edges at the boundary of regions
    boost::vector_property_map<K::Vector_3> normal_map; // normal of the supporting planes of the regions detected

    // detect planar regions in the mesh
    std::size_t nb_regions =
      PMP::region_growing_of_planes_on_faces(result_sm,
                                             CGAL::make_random_access_property_map(region_ids),
                                             CGAL::parameters::maximum_distance(threshold_lgth)
                                                              .region_primitive_map(normal_map));

    // detect corner vertices on the boundary of planar regions
    std::size_t nb_corners =
      PMP::detect_corners_of_regions(result_sm,
                                     CGAL::make_random_access_property_map(region_ids),
                                     nb_regions,
                                     CGAL::make_random_access_property_map(corner_id_map),
                                     CGAL::parameters::maximum_distance(threshold_lgth).
                                                       edge_is_constrained_map(CGAL::make_random_access_property_map(ecm)));

    Mesh remeshed_sm;
    PMP::remesh_almost_planar_patches(result_sm, remeshed_sm,
                                      nb_regions, nb_corners,
                                      CGAL::make_random_access_property_map(region_ids),
                                      CGAL::make_random_access_property_map(corner_id_map),
                                      CGAL::make_random_access_property_map(ecm),
                                      CGAL::parameters::patch_normal_map(normal_map),
                                      CGAL::parameters::do_not_triangulate_faces(true));

    std::cout << "Result mesh @ " << save_time << " : " << vertices(remeshed_sm).size() << " NV " << faces(remeshed_sm).size() << " NF (remeshed)" << std::endl;

    std::stringstream out_rss;
    out_rss << "remeshed_result_" << save_time << ".off";

    if (!CGAL::IO::write_OFF(out_rss.str(), remeshed_sm, CGAL::parameters::stream_precision(17))) {
      std::cerr << "Error: failed to write result " << out_rss.str() << std::endl;
      return EXIT_FAILURE;
    }
  }

  std::cout << "Done" << std::endl;
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
