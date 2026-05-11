#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Straight_skeleton_3/create_offset_polyhedra.h>

#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/remesh_planar_patches.h>

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
    CGAL::parameters::face_weight_map(face_weight_map),
    CGAL::parameters::face_to_face_map(boost::make_assoc_property_map(offset_faces_to_input_faces[0])),
    CGAL::parameters::face_to_face_map(boost::make_assoc_property_map(offset_faces_to_input_faces[1])),
    CGAL::parameters::face_to_face_map(boost::make_assoc_property_map(offset_faces_to_input_faces[2])));

  timer.stop();
  std::cout << "Elapsed: " << timer.time() << std::endl;

  // save the results
  for (std::size_t i=0; i<results.size(); ++i) {
    const FT save_time = save_times[i];
    const Mesh& result_sm = results[i];

    std::cout << "Result mesh @ " << save_time << " : " << vertices(result_sm).size() << " NV " << faces(result_sm).size() << " NF (triangulated)" << std::endl;

    // Create a color property map for this result using the F2F map
    auto result_color_map = CGAL::get(CGAL::dynamic_face_property_t<CGAL::IO::Color>(), result_sm);
    for (auto f : faces(result_sm))
      put(result_color_map, f, get(input_color_map, offset_faces_to_input_faces[i][f]));

    std::stringstream out_ss;
    out_ss << "result_" << save_time << ".off";

    if (!CGAL::IO::write_OFF(out_ss.str(), result_sm,
                             CGAL::parameters::face_color_map(result_color_map)
                                              .stream_precision(17))) {
      std::cerr << "Error: failed to write result " << out_ss.str() << std::endl;
      return EXIT_FAILURE;
    }

    // Now, let's simplify the mesh to go back to the original (non-triangulated) faces

    // corners are vertices that have 3 unique colors amongst incident faces

    std::vector<std::size_t> region_ids(num_faces(result_sm));
    std::vector<std::size_t> corner_ids(num_vertices(result_sm), -1); // corner status of vertices
    std::vector<bool> patch_borders(num_edges(result_sm), false); // mark edges at the boundary of regions

    // Build maps from input face ID to region ID and color

    // We proceed by flooding the output mesh, assigning the same region ID to faces that
    // are edge-adjacent and that have the same input face as ancestor (using the F2F map).
    // This way, we can assign a region ID to each output face, and then identify corners and patch borders based on these region IDs.
    // Note that an input face can thus correspond to different regions in the output
    std::vector<bool> visited(num_faces(result_sm), false);
    std::vector<face_descriptor> region_to_input_face;
    std::size_t current_region_id = 0;
    for (auto f : faces(result_sm)) {
      if (visited[f])
        continue;
      face_descriptor input_face = offset_faces_to_input_faces[i][f];
      std::vector<face_descriptor> stack = { f };
      region_to_input_face.push_back(input_face);
      std::cout << "region " << current_region_id << " corresponds to input face " << input_face << std::endl;
      std::cout << "  with color " << get(input_color_map, input_face) << std::endl;
      while (!stack.empty()) {
        face_descriptor current = stack.back();
        stack.pop_back();
        if (visited[current])
          continue;
        visited[current] = true;
        region_ids[current] = current_region_id;
        for (halfedge_descriptor h : halfedges_around_face(halfedge(current, result_sm), result_sm)) {
          face_descriptor neighbor = face(opposite(h, result_sm), result_sm);
          if (neighbor == Mesh::null_face() || visited[neighbor]) continue;
          if (offset_faces_to_input_faces[i][neighbor] == input_face)
            stack.push_back(neighbor);
        }
      }
      ++current_region_id;
    }

    CGAL_assertion(region_to_input_face.size() == current_region_id);

   // Identify corners: vertices with exactly 3 unique regions among incident faces
    std::size_t corner_count = 0;
    std::ofstream corner_out("corners.xyz");
    corner_out.precision(17);
    for (auto v : vertices(result_sm)) {
      std::set<std::size_t> incident_regions;
      for (auto f : faces_around_target(halfedge(v, result_sm), result_sm)) {
        incident_regions.insert(region_ids[f]);
      }

      // std::cout << incident_regions.size() << " incident regions" << std::endl;
      if (incident_regions.size() == 3) {
        corner_ids[v] = corner_count++;
        corner_out << result_sm.point(v) << "\n";
      }
    }

    // Identify patch borders: edges whose incident faces have different regions
    std::ofstream pb_out("patch_borders.polylines.txt");
    pb_out.precision(17);
    for (edge_descriptor e : edges(result_sm)) {
      CGAL_assertion(!is_border(e, result_sm));
      halfedge_descriptor h1 = halfedge(e, result_sm);
      halfedge_descriptor h2 = opposite(h1, result_sm);
      patch_borders[e] = (region_ids[face(h1, result_sm)] != region_ids[face(h2, result_sm)]);
      if (patch_borders[e]) {
        // std::cout << source(h1, result_sm) << " " << target(h1, result_sm) << " is a patch border" << std::endl;
        pb_out << "2 " << result_sm.point(source(h1, result_sm)) << " " << result_sm.point(target(h1, result_sm)) << "\n";
      }
    }

    std::size_t nb_regions = current_region_id;
    std::size_t nb_corners = corner_count;
    std::cout << "Identified " << nb_regions << " regions and " << nb_corners << " corners" << std::endl;

    // Now, let's remesh
    Mesh remeshed_sm;
    std::vector<std::size_t> remesh_face_patches(num_faces(result_sm));

    PMP::remesh_almost_planar_patches(result_sm, remeshed_sm,
                                      nb_regions, nb_corners,
                                      CGAL::make_random_access_property_map(region_ids),
                                      CGAL::make_random_access_property_map(corner_ids),
                                      CGAL::make_random_access_property_map(patch_borders),
                                      CGAL::parameters::default_values(),
                                      CGAL::parameters::do_not_triangulate_faces(true)
                                                       .face_patch_map(CGAL::make_random_access_property_map(remesh_face_patches)));

    std::cout << "Remeshed mesh @ " << save_time << " : " << vertices(remeshed_sm).size() << " NV " << faces(remeshed_sm).size() << " NF" << std::endl;

    // Reassign colors to remeshed faces based on their patch IDs
    for (face_descriptor f : faces(remeshed_sm)) {
      std::size_t patch_id = remesh_face_patches[f];
      face_descriptor input_face_id = region_to_input_face[patch_id];
      put(result_color_map, f, get(input_color_map, input_face_id));
    }

    std::stringstream out_rss;
    out_rss << "remeshed_result_" << save_time << ".off";

    if (!CGAL::IO::write_OFF(out_rss.str(), remeshed_sm,
                             CGAL::parameters::face_color_map(result_color_map)
                                              .stream_precision(17))) {
      std::cerr << "Error: failed to write result " << out_rss.str() << std::endl;
      return EXIT_FAILURE;
    }
  }

  std::cout << "Done" << std::endl;
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
