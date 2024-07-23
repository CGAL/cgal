#include "hexmeshing.h"
#include "utils.h"

#include <CGAL/Combinatorial_map_save_load.h>
#include <CGAL/config.h>
#include <CGAL/draw_polyhedron.h>
#include <CGAL/Graphics_scene_options.h>
#include <cstdlib>
#include <filesystem>
#include <iostream>

// int pattern_substitution()
// {
//   LCC lcc;
//   auto d1=
//     lcc.make_hexahedron(Point(0,0,0), Point(5,0,0),
//                         Point(5,5,0), Point(0,5,0),
//                         Point(0,5,5), Point(0,0,5),
//                         Point(5,0,5), Point(5,5,5));

//   auto d2=
//     lcc.make_hexahedron(Point(5,0,0), Point(10,0,0),
//                         Point(10,5,0), Point(5,5,0),
//                         Point(5,5,5), Point(5,0,5),
//                         Point(10,0,5), Point(10,5,5));

//   auto d3=
//     lcc.make_hexahedron(Point(0,0,-5), Point(5,0,-5),
//                         Point(5,5,-5), Point(0,5,-5),
//                         Point(0,5,0), Point(0,0,0),
//                         Point(5,0,0), Point(5,5,0));

//   lcc.sew<3>(lcc.beta(d1, 1, 1, 2), lcc.beta(d2, 2));
//   lcc.sew<3>(lcc.beta(d1, 0, 2), lcc.beta(d3, 1, 2));

//   lcc.display_characteristics(std::cout);

//   auto toprocess_mark = lcc.get_new_mark();
//   debug_edge_mark = lcc.get_new_mark();

//   lcc.mark_cell<0>(d1, toprocess_mark);
//   auto md1 = lcc.beta(d1, 1);
//   lcc.mark_cell<0>(md1, toprocess_mark);
//   auto md2 = lcc.beta(d2, 1);
//   lcc.mark_cell<0>(md2, toprocess_mark);
//   auto md3 = lcc.beta(d2, 1, 1);
//   lcc.mark_cell<0>(md3, toprocess_mark);
//   auto md4 = d3;
//   lcc.mark_cell<0>(md4, toprocess_mark);

//   // lcc.mark_cell<0>(md3, toprocess_mark);

//   std::vector<Dart_handle> marked_cells;

//   // We assume we already filled out marked_cells earlier

//   marked_cells.push_back(d1);
//   marked_cells.push_back(md1);
//   marked_cells.push_back(md2);
//   marked_cells.push_back(md3);
//   marked_cells.push_back(md4);

//   CGAL::HexRefinement::tworefinement_one_step_from_marked(lcc, marked_cells, toprocess_mark);

//   render(lcc, toprocess_mark, debug_edge_mark);
//   return EXIT_SUCCESS;
// }

// int grid_test() {
//   using namespace CGAL::HexRefinement::TwoRefinement;

//   LCC lcc;
//   auto grid = generate_grid(lcc, Point(0,0,0), 65.059600000000003, 10);

//   lcc.display_characteristics(std::cout);

//   auto m1 = lcc.get_new_mark();
//   auto m2 = lcc.get_new_mark();

//   mark_nodes_in_cell_pair(lcc, grid, Plane::XY, m1, m1);

//   // auto a =  lcc.beta(lcc.first_dart(), 0, 2, 0);
//   // mark_face(lcc, a, m2);
//   // lcc.mark_cell<0>(a, m1);

//   render<LCC>(lcc, m1, m2);
//   lcc.display_characteristics(std::cout);

//   return EXIT_SUCCESS;
// }


int surface_test() {
  using namespace CGAL::HexRefinement;

  LCC lcc = two_refinement(CGAL::data_file_path("query_replace/mesh2.off"), 10, TwoRefinement::mark_intersecting_volume_with_poly);
  lcc.display_characteristics(std::cout);

  // render<LCC>(lcc, debug_node_mark, debug_edge_mark );
  render_two_refinement_result(lcc);

  return EXIT_SUCCESS;
}

int two_stacked_3_template_test() {
  using namespace CGAL::HexRefinement::TwoRefinement;

  HexMeshingData hdata;
  LCC& lcc = hdata.lcc;

  hdata.grid = generate_grid(hdata.lcc, Point(0,0,0), 5, 4);
  hdata.identified_mark = lcc.get_new_mark();
  hdata.template_mark = lcc.get_new_mark();
  hdata.corner_mark = lcc.get_new_mark();
  hdata.propagation_face_mark = lcc.get_new_mark();

  lcc.negate_mark(hdata.corner_mark);

  load_patterns(hdata.regular_templates, hdata.partial_templates);

  RefinementData plane;
  Dart_handle first = lcc.first_dart();
  first = lcc.beta(first, 0, 0, 2, 3, 2);

  // TODO mark 2 volumes to shape two 3 templates using get_or_create_attr<3>

  lcc.mark_cell<0>(lcc.beta(first, 1, 2), hdata.identified_mark);
  lcc.mark_cell<0>(lcc.beta(first, 1, 2, 1), hdata.identified_mark);
  lcc.mark_cell<0>(lcc.beta(first, 1, 2, 0), hdata.identified_mark);

  setup_initial_planes(hdata);
  extract_darts_from_even_planes(hdata, plane, Plane::XY);
  create_vertices_for_templates(hdata, plane);
  refine_marked_hexes(hdata, plane);

  lcc.unmark_all(debug_node_mark);
  lcc.unmark_all(debug_edge_mark);
  Dart_handle d = lcc.beta(first, 1, 2, 3, 2, 1, 2, 3, 0, 0, 2, 3, 2, 1, 1, 2, 3, 2,  1, 2, 3, 2, 1,       2, 3, 2    );
  mark_edge(lcc, d, debug_node_mark);
  mark_face(lcc, d, debug_edge_mark);

  CGAL::save_combinatorial_map(lcc, "work.3map");
  render<LCC>(lcc, debug_node_mark, debug_edge_mark);

  return 0;
}

int propagation_face(){
  using namespace CGAL::HexRefinement::TwoRefinement;

  LCC lcc;
  auto d1=
    lcc.make_hexahedron(Point(0,0,0), Point(5,0,0),
                        Point(5,5,0), Point(0,5,0),
                        Point(0,5,5), Point(0,0,5),
                        Point(5,0,5), Point(5,5,5));

  auto d2=
    lcc.make_hexahedron(Point(5,0,0), Point(10,0,0),
                        Point(10,5,0), Point(5,5,0),
                        Point(5,5,5), Point(5,0,5),
                        Point(10,0,5), Point(10,5,5));

  auto d3=
    lcc.make_hexahedron(Point(0,0,-5), Point(5,0,-5),
                        Point(5,5,-5), Point(0,5,-5),
                        Point(0,5,0), Point(0,0,0),
                        Point(5,0,0), Point(5,5,0));

  lcc.sew<3>(lcc.beta(d1, 1, 1, 2), lcc.beta(d2, 2));
  lcc.sew<3>(lcc.beta(d1, 0, 2), lcc.beta(d3, 1, 2));

  auto m1 = lcc.get_new_mark();
  auto m2 = lcc.get_new_mark();

  auto a = lcc.beta(lcc.first_dart(), 1, 1, 2);
  mark_half_face_unchecked(lcc, a, m1);
  assert(!is_half_face_marked(lcc, lcc.beta(a,3), m1));

  render(lcc, m1, m2);

  return 1;
}


int create_validation_files(){
  namespace fs = std::filesystem;
  const std::string dir = CGAL::data_file_path("hexmeshing");

  for (const auto& entry : fs::directory_iterator(dir + "/mesh")){
    if (!entry.is_regular_file() or entry.path().extension() != ".off") continue;

    LCC lcc = CGAL::HexRefinement::two_refinement(entry.path().string(), 20, CGAL::HexRefinement::TwoRefinement::mark_intersecting_volume_with_poly);

    std::string out = dir + "/validation/" + entry.path().filename().replace_extension().string() + ".3map";

    if (!CGAL::save_combinatorial_map(lcc, out.c_str())){
      std::cerr << "Error while saving result of tworefinement for " << entry.path().filename() << std::endl;
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}

int validation_test(){
  namespace fs = std::filesystem;
  const std::string dir = CGAL::data_file_path("hexmeshing");

  for (const auto& entry : fs::directory_iterator(dir + "/validation")){
    if (!entry.is_regular_file() or entry.path().extension() != ".3map") continue;

    fs::directory_entry off_file(dir + "/mesh/" + entry.path().filename().replace_extension().string() + ".off");

    if (!off_file.exists()) return EXIT_FAILURE;

    LCC result = CGAL::HexRefinement::two_refinement(off_file.path().string(), 20, CGAL::HexRefinement::TwoRefinement::mark_intersecting_volume_with_poly);
    LCC validation;

    if (!CGAL::load_combinatorial_map(entry.path().string().c_str(), validation)) {
      std::cerr << "Could not load combinatorial map from "
                << entry.path().filename()
                << std::endl;
      return EXIT_FAILURE;
    }

    if (!result.is_isomorphic_to(validation, true, false, true)) {
      std::cerr << "Resulting tworefinement from polyhedron : "
                << entry.path().filename()
                << " is not isomorphic to its validation file "
                << std::endl;
      return EXIT_FAILURE;
    }
  }

  std::cout << "All meshes have been validated successfully" << std::endl;

  return EXIT_SUCCESS;
}

int main(){
  // /!\ Run validations only in Release mode, some assertions in Debug mode are very costly
  return surface_test();
}
