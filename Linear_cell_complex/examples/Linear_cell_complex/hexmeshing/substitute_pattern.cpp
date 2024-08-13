#include "hexmeshing.h"
#include "utils.h"

#include <CGAL/Combinatorial_map_save_load.h>
#include <CGAL/config.h>
#include <CGAL/draw_polyhedron.h>
#include <CGAL/Graphics_scene_options.h>
#include <cstdlib>
#include <filesystem>
#include <iostream>

auto default_two_refinement(const std::string& path, int cube_cells_per_dim, int nb_levels = 1){
  using namespace CGAL::HexRefinement;
  using namespace CGAL::HexRefinement::TwoRefinement;
  Tree aabb;
  // Meh
  return two_refinement_mt(path, 2,
    [&](Polyhedron& poly){
      simple_voxelisation_setup(aabb)(poly);
      auto bbox = aabb.bbox();

      Point center = {bbox.xmin() + (bbox.x_span()/2),
                      bbox.ymin() + (bbox.y_span()/2),
                      bbox.zmin() + (bbox.z_span()/2)};

      double max_size = std::max(std::max(bbox.x_span(), bbox.y_span()), bbox.z_span());
      return Grid::make_centered_cube(center, max_size / cube_cells_per_dim, cube_cells_per_dim);
    },
    mark_intersecting_volume_with_poly(aabb),
    nb_levels
  );
}

// auto default_two_refinement(const std::string& path, int cube_cells_per_dim, int nb_levels = 1){
//   using namespace CGAL::HexRefinement;
//   using namespace CGAL::HexRefinement::TwoRefinement;
//   Tree aabb;
//   // Meh
//   return two_refinement(path,
//     [&](Polyhedron& poly){
//       simple_voxelisation_setup(aabb)(poly);
//       auto bbox = aabb.bbox();

//       Point center = {bbox.xmin() + (bbox.x_span()/2),
//                       bbox.ymin() + (bbox.y_span()/2),
//                       bbox.zmin() + (bbox.z_span()/2)};

//       double max_size = std::max(std::max(bbox.x_span(), bbox.y_span()), bbox.z_span());
//       return Grid::make_centered_cube(center, max_size / cube_cells_per_dim, cube_cells_per_dim);
//     },
//     mark_intersecting_volume_with_poly(aabb),
//     nb_levels
//   );
// }

int surface_test() {
  using namespace CGAL::HexRefinement;

  auto [lcc, _]= default_two_refinement(
    CGAL::data_file_path("query_replace/mesh1.off"), 10);

  lcc.display_characteristics(std::cout);

  // render_two_refinement_result(lcc, aabb);
  CGAL::draw(lcc);

  return EXIT_SUCCESS;
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

  HexMeshingData hdata;
  RefinementData rdata;

  // initial_setup(hdata, cellIdentifier);

  render(lcc, m1, m2);

  return 1;
}

using InputFile = const std::string;
using OutputFile = const std::string;
using GeneratorFunction = std::function<LCC(const std::string&)>;
using ValidationData = std::tuple<InputFile, OutputFile, GeneratorFunction>;

ValidationData validation_data[] = {
  {"mesh1.off", "mesh1-1r-20x20x20", [](const std::string& path){
    auto tuple = default_two_refinement(path, 20);
    return std::get<LCC>(tuple);
  }},
  {"mesh2.off", "mesh2-1r-20x20x20", [](const std::string& path){
    auto tuple = default_two_refinement(path, 20);
    return std::get<LCC>(tuple);
  }},
  {"mesh3.off", "mesh3-1r-20x20x20", [](const std::string& path){
    auto tuple = default_two_refinement(path, 20);
    return std::get<LCC>(tuple);
  }},

  // 2 Levels of refinement
  {"mesh1.off", "mesh1-2r-20x20x20", [](const std::string& path){
    auto tuple = default_two_refinement(path, 20, 2);
    return std::get<LCC>(tuple);
  }},

  // 3 Levels of refinement
  // Do not try, it just takes ages to save the combinatorial map
  // {"mesh1.off", "mesh1-3r-20x20x20", [](const std::string& path){
  //   auto tuple = CGAL::HexRefinement::two_refinement(path, 20, CGAL::HexRefinement::TwoRefinement::mark_intersecting_volume_with_poly, 3);
  //   return std::get<LCC>(tuple);
  // }},
};

int create_validation_files(){
  namespace fs = std::filesystem;
  const fs::path dir = CGAL::data_file_path("hexmeshing");
  const fs::path meshDir = dir  / "mesh";
  const fs::path validationDir = dir  / "validation";

  for (auto& [in_file, out_file, generator] : validation_data){
    fs::directory_entry entry(meshDir / in_file);
    if (!entry.is_regular_file()) {
      std::cerr << "File not found " << in_file << ", search path : " << dir / in_file;
      return EXIT_FAILURE;
    };

    LCC lcc = generator(entry.path().string());

    fs::path out = validationDir / (out_file + ".3map");

    if (!CGAL::save_combinatorial_map(lcc, out.string().c_str())){
      std::cerr << "Error while saving result of tworefinement for " << out.string() << std::endl;
      return EXIT_FAILURE;
    }
  }

  std::cout << "All validation files have been created" << std::endl;

  return EXIT_SUCCESS;
}

int validation_test(){
  namespace fs = std::filesystem;
  const fs::path dir = CGAL::data_file_path("hexmeshing");
  const fs::path meshDir = dir  / "mesh";
  const fs::path validationDir = dir  / "validation";

  for (auto& [in_file, out_file, generator] : validation_data){
    fs::directory_entry validationEntry(validationDir / (out_file + ".3map"));
    fs::directory_entry meshEntry(meshDir / in_file);

    if (!validationEntry.is_regular_file()) {
      std::cerr << "File not found " << out_file << ".3map, search path : " << dir / (out_file + ".3map") ;
      return EXIT_FAILURE;
    }
    if (!meshEntry.is_regular_file()) {
      std::cerr << "File not found " << in_file << ", search path : " << dir / in_file;
      return EXIT_FAILURE;
    }

    LCC result = generator(meshEntry.path().string());
    LCC validation;

    if (!CGAL::load_combinatorial_map(validationEntry.path().string().c_str(), validation)) {
      std::cerr << "Could not load combinatorial map from "
                << validationEntry.path().filename()
                << std::endl;
      return EXIT_FAILURE;
    }

    if (!result.is_isomorphic_to(validation, false, false, false)) {
      std::cerr << "Resulting tworefinement from polyhedron : "
                << meshEntry.path().filename()
                << " is not isomorphic to its validation file "
                << validationEntry.path().filename()
                << std::endl;
      return EXIT_FAILURE;
    }
  }

  std::cout << "All meshes have been validated successfully" << std::endl;

  return EXIT_SUCCESS;
}

// Multi threading : Test if two decomposition leads to the same result
int test_deterministic_mt_refinement(){
  return false;
}

int main(){
  // /!\ Run validations only in Release mode, some assertions in Debug mode are very costly
  // /!\ Github cannot store >50Mb files, and thus we cannot track the .3map on github,
  // To run validations, first run create_validation_files in a older commit, also check if the
  // Generator function and data is still the same. Then fast forward to latest commit and run validation_test()
  return validation_test();
}
