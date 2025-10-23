#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/variational_medial_axis_sampling.h>
#include <CGAL/extract_variational_medial_skeleton.h>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::Surface_mesh<Point> Mesh;
typedef CGAL::Variational_medial_axis<Mesh> VMAS;
typedef CGAL::Medial_skeleton<Mesh> Skeleton;

// Test parameters structure
struct TestParams
{
  std::string mesh_name;
  double lambda;
  int num_spheres;
};

// Mimic the "check_value_equal" function in Surface Mesh Skeletonization
template <class T> bool check_value_close(T a, T b, T tolerance = T(1e-6), const std::string& description = "") {
  if(std::abs(a - b) > tolerance) {
    std::cerr << "Value not close enough! " << description << "\n"
              << " Expected: " << b << ", Got: " << a << "\n"
              << " line " << __LINE__ << "\n"
              << " file " << __FILE__ << "\n"
              << ", Difference: " << std::abs(a - b) << "\n";
    return false;
  }
  return true;
}

bool compare_skeletons(const Skeleton& skeleton1, const Skeleton& skeleton2, double tolerance = 1e-6) {
  // Compare vertex count
  if(skeleton1.number_of_vertices() != skeleton2.number_of_vertices()) {
    std::cerr << "Vertex count mismatch: " << skeleton1.number_of_vertices() << " vs " << skeleton2.number_of_vertices()
              << std::endl;
    return false;
  }

  // Compare edge count
  if(skeleton1.number_of_edges() != skeleton2.number_of_edges()) {
    std::cerr << "Edge count mismatch: " << skeleton1.number_of_edges() << " vs " << skeleton2.number_of_edges()
              << std::endl;
    return false;
  }

  // Compare face count
  if(skeleton1.number_of_faces() != skeleton2.number_of_faces()) {
    std::cerr << "Face count mismatch: " << skeleton1.number_of_faces() << " vs " << skeleton2.number_of_faces()
              << std::endl;
    return false;
  }

  // Compare vertices
  const auto& vertices1 = skeleton1.vertices();
  const auto& vertices2 = skeleton2.vertices();

  for(std::size_t i = 0; i < vertices1.size(); ++i) {
    const auto& sphere1 = vertices1[i];
    const auto& sphere2 = vertices2[i];

    // Compare sphere centers
    auto center1 = sphere1.center();
    auto center2 = sphere2.center();

    if(!check_value_close(center1.x(), center2.x(), tolerance, "sphere center x") ||
       !check_value_close(center1.y(), center2.y(), tolerance, "sphere center y") ||
       !check_value_close(center1.z(), center2.z(), tolerance, "sphere center z"))
    {
      return false;
    }

    // Compare sphere radii
    auto radius1 = CGAL::sqrt(sphere1.squared_radius());
    auto radius2 = CGAL::sqrt(sphere2.squared_radius());

    if(!check_value_close(radius1, radius2, tolerance, "sphere radius")) {
      return false;
    }
  }

  // Compare edges
  const auto& edges1 = skeleton1.edges();
  const auto& edges2 = skeleton2.edges();

  for(std::size_t i = 0; i < edges1.size(); ++i) {
    const auto& edge1 = edges1[i];
    const auto& edge2 = edges2[i];
    bool same_edge = (edge1.first == edge2.first && edge1.second == edge2.second) ||
                     (edge1.first == edge2.second && edge1.second == edge2.first);
    if(!same_edge) {
      std::cerr << "Edge " << i << " mismatch: (" << edge1.first << "," << edge1.second << ") vs (" << edge2.first
                << "," << edge2.second << ")" << std::endl;
      return false;
    }
  }

  // Compare faces
  const auto& faces1 = skeleton1.faces();
  const auto& faces2 = skeleton2.faces();

  for(std::size_t i = 0; i < faces1.size(); ++i) {
    std::array<std::size_t, 3> sorted_face1 = {faces1[i][0], faces1[i][1], faces1[i][2]};
    std::array<std::size_t, 3> sorted_face2 = {faces2[i][0], faces2[i][1], faces2[i][2]};

    std::sort(sorted_face1.begin(), sorted_face1.end());
    std::sort(sorted_face2.begin(), sorted_face2.end());

    if(sorted_face1 != sorted_face2) {
      std::cerr << "Face " << i << " mismatch: "
                << "[" << sorted_face1[0] << "," << sorted_face1[1] << "," << sorted_face1[2] << "] vs "
                << "[" << sorted_face2[0] << "," << sorted_face2[1] << "," << sorted_face2[2] << "]" << std::endl;
      return false;
    }
  }

  return true;
}

bool test_correcteness() {

  Mesh mesh;
  const std::string filename(CGAL::data_file_path("meshes/elephant.off"));

  if(!CGAL::IO::read_polygon_mesh(filename, mesh)) {
    std::cerr << "Cannot open" << filename << std::endl;
    return false;
  }

  auto skeleton_bvh = CGAL::extract_variational_medial_skeleton(
      mesh, CGAL::parameters::lambda(0.2).random_seed(42).number_of_spheres(300).acceleration_structure(CGAL::BVH_tag{}));
  auto skeleton_kd_tree = CGAL::extract_variational_medial_skeleton(
      mesh, CGAL::parameters::lambda(0.2).random_seed(42).number_of_spheres(300).acceleration_structure(CGAL::KD_tree_tag{}));

  std::string saved_file_bvh = CGAL::data_file_path("meshes/elephant_0.2_300_seed_42_BVH.ply");
  std::string saved_file_kd_tree = CGAL::data_file_path("meshes/elephant_0.2_300_seed_42_KD_tree.ply");

  Skeleton saved_skeleton_bvh,saved_skeleton_kd_tree;

  if(saved_skeleton_bvh.read_skeleton_from_PLY(saved_file_bvh)) {
    if(compare_skeletons(skeleton_bvh, saved_skeleton_bvh)) {
      std::cout << "Skeletons match after reading from file." << std::endl;
    } else {
      std::cerr << " Skeletons do not match after reading from file." << std::endl;
      return false;
    }
  }

  if(saved_skeleton_kd_tree.read_skeleton_from_PLY(saved_file_kd_tree)) {
    if(compare_skeletons(skeleton_kd_tree, saved_skeleton_kd_tree)) {
      std::cout << "KD-tree Skeletons match after reading from file." << std::endl;
    } else {
      std::cerr << "KD-tree  Skeletons do not match after reading from file." << std::endl;
      return false;
    }
  }
  return true;
}

bool test_API() {
  Mesh mesh;
  const std::string filename(CGAL::data_file_path("meshes/chair.off"));

  if(!CGAL::IO::read_polygon_mesh(filename, mesh)) {
    std::cerr << "Cannot open" << filename << std::endl;
    return false;
  }
  VMAS vmas(mesh);
  bool success = vmas.sample(
      CGAL::parameters::lambda(0.2).number_of_spheres(200).concurrency_tag(CGAL::Sequential_tag{}));
  if(!success)
  {
    std::cerr << "compute_variational_medial_axis_sampling failed." << std::endl;
    return false;
  } else {
    std::cout << "compute_variational_medial_axis_sampling pass." << std::endl;
  }
  auto skeleton = vmas.export_skeleton();
  success = vmas.add_spheres(2);
  auto skeleton2 = vmas.export_skeleton();
  if(!success || skeleton2.number_of_vertices() != skeleton.number_of_vertices() + 2) {
    std::cerr << "Adding spheres failed: expected " << skeleton.number_of_vertices() + 2 << " vertices, got "
              << skeleton2.number_of_vertices() << std::endl;
    return false;
  } else {
    std::cout << "Adding spheres succeeded: new vertex count is " << skeleton2.number_of_vertices() << std::endl;
  }
  vmas.add_sphere_by_id(5);
  auto skeleton3 = vmas.skeleton();
  if(skeleton3.number_of_vertices() != skeleton.number_of_vertices() + 3) {
    std::cerr << "Adding spheres failed: expected " << skeleton.number_of_vertices() + 3 << " vertices, got "
              << skeleton3.number_of_vertices() << std::endl;
    return false;
  } else {
    std::cout << "Adding spheres succeeded: new vertex count is " << skeleton3.number_of_vertices() << std::endl;
  }
  vmas.remove_sphere_by_id(0);
  vmas.remove_sphere_by_id(1);
  vmas.remove_sphere_by_id(2);
  auto skeleton4 = vmas.skeleton();
  if(skeleton4.number_of_vertices() != skeleton.number_of_vertices()) {
    std::cerr << "Removing sphere failed: expected " << skeleton.number_of_vertices() << " vertices, got "
              << skeleton4.number_of_vertices() << std::endl;
    return false;
  } else {
    std::cout << "Removing sphere succeeded: new vertex count is " << skeleton4.number_of_vertices() << std::endl;
  }
  return true;
}

// Test determinism (same parameters should produce same results)
bool test_determinism(const TestParams& params, const std::string& mesh_file_path) {
  std::cout << "Testing determinism for " << params.mesh_name << "..." << std::endl;

  Mesh mesh;
  if(!CGAL::IO::read_polygon_mesh(mesh_file_path, mesh)) {
    std::cerr << "Cannot read mesh file: " << mesh_file_path << std::endl;
    return false;
  }
  const int num_runs = 3;
  std::vector<Skeleton> skeletons;
  skeletons.reserve(num_runs);

  // Run algorithm multiple times with same parameters
  for(int i = 0; i < num_runs; ++i) {
    VMAS vmas(mesh);
    vmas.sample(CGAL::parameters::lambda(params.lambda)
                .number_of_spheres(params.num_spheres)
                //.verbose(true)
                .concurrency_tag(CGAL::Sequential_tag{}));

    skeletons.push_back(vmas.skeleton());
  }

  // Compare all results with first one
  for(int i = 1; i < num_runs; ++i) {
    if(!compare_skeletons(skeletons[0], skeletons[i])) {
      std::cerr << "Determinism test failed: Run " << i << " differs from run 0" << std::endl;
      return false;
    }
  }
  return true;
}

int main() {
  std::vector<TestParams> test_cases = {{"chair", 0.1, 100}, {"chair", 0.2, 150}, {"bug", 0.1, 100}, {"bug", 0.2, 150}};

  std::cout << "=== Variational Medial Axis Sampling Tests ===" << std::endl;

  std::cout << "\n--- Test 1: Correcteness ---" << std::endl;
  if(test_correcteness()) {
    std::cout << "Correcteness test passed for elephant" << std::endl;
  } else {
    std::cerr << "Correcteness test failed for elephant" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "\n--- Test 2: API---" << std::endl;
  if(test_API()) {
    std::cout << "API test passed for chair" << std::endl;
  } else {
    std::cerr << "API test failed for chair" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "\n--- Test 3: Determinism ---" << std::endl;
  for(const auto& params : test_cases) {
    std::string mesh_file = CGAL::data_file_path("meshes/" + params.mesh_name + ".off");

    // Check if mesh file exists
    std::ifstream test_file(mesh_file);
    if(!test_file) {
      std::cout << "Skipping tests for " << mesh_file << ": mesh file not found" << std::endl;
      continue;
    }
    test_file.close();

    std::cout << "\n--- Test: Determinism for " << params.mesh_name << ", lambda: " << params.lambda
              << ", num_spheres: " << params.num_spheres << " ---" << std::endl;
    if(test_determinism(params, mesh_file)) {
      std::cout << "Determinism test passed for " << params.mesh_name << ", lambda: " << params.lambda
                << ", num_spheres: " << params.num_spheres << std::endl;
    } else {
      std::cout << "Determinism test failed for " << params.mesh_name << ", lambda: " << params.lambda
                << ", num_spheres: " << params.num_spheres << std::endl;
      return EXIT_FAILURE;
    }
  }
  std::cout << "All tests passed!" << std::endl;
  return EXIT_SUCCESS;
}
