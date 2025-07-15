#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/variational_medial_axis_sampling.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/IO/polygon_mesh_io.h>

#include <fstream>
#include <iostream>
#include <filesystem>
#include <cmath>

typedef CGAL::Exact_predicates_inexact_constructions_kernel               Kernel;
typedef Kernel::Point_3                                                   Point;
typedef Kernel::Vector_3                                                  Vector;
typedef CGAL::Surface_mesh<Point>                                         Mesh;
typedef CGAL::Variational_medial_axis<Mesh, Kernel>                       VMAS;
typedef CGAL::Skeleton<Mesh, Kernel>                                      Skeleton;

// Test parameters structure
struct TestParams {
    std::string mesh_name;
    double lambda;
    int num_spheres;
};

void normalize_mesh(Mesh& mesh) {
  CGAL::Bbox_3 bbox;
  for(auto v : mesh.vertices())
    bbox = bbox + mesh.point(v).bbox();

  double cx = (bbox.xmin() + bbox.xmax()) / 2.0;
  double cy = (bbox.ymin() + bbox.ymax()) / 2.0;
  double cz = (bbox.zmin() + bbox.zmax()) / 2.0;
  double max_dim = std::max({bbox.xmax() - bbox.xmin(), bbox.ymax() - bbox.ymin(), bbox.zmax() - bbox.zmin()});

  for(auto v : mesh.vertices()) {
    Point p = mesh.point(v);
    double nx = (p.x() - bbox.xmin()) / max_dim;
    double ny = (p.y() - bbox.ymin()) / max_dim;
    double nz = (p.z() - bbox.zmin()) / max_dim;
    mesh.point(v) = Point(nx, ny, nz);
  }
}

// Mimic the "check_value_equal" function in Surface Mesh Skeletonization
template <class T> 
bool check_value_close(T a, T b, T tolerance = T(1e-6), const std::string& description = "") {
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

bool compare_skeletons(const Skeleton& skeleton1, const Skeleton& skeleton2, double tolerance = 1e-6)
{
    // Compare vertex count
    if (skeleton1.number_of_vertices() != skeleton2.number_of_vertices()) {
        std::cerr << "Vertex count mismatch: " << skeleton1.number_of_vertices() 
                  << " vs " << skeleton2.number_of_vertices() << std::endl;
        return false;
    }
    
    // Compare edge count  
    if (skeleton1.number_of_edges() != skeleton2.number_of_edges()) {
        std::cerr << "Edge count mismatch: " << skeleton1.number_of_edges() 
                  << " vs " << skeleton2.number_of_edges() << std::endl;
        return false;
    }
    
    // Compare face count
    if (skeleton1.number_of_faces() != skeleton2.number_of_faces()) {
        std::cerr << "Face count mismatch: " << skeleton1.number_of_faces() 
                  << " vs " << skeleton2.number_of_faces() << std::endl;
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
    
    for (std::size_t i = 0; i < edges1.size(); ++i) {
        const auto& edge1 = edges1[i];
        const auto& edge2 = edges2[i];
        bool same_edge = (edge1.first == edge2.first && edge1.second == edge2.second) ||
                         (edge1.first == edge2.second && edge1.second == edge2.first);
        if (!same_edge) {
            std::cerr << "Edge " << i << " mismatch: (" << edge1.first 
                      << "," << edge1.second << ") vs (" 
                      << edge2.first << "," << edge2.second << ")" << std::endl;
            return false;
        }
    }
    
    // Compare faces
    const auto& faces1 = skeleton1.faces();
    const auto& faces2 = skeleton2.faces();
    
    for (std::size_t i = 0; i < faces1.size(); ++i) {
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

bool test_correcteness()
{
    
    Mesh mesh;
    const std::string filename(CGAL::data_file_path("meshes/elephant_dense.off"));
    
    if(!CGAL::IO::read_polygon_mesh(filename, mesh)) {
        std::cerr << "Cannot open" << filename << std::endl;
        return false;
    }
    normalize_mesh(mesh);
    VMAS vmas(mesh);
    vmas.init();
    vmas.compute(CGAL::parameters::lambda(0.2)
                                .number_of_spheres(300)
                                .concurrency_tag(CGAL::Sequential_tag{}));
   
    auto skeleton = vmas.export_skeleton(); 
    std::string result_file = CGAL::data_file_path("meshes/elephant_dense_0.2_300.ply");
    auto skeleton2 = vmas.read_skeleton_from_ply(result_file);
    if(compare_skeletons(skeleton, skeleton2)) {
      std::cout << "Skeletons match after reading from file." << std::endl;
    } else {
      std::cerr << " Skeletons do not match after reading from file." << std::endl;
      return false;
    }
    vmas.add_sphere(0);
    auto skeleton3 = vmas.export_skeleton();
    if(skeleton3.number_of_vertices() != skeleton.number_of_vertices() + 1) {
      std::cerr << "Adding sphere failed: expected " << skeleton.number_of_vertices() + 1 << " vertices, got "
                << skeleton3.number_of_vertices() << std::endl;
      return false;
    }
    vmas.remove_sphere(0);
    auto skeleton4 = vmas.export_skeleton();
    if(skeleton4.number_of_vertices() != skeleton.number_of_vertices()) {
      std::cerr << "Removing sphere failed: expected " << skeleton.number_of_vertices() << " vertices, got "
                << skeleton4.number_of_vertices() << std::endl;
      return false;
    }
    return true;
}

// Test determinism (same parameters should produce same results)
bool test_determinism(const TestParams& params, const std::string& mesh_file_path)
{
    std::cout << "Testing determinism for " << params.mesh_name << "..." << std::endl;
    
    Mesh mesh;
    if (!CGAL::IO::read_polygon_mesh(mesh_file_path, mesh)) {
        std::cerr << "Cannot read mesh file: " << mesh_file_path << std::endl;
        return false;
    }
    normalize_mesh(mesh);
    const int num_runs = 3;
    std::vector<Skeleton> skeletons;
    skeletons.reserve(num_runs);
    
    // Run algorithm multiple times with same parameters
    for (int i = 0; i < num_runs; ++i) {
        VMAS vmas(mesh);
        vmas.compute(CGAL::parameters::lambda(params.lambda)
                                   .number_of_spheres(params.num_spheres)
                                   .concurrency_tag(CGAL::Sequential_tag{}));
        
        skeletons.push_back(vmas.export_skeleton());
    }
    
    // Compare all results with first one
    for (int i = 1; i < num_runs; ++i) {
        if (!compare_skeletons(skeletons[0], skeletons[i])) {
            std::cerr << "Determinism test failed: Run " << i << " differs from run 0" << std::endl;
            return false;
        }
    }
    return true;
}

int main()
{
  std::vector<TestParams> test_cases = {
      {"chair", 0.1, 100}, {"chair", 0.2, 150}, {"bug", 0.1, 100}, {"bug", 0.2, 150}};
    
    int total_tests = 0;
    int passed_tests = 0;
    
    std::cout << "=== Variational Medial Axis Sampling Tests ===" << std::endl;
    
    std::cout << "\n--- Test 1: Correcteness ---" << std::endl;
    total_tests++;
    if(test_correcteness()) {
      std::cout << "Correcteness test passed for elephant_dense" << std::endl;
    } else {
        std::cout << "Correcteness test failed for elephant_dense" << std::endl;
        return EXIT_FAILURE;
    }
    for (const auto& params : test_cases) {
      std::string mesh_file = CGAL::data_file_path("meshes/"+ params.mesh_name + ".off");
        
        // Check if mesh file exists
        std::ifstream test_file(mesh_file);
        if (!test_file) {
            std::cout << "Skipping tests for " <<mesh_file << ": mesh file not found" << std::endl;
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
    return EXIT_SUCCESS;
   
}
