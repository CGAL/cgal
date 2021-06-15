#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Barycentric_coordinates_3/Wachspress_coordinates_3.h>

//Typedefs
using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;

using FT = Kernel::FT;
using Point_3 = Kernel::Point_3;
using Mesh =  CGAL::Surface_mesh<Point_3>;

// REVIEW: Look here: https://github.com/danston/cgal/blob/Weights-new_package-danston/Weights/test/Weights/test_wachspress_weights.cpp
// and test all three kernels: SCKER, EPICK, and EPECK. Do not copy the test, just how the kernels are tested.

int main(){
  // Tetrahedron Mesh
  Mesh ms;

  // Set cout precision
  std::cout.precision(20);

  // Construct tetrahedron
  const Point_3 p0(0.0, 0.0, 0.0);
  const Point_3 p1(1.0, 0.0, 0.0);
  const Point_3 p2(0.0, 1.0, 0.0);
  const Point_3 p3(0.0, 0.0, 1.0);

  CGAL::make_tetrahedron(p0, p1, p2, p3, ms);

  // Instantiate some interior, boundary, and exterior query points for which we compute coordinates
  const std::vector<Point_3> queries = {Point_3(0.25f, 0.25f, 0.25f), Point_3(0.3f, 0.2f, 0.3f),
                                        Point_3(0.1f, 0.1f, 0.1f), Point_3(0.2f, 0.5f, 0.3f),
                                        Point_3(0.5f, 0.5f, 0.5f), Point_3(-1.0f, -1.0f, 1.0f),
                                        Point_3(0.5f, 0.5f, -2.0f)
                                      };

  CGAL::Barycentric_coordinates::Wachspress_coordinates_3<Mesh, Kernel> ws(ms);
  std::vector<FT> coordinates;

  for(const auto q:queries){
    // Store results
    coordinates.clear();
    ws(q, std::back_inserter(coordinates));

    std::cout << "Coordinates: " << "(" << q << ") ---->>> " << coordinates[0] << " " <<
    coordinates[1] << " " << coordinates[2] << " " << coordinates[3] << std::endl;
  }


  // REVIEW: Actually, I do not see any tests here! The test must be something like here:
  // https://github.com/danston/cgal/blob/Barycentric_coordinates_2-danston/Barycentric_coordinates_2/test/Barycentric_coordinates_2/test_wp_triangle.cpp
  // that is it must compute both Wachspress and tetrahedron coordinates and then use assert(wp_coord == tetra_coord)!
  // This looks more like an example rather than a test.

  return EXIT_SUCCESS;
}
