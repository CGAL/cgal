#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/boundary_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/Discrete_harmonic_coordinates_2.h>

// Typedefs.
using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;

using FT      = Kernel::FT;
using Point_2 = Kernel::Point_2;

struct Info {

  Info(const std::string _name) :
  name(_name) { }
  std::string name;
};

using Vertex       = std::pair<Point_2, Info>;
using Point_map    = CGAL::First_of_pair_property_map<Vertex>;
using Vertex_range = std::vector<Vertex>;

using Discrete_harmonic_coordinates_2 =
  CGAL::Barycentric_coordinates::Discrete_harmonic_coordinates_2<Vertex_range, Kernel, Point_map>;
using Policy =
  CGAL::Barycentric_coordinates::Computation_policy_2;

int main() {

  Kernel kernel;
  Point_map point_map;

  // Construct a unit square.
  const std::vector<Vertex> square = {
    std::make_pair(Point_2(0, 0), Info("1")), std::make_pair(Point_2(1, 0), Info("2")),
    std::make_pair(Point_2(1, 1), Info("3")), std::make_pair(Point_2(0, 1), Info("4"))
  };

  // Construct the class with discrete harmonic weights.
  // We do not check for edge cases since we know the exact positions
  // of all our points. We speed up the computation by using the O(n) algorithm.
  const Policy policy = Policy::FAST;
  Discrete_harmonic_coordinates_2 discrete_harmonic_2(square, policy, kernel, point_map);

  // Construct the center point of the unit square.
  const Point_2 center(FT(1) / FT(2), FT(1) / FT(2));

  // Compute discrete harmonic weights for the center point.
  std::list<FT> weights;
  discrete_harmonic_2.weights(center, std::back_inserter(weights));

  std::cout << std::endl << "discrete harmonic weights (center): ";
  for (const FT& weight : weights) {
    std::cout << weight << " ";
  }
  std::cout << std::endl;

  // Compute discrete harmonic coordinates for the center point.
  std::list<FT> coordinates;
  discrete_harmonic_2(center, std::back_inserter(coordinates));

  std::cout << std::endl << "discrete harmonic coordinates (center): ";
  for (const FT& coordinate : coordinates) {
    std::cout << coordinate << " ";
  }
  std::cout << std::endl;

  // Construct several interior points.
  const std::vector<Point_2> interior_points = {
    Point_2(FT(1) / FT(5), FT(1) / FT(5)),
    Point_2(FT(4) / FT(5), FT(1) / FT(5)),
    Point_2(FT(4) / FT(5), FT(4) / FT(5)),
    Point_2(FT(1) / FT(5), FT(4) / FT(5)) };

  // Compute discrete harmonic weights for all interior points.
  std::cout << std::endl << "discrete harmonic weights (interior): " << std::endl << std::endl;

  std::vector<FT> ws;
  for (const auto& query : interior_points) {
    ws.clear();
    discrete_harmonic_2.weights(query, std::back_inserter(ws));
    for (std::size_t i = 0; i < ws.size() - 1; ++i) {
      std::cout << ws[i] << ", ";
    }
    std::cout << ws[ws.size() - 1] << std::endl;
  }

  // Compute discrete harmonic coordinates for all interior point.
  std::cout << std::endl << "discrete harmonic coordinates (interior): " << std::endl << std::endl;

  std::vector<FT> bs;
  for (const auto& query : interior_points) {
    bs.clear();
    discrete_harmonic_2(query, std::back_inserter(bs));
    for (std::size_t i = 0; i < bs.size() - 1; ++i) {
      std::cout << bs[i] << ", ";
    }
    std::cout << bs[bs.size() - 1] << std::endl;
  }

  // Construct 2 boundary points on the second and fourth edges.
  const Point_2 e2(1, FT(4) / FT(5));
  const Point_2 e4(0, FT(4) / FT(5));

  // Compute discrete harmonic coordinates = boundary coordinates
  // for these 2 points one by one.
  coordinates.clear();
  CGAL::Barycentric_coordinates::boundary_coordinates_2(
    square, e2, std::back_inserter(coordinates), kernel, point_map);
  CGAL::Barycentric_coordinates::boundary_coordinates_2(
    square, e4, std::back_inserter(coordinates), kernel, point_map);

  std::cout << std::endl << "boundary coordinates (edge 2 and edge 4): ";
  for (const FT& coordinate : coordinates) {
    std::cout << coordinate << " ";
  }
  std::cout << std::endl;

  // Construct 6 other boundary points: 2 on the first and third edges respectively
  // and 4 at the vertices.
  const std::vector<Point_2> es13 = {
    Point_2(FT(1) / FT(2), 0), // edges
    Point_2(FT(1) / FT(2), 1),

    // vertices
    Point_2(0, 0), Point_2(1, 0),
    Point_2(1, 1), Point_2(0, 1)
  };

  // Compute discrete harmonic coordinates = boundary coordinates for all 6 points.
  std::cout << std::endl << "boundary coordinates (edge 1, edge 3, and vertices): " << std::endl << std::endl;

  for (const auto& query : es13) {
    bs.clear();
    CGAL::Barycentric_coordinates::boundary_coordinates_2(
      square, query, std::back_inserter(bs), point_map); // we can skip kernel here
    for (std::size_t i = 0; i < bs.size() - 1; ++i) {
      std::cout << bs[i] << ", ";
    }
    std::cout << bs[bs.size() - 1] << std::endl;
  }

  // Construct 2 points outside the unit square - one from the left and one from the right.
  // Even if discrete harmonic coordinates may not be valid for some exterior points,
  // we can still do it.
  const Point_2 l(FT(-1) / FT(2), FT(1) / FT(2));
  const Point_2 r(FT(3)  / FT(2), FT(1) / FT(2));

  // Compute discrete harmonic coordinates for all exterior points.
  coordinates.clear();
  discrete_harmonic_2(l, std::back_inserter(coordinates));
  discrete_harmonic_2(r, std::back_inserter(coordinates));

  std::cout << std::endl << "discrete harmonic coordinates (exterior): ";
  for (const FT& coordinate : coordinates) {
    std::cout << coordinate << " ";
  }
  std::cout << std::endl << std::endl;

  return EXIT_SUCCESS;
}
