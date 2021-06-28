#ifndef CGAL_BARYCENTRIC_TESTS_UTILS_H
#define CGAL_BARYCENTRIC_TESTS_UTILS_H

// STL includes
#include <vector>
#include <array>

// CGAL includes
#include <CGAL/boost/graph/helpers.h>

namespace tests{

  template<typename FT>
  FT get_tolerance() {
    return FT(1) / FT(10000000000);
  }

  template<typename Kernel, typename Mesh>
  std::pair<Mesh, std::vector<typename Kernel::Point_3>> get_irregular_tetrahedron(){

    using Point_3 = typename Kernel::Point_3;

    Mesh tetrahedron0;
    std::vector<Point_3> coords = {Point_3(0.0, 0.0, 0.0), Point_3(1.0, 0.0, 0.0),
                                   Point_3(0.0, 1.0, 0.0), Point_3(0.0, 0.0, 1.0)};

    CGAL::make_tetrahedron(coords[0], coords[1],
                           coords[2], coords[3], tetrahedron0);

    return {tetrahedron0, coords};
  }

  template<typename Kernel, typename Mesh>
  std::pair<Mesh, std::vector<typename Kernel::Point_3>> get_regular_tetrahedron(){

    using Point_3 = typename Kernel::Point_3;

    Mesh tetrahedron0;
    std::vector<Point_3> coords = {Point_3(1.0, 1.0, 1.0), Point_3(-1.0, 1.0, -1.0),
                                   Point_3(1.0, -1.0, -1.0), Point_3(-1.0, -1.0, 1.0)};

    CGAL::make_tetrahedron(coords[0], coords[1],
                           coords[2], coords[3], tetrahedron0);

    return {tetrahedron0, coords};
  }

  template<typename Kernel, typename Mesh>
  std::pair<Mesh, std::vector<typename Kernel::Point_3>> get_hexahedron(){

    using Point_3 = typename Kernel::Point_3;

    Mesh hexahedron0;
    std::vector<Point_3> coords = {Point_3(1.0, 0.0, 0.0), Point_3(1.0, 1.0, 0.0),
                              Point_3(0.0, 1.0, 0.0), Point_3(0.0, 0.0, 0.0),
                              Point_3(0.0, 0.0, 1.0), Point_3(1.0, 0.0, 1.0),
                              Point_3(1.0, 1.0, 1.0), Point_3(0.0, 1.0, 1.0)};

    CGAL::make_hexahedron(coords[0], coords[1], coords[2], coords[3], coords[4], coords[5],
                          coords[6], coords[7], hexahedron0);

    return {hexahedron0, coords};
  }

  template<typename Kernel>
  void test_partition_of_unity(std::vector<typename Kernel::FT>& coords){

    using FT = typename Kernel::FT;

    FT tol = get_tolerance<FT>();

    FT sum = FT(0);
    for(auto& coord : coords)
      sum += coord;

    assert(CGAL::abs(FT(1) - sum) < tol);
  }

  template<typename Kernel>
  void test_barycenter(std::vector<typename Kernel::FT>& coords){

    using FT = typename Kernel::FT;

    std::size_t num_coords = coords.size();
    assert(num_coords != 0);

    for(auto& coord : coords)
      assert(coord == FT(1)/FT(num_coords));
  }

  template<typename Kernel>
  void test_linear_precision(
    std::vector<typename Kernel::FT>& coords,
    const std::vector<typename Kernel::Point_3>& vertices,
    const typename Kernel::Point_3& query){

    using FT = typename Kernel::FT;

    FT tol = get_tolerance<FT>();

    std::size_t num_coords = coords.size();
    FT x_linear_comb = FT(0);
    FT y_linear_comb = FT(0);
    FT z_linear_comb = FT(0);

    for(std::size_t i = 0; i < num_coords; i++){

      x_linear_comb += vertices[i].x() * coords[i];
      y_linear_comb += vertices[i].y() * coords[i];
      z_linear_comb += vertices[i].z() * coords[i];
    }

    assert(CGAL::abs(x_linear_comb - query.x()) < tol &&
           CGAL::abs(y_linear_comb - query.y()) < tol &&
           CGAL::abs(z_linear_comb - query.z()) < tol);
  }

  template<typename Kernel>
  void test_positivity(std::vector<typename Kernel::FT>& coords){

    using FT = typename Kernel::FT;

    for(auto& coord : coords)
      assert(coord >= FT(0) && coord <= FT(1));
  }
}

#endif
