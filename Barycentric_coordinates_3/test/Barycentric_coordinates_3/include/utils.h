#ifndef CGAL_BARYCENTRIC_TESTS_UTILS_H
#define CGAL_BARYCENTRIC_TESTS_UTILS_H

// STL includes
#include <vector>
#include <array>

// CGAL includes
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/point_generators_3.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>

namespace tests{

  template<typename FT>
  FT get_tolerance() {
    return FT(1.0) / FT(10000000000.0);
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
  std::pair<Mesh, std::vector<typename Kernel::Point_3>> get_regular_tetrahedron(
    typename Kernel::FT scale){

    using Point_3 = typename Kernel::Point_3;

    Mesh tetrahedron0;
    std::vector<Point_3> coords = {Point_3(scale, scale, scale), Point_3(-scale, scale, -scale),
                                   Point_3(scale, -scale, -scale), Point_3(-scale, -scale, scale)};

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

    unsigned int num_coords = coords.size();
    assert(num_coords != 0);

    for(auto& coord : coords)
      assert(coord == FT(1.0)/FT(num_coords));
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

  template<typename Kernel, typename OutIterator>
  OutIterator random_points_tetrahedron(
    std::vector<typename Kernel::Point_3>& coords,
    OutIterator out, int n){

    using Tetrahedron_3 = typename Kernel::Tetrahedron_3;
    using Mesh = typename CGAL::Surface_mesh<typename Kernel::Point_3>;

    CGAL_assertion(coords.size() == 4);
    Tetrahedron_3 tetra_inside(coords[0], coords[1], coords[2], coords[3]);
    Mesh tetra_surf;
    CGAL::make_tetrahedron(coords[0], coords[1],
                           coords[2], coords[3], tetra_surf);

    CGAL::Random_points_in_tetrahedron_3<typename Kernel::Point_3> gen_in(tetra_inside);
    CGAL::Random_points_in_triangle_mesh_3<Mesh> gen_surf(tetra_surf);
    out = std::copy_n(gen_in, n/2, out);
    out = std::copy_n(gen_surf, n/2, out);

    return out;
  }
}

#endif
