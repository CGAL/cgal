#ifndef CGAL_BARYCENTRIC_TESTS_UTILS_H
#define CGAL_BARYCENTRIC_TESTS_UTILS_H

// STL includes
#include <vector>
#include <array>

// CGAL includes
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>

namespace tests{

  template<typename FT>
  FT get_tolerance() {
    return FT(1) / FT(10000000000);
  }

  template<typename Kernel, typename Mesh>
  std::pair<Mesh, std::vector<typename Kernel::Point_3>> get_irregular_tetrahedron(){

    using Point_3 = typename Kernel::Point_3;

    Mesh tetrahedron0;
    std::vector<Point_3> coords = {Point_3(1.0, 0.0, 0.0), Point_3(0.0, 1.0, 0.0),
                                   Point_3(0.0, 0.0, 0.0), Point_3(0.0, 0.0, 1.0)};

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

  template<typename Kernel, typename Mesh>
  Mesh get_regular_prism(){

    using Point_3 = typename Kernel::Point_3;
    using FT = typename Kernel::FT;

    Mesh prism0;
    CGAL::make_regular_prism(3, prism0, Point_3(0.0, 0.0, 0.0), FT(1.0), FT(1.0), true);

    return prism0;
  }

  template<
    typename Point_3,
    typename OutputIterator>
  OutputIterator sample_random_inside_tetrahedron(
    Point_3 p0,
    Point_3 p1,
    Point_3 p2,
    Point_3 p3,
    OutputIterator points,
    int num_points){

    using Kernel  = typename CGAL::Kernel_traits<Point_3>::Kernel;
    using Tetrahedron_3 = typename Kernel::Tetrahedron_3;
    using Point_generator = CGAL::Random_points_in_tetrahedron_3<Point_3> ;

    Tetrahedron_3 tr(p0, p1, p2, p3);
    Point_generator point_gen(tr);

    std::copy_n(point_gen, num_points, points);
    return points;
  }

  /*
  template<typename Kernel, typename Mesh, typename VertexToPointMap, typename OutputIterator>
  void sample_mesh_interior_points(Mesh& mesh, VertexToPointMap& vertex_to_point_map){

    using Point_3 = typename Kernel::Point_3;
    typedef CGAL::Delaunay_triangulation_3< Kernel > DT;
    typedef typename DT::Segment_simplex_iterator Segment_simplex_iterator;
    typedef CGAL::Triangulation_simplex_3<typename DT::Triangulation_data_structure> Simplex;

    std::vector<Point_3> points;

    const auto vd = vertices(mesh);
    for(auto& vertex : vd){

      Point_3 vertex_value = get(vertex_to_point_map, vertex);
      points.push_back(vertex_value);
    }

    DT dt(points.begin(), points.end());
    assert(dt.is_valid());
  }
  */

}

#endif
