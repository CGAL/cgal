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
  Mesh get_irregular_tetrahedron(){

    using Point_3 = typename Kernel::Point_3;

    Mesh tetrahedron0;
    CGAL::make_tetrahedron(Point_3(0.0, 0.0, 0.0), Point_3(1.0, 0.0, 0.0),
                     Point_3(0.0, 1.0, 0.0), Point_3(0.0, 0.0, 1.0), tetrahedron0);

    return tetrahedron0;
  }

  template<typename Kernel, typename Mesh>
  Mesh get_hexahedron(){

    using Point_3 = typename Kernel::Point_3;

    Mesh hexahedron0;
    CGAL::make_hexahedron(Point_3(1.0, 0.0, 0.0), Point_3(1.0, 1.0, 0.0), Point_3(0.0, 1.0, 0.0),
                          Point_3(0.0, 0.0, 0.0), Point_3(0.0, 0.0, 1.0), Point_3(1.0, 0.0, 1.0),
                          Point_3(1.0, 1.0, 1.0), Point_3(0.0, 1.0, 1.0), hexahedron0);

    return hexahedron0;
  }

  template<typename Kernel, typename Mesh>
  Mesh get_regular_prism(){

    using Point_3 = typename Kernel::Point_3;
    using FT = typename Kernel::FT;

    Mesh prism0;
    CGAL::make_regular_prism(3, prism0, Point_3(0.0, 0.0, 0.0), FT(1.0), FT(1.0), true);

    return prism0;
  }

}

#endif
