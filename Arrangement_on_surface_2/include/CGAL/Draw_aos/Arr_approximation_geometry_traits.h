#ifndef CGAL_DRAW_AOS_ARR_APPROXIMATION_GEOMETRY_TRAITS_H
#define CGAL_DRAW_AOS_ARR_APPROXIMATION_GEOMETRY_TRAITS_H

#include "CGAL/Simple_cartesian.h"

namespace CGAL {
class Arr_approximation_geometry_traits
{
public:
  using Approximation_kernel = Simple_cartesian<double>;
  using Approx_point = Approximation_kernel::Point_2;
  using FT = double;
  using Point_geom = Approx_point;
  using Apporx_point_vec = std::vector<Point_geom>;
  using Polyline_geom = Apporx_point_vec;
  using Triangle = std::array<std::size_t, 3>;
  using Triangle_vec = std::vector<Triangle>;
  struct Triangulated_face
  {
    Apporx_point_vec points;
    Triangle_vec triangles;
  };
};

} // namespace CGAL
#endif // CGAL_DRAW_AOS_ARR_APPROXIMATION_GEOMETRY_TRAITS_H