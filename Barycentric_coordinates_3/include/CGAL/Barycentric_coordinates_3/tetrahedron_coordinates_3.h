#ifndef CGAL_BARYCENTRIC_TRIANGLE_COORDINATES_3_H
#define CGAL_BARYCENTRIC_TRIANGLE_COORDINATES_3_H

#include <CGAL/Barycentric_coordinates_3/internal/utils_3.h>

namespace CGAL{
namespace Barycentric_coordinates{

  //return iterator (infer from traits)
  template<
  typename OutIterator,
  typename GeomTraits>
  OutIterator triangle_coordinates_3(
    const typename GeomTraits::Point_3& p0,
    const typename GeomTraits::Point_3& p1,
    const typename GeomTraits::Point_3& p2,
    const typename GeomTraits::Point_3& p3,
    const typename GeomTraits::Point_2& query,
    OutIterator c_begin,
    const GeomTraits& traits) {

    return internal::planar_coordinates_3(
      p0, p1, p2, p3, query, c_begin, traits);
  }

  //return iterator(infer from point_3)
  template<
  typename Point_3,
  typename OutIterator>
  OutIterator triangle_coordinates_3(
    const Point_3& p0,
    const Point_3& p1,
    const Point_3& p2,
    const Point_3& p3,
    const Point_3& query,
    OutIterator c_begin) {

    using GeomTraits = typename Kernel_traits<Point_3>::Kernel;
    const GeomTraits traits;
    return triangle_coordinates_3(
      p0, p1, p2, p3, query, c_begin, traits);
  }

  //return tuple (infer from traits)
  template<typename GeomTraits>
  std::tuple<
  typename GeomTraits::FT,
  typename GeomTraits::FT,
  typename GeomTraits::FT,
  typename GeomTraits::FT>
  triangle_coordinates_in_tuple_3(
    const typename GeomTraits::Point_3& p0,
    const typename GeomTraits::Point_3& p1,
    const typename GeomTraits::Point_3& p2,
    const typename GeomTraits::Point_3& p3,
    const typename GeomTraits::Point_3& query,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    std::vector<FT> coordinates;
    coordinates.reserve(4);
    internal::planar_coordinates_3(
      p0, p1, p2, p3, query, std::back_inserter(coordinates), traits);
    CGAL_assertion(coordinates.size() == 4);
    return std::make_tuple(coordinates[0], coordinates[1], coordinates[2], coordinates[3]);
  }

  //return tuple (infer from point_3)
  template<typename Point_3>
  std::tuple<
  typename Kernel_traits<Point_3>::Kernel::FT,
  typename Kernel_traits<Point_3>::Kernel::FT,
  typename Kernel_traits<Point_3>::Kernel::FT,
  typename Kernel_traits<Point_3>::Kernel::FT>
  triangle_coordinates_in_tuple_3(
    const Point_3& p0,
    const Point_3& p1,
    const Point_3& p2,
    const Point_3& p3,
    const Point_3& query) {

    using GeomTraits = typename Kernel_traits<Point_3>::Kernel;
    const GeomTraits traits;
    return triangle_coordinates_in_tuple_3(
      p0, p1, p2, p3, query, traits);
  }



}
}


#endif