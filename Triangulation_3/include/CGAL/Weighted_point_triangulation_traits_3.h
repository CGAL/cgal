
#ifndef CGAL_WEIGHTED_POINT_TRIANGULATION_TRAITS_3_H
#define CGAL_WEIGHTED_POINT_TRIANGULATION_TRAITS_3_H

namespace CGAL {

  template <typename K>
  struct Weighted_point_triangulation_traits_3 : public K {

    Weighted_point_triangulation_traits_3(const K& k)
      : K(k)
    {}

    typedef typename K::Weighted_point_3 Point_3;
    typedef typename K::Point_3 Bare_point;
  };
  
} // namespace CGAL

#endif // CGAL_WEIGHTED_POINT_TRIANGULATION_TRAITS_3_H
