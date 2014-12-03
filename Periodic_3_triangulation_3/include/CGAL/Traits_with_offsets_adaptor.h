#ifndef CGAL_TRAITS_WITH_OFFSETS_ADAPTOR_H
#define CGAL_TRAITS_WITH_OFFSETS_ADAPTOR_H

namespace CGAL {

template < class K, class Functor_ >
  class Traits_with_offsets_adaptor {
  typedef K Kernel;
  typedef Functor_ Functor;

  typedef typename Kernel::Point_3       Point;
  typedef typename Kernel::Offset        Offset;

public:
  typedef typename Kernel::Iso_cuboid_3  Iso_cuboid_3;
  typedef typename Kernel::Construct_point_3 Construct_point_3;
  typedef typename Functor::result_type result_type;

  Traits_with_offsets_adaptor(const Iso_cuboid_3 * dom) : _domain(dom) { }

  result_type operator()(const Point& p0, const Point& p1,
      const Offset& o0, const Offset& o1) const {
    return Functor()(pp(p0,o0),pp(p1,o1));
  }
  result_type operator()(const Point& p0, const Point& p1, const Point& p2,
      const Offset& o0, const Offset& o1, const Offset& o2) const {
    return Functor()(pp(p0,o0),pp(p1,o1),pp(p2,o2));
  }
  result_type operator()(const Point& p0, const Point& p1,
      const Point& p2, const Point& p3,
      const Offset& o0, const Offset& o1,
      const Offset& o2, const Offset& o3) const {
    return Functor()(pp(p0,o0),pp(p1,o1),pp(p2,o2),pp(p3,o3));
  }
  result_type operator()(const Point& p0, const Point& p1,
      const Point& p2, const Point& p3, const Point& p4,
      const Offset& o0, const Offset& o1, const Offset& o2,
      const Offset& o3, const Offset& o4) const {
    return Functor()(pp(p0,o0),pp(p1,o1),pp(p2,o2),
	pp(p3,o3),pp(p4,o4));
  }

  result_type operator()(const Point& p0, const Point& p1) const {
    return Functor()(p0, p1);
  }
  result_type operator()(const Point& p0, const Point& p1,
      const Point& p2) const {
    return Functor()(p0, p1, p2);
  }
  result_type operator()(const Point& p0, const Point& p1,
      const Point& p2, const Point& p3) const {
    return Functor()(p0, p1, p2, p3);
  }
  result_type operator()(const Point& p0, const Point& p1,
      const Point& p2, const Point& p3, const Point& p4) const {
    return Functor()(p0, p1, p2, p3, p4);
  }

protected:
  Point pp(const Point &p, const Offset &o) const {
    return Construct_point_3(*_domain)(p,o);
  }
 public:
  const Iso_cuboid_3* _domain;
};
}  // namespace CGAL

#endif /* CGAL_TRAITS_WITH_OFFSETS_ADAPTOR_H */
