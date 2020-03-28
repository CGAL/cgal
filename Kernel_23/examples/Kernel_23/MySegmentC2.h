#ifndef MY_SEGMENT_2_H
#define MY_SEGMENT_2_H

#include <boost/config.hpp>

template < class R_ >
class MySegmentC2
{
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_2              Point_2;
  typedef typename R_::Vector_2             Vector_2;
  typedef typename R_::Direction_2          Direction_2;
  typedef typename R_::Line_2               Line_2;
  typedef typename R_::Segment_2            Segment_2;
  typedef typename R_::Aff_transformation_2 Aff_transformation_2;

  Point_2 sp_, tp_;
public:
  typedef R_                                     R;

  MySegmentC2() {}

  MySegmentC2(const Point_2 &sp, const Point_2 &tp)
    : sp_(sp), tp_(tp) {}

  bool        is_horizontal() const;
  bool        is_vertical() const;
  bool        has_on(const Point_2 &p) const;
  bool        collinear_has_on(const Point_2 &p) const;

  bool        operator==(const MySegmentC2 &s) const;
  bool        operator!=(const MySegmentC2 &s) const;

  const Point_2 &   source() const
  {
      return sp_;
  }

  const Point_2 &   target() const
  {
      return tp_;
  }
  const Point_2 &    start() const;
  const Point_2 &    end() const;

  const Point_2 &   min BOOST_PREVENT_MACRO_SUBSTITUTION () const;
  const Point_2 &   max BOOST_PREVENT_MACRO_SUBSTITUTION () const;
  const Point_2 &   vertex(int i) const;
  const Point_2 &   point(int i) const;
  const Point_2 &   operator[](int i) const;

  FT          squared_length() const;

  Direction_2 direction() const;
  Vector_2    to_vector() const;
  Line_2      supporting_line() const;
  Segment_2   opposite() const;

  Segment_2   transform(const Aff_transformation_2 &t) const
  {
    return Segment_2(t.transform(source()), t.transform(target()));
  }

  bool        is_degenerate() const;
  CGAL::Bbox_2      bbox() const;
};

template < class R >
inline
bool
MySegmentC2<R>::operator==(const MySegmentC2<R> &s) const
{
  return source() == s.source() && target() == s.target();
}

template < class R >
inline
bool
MySegmentC2<R>::operator!=(const MySegmentC2<R> &s) const
{
  return !(*this == s);
}



template < class R >
CGAL_KERNEL_INLINE
const typename MySegmentC2<R>::Point_2 &
MySegmentC2<R>::min BOOST_PREVENT_MACRO_SUBSTITUTION () const
{
  typename R::Less_xy_2 less_xy;
  return less_xy(source(),target()) ? source() : target();
}

template < class R >
CGAL_KERNEL_INLINE
const typename MySegmentC2<R>::Point_2 &
MySegmentC2<R>::max BOOST_PREVENT_MACRO_SUBSTITUTION() const
{
  typename R::Less_xy_2 less_xy;
  return less_xy(source(),target()) ? target() : source();
}

template < class R >
CGAL_KERNEL_INLINE
const typename MySegmentC2<R>::Point_2 &
MySegmentC2<R>::vertex(int i) const
{
  return (i%2 == 0) ? source() : target();
}

template < class R >
CGAL_KERNEL_INLINE
const typename MySegmentC2<R>::Point_2 &
MySegmentC2<R>::point(int i) const
{
  return (i%2 == 0) ? source() : target();
}

template < class R >
inline
const typename MySegmentC2<R>::Point_2 &
MySegmentC2<R>::operator[](int i) const
{
  return vertex(i);
}

template < class R >
CGAL_KERNEL_INLINE
typename MySegmentC2<R>::FT
MySegmentC2<R>::squared_length() const
{
  typename R::Compute_squared_distance_2 squared_distance;
  return squared_distance(source(), target());
}

template < class R >
CGAL_KERNEL_INLINE
typename MySegmentC2<R>::Direction_2
MySegmentC2<R>::direction() const
{
  typename R::Construct_vector_2 construct_vector;
  return Direction_2( construct_vector( source(), target()));
}

template < class R >
CGAL_KERNEL_INLINE
typename MySegmentC2<R>::Vector_2
MySegmentC2<R>::to_vector() const
{
  typename R::Construct_vector_2 construct_vector;
  return construct_vector( source(), target());
}

template < class R >
inline
typename MySegmentC2<R>::Line_2
MySegmentC2<R>::supporting_line() const
{
  typename R::Construct_line_2 construct_line;

  return construct_line(*this);
}

template < class R >
inline
typename MySegmentC2<R>::Segment_2
MySegmentC2<R>::opposite() const
{
  return MySegmentC2<R>(target(), source());
}

template < class R >
CGAL_KERNEL_INLINE
CGAL::Bbox_2
MySegmentC2<R>::bbox() const
{
  return source().bbox() + target().bbox();
}

template < class R >
inline
bool
MySegmentC2<R>::is_degenerate() const
{
  return R().equal_y_2_object()(source(), target());
}

template < class R >
CGAL_KERNEL_INLINE
bool
MySegmentC2<R>::is_horizontal() const
{
  return R().equal_y_2_object()(source(), target());
}

template < class R >
CGAL_KERNEL_INLINE
bool
MySegmentC2<R>::is_vertical() const
{
  return R().equal_x_2_object()(source(), target());
}

template < class R >
CGAL_KERNEL_INLINE
bool
MySegmentC2<R>::
has_on(const typename MySegmentC2<R>::Point_2 &p) const
{
  return R().collinear_are_ordered_along_line_2_object()(source(), p, target());
}

template < class R >
inline
bool
MySegmentC2<R>::
collinear_has_on(const typename MySegmentC2<R>::Point_2 &p) const
{
  return R().collinear_has_on_2_object()(*this, p);
}

template < class R >
std::ostream &
operator<<(std::ostream &os, const MySegmentC2<R> &s)
{
    switch(CGAL::get_mode(os)) {
    case CGAL::IO::ASCII :
        return os << s.source() << ' ' << s.target();
    case CGAL::IO::BINARY :
        return os << s.source() << s.target();
    default:
        return os << "MySegmentC2(" << s.source() <<  ", " << s.target() << ")";
    }
}



template < class R >
std::istream &
operator>>(std::istream &is, MySegmentC2<R> &s)
{
    typename R::Point_2 p, q;

    is >> p >> q;

    if (is)
        s = MySegmentC2<R>(p, q);
    return is;
}



#endif // CARTESIAN_SEGMENT_2_H
