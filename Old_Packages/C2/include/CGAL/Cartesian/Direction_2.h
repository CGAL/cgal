// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Herve Bronnimann

#ifndef CGAL_CARTESIAN_DIRECTION_2_H
#define CGAL_CARTESIAN_DIRECTION_2_H

#include <CGAL/Cartesian/redefine_names_2.h>
#include <CGAL/Twotuple.h>

CGAL_BEGIN_NAMESPACE

template < class _R >
class DirectionC2
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
// This is a partial specialization
<_R,Cartesian_tag>
#endif
  : public Handle
{
public:
  typedef _R                                    R;
  typedef typename R::FT                        FT;
  typedef typename R::RT                        RT;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef DirectionC2<R,Cartesian_tag>          Self;
  typedef typename R::Point_2                   Point_2;
  typedef typename R::Vector_2                  Vector_2;
  typedef typename R::Line_2                    Line_2;
  typedef typename R::Ray_2                     Ray_2;
  typedef typename R::Triangle_2                Triangle_2;
  typedef typename R::Segment_2                 Segment_2;
  typedef typename R::Iso_rectangle_2           Iso_rectangle_2;
  typedef typename R::Aff_transformation_2      Aff_transformation_2;
  typedef typename R::Circle_2                  Circle_2;
#else
  typedef DirectionC2<R>                        Self;
  typedef typename R::Point_2_base              Point_2;
  typedef typename R::Vector_2_base             Vector_2;
  typedef typename R::Line_2_base               Line_2;
  typedef typename R::Ray_2_base                Ray_2;
  typedef typename R::Triangle_2_base           Triangle_2;
  typedef typename R::Segment_2_base            Segment_2;
  typedef typename R::Iso_rectangle_2_base      Iso_rectangle_2;
  typedef typename R::Aff_transformation_2_base Aff_transformation_2;
  typedef typename R::Circle_2_base             Circle_2;
#endif

  DirectionC2();
  DirectionC2(const Self &d);
  DirectionC2(const Vector_2 &v);
  DirectionC2(const FT &x, const FT &y);
  ~DirectionC2();

  Self     &operator=(const Self &d);

  bool     operator==(const Self &d) const;
  bool     operator!=(const Self &d) const;
  bool     operator>=(const Self &d) const;
  bool     operator<=(const Self &d) const;
  bool     operator>(const Self &d) const;
  bool     operator<(const Self &d) const;
  bool     counterclockwise_in_between( const Self &d1, const Self &d2) const;
  int      id() const;

  Vector_2 to_vector() const;

  Self     perpendicular(const Orientation &o) const;
  Self     transform(const Aff_transformation_2 &t) const;

  Self     operator-() const;

  FT       delta(int i) const;
  FT       dx() const;
  FT       dy() const;

private:
  _Twotuple<FT>*   ptr() const;
};

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_CLASS_DEFINED
#include <CGAL/Cartesian/Direction_2.C>
#endif

#endif // CGAL_CARTESIAN_DIRECTION_2_H
