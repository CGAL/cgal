// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Herve Bronnimann

#ifndef CGAL_CARTESIAN_RAY_2_H
#define CGAL_CARTESIAN_RAY_2_H

#include <CGAL/Cartesian/redefine_names_2.h>
#include <CGAL/Twotuple.h>

CGAL_BEGIN_NAMESPACE

template < class _R >
class RayC2
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
  typedef RayC2<R,Cartesian_tag>                Self;
  typedef typename R::Point_2                   Point_2;
  typedef typename R::Vector_2                  Vector_2;
  typedef typename R::Direction_2               Direction_2;
  typedef typename R::Line_2                    Line_2;
  typedef typename R::Triangle_2                Triangle_2;
  typedef typename R::Segment_2                 Segment_2;
  typedef typename R::Iso_rectangle_2           Iso_rectangle_2;
  typedef typename R::Aff_transformation_2      Aff_transformation_2;
  typedef typename R::Circle_2                  Circle_2;
#else
  typedef RayC2<R>                              Self;
  typedef typename R::Point_2_base              Point_2;
  typedef typename R::Vector_2_base             Vector_2;
  typedef typename R::Direction_2_base          Direction_2;
  typedef typename R::Line_2_base               Line_2;
  typedef typename R::Triangle_2_base           Triangle_2;
  typedef typename R::Segment_2_base            Segment_2;
  typedef typename R::Iso_rectangle_2_base      Iso_rectangle_2;
  typedef typename R::Aff_transformation_2_base Aff_transformation_2;
  typedef typename R::Circle_2_base             Circle_2;
#endif

  RayC2();
  RayC2(const Self &r);
  RayC2(const Point_2 &sp, const Point_2 &secondp);
  RayC2(const Point_2 &sp, const Direction_2 &d);
  ~RayC2();

  Self        &operator=(const Self &r);

  bool        operator==(const Self &r) const;
  bool        operator!=(const Self &r) const;
  int         id() const;

  Point_2     start() const;
  Point_2     source() const;
  Point_2     point(int i) const;
  Point_2     second_point() const;

  Direction_2 direction() const;
  Line_2      supporting_line() const;
  Self        opposite() const;

  Self        transform(const Aff_transformation_2 &t) const;

  bool        is_horizontal() const;
  bool        is_vertical() const;
  bool        is_degenerate() const;
  bool        has_on(const Point_2 &p) const;
  bool        collinear_has_on(const Point_2 &p) const;

private:
  _Twotuple< Point_2 > *ptr() const;
};

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_CLASS_DEFINED
#include <CGAL/Cartesian/Ray_2.C>
#endif 

#endif // CGAL_CARTESIAN_RAY_2_H
