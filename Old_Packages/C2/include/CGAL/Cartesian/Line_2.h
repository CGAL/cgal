// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Herve Bronnimann

#ifndef CGAL_CARTESIAN_LINE_2_H
#define CGAL_CARTESIAN_LINE_2_H

#include <CGAL/Cartesian/redefine_names_2.h>
#include <CGAL/Threetuple.h>

CGAL_BEGIN_NAMESPACE

template < class _R >
class LineC2
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
  typedef LineC2<R,Cartesian_tag>               Self;
  typedef typename R::Point_2                   Point_2;
  typedef typename R::Vector_2                  Vector_2;
  typedef typename R::Direction_2               Direction_2;
  typedef typename R::Ray_2                     Ray_2;
  typedef typename R::Triangle_2                Triangle_2;
  typedef typename R::Segment_2                 Segment_2;
  typedef typename R::Iso_rectangle_2           Iso_rectangle_2;
  typedef typename R::Aff_transformation_2      Aff_transformation_2;
  typedef typename R::Circle_2                  Circle_2;
#else
  typedef LineC2<R>                             Self;
  typedef typename R::Point_2_base              Point_2;
  typedef typename R::Vector_2_base             Vector_2;
  typedef typename R::Direction_2_base          Direction_2;
  typedef typename R::Ray_2_base                Ray_2;
  typedef typename R::Triangle_2_base           Triangle_2;
  typedef typename R::Segment_2_base            Segment_2;
  typedef typename R::Iso_rectangle_2_base      Iso_rectangle_2;
  typedef typename R::Aff_transformation_2_base Aff_transformation_2;
  typedef typename R::Circle_2_base             Circle_2;
#endif

  LineC2();
  LineC2(const Self  &l);
  LineC2(const Point_2 &p, const Point_2 &q);
  LineC2(const FT &a, const FT &b, const FT &c);
  LineC2(const Segment_2 &s);
  LineC2(const Ray_2 &r);
  LineC2(const Point_2 &p, const Direction_2 &d);
  ~LineC2();

  Self            &operator=(const Self &l);

  bool            operator==(const Self &l) const;
  bool            operator!=(const Self &l) const;
  int             id() const;

  FT              a() const;
  FT              b() const;
  FT              c() const;

  FT              x_at_y(const FT &y) const;
  FT              y_at_x(const FT &x) const;

  Self            perpendicular(const Point_2 &p) const;
  Self            opposite() const;
  Point_2         point(int i) const;

  Point_2         point() const;
  Point_2         projection(const Point_2 &p) const;

  Direction_2     direction() const;

  Oriented_side   oriented_side(const Point_2 &p) const;
  bool            has_on_boundary(const Point_2 &p) const;
  bool            has_on_positive_side(const Point_2 &p) const;
  bool            has_on_negative_side(const Point_2 &p) const;
  bool            has_on(const Point_2 &p) const { return has_on_boundary(p); }

  bool            is_horizontal() const;
  bool            is_vertical() const;
  bool            is_degenerate() const;

  Self            transform(const Aff_transformation_2 &t) const;

private:
  _Threetuple<FT>* ptr() const;
  void            new_rep(const Point_2 &p, const Point_2 &q);
  void            new_rep(const FT &a, const FT &b, const FT &c);
};

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_CLASS_DEFINED
#include <CGAL/Cartesian/Line_2.C>
#endif 

#endif // CGAL_CARTESIAN_LINE_2_H
