// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Herve Bronnimann

#ifndef CGAL_CARTESIAN_CIRCLE_REP_2_H
#define CGAL_CARTESIAN_CIRCLE_REP_2_H

CGAL_BEGIN_NAMESPACE

template < class R >
class Circle_repC2 : public Rep
{
public:
  typedef typename R::FT                        FT;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef typename R::Point_2                   Point_2;
#else
  typedef typename R::Point_2_base              Point_2;
#endif

  Point_2      center;
  FT           squared_radius;
  Orientation  orient;


  Circle_repC2() {}

  Circle_repC2(const Point_2 & c, const FT & r, const Orientation &o)
    : center(c), squared_radius(r), orient(o) {}

  ~Circle_repC2() {}
};

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_CIRCLE_REP_2_H
