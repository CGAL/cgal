// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri

#ifndef CGAL_CARTESIAN_SPHERE_REP_3_H
#define CGAL_CARTESIAN_SPHERE_REP_3_H

CGAL_BEGIN_NAMESPACE

template < class R >
class Sphere_repC3 : public Rep
{
public:
  typedef typename R::FT                        FT;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef typename R::Point_3                   Point_3;
#else
  typedef typename R::Point_3_base              Point_3;
#endif

  Point_3      center;
  FT           squared_radius;
  Orientation  orient;


  Sphere_repC3() {}

  Sphere_repC3(const Point_3 & c, const FT & r, const Orientation &o)
    : center(c), squared_radius(r), orient(o) {}

  ~Sphere_repC3() {}
};

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_SPHERE_REP_3_H
