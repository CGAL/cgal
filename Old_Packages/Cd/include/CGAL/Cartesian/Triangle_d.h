// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr


#ifndef CGAL_CARTESIAN_TRIANGLE_D_H
#define CGAL_CARTESIAN_TRIANGLE_D_H

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_D_H
#include <CGAL/Cartesian/redefine_names_d.h>
#endif

#ifndef CGAL_CARTESIAN_THREETUPLE_H
#include <CGAL/Threetuple.h>
#endif // CGAL_CARTESIAN_THREETUPLE_H

CGAL_BEGIN_NAMESPACE

template <class _R>
class TriangleCd
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
// This is a partial specialization
<_R,Cartesian_tag>
#endif
  : public Handle
{
public:
  typedef _R                               R;
  typedef typename R::FT                   FT;
  typedef typename R::RT                   RT;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef TriangleCd<R CGAL_CTAG>          Self;
  typedef typename R::Point_d              Point_d;
  typedef typename R::Vector_d             Vector_d;
  typedef typename R::Plane_d              Plane_d;
  typedef typename R::Aff_transformation_d Aff_transformation_d;
#else
  typedef TriangleCd<R>                         Self;
  typedef typename R::Point_d_base              Point_d;
  typedef typename R::Vector_d_base             Vector_d;
  typedef typename R::Plane_d_base              Plane_d;
  typedef typename R::Aff_transformation_d_base Aff_transformation_d;
#endif

  TriangleCd();
  TriangleCd(const Self &t);
  TriangleCd(const Point_d &p, const Point_d &q, const Point_d &r);
  ~TriangleCd();

  Self       &operator=(const Self &t);

  Point_d    vertex(int i) const;
  Point_d    operator[](int i) const;

  bool       operator==(const Self &t) const;
  bool       operator!=(const Self &t) const;
  long       id() const;

  Self       transform(const Aff_transformation_d &t) const;

  bool       has_on(const Point_d &p) const;
  bool       is_degenerate() const;

  Bbox_d     bbox() const;

private:
  _Threetuple< Point_d >*   ptr() const;
};

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_CLASS_DEFINED
#ifndef CGAL_CARTESIAN_TRIANGLE_D_C
#include <CGAL/Cartesian/Triangle_d.C>
#endif // CGAL_CARTESIAN_TRIANGLE_D_C
#endif 

#endif // CGAL_CARTESIAN_TRIANGLE_D_H
