
// release       :
// release_date  :
//
// file          : include/CGAL/Cartesian/Tetrahedron_d.h
// source        : include/CGAL/Cartesian/Tetrahedron_d.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr
//
// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann@sophia.inria.fr)
//
// ==========================================================================


#ifndef CGAL_CARTESIAN_TETRAHEDRON_D_H
#define CGAL_CARTESIAN_TETRAHEDRON_D_H

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_D_H
#include <CGAL/Cartesian/redefine_names_d.h>
#endif

#ifndef CGAL_CARTESIAN_FOURTUPLE_H
#include <CGAL/Fourtuple.h>
#endif // CGAL_CARTESIAN_FOURTUPLE_H

CGAL_BEGIN_NAMESPACE

template <class _R>
class TetrahedronCd
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
  typedef TetrahedronCd<R CGAL_CTAG>       Self;
  typedef typename R::Point_d              Point_d;
  typedef typename R::Plane_d              Plane_d;
  typedef typename R::Aff_transformation_d Aff_transformation_d;
#else
  typedef TetrahedronCd<R>                      Self;
  typedef typename R::Point_d_base              Point_d;
  typedef typename R::Plane_d_base              Plane_d;
  typedef typename R::Aff_transformation_d_base Aff_transformation_d;
#endif

  TetrahedronCd();
  TetrahedronCd(const Self &t);
  TetrahedronCd(const Point_d &p,
                const Point_d &q,
                const Point_d &r,
                const Point_d &s);
  ~TetrahedronCd();

  Self &operator=(const Self &t);

  Point_d    vertex(int i) const;
  Point_d    operator[](int i) const;

  bool       operator==(const Self &t) const;
  bool       operator!=(const Self &t) const;
  long       id() const;

  Bbox_d     bbox() const;

  Self       transform(const Aff_transformation_d &t) const;

  bool       has_on(const Point_d &p) const;
  bool       is_degenerate() const;

private:
  _Fourtuple< Point_d >*   ptr() const;
};

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_CLASS_DEFINED
#ifndef CGAL_CARTESIAN_TETRAHEDRON_D_C
#include <CGAL/Cartesian/Tetrahedron_d.C>
#endif
#endif 

#endif
