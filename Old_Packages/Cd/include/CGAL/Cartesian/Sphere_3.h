// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr


#ifndef CGAL_CARTESIAN_SPHERE_D_H
#define CGAL_CARTESIAN_SPHERE_D_H

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_D_H
#include <CGAL/Cartesian/redefine_names_d.h>
#endif
#ifndef CGAL_CARTESIAN_SPHERE_REP_D_H
#include <CGAL/Cartesian/Sphere_rep_d.h>
#endif // CGAL_CARTESIAN_SPHERE_REP_D_H

CGAL_BEGIN_NAMESPACE

template <class _R>
class SphereC3
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
// This is a partial specialization
<_R,Cartesian_tag>
#endif
{
public:
  typedef _R                                    R;
  typedef typename R::FT                        FT;
  typedef typename R::RT                        RT;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef SphereC3<R CGAL_CTAG>                 Self;
  typedef typename R::Point_d                   Point_d;
  // typedef typename R::Aff_transformation_d      Aff_transformation_d;
#else
  typedef SphereC3<R>                           Self;
  typedef typename R::Point_D_base              Point_d;
  // typedef typename R::Aff_transformation_D_base Aff_transformation_d;
#endif

  Sphere_d(const Point_d &p, const FT &s,
           const Orientation &o = COUNTERCLOCKWISE);
  // Sphere with center p, squared radius s, orientation o
  Sphere_d(const R::Point_d &p, const R::Point_d &q,
           const R::Point_d &r, const R::Point_d &u);
  // Sphere passing through p,q,r,u, oriented by p, q, r, u
  Sphere_d(const R::Point_d &p, const R::Point_d &q, const R::Point_d &r,
	   const Orientation &o = COUNTERCLOCKWISE);
  // Sphere with great circle passing through p,q,r, oriented by o
  Sphere_d(const Point_d & p, const Point_d & q,
           const Orientation &o = COUNTERCLOCKWISE);
  // Sphere with diameter pq and orientation o
  Sphere_d(const Point_d & p,
           const Orientation& o = COUNTERCLOCKWISE);
  // Sphere centered at p, radius 0, orientation o

  Point  center() const;
  // Returns the center of c
  FT     squared_radius() const;
  // Returns the square of the radius (instead of the radius itself,
  // which would require square roots)
  Orientation orientation() const;
  // Returns the orientation of c

  Self   orthogonal_transform(const Aff_transformation_d &t) const;
  //! precond: t.is_orthogonal() (*UNDEFINED*)
  // Returns the image of c by t. Since t is orthogonal, the image is
  // always a circle

  bool   is_degenerate() const;
  // A circle is degenerate if its (squared) radius is null or negative
  Self   opposite() const;
  // Returns a circle with opposite orientation

  Oriented_side  oriented_side(const Point_d &p) const;
  //! precond: ! x.is_degenerate() (when available)
  // Returns R::ON_POSITIVE_SIDE, R::ON_ORIENTED_BOUNDARY or
  // R::ON_NEGATIVE_SIDE
  bool   has_on_boundary(const Point_d &p) const
  { return oriented_side(p)==ON_ORIENTED_BOUNDARY; }
  bool   has_on_positive_side(const Point_d &p) const;
  { return oriented_side(p)==ON_POSITIVE_SIDE; }
  bool   has_on_negative_side(const Point_d &p) const;
  { return oriented_side(p)==ON_NEGATIVE_SIDE; }

  Bounded_side   bounded_side(const Point_d &p) const;
  //! precond: ! x.is_degenerate() (when available)
  // Returns R::ON_BOUNDED_SIDE, R::ON_BOUNDARY or R::ON_UNBOUNDED_SIDE
  bool   has_on_bounded_side(const Point_d &p) const;
  { return bounded_side(p)==ON_BOUNDED_SIDE; }
  bool   has_on_unbounded_side(const Point_d &p) const;
  { return bounded_side(p)==ON_UNBOUNDED_SIDE; }

protected:
  _Sphere_repC3<R> *ptr();
};

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_SPHERE_D_H
