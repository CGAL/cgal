// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr


#ifndef CGAL_CARTESIAN_AFF_TRANSFORMATION_D_H
#define CGAL_CARTESIAN_AFF_TRANSFORMATION_D_H

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_D_H
#include <CGAL/Cartesian/redefine_names_d.h>
#endif

#include <cmath>

CGAL_BEGIN_NAMESPACE

class Identity_transformation;
template <class R> class Aff_transformation_rep_baseCd;
template <class R> class Aff_transformation_repCd;
// template <class R> class Translation_repCd;
// template <class R> class Scaling_repCd;
// template <class R> class Homothecy_repCd;

CGAL_END_NAMESPACE

#include <CGAL/Cartesian/Aff_transformation_rep_d.h>
// #include <CGAL/Cartesian/Translation_rep_d.h>
// #include <CGAL/Cartesian/Homothecy_rep_d.h>
// #include <CGAL/Cartesian/Scaling_rep_d.h>

CGAL_BEGIN_NAMESPACE

template < class _R >
class Aff_transformationCd
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
// This is a partial specialization
<_R,Cartesian_tag>
#endif
  : public Handle
{
  friend class PlaneCd<_R CGAL_CTAG>;

public:
  typedef _R                               R;
  typedef typename R::FT                   FT;
  typedef typename R::FT                   RT;
  typedef Aff_transformation_rep_baseCd<R> Aff_t_base;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef Aff_transformationCd<R,Cartesian_tag> Self;
  typedef typename R::Point_d              Point_d;
  typedef typename R::Vector_d             Vector_d;
  typedef typename R::Direction_d          Direction_d;
  typedef typename R::Plane_d              Plane_d;
#else
  typedef Aff_transformationCd<R>          Self;
  typedef typename R::Point_d_base         Point_d;
  typedef typename R::Vector_d_base        Vector_d;
  typedef typename R::Direction_d_base     Direction_d;
  typedef typename R::Plane_d_base         Plane_d;
#endif

  Aff_transformationCd();
  // Aff_transformationCd(const Self &t); // Provided by default

  // Identity constructor:
  Aff_transformationCd(const Identity_transformation &);

  // Translation:
  Aff_transformationCd(const Translation,
                       const Vector_d &v);

  // Scaling:
  Aff_transformationCd(const Scaling,
                       const FT &s,
                       const FT &w = FT(1));
  template < class InputIterator >
  Aff_transformationCd(const Scaling,
                       const InputIterator &begin, const InputIterator &last);

  // General form: without translation
  template < class InputIterator >
  Aff_transformationCd(int d,
                       const InputIterator &begin, const InputIterator &last,
                       const FT& w= FT(1));

  // General form: with translation
  Aff_transformationCd(int d,
                       const InputIterator &begin, const InputIterator &last,
                       const InputIterator &translation_begin,
		       const InputIterator &translation_last,
                       const FT& w = FT(1));

  ~Aff_transformationCd();

  Point_d     transform(const Point_d &p) const { return ptr()->transform(p); }
  Point_d     operator()(const Point_d &p) const { return transform(p); }

  Vector_d    transform(const Vector_d &v) const
                                           { return ptr()->transform(v); }
  Vector_d    operator()(const Vector_d &v) const { return transform(v); }

  Direction_d transform(const Direction_d &d) const
                                              { return ptr()->transform(d); }
  Direction_d operator()(const Direction_d &d) const { return transform(d); }

  Plane_d     transform(const Plane_d& p) const { return p.transform(*this); }
  Plane_d     operator()(const Plane_d& p) const { return transform(l); }

  Self        inverse() const { return ptr()->inverse(); }
  
  bool        is_even() const { return ptr()->is_even(); }
  bool        is_odd() const { return  ! (ptr()->is_even()); }
  
  FT          cartesian(int i, int j) const { return ptr()->cartesian(i,j); }
  FT          homogeneous(int i, int j) const { return cartesian(i,j); }
  FT          m(int i, int j) const { return cartesian(i,j); }
  FT          hm(int i, int j) const { return cartesian(i,j); }

  Self        operator*(const Self &t) const { return (*ptr()) * (*t.ptr()); }

protected:
  Self        transpose() const { return ptr()->transpose(); }

private:
  Aff_t_base* ptr() const { return  (Aff_t_base*)PTR; }
};

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_CLASS_DEFINED
#ifndef CGAL_CARTESIAN_AFF_TRANSFORMATION_d_C
#include <CGAL/Cartesian/Aff_transformation_d.C>
#endif // CGAL_CARTESIAN_AFF_TRANSFORMATION_d_C
#endif 

#endif // CGAL_CARTESIAN_AFF_TRANSFORMATION_D_H
