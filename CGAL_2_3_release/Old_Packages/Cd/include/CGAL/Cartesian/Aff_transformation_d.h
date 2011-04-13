// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Cartesian/Aff_transformation_d.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_AFF_TRANSFORMATION_D_H
#define CGAL_CARTESIAN_AFF_TRANSFORMATION_D_H

#include <CGAL/Cartesian/redefine_names_d.h>
#include <cmath>

CGAL_BEGIN_NAMESPACE

class Identity_transformation;
class Reflexion;
template <class R> class Aff_transformation_rep_baseCd;
template <class R> class Aff_transformation_repCd;
//TODO template <class R> class Translation_repCd;
//TODO template <class R> class Homothecy_repCd;
//TODO template <class R> class Scaling_repCd;
//TODO template <class R> class Reflexion_repCd;

CGAL_END_NAMESPACE

#include <CGAL/Cartesian/Aff_transformation_rep_d.h>
//TODO #include <CGAL/Cartesian/Translation_rep_d.h>
//TODO #include <CGAL/Cartesian/Homothecy_rep_d.h>
//TODO #include <CGAL/Cartesian/Scaling_rep_d.h>
//TODO #include <CGAL/Cartesian/Reflexion_rep_d.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class Aff_transformationCd CGAL_ADVANCED_KERNEL_PARTIAL_SPEC
  : public Handle
{
  friend class PlaneCd<R_ CGAL_CTAG>;

public:
  typedef R_                               R;
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

  Aff_transformationCd()
    {
      PTR = new Aff_transformation_repCd<R>();
    }

  // Identity constructor:
  Aff_transformationCd(const Identity_transformation &, int d)
    {
      typename R::LA::Matrix M(d);
      std::fill(M.begin(),M.end(),FT(0));
      int i; for (i=0; i<d; ++i) M[i][i] = FT(1);
      PTR = new Aff_transformation_repCd<R>(d,M.begin(),M.end(),FT(1));
    }

  // Translation:
  Aff_transformationCd(const Translation,
                       const Vector_d &v)
    {
      typename R::LA::Matrix M(v.dimension());
      std::fill(M.begin(),M.end(),FT(0));
      int i; for (i=0; i<v.dimension(); ++i) M[i][i] = FT(1);
      PTR = new Aff_transformation_repCd<R>(v.dimension(),
	           M.begin(),M.end(),v.begin(),v.end(),FT(1));
    }

  // Homothecy:
  Aff_transformationCd(const Scaling, int d,
                       const FT &s,
                       const FT &w = FT(1))
    {
      typename R::LA::Matrix M(d);
      std::fill(M.begin(),M.end(),FT(0));
      int i; for (i=0; i<d; ++i) M[i][i] = s/w;
      PTR = new Aff_transformation_repCd<R>(d,M.begin(),M.end(),FT(1));
    }

  // Scaling:
  template < class InputIterator >
  Aff_transformationCd(const Scaling, int d,
                       const InputIterator &first, const InputIterator &last)
    {
      CGAL_kernel_precondion( last-first==d );
      typename R::LA::Matrix M(d);
      std::fill(M.begin(),M.end(),FT(0));
      int i; InputIterator it;
      for (i=0,it=first; it!=last; ++i,++it) M[i][i] = *it;
      PTR = new Aff_transformation_repCd<R>(d,M.begin(),M.end(),FT(1));
    }

  // Reflexion
  Aff_transformationCd(const Reflexion,
                       const Plane_d &p)
    {
      typename R::LA::Matrix M(d);
      typename R::LA::Vector w(d);
      //TODO aff_transformation_reflexion_from_plane(p,M.begin(),w.begin());
      PTR = new Aff_transformation_repCd<R>(d,
	           M.begin(),M.end(),w.begin(),w.end(),FT(1));
    }

  // General form: without translation
  template < class InputIterator >
  Aff_transformationCd(int d,
                       const InputIterator &first, const InputIterator &last,
                       const FT& w= FT(1))
    {
      PTR = new Aff_transformation_repCd<R>(d, first, last, w);
    }

  // General form: with translation
  template < class InputIterator1, class InputIterator2 >
  Aff_transformationCd(int d,
                       const InputIterator1 &first, const InputIterator1 &last,
                       const InputIterator2 &translation_first,
		       const InputIterator2 &translation_last,
                       const FT& w = FT(1))
    {
      PTR = new Aff_transformation_repCd<R>(d,
	           first, last, translation_first, translation_last, w);
    }

  ~Aff_transformationCd() {}

  Point_d     transform(const Point_d &p) const { return ptr()->transform(p); }
  Point_d     operator()(const Point_d &p) const { return transform(p); }

  Vector_d    transform(const Vector_d &v) const { return ptr()->transform(v); }
  Vector_d    operator()(const Vector_d &v) const { return transform(v); }

  Direction_d transform(const Direction_d &d) const
                                              { return ptr()->transform(d); }
  Direction_d operator()(const Direction_d &d) const { return transform(d); }

  Plane_d     transform(const Plane_d& h) const { return ptr()->transform(h); }
  Plane_d     operator()(const Plane_d& h) const { return transform(p); }

  Self        inverse() const { return ptr()->inverse(); }
  
  int         dimension() const { return ptr()->dimension(); }
  bool        is_even() const { return ptr()->is_even(); }
  bool        is_odd() const { return  ! (ptr()->is_even()); }
  
  FT          cartesian(int i, int j) const { return ptr()->cartesian(i,j); }
  FT          homogeneous(int i, int j) const
                {
		  if (i==dimension()) 
		     return (j==dimension()) ? FT(1) : FT(0);
		  return cartesian(i,j); 
	        }
  FT          m(int i, int j) const { return cartesian(i,j); }
  FT          hm(int i, int j) const { return homogeneous(i,j); }

  Self        operator*(const Self &t) const { return (*ptr()) * (*t.ptr()); }

protected:
  Self        transpose() const { return ptr()->transpose(); }

private:
  Aff_t_base* ptr() const { return  (Aff_t_base*)PTR; }
};

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_AFF_TRANSFORMATION_D_H
