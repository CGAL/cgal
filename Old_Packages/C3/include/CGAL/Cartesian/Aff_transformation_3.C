#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_3_H
#define CGAL_CTAG
#endif

#ifndef CGAL_CARTESIAN_AFF_TRANSFORMATION_3_C
#define CGAL_CARTESIAN_AFF_TRANSFORMATION_3_C

#ifndef CGAL_CARTESIAN_AFF_TRANSFORMATION_REP_3_C
#include <CGAL/Cartesian/Aff_transformation_rep_3.C>
#endif // CGAL_CARTESIAN_AFF_TRANSFORMATION_REP_3_C

#ifndef CGAL_CARTESIAN_DETERMINANT_H
#include <CGAL/determinant.h>
#endif // CGAL_CARTESIAN_DETERMINANT_H

CGAL_BEGIN_NAMESPACE

template < class R >
Aff_transformationC3<R CGAL_CTAG>::
Aff_transformationC3()
{
  FT ft1(1), ft0(0);
  PTR = new Aff_transformation_repC3<R>(ft1, ft0, ft0,
                                        ft0, ft1, ft0,
                                        ft0, ft0, ft1);
}

template < class R >
Aff_transformationC3<R CGAL_CTAG>::
Aff_transformationC3(const Identity &)
{
  FT ft1(1), ft0(0);
  PTR = new Aff_transformation_repC3<R>(ft1, ft0, ft0,
                                        ft0, ft1, ft0,
                                        ft0, ft0, ft1);
}

template < class R >
Aff_transformationC3<R CGAL_CTAG>::
Aff_transformationC3(const Translation,
                     const typename Aff_transformationC3<R CGAL_CTAG>::Vector_3 &v)
{
  PTR = new Translation_repC3<R>(v);
}

template < class R >
Aff_transformationC3<R CGAL_CTAG>::
Aff_transformationC3(const Scaling,
                     const typename R::FT &s, const typename R::FT &w)
{
  if (w != FT(1))
    PTR = new Scaling_repC3<R>(s/w);
  else
    PTR = new Scaling_repC3<R>(s);
}


template < class R >
Aff_transformationC3<R CGAL_CTAG>::
Aff_transformationC3(
  const typename R::FT& m11, const typename R::FT& m12,
  const typename R::FT& m13, const typename R::FT& m14,
  const typename R::FT& m21, const typename R::FT& m22,
  const typename R::FT& m23, const typename R::FT& m24,
  const typename R::FT& m31, const typename R::FT& m32,
  const typename R::FT& m33, const typename R::FT& m34,
  const typename R::FT &w)
{
  if (w != FT(1))
    PTR = new Aff_transformation_repC3<R>(m11/w, m12/w, m13/w, m14/w,
                                          m21/w, m22/w, m23/w, m24/w,
                                          m31/w, m32/w, m33/w, m34/w);
  else
    PTR = new Aff_transformation_repC3<R>(m11, m12, m13, m14,
                                          m21, m22, m23, m24,
                                          m31, m32, m33, m34);
}

template < class R >
Aff_transformationC3<R CGAL_CTAG>::Aff_transformationC3(
  const typename R::FT& m11, const typename R::FT& m12, const typename R::FT& m13,
  const typename R::FT& m21, const typename R::FT& m22, const typename R::FT& m23,
  const typename R::FT& m31, const typename R::FT& m32, const typename R::FT& m33,
  const typename R::FT &w)
{
  if (w != FT(1))
    PTR = new Aff_transformation_repC3<R>(m11/w, m12/w, m13/w,
                                          m21/w, m22/w, m23/w,
                                          m31/w, m32/w, m33/w);
  else
    PTR = new Aff_transformation_repC3<R>(m11, m12, m13,
                                          m21, m22, m23,
                                          m31, m32, m33);
}

template < class R >
Aff_transformationC3<R CGAL_CTAG>::~Aff_transformationC3()
{}
template < class R >
typename Aff_transformationC3<R CGAL_CTAG>::Point_3
Aff_transformationC3<R CGAL_CTAG>::
transform(const typename Aff_transformationC3<R CGAL_CTAG>::Point_3 &p) const
{
  return ptr()->transform(p);
}

template < class R >
inline
typename Aff_transformationC3<R CGAL_CTAG>::Point_3
Aff_transformationC3<R CGAL_CTAG>::
operator()(const typename Aff_transformationC3<R CGAL_CTAG>::Point_3 &p) const
{
  return transform(p);
}

template < class R >
inline
typename Aff_transformationC3<R CGAL_CTAG>::Vector_3
Aff_transformationC3<R CGAL_CTAG>::
transform(const typename Aff_transformationC3<R CGAL_CTAG>::Vector_3 &v) const
{
  return ptr()->transform(v);
}

template < class R >
inline
typename Aff_transformationC3<R CGAL_CTAG>::Vector_3
Aff_transformationC3<R CGAL_CTAG>::
operator()(const typename Aff_transformationC3<R CGAL_CTAG>::Vector_3 &v) const
{
  return transform(v);
}

template < class R >
inline
typename Aff_transformationC3<R CGAL_CTAG>::Direction_3
Aff_transformationC3<R CGAL_CTAG>::
transform(const typename Aff_transformationC3<R CGAL_CTAG>::Direction_3 &d) const
{
  return ptr()->transform(d);
}

template < class R >
inline
typename Aff_transformationC3<R CGAL_CTAG>::Direction_3
Aff_transformationC3<R CGAL_CTAG>::
operator()(const typename Aff_transformationC3<R CGAL_CTAG>::Direction_3& d) const
{
  return transform(d);
}

template < class R >
inline
typename Aff_transformationC3<R CGAL_CTAG>::Plane_3
Aff_transformationC3<R CGAL_CTAG>::
transform(const typename Aff_transformationC3<R CGAL_CTAG>::Plane_3& p) const
{
  return p.transform(*this);
}

template < class R >
inline
typename Aff_transformationC3<R CGAL_CTAG>::Plane_3
Aff_transformationC3<R CGAL_CTAG>::
operator()(const typename Aff_transformationC3<R CGAL_CTAG>::Plane_3& p) const
{
  return transform(p);
}

template < class R >
Aff_transformationC3<R CGAL_CTAG>
Aff_transformationC3<R CGAL_CTAG>::inverse() const
{
  return ptr()->inverse();
}


template < class R >
Aff_transformationC3<R CGAL_CTAG>
Aff_transformationC3<R CGAL_CTAG>::general_form() const
{
  return ptr()->general_form();
}

template < class R >
Aff_transformationC3<R CGAL_CTAG>
Aff_transformationC3<R CGAL_CTAG>::transpose() const
{
  return ptr()->transpose();
}

template < class R >
bool  Aff_transformationC3<R CGAL_CTAG>::is_even() const
{
  return ptr()->is_even();
}

template < class R >
bool  Aff_transformationC3<R CGAL_CTAG>::is_odd() const
{
  return ! (ptr()->is_even());
}



template < class R >
Aff_transformationC3<R CGAL_CTAG>
_general_transformation_composition(const Aff_transformation_repC3<R> &l,
                                    const Aff_transformation_repC3<R> &r )
{
return Aff_transformationC3<R CGAL_CTAG>(
            l.t11*r.t11 + l.t12*r.t21 + l.t13*r.t31,
            l.t11*r.t12 + l.t12*r.t22 + l.t13*r.t32,
            l.t11*r.t13 + l.t12*r.t23 + l.t13*r.t33,
            l.t11*r.t14 + l.t12*r.t24 + l.t13*r.t34 + l.t14,

            l.t21*r.t11 + l.t22*r.t21 + l.t23*r.t31,
            l.t21*r.t12 + l.t22*r.t22 + l.t23*r.t32,
            l.t21*r.t13 + l.t22*r.t23 + l.t23*r.t33,
            l.t21*r.t14 + l.t22*r.t24 + l.t23*r.t34 + l.t24,

            l.t31*r.t11 + l.t32*r.t21 + l.t33*r.t31,
            l.t31*r.t12 + l.t32*r.t22 + l.t33*r.t32,
            l.t31*r.t13 + l.t32*r.t23 + l.t33*r.t33,
            l.t31*r.t14 + l.t32*r.t24 + l.t33*r.t34 + l.t34);
}


// this is really quick and dirty.  As in the 2D case the composition of
// translations or scalings should be a translation or a scaling
template < class R >
Aff_transformationC3<R CGAL_CTAG>
operator*(const Aff_transformationC3<R CGAL_CTAG> &a,
          const Aff_transformationC3<R CGAL_CTAG> &b)
{
  return _general_transformation_composition(
         *(Aff_transformation_repC3<R>*)(a.general_form().ptr()),
         *(Aff_transformation_repC3<R>*)(b.general_form().ptr()));
}


#ifndef CGAL_NO_OSTREAM_INSERT_AFF_TRANSFORMATIONC3
template < class R >
std::ostream &operator<<(std::ostream &os,
                         const Aff_transformationC3<R CGAL_CTAG> &t)
{
    t.print(os);
    return os;
}
#endif // CGAL_NO_OSTREAM_INSERT_AFF_TRANSFORMATIONC3

#ifndef CGAL_NO_ISTREAM_EXTRACT_AFF_TRANSFORMATIONC3
#endif // CGAL_NO_ISTREAM_EXTRACT_AFF_TRANSFORMATIONC3

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_AFF_TRANSFORMATION_3_C
