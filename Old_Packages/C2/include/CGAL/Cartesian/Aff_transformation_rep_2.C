// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Lutz Kettner

#ifndef CGAL_CARTESIAN_AFF_TRANSFORMATION_REP_2_C
#define CGAL_CARTESIAN_AFF_TRANSFORMATION_REP_2_C

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

CGAL_BEGIN_NAMESPACE

template < class R >
CGAL_KERNEL_LARGE_INLINE
Aff_transformation_repC2<R>::Aff_transformation_2
Aff_transformation_repC2<R>::
inverse() const
{
  FT det = FT(1) / (t11 * t22 - t12 * t21);
  return Aff_transformation_2(
    det * t22,    det * (-t12), det * (t12*t23-t13*t22),
    det * (-t21), det * t11 ,   det * (t13*t21-t11*t23));
}

template < class R >
CGAL_KERNEL_INLINE
typename Aff_transformation_repC2<R>::Aff_transformation_2
Aff_transformation_repC2<R>::
operator*(const Aff_transformation_rep_baseC2<R> &t) const
{
  return t.compose(*this);
}
 
template < class R >
CGAL_KERNEL_LARGE_INLINE
typename Aff_transformation_repC2<R>::Aff_transformation_2
Aff_transformation_repC2<R>::
compose(const Aff_transformation_repC2<R> &t) const
{
  return Aff_transformation_2(t.t11*t11 + t.t12*t21,
                              t.t11*t12 + t.t12*t22,
                              t.t11*t13 + t.t12*t23 + t.t13,
                              t.t21*t11 + t.t22*t21,
                              t.t21*t12 + t.t22*t22,
                              t.t21*t13 + t.t22*t23 + t.t23 );
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
typename Aff_transformation_repC2<R>::Aff_transformation_2
Aff_transformation_repC2<R>::
compose(const Translation_repC2<R> &t) const
{
  return Aff_transformation_2(t11,
                              t12,
                              t13 + t._translationvector.x(),
                              t21,
                              t22,
                              t23 + t._translationvector.y());
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
typename Aff_transformation_repC2<R>::Aff_transformation_2
Aff_transformation_repC2<R>::
compose(const Rotation_repC2<R> &t) const
{
  return Aff_transformation_2(t._cosinus*t11 - t._sinus*t21,
                              t._cosinus*t12 - t._sinus*t22,
                              t._cosinus*t13 - t._sinus*t23,
                              t._sinus*t11 + t._cosinus*t21,
                              t._sinus*t12 + t._cosinus*t22,
                              t._sinus*t13 + t._cosinus*t23);
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
typename Aff_transformation_repC2<R>::Aff_transformation_2
Aff_transformation_repC2<R>::
compose(const Scaling_repC2<R> &t) const
{
   return Aff_transformation_2(t._scalefactor * t11,
                               t._scalefactor * t12,
                               t._scalefactor * t13,
                               t._scalefactor * t21,
                               t._scalefactor * t22,
                               t._scalefactor * t23);
}

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

#endif // CGAL_CARTESIAN_AFF_TRANSFORMATION_REP_2_C
