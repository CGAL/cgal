// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann

#ifndef CGAL_CARTESIAN_AFF_TRANSFORMATIONR_EP_3_C
#define CGAL_CARTESIAN_AFF_TRANSFORMATIONR_EP_3_C

#ifndef CGAL_CARTESIANR_EDEFINE_NAMES_3_H
#define CGAL_CTAG
#endif

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

CGAL_BEGIN_NAMESPACE

template < class R >
CGAL_KERNEL_LARGE_INLINE
typename Aff_transformation_repC3<R>::Aff_transformation_3
Aff_transformation_repC3<R>::inverse() const
{
  return Aff_transformation_3(
      det2x2_by_formula( t22, t23, t32, t33),         // i 11
     -det2x2_by_formula( t12, t13, t32, t33),         // i 12
      det2x2_by_formula( t12, t13, t22, t23),         // i 13
     -det3x3_by_formula( t12, t13, t14, t22, t23, t24, t32, t33, t34 ),

     -det2x2_by_formula( t21, t23, t31, t33),         // i 21 
      det2x2_by_formula( t11, t13, t31, t33),         // i 22
     -det2x2_by_formula( t11, t13, t21, t23),         // i 23
      det3x3_by_formula( t11, t13, t14, t21, t23, t24, t31, t33, t34 ),

      det2x2_by_formula( t21, t22, t31, t32),         // i 31
     -det2x2_by_formula( t11, t12, t31, t32),         // i 32
      det2x2_by_formula( t11, t12, t21, t22),         // i 33
     -det3x3_by_formula( t11, t12, t14, t21, t22, t24, t31, t32, t34 ),

      det3x3_by_formula( t11, t12, t13, t21, t22, t23, t31, t32, t33 ));
}

template < class R >
CGAL_KERNEL_INLINE
typename Aff_transformation_repC3<R>::Aff_transformation_3
Aff_transformation_repC3<R>::
operator*(const Aff_transformation_rep_baseC3<R> &t) const
{
  return t.compose(*this);
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
typename Aff_transformation_repC3<R>::Aff_transformation_3
Aff_transformation_repC3<R>::
compose(const Aff_transformation_repC3<R> &t) const
{
  return Aff_transformation_3(t.t11*t11 + t.t12*t21 + t.t13*t31,
                              t.t11*t12 + t.t12*t22 + t.t13*t32,
                              t.t11*t13 + t.t12*t23 + t.t13*t33,
                              t.t11*t14 + t.t12*t24 + t.t13*t34 + t.t14,

                              t.t21*t11 + t.t22*t21 + t.t23*t31,
                              t.t21*t12 + t.t22*t22 + t.t23*t32,
                              t.t21*t13 + t.t22*t23 + t.t23*t33,
                              t.t21*t14 + t.t22*t24 + t.t23*t34 + t.t24,

                              t.t31*t11 + t.t32*t21 + t.t33*t31,
                              t.t31*t12 + t.t32*t22 + t.t33*t32,
                              t.t31*t13 + t.t32*t23 + t.t33*t33,
                              t.t31*t14 + t.t32*t24 + t.t33*t34 + t.t34);
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
typename Aff_transformation_repC3<R>::Aff_transformation_3
Aff_transformation_repC3<R>::
compose(const Translation_repC3<R> &t) const
{
  return Aff_transformation_3(t11,
                              t12,
                              t13,
                              t14 + t._translationvector.x(),

                              t21,
                              t22,
                              t23,
                              t24 + t._translationvector.y(),

                              t31,
                              t32,
                              t33,
                              t34 + t._translationvector.z());
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
typename Aff_transformation_repC3<R>::Aff_transformation_3
Aff_transformation_repC3<R>::
compose(const Scaling_repC3<R> &t) const
{
  return Aff_transformation_3(t._scalefactor * t11,
                              t._scalefactor * t12,
                              t._scalefactor * t13,
                              t._scalefactor * t14,
			      
                              t._scalefactor * t21,
                              t._scalefactor * t22,
                              t._scalefactor * t23,
                              t._scalefactor * t24,
                              
			      t._scalefactor * t31,
                              t._scalefactor * t32,
                              t._scalefactor * t33,
                              t._scalefactor * t34);
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
typename Aff_transformation_repC3<R>::Aff_transformation_3
Aff_transformation_repC3<R>::transpose() const
{
  FT ft0(0);
  return Aff_transformation_3( t11, t21, t31, t14,
                               t12, t22, t32, t24,
                               t13, t23, t33, t34);
}

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

#endif // CGAL_CARTESIAN_AFF_TRANSFORMATIONR_EP_3_H
