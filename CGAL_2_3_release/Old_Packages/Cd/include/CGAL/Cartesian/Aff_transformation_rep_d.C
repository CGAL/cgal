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
// file          : include/CGAL/Cartesian/Aff_transformation_rep_d.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================


#ifndef CGAL_CARTESIAN_AFF_TRANSFORMATION_REP_D_C
#define CGAL_CARTESIAN_AFF_TRANSFORMATION_REP_D_C

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_D_H
#define CGAL_CTAG
#endif

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

CGAL_BEGIN_NAMESPACE

template < class R >
CGAL_KERNEL_LARGE_INLINE
typename Aff_transformation_repCd<R>::Aff_transformation_d
Aff_transformation_repCd<R>::inverse() const
{
  FT D;
  Matrix M( LA().inverse(_linear_transformation,D) );
  Vector v( - M * _translation_vector );
  return Aff_transformation_d(dimension(), M.begin(), M.end(),
                              v.begin(), v.end(), D);
}

template < class R >
CGAL_KERNEL_INLINE
typename Aff_transformation_repCd<R>::Aff_transformation_d
Aff_transformation_repCd<R>::
operator*(const Aff_transformation_rep_baseCd<R> &t) const
{
  return t.compose(*this);
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
typename Aff_transformation_repCd<R>::Aff_transformation_d
Aff_transformation_repCd<R>::
compose(const Aff_transformation_repCd<R> &t) const
{
  Matrix M( _linear_transformation * t._linear_transformation );
  Vector v( _translation_vector+_linear_transformation*t._translation_vector );
  return Aff_transformation_d(dimension(), M.begin(), M.end(),
                              v.begin(), v.end(), FT(1));
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
typename Aff_transformation_repCd<R>::Aff_transformation_d
Aff_transformation_repCd<R>::transpose() const
{
  Matrix M( LA().transpose(_linear_transformation) );
  return Aff_transformation_d(dimension(), M.begin(), M.end(),
            _translation_vector.begin(), _translation_vector.end(), FT(1));
}

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

#endif // CGAL_CARTESIAN_AFF_TRANSFORMATION_REP_D_H
