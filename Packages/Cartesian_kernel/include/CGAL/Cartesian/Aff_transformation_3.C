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
// file          : include/CGAL/Cartesian/Aff_transformation_3.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_AFF_TRANSFORMATION_3_C
#define CGAL_CARTESIAN_AFF_TRANSFORMATION_3_C

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

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
Aff_transformationC3(const Identity_transformation &)
{
  FT ft1(1), ft0(0);
  PTR = new Aff_transformation_repC3<R>(ft1, ft0, ft0,
                                        ft0, ft1, ft0,
                                        ft0, ft0, ft1);
}

template < class R >
Aff_transformationC3<R CGAL_CTAG>::
Aff_transformationC3(
	const Translation,
        const typename Aff_transformationC3<R CGAL_CTAG>::Vector_3 &v)
{
  PTR = new Translation_repC3<R>(v);
}

template < class R >
Aff_transformationC3<R CGAL_CTAG>::
Aff_transformationC3(const Scaling,
                     const typename Aff_transformationC3<R CGAL_CTAG>::FT &s,
                     const typename Aff_transformationC3<R CGAL_CTAG>::FT &w)
{
  if (w != FT(1))
    PTR = new Scaling_repC3<R>(s/w);
  else
    PTR = new Scaling_repC3<R>(s);
}

template < class R >
Aff_transformationC3<R CGAL_CTAG>::
Aff_transformationC3(const typename Aff_transformationC3<R CGAL_CTAG>::FT& m11,
                     const typename Aff_transformationC3<R CGAL_CTAG>::FT& m12,
                     const typename Aff_transformationC3<R CGAL_CTAG>::FT& m13,
                     const typename Aff_transformationC3<R CGAL_CTAG>::FT& m14,
                     const typename Aff_transformationC3<R CGAL_CTAG>::FT& m21,
                     const typename Aff_transformationC3<R CGAL_CTAG>::FT& m22,
                     const typename Aff_transformationC3<R CGAL_CTAG>::FT& m23,
                     const typename Aff_transformationC3<R CGAL_CTAG>::FT& m24,
                     const typename Aff_transformationC3<R CGAL_CTAG>::FT& m31,
                     const typename Aff_transformationC3<R CGAL_CTAG>::FT& m32,
                     const typename Aff_transformationC3<R CGAL_CTAG>::FT& m33,
                     const typename Aff_transformationC3<R CGAL_CTAG>::FT& m34,
                     const typename Aff_transformationC3<R CGAL_CTAG>::FT &w)
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
Aff_transformationC3<R CGAL_CTAG>::
Aff_transformationC3(const typename Aff_transformationC3<R CGAL_CTAG>::FT& m11,
                     const typename Aff_transformationC3<R CGAL_CTAG>::FT& m12,
                     const typename Aff_transformationC3<R CGAL_CTAG>::FT& m13,
                     const typename Aff_transformationC3<R CGAL_CTAG>::FT& m21,
                     const typename Aff_transformationC3<R CGAL_CTAG>::FT& m22,
                     const typename Aff_transformationC3<R CGAL_CTAG>::FT& m23,
                     const typename Aff_transformationC3<R CGAL_CTAG>::FT& m31,
                     const typename Aff_transformationC3<R CGAL_CTAG>::FT& m32,
                     const typename Aff_transformationC3<R CGAL_CTAG>::FT& m33,
                     const typename Aff_transformationC3<R CGAL_CTAG>::FT &w)
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

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

#endif // CGAL_CARTESIAN_AFF_TRANSFORMATION_3_C
