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
// file          : include/CGAL/Cartesian/Aff_transformation_d.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Brönnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_AFF_TRANSFORMATION_D_C
#define CGAL_CARTESIAN_AFF_TRANSFORMATION_D_C

#include <CGAL/Cartesian/Aff_transformation_rep_d.C>
#include <iostream>

#ifndef CGAL_CTAG
#define CGAL_CTAG
#endif

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

CGAL_BEGIN_NAMESPACE

#ifndef CGAL_NO_OSTREAM_INSERT_AFF_TRANSFORMATIONCD
template < class R >
std::ostream &
operator<<(std::ostream &os,
           const Aff_transformationCd<R CGAL_CTAG> &t)
{
  t.print(os);
  return os;
}
#endif // CGAL_NO_OSTREAM_INSERT_AFF_TRANSFORMATIONCD

#ifndef CGAL_NO_ISTREAM_EXTRACT_AFF_TRANSFORMATIONCD
template < class R >
std::istream &
operator>>(std::istream &is,
           Aff_transformationCd<R CGAL_CTAG> &t)
{
  // FIXME : TODO
  return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_AFF_TRANSFORMATIONCD

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

#endif // CGAL_CARTESIAN_AFF_TRANSFORMATION_D_C
