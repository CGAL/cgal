// ==========================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// --------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Cartesian/Aff_transformation_rep_3.h
// source        : include/CGAL/Cartesian/Aff_transformation_rep_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr
//
// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann@sophia.inria.fr)
//
// ==========================================================================

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_3_H
#define CGAL_CTAG
#endif

#ifdef _MSC_VER
#define typename
#endif

#ifndef CGAL_CARTESIAN_AFF_TRANSFORMATION_REP_3_C
#define CGAL_CARTESIAN_AFF_TRANSFORMATION_REP_3_C

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
CGAL_KERNEL_LARGE_INLINE
typename Aff_transformation_repC3<R>::Aff_transformation_3
Aff_transformation_repC3<R>::general_form() const
  {
    return Aff_transformation_3(t11, t12, t13, t14,
                                t21, t22, t23, t24,
                                t31, t32, t33, t34);
  }

template < class R >
CGAL_KERNEL_LARGE_INLINE
typename Aff_transformation_repC3<R>::Aff_transformation_3
Aff_transformation_repC3<R>::transpose() const
  {
    FT ft0(0);
    return Aff_transformation_3( t11, t21, t31, ft0,
                                 t12, t22, t32, ft0,
                                 t13, t23, t33, ft0);
  }

CGAL_END_NAMESPACE

#ifdef _MSC_VER
#undef typename
#endif

#endif // CGAL_CARTESIAN_AFF_TRANSFORMATION_REP_3_H
