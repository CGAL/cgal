// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Cartesian_dynamic_d.h
// source        : include/CGAL/Cartesian_dynamic_d.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr
//
// coordinator   : INRIA Sophia-Antipolis
//
// ============================================================================


#ifndef CGAL_CARTESIAN_DYNAMIC_D_H
#define CGAL_CARTESIAN_DYNAMIC_D_H

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif // CGAL_BASIC_H

#ifndef CGAL_CARTESIAN_CLASSES_H
#include <CGAL/cartesian_classes.h>
#endif // CGAL_CARTESIAN_CLASSES_H

#ifdef CGAL_CFG_NO_ADVANCED_KERNEL
  // Because we cannot use Michael's scheme, we need the wrapper classes
  // We include them (they are common to Cartesian and Homogeneous)
  #ifndef CGAL_USER_CLASSES_H
  #include <CGAL/user_classes.h>
  #endif // CGAL_USER_CLASSES_H
#endif // CGAL_CFG_NO_ADVANCED_KERNEL

#define CGAL_REP_CLASS_DEFINED
#define CGAL_CARTESIAN_CLASS_DEFINED

CGAL_BEGIN_NAMESPACE

template< class R, class _FT >
struct Cartesian_base_dynamic_d
{
    // Number types and representation tag
    typedef _FT                                 RT;
    typedef _FT                                 FT;
    typedef Cartesian_tag                       Rep_tag;

#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
    // Because of partial specialization, CGAL::Point_d<R> is selected as
    // CGAL::Point_d<Cartesian_d<FT>,Cartesian_tag>
    typedef PointCd<R>                          Point_d;
#else
    typedef PointCd<R>                          Point_d;
#endif // CGAL_CFG_NO_ADVANCED_KERNEL
};

#include <CGAL/PointCd.h>

// This class is a restricted dD geometric kernel
// It is useful only if you do not need the 3D kernel
// If you need both, you should be using Cartesian<FT>

template< class _FT >
struct Cartesian_dynamic_d
  : public Cartesian_base_dynamic_d< Cartesian_dynamic_d<_FT>, _FT >
{
    typedef Cartesian_dynamic_d<_FT>            Self;
    typedef Cartesian_base_dynamic_d<Self,_FT>  Kernel_base;

#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
    // The lasses are inherited and because of partial specialization,
    // Cartesian_d<FT>::Point_d is exactly CGAL::Point_d< Cartesian_d<FT> >
    // We still need to inherit explicitly, see Cartesian.h for explanation

    typedef Kernel_base::Point_d                Point_d;

#else
    // Now CGAL::Point_d<R> is only a wrapper around CGAL::PointCd<R>
    // It is necessary to redefine here the classes to ensure that
    // Cartesian_d<FT>::Point_d is exactly CGAL::Point_d< Cartesian_d<FT> >

    // Cartesian_d<FT>::Base is needed so that CGAL::Point_d< Cartesian_d<FT> >
    // can inherit from Cartesian_d<FT>::Kernel_base::Point_d

    typedef Kernel_base::Point_d                Point_d_base;

    // Note: necessary to qualify Point_d by ::CGAL:: to disambiguate between
    // Point_d in the current namespace (nested within CGAL) and
    // CGAL::Point_d< Cartesian_d<FT> > (which is in the CGAL namespace)

    typedef CGAL::Point_d<Self>                 Point_d;

    // TODO: cleanup
    static   FT make_FT(const RT & num, const RT& denom) { return num/denom;}
    static   FT make_FT(const RT & num)                  { return num;}
    static   RT FT_numerator(const FT &r)                { return r;}
    static   RT FT_denominator(const FT &)               { return RT(1);}

#endif // CGAL_CFG_NO_ADVANCED_KERNEL
};

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_DYNAMIC_D_H
