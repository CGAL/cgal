// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
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
// file          : include/CGAL/Cartesian_3.h
// source        : include/CGAL/Cartesian_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr
//
// coordinator   : INRIA Sophia-Antipolis
//
// ============================================================================


#ifndef CGAL_CARTESIAN_3_H
#define CGAL_CARTESIAN_3_H

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
struct Cartesian_base_3
{
    // Number types and representation tag
    typedef _FT                                   RT;
    typedef _FT                                   FT;
    typedef Cartesian_tag                         Rep_tag;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
    typedef CGAL::Point_3<R,Rep_tag>              Point_3;
    typedef CGAL::Vector_3<R,Rep_tag>             Vector_3;
    typedef CGAL::Direction_3<R,Rep_tag>          Direction_3;
    typedef CGAL::Line_3<R,Rep_tag>               Line_3;
    typedef CGAL::Plane_3<R,Rep_tag>              Plane_3;
    typedef CGAL::Ray_3<R,Rep_tag>                Ray_3;
    typedef CGAL::Segment_3<R,Rep_tag>            Segment_3;
    typedef CGAL::Triangle_3<R,Rep_tag>           Triangle_3;
    typedef CGAL::Tetrahedron_3<R,Rep_tag>        Tetrahedron_3;
    typedef CGAL::Aff_transformation_3<R,Rep_tag> Aff_transformation_3;
#else
    typedef PointC3<R>                          Point_3;
    typedef VectorC3<R>                         Vector_3;
    typedef DirectionC3<R>                      Direction_3;
    typedef LineC3<R>                           Line_3;
    typedef PlaneC3<R>                          Plane_3;
    typedef RayC3<R>                            Ray_3;
    typedef SegmentC3<R>                        Segment_3;
    typedef TriangleC3<R>                       Triangle_3;
    typedef TetrahedronC3<R>                    Tetrahedron_3;
    typedef Aff_transformationC3<R>             Aff_transformation_3;
#endif // CGAL_CFG_NO_ADVANCED_KERNEL
};

CGAL_END_NAMESPACE
 
#include <CGAL/Cartesian/Point_3.h>
#include <CGAL/Cartesian/Vector_3.h>
#include <CGAL/Cartesian/Direction_3.h>
#include <CGAL/Cartesian/Line_3.h>
#include <CGAL/Cartesian/Plane_3.h>
#include <CGAL/Cartesian/Ray_3.h>
#include <CGAL/Cartesian/Segment_3.h>
#include <CGAL/Cartesian/Triangle_3.h>
#include <CGAL/Cartesian/Tetrahedron_3.h>
#include <CGAL/Cartesian/Aff_transformation_3.h>

#include <CGAL/Cartesian/global_operators_3.h>
#include <CGAL/Cartesian/predicates_on_points_3.h>
#include <CGAL/Cartesian/distance_predicates_3.h>
#include <CGAL/Cartesian/basic_constructions_3.h>

#include <CGAL/Cartesian/Point_3.C>
#include <CGAL/Cartesian/Vector_3.C>
#include <CGAL/Cartesian/Direction_3.C>
#include <CGAL/Cartesian/Line_3.C>
#include <CGAL/Cartesian/Plane_3.C>
#include <CGAL/Cartesian/Ray_3.C>
#include <CGAL/Cartesian/Segment_3.C>
#include <CGAL/Cartesian/Triangle_3.C>
#include <CGAL/Cartesian/Tetrahedron_3.C>
#include <CGAL/Cartesian/Aff_transformation_3.C>

CGAL_BEGIN_NAMESPACE

// This class is a restricted 3D geometric kernel
// It is useful only if you do not need the 2D kernel
// If you need both, you should be using Cartesian<FT>

template< class _FT >
struct Cartesian_3 :
  public Cartesian_base_3< Cartesian_3<_FT>, _FT >
{
    // Number types and representation tag
    typedef _FT                                 RT;
    typedef _FT                                 FT;
    typedef Cartesian_tag                       Rep_tag;

    typedef Cartesian_3<_FT>                    Self;
    typedef Cartesian_base_3<Self,_FT>          Base;

#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
    // The other classes are inherited and because of partial specialization,
    // Cartesian_3<FT>::Point_3 is exactly CGAL::Point_3< Cartesian_3<FT> >
    // We still need to inherit explicitly, see Cartesian.h for explanation

    typedef Base::Point_3                       Point_3;
    typedef Base::Vector_3                      Vector_3;
    typedef Base::Direction_3                   Direction_3;
    typedef Base::Line_3                        Line_3;
    typedef Base::Plane_3                       Plane_3;
    typedef Base::Ray_3                         Ray_3;
    typedef Base::Segment_3                     Segment_3;
    typedef Base::Triangle_3                    Triangle_3;
    typedef Base::Tetrahedron_3                 Tetrahedron_3;
    typedef Base::Aff_transformation_3          Aff_transformation_3;

 #else
    // Now CGAL::Point_3<R> is only a wrapper around CGAL::PointC3<R>
    // It is necessary to redefine here the classes to ensure that
    // Cartesian_3<FT>::Point_3 is exactly CGAL::Point_3< Cartesian_3<FT> >

    // Cartesian_3<FT>::Base is needed so that CGAL::Point_3< Cartesian_3<FT> >
    // can inherit from Cartesian_3<FT>::Base::Point_3

    typedef Base::Point_3                       Point_3_base;
    typedef Base::Vector_3                      Vector_3_base;
    typedef Base::Direction_3                   Direction_3_base;
    typedef Base::Line_3                        Line_3_base;
    typedef Base::Plane_3                       Plane_3_base;
    typedef Base::Ray_3                         Ray_3_base;
    typedef Base::Segment_3                     Segment_3_base;
    typedef Base::Triangle_3                    Triangle_3_base;
    typedef Base::Tetrahedron_3                 Tetrahedron_3_base;
    typedef Base::Aff_transformation_3          Aff_transformation_3_base;

    // Note: necessary to qualify Point_3 by CGAL:: to disambiguate between
    // Point_2 in the current namespace (nested within CGAL) and
    // CGAL::Point_3< Cartesian_3<FT> > (which is in the CGAL namespace)

    typedef CGAL::Point_3<Self>                 Point_3;
    typedef CGAL::Vector_3<Self>                Vector_3;
    typedef CGAL::Direction_3<Self>             Direction_3;
    typedef CGAL::Line_3<Self>                  Line_3;
    typedef CGAL::Plane_3<Self>                 Plane_3;
    typedef CGAL::Ray_3<Self>                   Ray_3;
    typedef CGAL::Segment_3<Self>               Segment_3;
    typedef CGAL::Triangle_3<Self>              Triangle_3;
    typedef CGAL::Tetrahedron_3<Self>           Tetrahedron_3;
    typedef CGAL::Aff_transformation_3<Self>    Aff_transformation_3;

    // TODO: cleanup
    static   FT make_FT(const RT & num, const RT& denom) { return num/denom;}
    static   FT make_FT(const RT & num)                  { return num;}
    static   RT FT_numerator(const FT &r)                { return r;}
    static   RT FT_denominator(const FT &)               { return RT(1);}

#endif // CGAL_CFG_NO_ADVANCED_KERNEL
};

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_H
