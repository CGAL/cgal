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
// file          : include/CGAL/Cartesian_2.h
// source        : include/CGAL/Cartesian_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr
//
// coordinator   : INRIA Sophia-Antipolis
//
// ============================================================================


#ifndef CGAL_CARTESIAN_2_H
#define CGAL_CARTESIAN_2_H

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
  // Forgotten in user_classes, for CGAL-2.0!!!
  // #warning "Forgotten classes in <CGAL/user_class.h> "
  // template <class R > class Circle_2;
#endif // CGAL_CFG_NO_ADVANCED_KERNEL

#define CGAL_REP_CLASS_DEFINED
#define CGAL_CARTESIAN_CLASS_DEFINED

CGAL_BEGIN_NAMESPACE

template< class R, class _FT >
struct Cartesian_base_2
{
    // Number types and representation tag
    typedef _FT                                    RT;
    typedef _FT                                    FT;
    typedef Cartesian_tag                          Rep_tag;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
    typedef CGAL::Point_2<R,Rep_tag>             Point_2;
    typedef CGAL::Vector_2<R,Rep_tag>            Vector_2;
    typedef CGAL::Direction_2<R,Rep_tag>         Direction_2;
    typedef CGAL::Segment_2<R,Rep_tag>           Segment_2;
    typedef CGAL::Line_2<R,Rep_tag>              Line_2;
    typedef CGAL::Ray_2<R,Rep_tag>               Ray_2;
    typedef CGAL::Triangle_2<R,Rep_tag>          Triangle_2;
    typedef CGAL::Circle_2<R,Rep_tag>            Circle_2;
    typedef CGAL::Iso_rectangle_2<R,Rep_tag>     Iso_rectangle_2;
    typedef CGAL::Aff_transformation_2<R,Rep_tag> Aff_transformation_2;
    typedef CGAL::Data_accessor_2<R,Rep_tag>     Data_accessor_2;
    typedef CGAL::ConicCPA2<Point_2,Data_accessor_2> Conic_2;
#else
    typedef PointC2<R>                             Point_2;
    typedef VectorC2<R>                            Vector_2;
    typedef DirectionC2<R>                         Direction_2;
    typedef SegmentC2<R>                           Segment_2;
    typedef LineC2<R>                              Line_2;
    typedef RayC2<R>                               Ray_2;
    typedef TriangleC2<R>                          Triangle_2;
    typedef CircleC2<R>                            Circle_2;
    typedef Iso_rectangleC2<R>                     Iso_rectangle_2;
    typedef Aff_transformationC2<R>                Aff_transformation_2;
    typedef Data_accessorC2<R>                     Data_accessor_2;
    typedef ConicCPA2<Point_2,Data_accessor_2>     Conic_2;
#endif // CGAL_CFG_NO_ADVANCED_KERNEL
};

CGAL_END_NAMESPACE

#include <CGAL/Cartesian/Point_2.h>
#include <CGAL/Cartesian/Vector_2.h>
#include <CGAL/Cartesian/Direction_2.h>
#include <CGAL/Cartesian/Line_2.h>
#include <CGAL/Cartesian/Ray_2.h>
#include <CGAL/Cartesian/Segment_2.h>
#include <CGAL/Cartesian/Triangle_2.h>
#include <CGAL/Cartesian/Circle_2.h>
#include <CGAL/Cartesian/Iso_rectangle_2.h>
#include <CGAL/Cartesian/Aff_transformation_2.h>
#include <CGAL/Cartesian/Data_accessor_2.h>

#include <CGAL/Cartesian/global_operators_2.h>
#include <CGAL/Cartesian/predicates_on_points_2.h>
#include <CGAL/Cartesian/predicates_on_directions_2.h>
#include <CGAL/Cartesian/predicates_on_lines_2.h>
#include <CGAL/Cartesian/distance_predicates_2.h>
#include <CGAL/Cartesian/basic_constructions_2.h>

#include <CGAL/Cartesian/Point_2.C>
#include <CGAL/Cartesian/Vector_2.C>
#include <CGAL/Cartesian/Direction_2.C>
#include <CGAL/Cartesian/Line_2.C>
#include <CGAL/Cartesian/Ray_2.C>
#include <CGAL/Cartesian/Segment_2.C>
#include <CGAL/Cartesian/Triangle_2.C>
#include <CGAL/Cartesian/Circle_2.C>
#include <CGAL/Cartesian/Iso_rectangle_2.C>
#include <CGAL/Cartesian/Aff_transformation_2.C>

CGAL_BEGIN_NAMESPACE

// This class is a restricted 2D geometric kernel
// It is useful only if you do not need the 3D kernel
// If you need both, you should be using Cartesian<FT>

template< class _FT >
struct Cartesian_2 : public Cartesian_base_2< Cartesian_2<_FT>, _FT >
{
    // Number types and representation tag
    typedef _FT                                 RT;
    typedef _FT                                 FT;
    typedef Cartesian_tag                       Rep_tag;

    typedef Cartesian_2<_FT>                    Self;
    typedef Cartesian_base_2<Self,_FT>          Base;
    
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
    // The classes are inherited and because of partial specialization,
    // Cartesian_2<FT>::Point_2 is exactly CGAL::Point_2< Cartesian_2<FT> >
    // We still need to inherit explicitly, see Cartesian.h for explanation

    typedef Base::Point_2                       Point_2;
    typedef Base::Vector_2                      Vector_2;
    typedef Base::Direction_2                   Direction_2;
    typedef Base::Segment_2                     Segment_2;
    typedef Base::Line_2                        Line_2;
    typedef Base::Ray_2                         Ray_2;
    typedef Base::Triangle_2                    Triangle_2;
    typedef Base::Circle_2                      Circle_2;
    typedef Base::Iso_rectangle_2               Iso_rectangle_2;
    typedef Base::Aff_transformation_2          Aff_transformation_2;

#else
    // Now CGAL::Point_2<R> is only a wrapper around CGAL::PointC2<R>
    // It is necessary to redefine here the classes to ensure that
    // Cartesian_2<FT>::Point_2 is exactly CGAL::Point_2< Cartesian_2<FT> >

    // Cartesian_2<FT>::Base is needed so that CGAL::Point_2< Cartesian_2<FT> >
    // can inherit from Cartesian_2<FT>::Point_2_base

    typedef Base::Point_2                       Point_2_base;
    typedef Base::Vector_2                      Vector_2_base;
    typedef Base::Direction_2                   Direction_2_base;
    typedef Base::Segment_2                     Segment_2_base;
    typedef Base::Line_2                        Line_2_base;
    typedef Base::Ray_2                         Ray_2_base;
    typedef Base::Triangle_2                    Triangle_2_base;
    typedef Base::Circle_2                      Circle_2_base;
    typedef Base::Iso_rectangle_2               Iso_rectangle_2_base;
    typedef Base::Aff_transformation_2          Aff_transformation_2_base;

    // Note: necessary to qualify Point_2 by CGAL:: to disambiguate between
    // Point_2 in the current namespace (nested within CGAL)
    // CGAL::Point_2< Cartesian_2<FT> > (which is in the CGAL namespace)

    typedef CGAL::Point_2<Self>                 Point_2;
    typedef CGAL::Vector_2<Self>                Vector_2;
    typedef CGAL::Direction_2<Self>             Direction_2;
    typedef CGAL::Segment_2<Self>               Segment_2;
    typedef CGAL::Line_2<Self>                  Line_2;
    typedef CGAL::Ray_2<Self>                   Ray_2;
    typedef CGAL::Triangle_2<Self>              Triangle_2;
    typedef CGAL::Circle_2<Self>                Circle_2;
    typedef CGAL::Iso_rectangle_2<Self>         Iso_rectangle_2;
    typedef CGAL::Aff_transformation_2<Self>    Aff_transformation_2;

    typedef Data_accessorC2<Self>               Data_accessor_2;
    typedef ConicCPA2<Point_2,Data_accessor_2>  Conic_2;

    // TODO: cleanup
    static   FT make_FT(const RT & num, const RT& denom) { return num/denom;}
    static   FT make_FT(const RT & num)                  { return num;}
    static   RT FT_numerator(const FT &r)                { return r;}
    static   RT FT_denominator(const FT &)               { return RT(1);}

#endif // CGAL_CFG_NO_ADVANCED_KERNEL
};

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_2_H
