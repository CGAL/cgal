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
// file          : include/CGAL/Cartesian.h
// source        : include/CGAL/Cartesian.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr
//
// coordinator   : INRIA Sophia-Antipolis
//
// ============================================================================


#ifndef CGAL_CARTESIAN_H
#define CGAL_CARTESIAN_H

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif // CGAL_BASIC_H

#ifndef CGAL_CARTESIAN_CLASSES_H
#include <CGAL/cartesian_classes.h>
#endif // CGAL_CARTESIAN_CLASSES_H

#ifndef CGAL_CARTESIAN_2_H
#include <CGAL/Cartesian_2.h>
#endif // CGAL_CARTESIAN_2_H
#ifndef CGAL_CARTESIAN_3_H
#include <CGAL/Cartesian_3.h>
#endif // CGAL_CARTESIAN_3_H
#ifndef CGAL_CARTESIAN_DYNAMIC_D_H
// #include <CGAL/Cartesian_dynamic_d.h>
#endif // CGAL_CARTESIAN_DYNAMIC_D_H

#define CGAL_REP_CLASS_DEFINED
#define CGAL_CARTESIAN_CLASS_DEFINED

CGAL_BEGIN_NAMESPACE

template< class R, class _FT >
struct Cartesian_base :
    public Cartesian_base_2<R,_FT>
    , public Cartesian_base_3<R,_FT>
    // , public Cartesian_base_dynamic_d<R,_FT>
{
    // Number types and representation tag (to avoid ambiguity in
    // inheritance tree)
    typedef _FT                                 RT;
    typedef _FT                                 FT;
    typedef Cartesian_tag                       Rep_tag;

    // All the classes are inherited, but because we inherit from a
    // template parameter, we need to explicitly write the inheritance
    // (see mail from Michael Hoffmann of July 28th 1999 in cgal-develop)

    typedef Cartesian_base_2<R,_FT>             Base_2;
    typedef Cartesian_base_3<R,_FT>             Base_3;
    // typedef Cartesian_base_dynamic_d<R,_FT>     Base_d;

    typedef typename Base_2::Point_2                     Point_2;
    typedef typename Base_2::Vector_2                    Vector_2;
    typedef typename Base_2::Direction_2                 Direction_2;
    typedef typename Base_2::Segment_2                   Segment_2;
    typedef typename Base_2::Line_2                      Line_2;
    typedef typename Base_2::Ray_2                       Ray_2;
    typedef typename Base_2::Triangle_2                  Triangle_2;
    typedef typename Base_2::Circle_2                    Circle_2;
    typedef typename Base_2::Iso_rectangle_2             Iso_rectangle_2;
    typedef typename Base_2::Aff_transformation_2        Aff_transformation_2;

    typedef typename Base_2::Data_accessor_2             Data_accessor_2;
    typedef typename Base_2::Conic_2                     Conic_2;

    typedef typename Base_3::Point_3                     Point_3;
    typedef typename Base_3::Vector_3                    Vector_3;
    typedef typename Base_3::Direction_3                 Direction_3;
    typedef typename Base_3::Line_3                      Line_3;
    typedef typename Base_3::Plane_3                     Plane_3;
    typedef typename Base_3::Ray_3                       Ray_3;
    typedef typename Base_3::Segment_3                   Segment_3;
    typedef typename Base_3::Triangle_3                  Triangle_3;
    typedef typename Base_3::Tetrahedron_3               Tetrahedron_3;
    typedef typename Base_3::Aff_transformation_3        Aff_transformation_3;

    // typedef Base_d::Point_d                     Point_d;
};

CGAL_END_NAMESPACE

#include <CGAL/Kernel/Construction_objects.h>
#include <CGAL/Kernel/Predicate_objects.h>

CGAL_BEGIN_NAMESPACE

template< class _FT >
struct Cartesian : public Cartesian_base< Cartesian<_FT>, _FT >
{
    // Number types and representation tag (to avoid ambiguity)
    typedef _FT                                 RT;
    typedef _FT                                 FT;
    typedef Cartesian_tag                       Rep_tag;

    typedef Cartesian<FT>                       Self;
    typedef Cartesian_base<Self,FT>             Base;

#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
    // The other classes are inherited and because of partial specialization,
    // Cartesian<FT>::Point_2 is exactly CGAL::Point_2< Cartesian<FT> >
    // As above, we still need to write down the inheritance explicitly

    typedef typename Base::Point_2                       Point_2;
    typedef typename Base::Vector_2                      Vector_2;
    typedef typename Base::Direction_2                   Direction_2;
    typedef typename Base::Segment_2                     Segment_2;
    typedef typename Base::Line_2                        Line_2;
    typedef typename Base::Ray_2                         Ray_2;
    typedef typename Base::Triangle_2                    Triangle_2;
    typedef typename Base::Circle_2                      Circle_2;
    typedef typename Base::Iso_rectangle_2               Iso_rectangle_2;
    typedef typename Base::Aff_transformation_2          Aff_transformation_2;

    typedef typename Base::Data_accessor_2               Data_accessor_2;
    typedef typename Base::Conic_2                       Conic_2;

    typedef typename Base::Point_3                       Point_3;
    typedef typename Base::Vector_3                      Vector_3;
    typedef typename Base::Direction_3                   Direction_3;
    typedef typename Base::Line_3                        Line_3;
    typedef typename Base::Plane_3                       Plane_3;
    typedef typename Base::Ray_3                         Ray_3;
    typedef typename Base::Segment_3                     Segment_3;
    typedef typename Base::Triangle_3                    Triangle_3;
    typedef typename Base::Tetrahedron_3                 Tetrahedron_3;
    typedef typename Base::Aff_transformation_3          Aff_transformation_3;

    // typedef typename Base::Point_d                       Point_d;

#else
    // Now CGAL::Point_2<R> is only a wrapper around CGAL::PointC2<R>
    // It is necessary to redefine here the classes to ensure that
    // Cartesian<FT>::Point_2 is exactly CGAL::Point_2< Cartesian<FT> >

    // Cartesian<FT>::Base is needed so that CGAL::Point_2< Cartesian<FT> >
    // can inherit from Cartesian<FT>::Point_2_base

    typedef typename Base::Point_2                    Point_2_base;
    typedef typename Base::Vector_2                   Vector_2_base;
    typedef typename Base::Direction_2                Direction_2_base;
    typedef typename Base::Segment_2                  Segment_2_base;
    typedef typename Base::Line_2                     Line_2_base;
    typedef typename Base::Ray_2                      Ray_2_base;
    typedef typename Base::Triangle_2                 Triangle_2_base;
    typedef typename Base::Circle_2                   Circle_2_base;
    typedef typename Base::Iso_rectangle_2            Iso_rectangle_2_base;
    typedef typename Base::Aff_transformation_2       Aff_transformation_2_base;

    typedef typename Base::Point_3                    Point_3_base;
    typedef typename Base::Vector_3                   Vector_3_base;
    typedef typename Base::Direction_3                Direction_3_base;
    typedef typename Base::Line_3                     Line_3_base;
    typedef typename Base::Plane_3                    Plane_3_base;
    typedef typename Base::Ray_3                      Ray_3_base;
    typedef typename Base::Segment_3                  Segment_3_base;
    typedef typename Base::Triangle_3                 Triangle_3_base;
    typedef typename Base::Tetrahedron_3              Tetrahedron_3_base;
    typedef typename Base::Aff_transformation_3       Aff_transformation_3_base;
  
    // Note: necessary to qualify Point_2 by CGAL:: to disambiguate between
    // Point_2 in the current namespace (nested within CGAL) and
    // CGAL::Point_2< Cartesian<FT> > (which is in the CGAL namespace)

    typedef CGAL::Point_2<Self>              Point_2;
    typedef CGAL::Vector_2<Self>             Vector_2;
    typedef CGAL::Direction_2<Self>          Direction_2;
    typedef CGAL::Segment_2<Self>            Segment_2;
    typedef CGAL::Line_2<Self>               Line_2;
    typedef CGAL::Ray_2<Self>                Ray_2;
    typedef CGAL::Triangle_2<Self>           Triangle_2;
    typedef CGAL::Circle_2<Self>             Circle_2;
    typedef CGAL::Iso_rectangle_2<Self>      Iso_rectangle_2;
    typedef CGAL::Aff_transformation_2<Self> Aff_transformation_2;

    typedef Data_accessorC2<Self>            Data_accessor_2;
    typedef ConicCPA2<Point_2,Data_accessor_2>  Conic_2;

    typedef CGAL::Point_3<Self>              Point_3;
    typedef CGAL::Vector_3<Self>             Vector_3;
    typedef CGAL::Direction_3<Self>          Direction_3;
    typedef CGAL::Line_3<Self>               Line_3;
    typedef CGAL::Plane_3<Self>              Plane_3;
    typedef CGAL::Ray_3<Self>                Ray_3;
    typedef CGAL::Segment_3<Self>            Segment_3;
    typedef CGAL::Triangle_3<Self>           Triangle_3;
    typedef CGAL::Tetrahedron_3<Self>        Tetrahedron_3;
    typedef CGAL::Aff_transformation_3<Self> Aff_transformation_3;

    // typedef CGAL::Point_d<Self>              Point_d;

#endif // CGAL_CFG_NO_ADVANCED_KERNEL

    // TODO: cleanup
    static   FT make_FT(const RT & num, const RT& denom) { return num/denom;}
    static   FT make_FT(const RT & num)                  { return num;}
    static   RT FT_numerator(const FT &r)                { return r;}
    static   RT FT_denominator(const FT &)               { return RT(1);}

    typedef Kernel_construction_objects<Self> Construction_objects;

    typedef Kernel_predicate_objects<Self>    Predicate_objects;

    // These macros are defined in Kernel_function_objects
    CGAL_UNPACK_KERNEL_CONSTRUCTION_OBJECTS(typename Construction_objects)
    CGAL_UNPACK_KERNEL_PREDICATE_OBJECTS(typename Predicate_objects)
};

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_H
