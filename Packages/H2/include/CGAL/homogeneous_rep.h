// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       : 
// release_date  : 
//
// file          : include/CGAL/homogeneous_rep.h
// package       : Kernel_basic
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CGAL_HOMOGENEOUS_REP_H
#define CGAL_HOMOGENEOUS_REP_H

#include <CGAL/Handle_for.h>
#include <CGAL/Handle.h>
#include <CGAL/Twotuple.h>
#include <CGAL/Threetuple.h>
#include <CGAL/Fourtuple.h>

#include <CGAL/RepH3.h>

#include <CGAL/Quotient.h>

#include <CGAL/homogeneous_classes.h>
#include <CGAL/representation_tags.h>
#include <CGAL/predicate_objects_on_points_2.h>
#include <CGAL/Kernel/function_objects.h>

CGAL_BEGIN_NAMESPACE

template <class R_, class RT_, class FT_>
class Homogeneous_base
{
  public:
    typedef RT_                                     RT;
    typedef FT_                                     FT;
    typedef CGAL::Object                            Object_2;
    // we have: template <class R> CGAL::Point_2 : public R::Point_2_base
    typedef CGAL::Point_2< R_ >                     Point_2;
    typedef CGAL::Vector_2< R_ >                    Vector_2;
    typedef CGAL::Direction_2< R_ >                 Direction_2;
    typedef CGAL::Segment_2< R_ >                   Segment_2;
    typedef CGAL::Line_2< R_ >                      Line_2;
    typedef CGAL::Ray_2< R_ >                       Ray_2;
    typedef CGAL::Circle_2< R_ >                    Circle_2;
    typedef CGAL::Triangle_2< R_ >                  Triangle_2;
    typedef CGAL::Iso_rectangle_2< R_ >             Iso_rectangle_2;
    typedef CGAL::Aff_transformation_2< R_ >        Aff_transformation_2;
    typedef CGAL::Object                            Object_3;
    typedef CGAL::Point_3< R_ >                     Point_3;
    typedef CGAL::Vector_3< R_ >                    Vector_3;
    typedef CGAL::Direction_3< R_ >                 Direction_3;
    typedef CGAL::Segment_3< R_ >                   Segment_3;
    typedef CGAL::Plane_3< R_ >                     Plane_3;
    typedef CGAL::Line_3< R_ >                      Line_3;
    typedef CGAL::Ray_3< R_ >                       Ray_3;
    typedef CGAL::Triangle_3< R_ >                  Triangle_3;
    typedef CGAL::Tetrahedron_3< R_ >               Tetrahedron_3;
    typedef CGAL::Iso_cuboid_3< R_ >                Iso_cuboid_3;
    typedef CGAL::Sphere_3< R_ >                    Sphere_3;
    typedef CGAL::Aff_transformation_3< R_ >        Aff_transformation_3;
    // we have: template <class R> CGAL::Point_d : public R::Point_d_base
    typedef CGAL::Point_d< R_ >                     Point_d;

    typedef CGAL::Handle_for< Threetuple<RT> >      Point_handle_2;
    typedef CGAL::Handle_for< Threetuple<RT> >      Vector_handle_2;
    typedef CGAL::Handle_for< Threetuple<RT> >      Direction_handle_2;
    typedef CGAL::Handle_for< Ray_repH2<R_> >       Ray_handle_2;
    typedef CGAL::Handle_for< Triangle_repH2< R_> > Triangle_handle_2;
    typedef CGAL::Handle_for< Circle_repH2<R_> >    Circle_handle_2;
    typedef CGAL::Handle_for< Twotuple<PointH2<R_> > > Iso_rectangle_handle_2;
    typedef CGAL::Handle_for< Threetuple<RT> >      Line_handle_2;
    typedef CGAL::Handle_for< Segment_repH2<R_> >   Segment_handle_2;
    typedef CGAL::Handle_for< Aff_transformation_rep_baseH2<R_>,
                       No_op_allocator< Aff_transformation_rep_baseH2<R_> > >
			                           Aff_transformation_handle_2;

    typedef CGAL::Handle_for< RepH3<RT> >          Point_handle_3;
    typedef CGAL::Handle_for< RepH3<RT> >          Vector_handle_3;
    typedef CGAL::Handle_for< RepH3<RT> >          Direction_handle_3;
    typedef CGAL::Handle_for< Fourtuple<RT> >      Plane_handle_3;
    typedef CGAL::Handle_for< Ray_repH3<R_> >      Ray_handle_3;
    typedef CGAL::Handle_for< Line_repH3<R_> >     Line_handle_3;
    typedef CGAL::Handle_for< Sphere_repH3<R_> >   Sphere_handle_3;
    typedef CGAL::Handle_for< Tetrahedron_repH3<R_> >   Tetrahedron_handle_3;
    typedef CGAL::Handle_for< Segment_repH3<R_> >  Segment_handle_3;
    typedef CGAL::Handle_for< Threetuple< PointH3<R_> > >  Triangle_handle_3;
    typedef CGAL::Handle_for< Twotuple<PointH3<R_> > > Iso_cuboid_handle_3;
    typedef CGAL::Handle                          Aff_transformation_handle_3;

};

template <class RT_, class FT_ = Quotient<RT_> >
class Homogeneous : public Homogeneous_base< Homogeneous<RT_,FT_>, RT_, FT_ >
{
  public:
    typedef Homogeneous< RT_, FT_ >                 R;
    typedef RT_                                     RT;
    typedef FT_                                     FT;
    typedef Homogeneous_tag                         Rep_tag;
    typedef PointH2< R >                            Point_2_base;
    typedef VectorH2< R >                           Vector_2_base;
    typedef DirectionH2< R >                        Direction_2_base;
    typedef SegmentH2< R >                          Segment_2_base;
    typedef LineH2< R >                             Line_2_base;
    typedef RayH2< R >                              Ray_2_base;
    typedef CircleH2< R >                           Circle_2_base;
    typedef TriangleH2< R >                         Triangle_2_base;
    typedef Iso_rectangleH2< R >                    Iso_rectangle_2_base;
    typedef Aff_transformationH2< R >               Aff_transformation_2_base;
    typedef Homogeneous_base< Homogeneous<RT_,FT_>, RT_, FT_ >   KernelBase;

    typedef typename KernelBase::Point_2               Point_2;
    typedef typename KernelBase::Vector_2              Vector_2;
    typedef typename KernelBase::Direction_2           Direction_2;
    typedef typename KernelBase::Line_2                Line_2;
    typedef typename KernelBase::Segment_2             Segment_2;
    typedef typename KernelBase::Ray_2                 Ray_2;
    typedef typename KernelBase::Circle_2              Circle_2;
    typedef typename KernelBase::Triangle_2            Triangle_2;
    typedef typename KernelBase::Iso_rectangle_2       Iso_rectangle_2;
    typedef typename KernelBase::Aff_transformation_2  Aff_transformation_2;
    typedef typename KernelBase::Object_2              Object_2;
    typedef typename KernelBase::Object_3              Object_3;

    typedef PointH3< R >                             Point_3_base;
    typedef VectorH3< R >                            Vector_3_base;
    typedef DirectionH3< R >                         Direction_3_base;
    typedef SegmentH3< R >                           Segment_3_base;
    typedef PlaneH3< R >                             Plane_3_base;
    typedef LineH3< R >                              Line_3_base;
    typedef RayH3< R >                               Ray_3_base;
    typedef TriangleH3< R >                          Triangle_3_base;
    typedef TetrahedronH3< R >                       Tetrahedron_3_base;
    typedef Iso_cuboidH3< R >                        Iso_cuboid_3_base;
    typedef SphereH3< R >                            Sphere_3_base;
    typedef Aff_transformationH3< R >                Aff_transformation_3_base;
    typedef typename KernelBase::Point_3             Point_3;
    typedef typename KernelBase::Vector_3            Vector_3;
    typedef typename KernelBase::Direction_3         Direction_3;
    typedef typename KernelBase::Plane_3             Plane_3;
    typedef typename KernelBase::Line_3              Line_3;
    typedef typename KernelBase::Segment_3           Segment_3;
    typedef typename KernelBase::Ray_3               Ray_3;
    typedef typename KernelBase::Sphere_3            Sphere_3;
    typedef typename KernelBase::Triangle_3          Triangle_3;
    typedef typename KernelBase::Tetrahedron_3       Tetrahedron_3;
    typedef typename KernelBase::Iso_cuboid_3        Iso_cuboid_3;
    typedef typename KernelBase::Aff_transformation_3
                                                     Aff_transformation_3;

    typedef PointHd< FT, RT>                         Point_d_base;
    typedef typename KernelBase::Point_d             Point_d;

    typedef Data_accessorH2<R>                       Data_accessor_2;
    typedef ConicHPA2<Point_2,Data_accessor_2>       Conic_2;

    typedef typename KernelBase::Point_handle_2        Point_handle_2;
    typedef typename KernelBase::Vector_handle_2       Vector_handle_2;
    typedef typename KernelBase::Direction_handle_2    Direction_handle_2;
    typedef typename KernelBase::Ray_handle_2          Ray_handle_2;
    typedef typename KernelBase::Triangle_handle_2     Triangle_handle_2;
    typedef typename KernelBase::Circle_handle_2       Circle_handle_2;
    typedef typename KernelBase::Iso_rectangle_handle_2 Iso_rectangle_handle_2;
    typedef typename KernelBase::Line_handle_2         Line_handle_2;
    typedef typename KernelBase::Segment_handle_2      Segment_handle_2;
    typedef typename KernelBase::Aff_transformation_handle_2
			                           Aff_transformation_handle_2;

    typedef typename KernelBase::Point_handle_3        Point_handle_3;
    typedef typename KernelBase::Vector_handle_3       Vector_handle_3;
    typedef typename KernelBase::Direction_handle_3    Direction_handle_3;
    typedef typename KernelBase::Plane_handle_3        Plane_handle_3;
    typedef typename KernelBase::Ray_handle_3          Ray_handle_3;
    typedef typename KernelBase::Line_handle_3         Line_handle_3;
    typedef typename KernelBase::Sphere_handle_3       Sphere_handle_3;
    typedef typename KernelBase::Tetrahedron_handle_3  Tetrahedron_handle_3;
    typedef typename KernelBase::Segment_handle_3      Segment_handle_3;
    typedef typename KernelBase::Triangle_handle_3     Triangle_handle_3;
    typedef typename KernelBase::Iso_cuboid_handle_3   Iso_cuboid_handle_3;
    typedef typename KernelBase::Aff_transformation_handle_3
			                           Aff_transformation_handle_3;

    static
    FT
    make_FT(const RT & num, const RT& denom)
    { return FT(num, denom); }

    static
    FT
    make_FT(const RT & num)
    { return FT(num); }

    static
    RT
    FT_numerator(const FT &r)
    { return r.numerator(); }

    static
    RT
    FT_denominator(const FT &r)
    { return r.denominator(); }

#include <CGAL/Kernel_traits_common.h>

};

CGAL_END_NAMESPACE

#endif // CGAL_HOMOGENEOUS_REP_H
