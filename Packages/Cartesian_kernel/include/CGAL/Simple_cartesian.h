// ======================================================================
//
// Copyright (c) 2000,2001 The CGAL Consortium
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
// file          : include/CGAL/Simple_cartesian.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sylvain Pion
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_SIMPLE_CARTESIAN_H
#define CGAL_SIMPLE_CARTESIAN_H

// Same Rep class as ref-counted Cartesian ?
#define CGAL_REP_CLASS_DEFINED
#define CGAL_CARTESIAN_CLASS_DEFINED

#include <CGAL/basic.h>
#include <CGAL/cartesian_classes.h>
#include <CGAL/Cartesian_dynamic_d.h>

#include <CGAL/Twotuple.h>
#include <CGAL/Threetuple.h>
#include <CGAL/Simple_Handle_for.h>

#ifdef CGAL_CFG_NO_ADVANCED_KERNEL
  // Because we cannot use Michael's scheme, we need the wrapper classes
  // We include them (they are common to Cartesian and Homogeneous)
#include <CGAL/user_classes.h>
#endif // CGAL_CFG_NO_ADVANCED_KERNEL

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

#include <CGAL/predicate_objects_on_points_2.h>
#include <CGAL/Kernel/function_objects.h>

#include <CGAL/Cartesian/Point_3.h>
#include <CGAL/Cartesian/Vector_3.h>
#include <CGAL/Cartesian/Direction_3.h>
#include <CGAL/Cartesian/Line_3.h>
#include <CGAL/Cartesian/Plane_3.h>
#include <CGAL/Cartesian/Ray_3.h>
#include <CGAL/Cartesian/Segment_3.h>
#include <CGAL/Cartesian/Triangle_3.h>
#include <CGAL/Cartesian/Tetrahedron_3.h>
#include <CGAL/Cartesian/Iso_cuboid_3.h>
#include <CGAL/Cartesian/Sphere_3.h>
#include <CGAL/Cartesian/Aff_transformation_3.h>

#include <CGAL/Cartesian/global_operators_3.h>
#include <CGAL/Cartesian/predicates_on_points_3.h>
#include <CGAL/Cartesian/predicates_on_planes_3.h>
#include <CGAL/Cartesian/distance_predicates_3.h>

#include <CGAL/Cartesian/basic_constructions_3.h>

CGAL_BEGIN_NAMESPACE

template< class R, class FT_ >
struct Simple_cartesian_base :
     public Cartesian_base_dynamic_d<R,FT_>
{
    typedef FT_                                         RT;
    typedef FT_                                         FT;
    typedef Cartesian_tag                               Rep_tag;
    typedef Cartesian_tag                               Kernel_tag;

    // Maybe Cd should be merged like C2 and C3 ?
    typedef Cartesian_base_dynamic_d<R,FT_>             Kernel_base_d;

    typedef CGAL::Object                                Object_2;
    typedef CGAL::Object                                Object_3;

#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
    typedef CGAL::Point_2<R,Rep_tag>                    Point_2;
    typedef CGAL::Vector_2<R,Rep_tag>                   Vector_2;
    typedef CGAL::Direction_2<R,Rep_tag>                Direction_2;
    typedef CGAL::Segment_2<R,Rep_tag>                  Segment_2;
    typedef CGAL::Line_2<R,Rep_tag>                     Line_2;
    typedef CGAL::Ray_2<R,Rep_tag>                      Ray_2;
    typedef CGAL::Triangle_2<R,Rep_tag>                 Triangle_2;
    typedef CGAL::Circle_2<R,Rep_tag>                   Circle_2;
    typedef CGAL::Iso_rectangle_2<R,Rep_tag>            Iso_rectangle_2;
    typedef CGAL::Aff_transformation_2<R,Rep_tag>       Aff_transformation_2;
    typedef CGAL::Data_accessor_2<R,Rep_tag>            Data_accessor_2;
    typedef CGAL::ConicCPA2<Point_2,Data_accessor_2>    Conic_2;

    typedef CGAL::Point_3<R,Rep_tag>                    Point_3;
    typedef CGAL::Vector_3<R,Rep_tag>                   Vector_3;
    typedef CGAL::Direction_3<R,Rep_tag>                Direction_3;
    typedef CGAL::Line_3<R,Rep_tag>                     Line_3;
    typedef CGAL::Plane_3<R,Rep_tag>                    Plane_3;
    typedef CGAL::Ray_3<R,Rep_tag>                      Ray_3;
    typedef CGAL::Segment_3<R,Rep_tag>                  Segment_3;
    typedef CGAL::Triangle_3<R,Rep_tag>                 Triangle_3;
    typedef CGAL::Tetrahedron_3<R,Rep_tag>              Tetrahedron_3;
    typedef CGAL::Iso_cuboid_3<R,Rep_tag>               Iso_cuboid_3;
    typedef CGAL::Sphere_3<R,Rep_tag>                   Sphere_3;
    typedef CGAL::Aff_transformation_3<R,Rep_tag>       Aff_transformation_3;
#else
    typedef PointC2<R>                                  Point_2;
    typedef VectorC2<R>                                 Vector_2;
    typedef DirectionC2<R>                              Direction_2;
    typedef SegmentC2<R>                                Segment_2;
    typedef LineC2<R>                                   Line_2;
    typedef RayC2<R>                                    Ray_2;
    typedef TriangleC2<R>                               Triangle_2;
    typedef CircleC2<R>                                 Circle_2;
    typedef Iso_rectangleC2<R>                          Iso_rectangle_2;
    typedef Aff_transformationC2<R>                     Aff_transformation_2;
    typedef Data_accessorC2<R>                          Data_accessor_2;
    typedef ConicCPA2<Point_2,Data_accessor_2>          Conic_2;

    typedef PointC3<R>                                  Point_3;
    typedef VectorC3<R>                                 Vector_3;
    typedef DirectionC3<R>                              Direction_3;
    typedef LineC3<R>                                   Line_3;
    typedef PlaneC3<R>                                  Plane_3;
    typedef RayC3<R>                                    Ray_3;
    typedef SegmentC3<R>                                Segment_3;
    typedef TriangleC3<R>                               Triangle_3;
    typedef TetrahedronC3<R>                            Tetrahedron_3;
    typedef Iso_cuboidC3<R>                             Iso_cuboid_3;
    typedef SphereC3<R>                                 Sphere_3;
    typedef Aff_transformationC3<R>                     Aff_transformation_3;
#endif // CGAL_CFG_NO_ADVANCED_KERNEL
  
    typedef typename Kernel_base_d::Point_d             Point_d;

    // The typedefs that allow to specify the handle of each type.

    typedef CGAL::Simple_Handle_for<CGAL::Simple_Twotuple<FT> >	
	                                                Point_handle_2;
    typedef CGAL::Simple_Handle_for<CGAL::Simple_Twotuple<FT> >
                                                        Vector_handle_2;
    typedef CGAL::Simple_Handle_for<CGAL::Simple_Twotuple<FT> >
                                                 	Direction_handle_2;
    typedef CGAL::Simple_Handle_for<CGAL::Simple_Threetuple<FT> >
                                                        Line_handle_2;
    typedef CGAL::Simple_Handle_for<CGAL::Simple_Twotuple<Point_2> >
                                                       	Ray_handle_2;
    typedef CGAL::Simple_Handle_for<CGAL::Simple_Twotuple<Point_2> >
                                                      	Segment_handle_2;
    typedef CGAL::Simple_Handle_for<CGAL::Simple_Circle_repC2<R> >
                                                      	Circle_handle_2;
    typedef CGAL::Simple_Handle_for<CGAL::Simple_Threetuple<Point_2> >
                                                       	Triangle_handle_2;
    typedef CGAL::Simple_Handle_for<CGAL::Simple_Twotuple<Point_2> >
                                                 	Iso_rectangle_handle_2;
    typedef CGAL::Handle			Aff_transformation_handle_2;


    typedef CGAL::Simple_Handle_for<CGAL::Simple_Threetuple<FT> >
                                                    	Point_handle_3;
    typedef CGAL::Simple_Handle_for<CGAL::Simple_Threetuple<FT> >
                                                   	Vector_handle_3;
    typedef CGAL::Simple_Handle_for<CGAL::Simple_Threetuple<FT> >
                                                  	Direction_handle_3;
    typedef CGAL::Simple_Handle_for<CGAL::Simple_Line_repC3<R> >
                                                     	Line_handle_3;
    typedef CGAL::Simple_Handle_for<CGAL::Simple_Fourtuple<FT> >
                                                   	Plane_handle_3;
    typedef CGAL::Simple_Handle_for<CGAL::Simple_Twotuple<Point_3> >
                                                  	Ray_handle_3;
    typedef CGAL::Simple_Handle_for<CGAL::Simple_Twotuple<Point_3> >
                                                   	Segment_handle_3;
    typedef CGAL::Simple_Handle_for<CGAL::Simple_Sphere_repC3<R> >
                                                    	Sphere_handle_3;
    typedef CGAL::Simple_Handle_for<CGAL::Simple_Threetuple<Point_3> >
                                                  	Triangle_handle_3;
    typedef CGAL::Simple_Handle_for<CGAL::Simple_Fourtuple<Point_3> >
                                                   	Tetrahedron_handle_3;
    typedef CGAL::Simple_Handle_for<CGAL::Simple_Twotuple<Point_3> >
                                                	Iso_cuboid_handle_3;
    typedef CGAL::Handle	Aff_transformation_handle_3;
};


template< class FT_ >
struct Simple_cartesian
  : public Simple_cartesian_base< Simple_cartesian<FT_>, FT_ >
{
    // Number types and representation tag (to avoid ambiguity)
    typedef FT_                                           RT;
    typedef FT_                                           FT;
    typedef Cartesian_tag                                 Rep_tag;
    typedef Cartesian_tag                                 Kernel_tag;

    typedef Simple_cartesian<FT>                       Self;
    typedef Simple_cartesian_base<Self,FT>             Kernel_base;

    typedef typename Kernel_base::Object_2                Object_2;
    typedef typename Kernel_base::Object_3                Object_3;          
  

#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
    // The other classes are inherited and because of partial specialization,
    // Cartesian<FT>::Point_2 is exactly CGAL::Point_2< Cartesian<FT> >

    // All the classes are inherited, but because we inherit from a
    // template parameter, we need to explicitly write the inheritance
    // (see mail from Michael Hoffmann of July 28th 1999 in cgal-develop)

    typedef typename Kernel_base::Point_2                 Point_2;
    typedef typename Kernel_base::Vector_2                Vector_2;
    typedef typename Kernel_base::Direction_2             Direction_2;
    typedef typename Kernel_base::Segment_2               Segment_2;
    typedef typename Kernel_base::Line_2                  Line_2;
    typedef typename Kernel_base::Ray_2                   Ray_2;
    typedef typename Kernel_base::Triangle_2              Triangle_2;
    typedef typename Kernel_base::Circle_2                Circle_2;
    typedef typename Kernel_base::Iso_rectangle_2         Iso_rectangle_2;
    typedef typename Kernel_base::Aff_transformation_2    Aff_transformation_2;

    typedef typename Kernel_base::Data_accessor_2         Data_accessor_2;
    typedef typename Kernel_base::Conic_2                 Conic_2;

    typedef typename Kernel_base::Point_3                 Point_3;
    typedef typename Kernel_base::Vector_3                Vector_3;
    typedef typename Kernel_base::Direction_3             Direction_3;
    typedef typename Kernel_base::Line_3                  Line_3;
    typedef typename Kernel_base::Plane_3                 Plane_3;
    typedef typename Kernel_base::Ray_3                   Ray_3;
    typedef typename Kernel_base::Segment_3               Segment_3;
    typedef typename Kernel_base::Triangle_3              Triangle_3;
    typedef typename Kernel_base::Tetrahedron_3           Tetrahedron_3;
    typedef typename Kernel_base::Sphere_3                Sphere_3;
    typedef typename Kernel_base::Iso_cuboid_3            Iso_cuboid_3;
    typedef typename Kernel_base::Aff_transformation_3    Aff_transformation_3;

    typedef typename Kernel_base::Point_d                 Point_d;

#else
    // Now CGAL::Point_2<R> is only a wrapper around CGAL::PointC2<R>
    // It is necessary to redefine here the classes to ensure that
    // Cartesian<FT>::Point_2 is exactly CGAL::Point_2< Cartesian<FT> >

    // Cartesian<FT>::Base is needed so that CGAL::Point_2< Cartesian<FT> >
    // can inherit from Cartesian<FT>::Point_2_base

    typedef typename Kernel_base::Point_2                 Point_2_base;
    typedef typename Kernel_base::Vector_2                Vector_2_base;
    typedef typename Kernel_base::Direction_2             Direction_2_base;
    typedef typename Kernel_base::Segment_2               Segment_2_base;
    typedef typename Kernel_base::Line_2                  Line_2_base;
    typedef typename Kernel_base::Ray_2                   Ray_2_base;
    typedef typename Kernel_base::Triangle_2              Triangle_2_base;
    typedef typename Kernel_base::Circle_2                Circle_2_base;
    typedef typename Kernel_base::Iso_rectangle_2         Iso_rectangle_2_base;
    typedef typename Kernel_base::Aff_transformation_2  
                                                     Aff_transformation_2_base;

    typedef typename Kernel_base::Point_3                 Point_3_base;
    typedef typename Kernel_base::Vector_3                Vector_3_base;
    typedef typename Kernel_base::Direction_3             Direction_3_base;
    typedef typename Kernel_base::Line_3                  Line_3_base;
    typedef typename Kernel_base::Plane_3                 Plane_3_base;
    typedef typename Kernel_base::Ray_3                   Ray_3_base;
    typedef typename Kernel_base::Segment_3               Segment_3_base;
    typedef typename Kernel_base::Triangle_3              Triangle_3_base;
    typedef typename Kernel_base::Tetrahedron_3           Tetrahedron_3_base;
    typedef typename Kernel_base::Sphere_3                Sphere_3_base;
    typedef typename Kernel_base::Iso_cuboid_3            Iso_cuboid_3_base;
    typedef typename Kernel_base::Aff_transformation_3    
                                                  Aff_transformation_3_base;
  
    typedef typename Kernel_base::Point_d                 Point_d_base;

    // Note: necessary to qualify Point_2 by CGAL:: to disambiguate between
    // Point_2 in the current namespace (nested within CGAL) and
    // CGAL::Point_2< Cartesian<FT> > (which is in the CGAL namespace)

    typedef CGAL::Point_2<Self>                           Point_2;
    typedef CGAL::Vector_2<Self>                          Vector_2;
    typedef CGAL::Direction_2<Self>                       Direction_2;
    typedef CGAL::Line_2<Self>                            Line_2;
    typedef CGAL::Ray_2<Self>                             Ray_2;
    typedef CGAL::Segment_2<Self>                         Segment_2;
    typedef CGAL::Triangle_2<Self>                        Triangle_2;
    typedef CGAL::Circle_2<Self>                          Circle_2;
    typedef CGAL::Iso_rectangle_2<Self>                   Iso_rectangle_2;
    typedef CGAL::Aff_transformation_2<Self>              Aff_transformation_2;

    typedef Data_accessorC2<Self>                         Data_accessor_2;
    typedef ConicCPA2<Point_2,Data_accessor_2>            Conic_2;

    typedef CGAL::Point_3<Self>                           Point_3;
    typedef CGAL::Vector_3<Self>                          Vector_3;
    typedef CGAL::Direction_3<Self>                       Direction_3;
    typedef CGAL::Line_3<Self>                            Line_3;
    typedef CGAL::Plane_3<Self>                           Plane_3;
    typedef CGAL::Ray_3<Self>                             Ray_3;
    typedef CGAL::Segment_3<Self>                         Segment_3;
    typedef CGAL::Triangle_3<Self>                        Triangle_3;
    typedef CGAL::Tetrahedron_3<Self>                     Tetrahedron_3;
    typedef CGAL::Sphere_3<Self>                          Sphere_3;
    typedef CGAL::Iso_cuboid_3<Self>                      Iso_cuboid_3;
    typedef CGAL::Aff_transformation_3<Self>              Aff_transformation_3;

    typedef CGAL::Point_d<Self>                           Point_d;

#endif // CGAL_CFG_NO_ADVANCED_KERNEL

    typedef typename Kernel_base::Point_handle_2	Point_handle_2;
    typedef typename Kernel_base::Vector_handle_2	Vector_handle_2;
    typedef typename Kernel_base::Direction_handle_2	Direction_handle_2;
    typedef typename Kernel_base::Line_handle_2		Line_handle_2;
    typedef typename Kernel_base::Ray_handle_2		Ray_handle_2;
    typedef typename Kernel_base::Segment_handle_2	Segment_handle_2;
    typedef typename Kernel_base::Circle_handle_2	Circle_handle_2;
    typedef typename Kernel_base::Triangle_handle_2	Triangle_handle_2;
    typedef typename Kernel_base::Iso_rectangle_handle_2
	                                                Iso_rectangle_handle_2;
    typedef typename Kernel_base::Aff_transformation_handle_2
	                                           Aff_transformation_handle_2;

    typedef typename Kernel_base::Point_handle_3	Point_handle_3;
    typedef typename Kernel_base::Vector_handle_3	Vector_handle_3;
    typedef typename Kernel_base::Direction_handle_3	Direction_handle_3;
    typedef typename Kernel_base::Line_handle_3		Line_handle_3;
    typedef typename Kernel_base::Plane_handle_3	Plane_handle_3;
    typedef typename Kernel_base::Ray_handle_3		Ray_handle_3;
    typedef typename Kernel_base::Segment_handle_3	Segment_handle_3;
    typedef typename Kernel_base::Sphere_handle_3	Sphere_handle_3;
    typedef typename Kernel_base::Triangle_handle_3	Triangle_handle_3;
    typedef typename Kernel_base::Tetrahedron_handle_3	Tetrahedron_handle_3;
    typedef typename Kernel_base::Iso_cuboid_handle_3
	                                                Iso_cuboid_handle_3;
    typedef typename Kernel_base::Aff_transformation_handle_3
	                                           Aff_transformation_handle_3;

    // TODO: cleanup
    static   FT make_FT(const RT & num, const RT& denom) { return num/denom;}
    static   FT make_FT(const RT & num)                  { return num;}
    static   RT FT_numerator(const FT &r)                { return r;}
    static   RT FT_denominator(const FT &)               { return RT(1);}

#include <CGAL/Kernel_traits_common.h>

};

CGAL_END_NAMESPACE

#include <CGAL/iterator_traits_pointer_specs_for_simple_cartesian_kernel.h>

#endif // CGAL_SIMPLE_CARTESIAN_H
