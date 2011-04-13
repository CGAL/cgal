#ifndef CGAL_CARTESIAN_DYNAMIC_D_H
#define CGAL_CARTESIAN_DYNAMIC_D_H

#include <CGAL/basic.h>
#include <CGAL/cartesian_classes.h>

#ifdef CGAL_CFG_NO_ADVANCED_KERNEL
  // Because we cannot use Michael's scheme, we need the wrapper classes
  // We include them (they are common to Cartesian and Homogeneous)
  #include <CGAL/user_classes.h>
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
    typedef CGAL::Linear_algebra_d<FT,Rep_tag>  LA;
    typedef CGAL::Point_d<R,Rep_tag>            Point_d;
    typedef CGAL::Vector_d<R,Rep_tag>           Vector_d;
    typedef CGAL::Direction_d<R,Rep_tag>        Direction_d;
    typedef CGAL::Plane_d<R,Rep_tag>            Plane_d;
    typedef CGAL::Line_d<R,Rep_tag>             Line_d;
    typedef CGAL::Ray_d<R,Rep_tag>              Ray_d;
    typedef CGAL::Segment_d<R,Rep_tag>          Segment_d;
    typedef CGAL::Triangle_d<R,Rep_tag>         Triangle_d;
    typedef CGAL::Tetrahedron_d<R,Rep_tag>      Tetrahedron_d;
    typedef CGAL::Simplex_d<R,Rep_tag>          Simplex_d;
    typedef CGAL::Aff_transformation_d<R,Rep_tag> Aff_transformation_d;
#else
    typedef Linear_algebraCd<FT>                LA;
    typedef PointCd<R>                          Point_d;
    typedef VectorCd<R>                         Vector_d;
    typedef DirectionCd<R>                      Direction_d;
    typedef PlaneCd<R>                          Plane_d;
    typedef LineCd<R>                           Line_d;
    typedef RayCd<R>                            Ray_d;
    typedef SegmentCd<R>                        Segment_d;
    typedef TriangleCd<R>                       Triangle_d;
    typedef TetrahedronCd<R>                    Tetrahedron_d;
    typedef SimplexCd<R>                        Simplex_d;
    typedef Aff_transformationCd<R>             Aff_transformation_d;
#endif // CGAL_CFG_NO_ADVANCED_KERNEL
};

CGAL_END_NAMESPACE

#include <CGAL/Cartesian/Linear_algebra_d.h>
#include <CGAL/Cartesian/Point_d.h>
#include <CGAL/Cartesian/Vector_d.h>
#include <CGAL/Cartesian/Direction_d.h>
#include <CGAL/Cartesian/Plane_d.h>
#include <CGAL/Cartesian/Line_d.h>
#include <CGAL/Cartesian/Ray_d.h>
#include <CGAL/Cartesian/Segment_d.h>
#include <CGAL/Cartesian/Triangle_d.h>
#include <CGAL/Cartesian/Tetrahedron_d.h>
#include <CGAL/Cartesian/Simplex_d.h>
#include <CGAL/Cartesian/Aff_transformation_d.h>

#include <CGAL/Cartesian/global_operators_d.h>
#include <CGAL/Cartesian/constructions_on_planes_d.h>
#include <CGAL/Cartesian/distance_computations_d.h>
#include <CGAL/Cartesian/predicates_on_points_d.h>
#include <CGAL/Cartesian/predicates_on_directions_d.h>
#include <CGAL/Cartesian/predicates_on_planes_d.h>
#include <CGAL/Cartesian/distance_predicates_d.h>

#include <CGAL/Cartesian/Linear_algebra_d.C>
#include <CGAL/Cartesian/Point_d.C>
#include <CGAL/Cartesian/Vector_d.C>
#include <CGAL/Cartesian/Direction_d.C>
#include <CGAL/Cartesian/Plane_d.C>
#include <CGAL/Cartesian/Line_d.C>
#include <CGAL/Cartesian/Ray_d.C>
#include <CGAL/Cartesian/Segment_d.C>
#include <CGAL/Cartesian/Triangle_d.C>
#include <CGAL/Cartesian/Tetrahedron_d.C>
#include <CGAL/Cartesian/Simplex_d.C>
#include <CGAL/Cartesian/Aff_transformation_d.C>

CGAL_BEGIN_NAMESPACE

// This class is a restricted dD geometric kernel
// It is useful only if you do not need the 3D kernel
// If you need both, you should be using Cartesian<FT>

template< class _FT >
struct Cartesian_dynamic_d
  : public Cartesian_base_dynamic_d< Cartesian_dynamic_d<_FT>, _FT >
{
    // Number types and representation tag
    typedef _FT                                 RT;
    typedef _FT                                 FT;
    typedef Cartesian_tag                       Rep_tag;

    typedef Cartesian_dynamic_d<_FT>            Self;
    typedef Cartesian_base_dynamic_d<Self,_FT>  Kernel_base;

#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
    // The lasses are inherited and because of partial specialization,
    // Cartesian_d<FT>::Point_d is exactly CGAL::Point_d< Cartesian_d<FT> >
    // We still need to inherit explicitly, see Cartesian.h for explanation

    typedef typename Kernel_base::Point_d       Point_d;
    typedef typename Kernel_base::Vector_d      Vector_d;
    typedef typename Kernel_base::Direction_d   Direction_d;
    typedef typename Kernel_base::Plane_d       Plane_d;
    typedef typename Kernel_base::Line_d        Line_d;
    typedef typename Kernel_base::Ray_d         Ray_d;
    typedef typename Kernel_base::Segment_d     Segment_d;
    typedef typename Kernel_base::Triangle_d    Triangle_d;
    typedef typename Kernel_base::Tetrahedron_d Tetrahedron_d;
    typedef typename Kernel_base::Simplex_d     Simplex_d;
    typedef typename Kernel_base::Aff_transformation_d
	                                        Aff_transformation_d;

#else
    // Now CGAL::Point_d<R> is only a wrapper around CGAL::PointCd<R>
    // It is necessary to redefine here the classes to ensure that
    // Cartesian_d<FT>::Point_d is exactly CGAL::Point_d< Cartesian_d<FT> >

    // Cartesian_d<FT>::Base is needed so that CGAL::Point_d< Cartesian_d<FT> >
    // can inherit from Cartesian_d<FT>::Kernel_base::Point_d

    typedef typename Kernel_base::Point_d       Point_d_base;
    typedef typename Kernel_base::Vector_d      Vector_d_base;
    typedef typename Kernel_base::Direction_d   Direction_d_base;
    typedef typename Kernel_base::Plane_d       Plane_d_base;
    typedef typename Kernel_base::Line_d        Line_d_base;
    typedef typename Kernel_base::Ray_d         Ray_d_base;
    typedef typename Kernel_base::Segment_d     Segment_d_base;
    typedef typename Kernel_base::Triangle_d    Triangle_d_base;
    typedef typename Kernel_base::Tetrahedron_d Tetrahedron_d_base;
    typedef typename Kernel_base::Simplex_d     Simplex_d_base;
    typedef typename Kernel_base::Aff_transformation_d
	                                        Aff_transformation_d_base;

    // Note: necessary to qualify Point_d by ::CGAL:: to disambiguate between
    // Point_d in the current namespace (nested within CGAL) and
    // CGAL::Point_d< Cartesian_d<FT> > (which is in the CGAL namespace)

    typedef CGAL::Point_d<Self>                 Point_d;
    typedef CGAL::Vector_d<Self>                Vector_d;
    typedef CGAL::Direction_d<Self>             Direction_d;
    typedef CGAL::Plane_d<Self>                 Plane_d;
    typedef CGAL::Line_d<Self>                  Line_d;
    typedef CGAL::Ray_d<Self>                   Ray_d;
    typedef CGAL::Segment_d<Self>               Segment_d;
    typedef CGAL::Triangle_d<Self>              Triangle_d;
    typedef CGAL::Tetrahedron_d<Self>           Tetrahedron_d;
    typedef CGAL::Simplex_d<Self>               Simplex_d;
    typedef CGAL::Aff_transformation_d<Self>    Aff_transformation_d;

    // TODO: cleanup
    static   FT make_FT(const RT & num, const RT& denom) { return num/denom;}
    static   FT make_FT(const RT & num)                  { return num;}
    static   RT FT_numerator(const FT &r)                { return r;}
    static   RT FT_denominator(const FT &)               { return RT(1);}

#endif // CGAL_CFG_NO_ADVANCED_KERNEL
};

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_DYNAMIC_D_H
