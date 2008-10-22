// Copyright (c) 2007-2008 Max-Planck-Institute Saarbruecken (Germany), 
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//                 Pavel Emeliyanenko <asm@mpi-sb.mpg.de>

#ifndef CGAL_CURVED_KERNEL_VIA_ANALYSIS_2L_FUNCTORS_H
#define CGAL_CURVED_KERNEL_VIA_ANALYSIS_2L_FUNCTORS_H

/*!\file include/CGAL/Curved_kernel_via_analysis_2l_functors.h
 * \brief Contains functors for \c Curved_kernel_via_analysis_2l
 *  
 * Functors for generic points and arcs on surfaces, lifted from 2D.
 */

#include <CGAL/config.h>

#include <CGAL/Curved_kernel_via_analysis_2/Curved_kernel_via_analysis_2_functors.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

#ifndef CERR
//#define CKvA_DEBUG_PRINT_CERR
#ifdef CKvA_DEBUG_PRINT_CERR
#define CERR(x) std::cerr << x
#else
#define CERR(x) static_cast<void>(0)
#endif
#endif

namespace Curved_kernel_via_analysis_2l_Functors {

#define CGAL_CKvA_2l_GRAB_BASE_FUNCTOR_TYPES \
    typedef typename Base::Curve_analysis_2 Curve_analysis_2; \
    typedef typename Base::Point_2 Point_2; \
    typedef typename Base::Arc_2 Arc_2; \
    typedef typename Base::X_coordinate_1 X_coordinate_1; \


template < class CurvedKernelViaAnalysis_2l >
class Construct_point_2l : public Curved_kernel_via_analysis_2_Functors::
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2l > {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;

    //! the base type
    typedef 
    Curved_kernel_via_analysis_2_Functors::
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2l >
    Base;
    
    CGAL_CKvA_2l_GRAB_BASE_FUNCTOR_TYPES;
    
    //! type of projected point
    typedef typename Point_2::Projected_point_2 Projected_point_2;
    
    //! type of surface
    typedef typename Point_2::Surface_3 Surface_3;

    //! the result type
    typedef Point_2 result_type;
    
    //! standard constructor
    Construct_point_2l(Curved_kernel_via_analysis_2l *kernel) :
        Base(kernel) {
    }
    
    Point_2 operator()(const Projected_point_2& xy, 
                       const Surface_3& surface,
                       int sheet) const {
        CGAL_precondition(sheet >= 0);
        Point_2 pt(xy, surface, sheet);
        return pt;
    }
};

//!\brief Functor to construct point on an arc
//! \c x on curve \c c with arc number \c arcno
//!
//! implies no boundary conditions in x/y
template < class CurvedKernelViaAnalysis_2l >
class Construct_point_on_arc_2 : public Curved_kernel_via_analysis_2_Functors::
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2l > {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;

    //! the base type
    typedef 
    Curved_kernel_via_analysis_2_Functors::
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2l >
    Base;
    
    CGAL_CKvA_2l_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef Point_2 result_type;
    
    //! standard constructor
    Construct_point_on_arc_2(Curved_kernel_via_analysis_2l *kernel) :
        Base(kernel) {
    }
    
    //! constructs points at x 
    Point_2 operator()(
            const typename Point_2::X_coordinate_1& x, 
            const typename Point_2::Curve_analysis_2& c, int arcno,
            const Arc_2& arc) {

        typename Curved_kernel_via_analysis_2l::Construct_projected_point_2 
            construct_projected_point = 
            Curved_kernel_via_analysis_2l::instance().
            construct_projected_point_2_object();
        typename Point_2::Projected_point_2 p_pt = 
            construct_projected_point(x, c, arcno);
        int sheet = arc.sheet();
        bool same_projection = false;
        if (arc.is_finite(CGAL::ARR_MIN_END)) {
            if (p_pt.compare_xy(
                        arc.curve_end(CGAL::ARR_MIN_END).projected_point()
                ) == CGAL::EQUAL) {
                same_projection = true;
                sheet = arc.sheet(CGAL::ARR_MIN_END);
            }
        }
        if (!same_projection && arc.is_finite(CGAL::ARR_MAX_END)) {
            if (p_pt.compare_xy(
                        arc.curve_end(CGAL::ARR_MAX_END).projected_point()
                ) == CGAL::EQUAL) {
                sheet = arc.sheet(CGAL::ARR_MAX_END);
            }
        }
        typename Curved_kernel_via_analysis_2l::Construct_point_2 
            construct_point_2 = 
            Curved_kernel_via_analysis_2l::instance().
            construct_point_2_object();
        
        Point_2 pt = construct_point_2(p_pt, arc.surface(), sheet);
        return pt;
    }
};


template < class CurvedKernelViaAnalysis_2l >
class Construct_arc_2l : public Curved_kernel_via_analysis_2_Functors::
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2l > {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;

    //! the base type
    typedef 
    Curved_kernel_via_analysis_2_Functors::
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2l >
    Base;
    
    CGAL_CKvA_2l_GRAB_BASE_FUNCTOR_TYPES;
    
    //! type of point on surface
    typedef typename Curved_kernel_via_analysis_2l::Point_2 Surface_point_2l;
    
    //! type of arc on surface
    typedef typename Curved_kernel_via_analysis_2l::Arc_2 Surface_arc_2l;
    
    //! type of projected point
    typedef typename Surface_point_2l::Projected_point_2 Projected_point_2;
    
    //! type of projceteda rc
    typedef typename Surface_arc_2l::Projected_arc_2 Projected_arc_2;

    //! type of surface
    typedef typename Surface_point_2l::Surface_3 Surface_3;

    //! the result type
    typedef Point_2 result_type;
     
    //! standard constructor
    Construct_arc_2l(Curved_kernel_via_analysis_2l *kernel) :
        Base(kernel) {
    }
    
    //!\name Constructors based on bounded planar arcs
    //!@{
    
    /*!\brief
     * Standard constructor for an bounded arc on xy-monotone part
     * of the surface.
     * It represents the arc on \c surface covertical to \c arc which
     * lies on \c sheet of the xy-monotone subsurface.
     *
     * \pre arc.curve_end(MIN) = p.projected_point()
     * \pre arc.curve_end(MAX) = q.projected_point()
     */
    Surface_arc_2l operator()(const Projected_arc_2& arc, 
                              const Surface_point_2l& p,
                              const Surface_point_2l& q,
                              const Surface_3& surface,
                              int sheet, int sheet_p, int sheet_q) {
        Surface_arc_2l surface_arc(arc, p, q, surface, sheet, sheet_p, sheet_q);
        return surface_arc;
    }
    
    /*!\brief
     * Standard constructor for an ray on a xy-monotone part
     * of the surface with a z-asymptotic behavior on one side.
     * It represents the arc on \c surface covertical to \c arc which
     * lies on \c sheet of the xy-monotone subsurface.
     *
     * \pre arc.curve_end(MIN) = p.projected_point() || 
     *      arc.curve_end(MAX) = p.projected_point()
     */
    Surface_arc_2l operator()(const Projected_arc_2& arc, 
                              const Surface_point_2l& p,
                              CGAL::Arr_curve_end z_inf_end_other,
                              const Surface_3& surface,
                              int sheet, int sheet_p) {
        Surface_arc_2l surface_arc(arc, p, z_inf_end_other, 
                                   surface, sheet, sheet_p);
        return surface_arc;
    }
    
    /*!\brief
     * Standard constructor for an unbounded arc on xy-monotone part
     * of the surface with z-asympotic behavior at both ends.
     * It represents the arc on \c surface covertical to \c arc which
     * lies on \c sheet of the xy-monotone subsurface.
     */
    Surface_arc_2l operator()(const Projected_arc_2& arc, 
                              CGAL::Arr_curve_end z_inf_end_p,
                              CGAL::Arr_curve_end z_inf_end_q,
                              const Surface_3& surface,
                              int sheet) {
        Surface_arc_2l surface_arc(arc, z_inf_end_p, z_inf_end_q, 
                                   surface, sheet);
        return surface_arc;
    }
    
    //!}
    
    //!\name Constructors based on planar rays
    //!@{
    
    /*!\brief
     * Standard constructor for a ray on xy-monotone part
     * of the surface.
     * It represents the arc on \c surface covertical to \c arc which
     * lies on \c sheet of the xy-monotone subsurface.
     *
     * \pre arc.curve_end(MIN) = p || arc.curve_end(MAX) == p
     */
    Surface_arc_2l operator()(const Projected_arc_2& arc, 
                              const Surface_point_2l& p,
                              const Surface_3& surface,
                              int sheet, int sheet_p) {
        Surface_arc_2l surface_arc(arc, p, surface, sheet, sheet_p);
        return surface_arc;    
    }
    
    /*!\brief
     * Standard constructor for a branch on a xy-monotone part
     * of the surface with a z-vertical asymptotic behaviour at the projected
     * finite end.
     * It represents the arc on \c surface covertical to \c arc which
     * lies on \c sheet of the xy-monotone subsurface
     */
    Surface_arc_2l operator()(const Projected_arc_2& arc, 
                              CGAL::Arr_curve_end z_inf_end,
                              const Surface_3& surface,
                              int sheet) {
        Surface_arc_2l surface_arc(arc, z_inf_end, surface, sheet);
        return surface_arc;
    }
    
    //!@}
    
    //!\name Constructors based on planar branches
    //!@{
    
    /*!\brief
     * Standard constructor for a branch on xy-monotone part
     * of the surface.
     * It represents the arc on \c surface covertical to \c arc which
     * lies on \c sheet of the xy-monotone subsurface.
     */
    Surface_arc_2l operator()(const Projected_arc_2& arc, 
                              const Surface_3& surface,
                              int sheet) {
        Surface_arc_2l surface_arc(arc, surface, sheet);
        return surface_arc;
    }
    
    //!@}

    //!\name Constructors for vertical arcs 
    //!@{
    
    // Remark for vertical arcs:
    // Their base is not an arc, i.e., the projection of the arc is a 
    // single point, so we have to deal with it throughout 
    // the whole class, i.e., 
    // This point must be used, as the default constructed Base is useless.
    
    //! represents a bounded vertical arc
    Surface_arc_2l operator()(const Surface_point_2l& p,
                              const Surface_point_2l& q,
                              const Surface_3& surface) {
        Surface_arc_2l surface_arc(p, q, surface);
        return surface_arc;
    }
    
    //! represents a vertical ray
    Surface_arc_2l operator()(const Surface_point_2l& p,
                              CGAL::Arr_curve_end z_inf_end,
                              const Surface_3& surface) {
        Surface_arc_2l surface_arc(p, z_inf_end, surface);
        return surface_arc;
    }

    //! represents a vertical branch
    Surface_arc_2l operator()(const Projected_point_2& p,
                              const Surface_3& surface) {
        Surface_arc_2l surface_arc(p, surface);
        return surface_arc;
    }
    
    //!@}
};


//!\brief Tests whether a point lies on a supporting curve
template < class CurvedKernelViaAnalysis_2l >
class Compare_xyz_3 : public Curved_kernel_via_analysis_2_Functors::
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2l > {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;

    //! the base type
    typedef 
    Curved_kernel_via_analysis_2_Functors::
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2l >
    Base;
    
    CGAL_CKvA_2l_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef CGAL::Comparison_result result_type;
    
    //! standard constructor
    Compare_xyz_3(Curved_kernel_via_analysis_2l *kernel) :
        Base(kernel) {
    }
    
    /*!
     *\brief Compares two points xyz-lexicographically
     * \param p1 The first point
     * \param p2 The second point
     * \return CGAL::SMALLER if p1 is lexicographically smaller than p2,
     *         CGAL::EQUAL if p1 is equal to p2, and
     *         CGAL::LARGER if p1 is lexicographically larger than p2,
     */
    result_type operator()(const Point_2& p1, const Point_2& p2,
                           bool equal_xy = false) const {

        CERR("\nqkvacompare_xyz: p1: " << p1 << "\n p2: " << p2 << "\n");

        result_type res = CGAL::EQUAL;

        if (!p1.is_identical(p2)) {
            res = (equal_xy ? 
                   CGAL::EQUAL : 
                   Curved_kernel_via_analysis_2l::instance().
                   compare_xy_2_object()(p1, p2)
            );
            if (res == CGAL::EQUAL) {
                
                if (p1.is_z_at_infinity()) {
                    if (p2.is_z_at_infinity()) {
                        CGAL::Arr_curve_end z1 = p1.z_infinity();
                        CGAL::Arr_curve_end z2 = p2.z_infinity();
                        if (z1 == z2) {
                            return EQUAL;
                        } else {
                            return (z1 == CGAL::ARR_MIN_END ?
                                    CGAL::SMALLER : CGAL::LARGER);
                        }
                    } else {
                        CGAL::Arr_curve_end z1 = p1.z_infinity();
                        return (z1 == CGAL::ARR_MIN_END ?
                                CGAL::SMALLER : CGAL::LARGER);
                    }
                } else {
                    if (p2.is_z_at_infinity()) {
                        CGAL::Arr_curve_end z2 = p2.z_infinity();
                        return (z2 == CGAL::ARR_MIN_END ?
                                CGAL::LARGER : CGAL::SMALLER);
                    }
                }
                
                // else both finite!
                
                typedef typename Curved_kernel_via_analysis_2l::Surface_pair_3
                    Surface_pair_3;
                typedef typename Surface_pair_3::Surface_3 Surface_3;
                const Surface_3& s1 = p1.surface();
                const Surface_3& s2 = p2.surface();
                int sheet1 = p1.sheet();
                int sheet2 = p2.sheet();
                
                // if same surface compare arc numbers, 
                if (s1.is_identical(s2)) {
                    return CGAL::sign(sheet1 - sheet2);
                }

                // otherwise use "z-stacke of surface pair 
                // (point location in 2d-map gives z-stacke 
                // 2d map consists of silhouettes-curves plus cut curve
                // slice will store sequence and also whether 
                // surfaces are vertical surface_info???
                
                Surface_pair_3 pair = 
                    Surface_pair_3::surface_pair_cache()(
                            std::make_pair(s1, s2)
                    );
                
                typedef typename 
                    Surface_pair_3::Restricted_cad_3::Z_stack Z_stack;
                
                Z_stack z_stack = pair.z_stack_for(
                        p1.projected_point()
                );
                
                int level1 = z_stack.z_level_of_sheet(s1, sheet1);
                int level2 = z_stack.z_level_of_sheet(s2, sheet2);
                
                return CGAL::sign(level1 - level2);
            }
        }
        
        CERR("result: " << res << "\n");
        
        return res;
    }
};


//!\brief Tests whether a point lies on a surface
template < class CurvedKernelViaAnalysis_2l >
class Is_on_3: public Curved_kernel_via_analysis_2_Functors::
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2l > {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;

    //! the base type
    typedef 
    Curved_kernel_via_analysis_2_Functors::
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2l >
    Base;
    
    CGAL_CKvA_2l_GRAB_BASE_FUNCTOR_TYPES;
    
    //! type of surface
    typedef typename Point_2::Surface_3 Surface_3;

    //! the result type
    typedef bool result_type;
    
    //! standard constructor
    Is_on_3(Curved_kernel_via_analysis_2l *kernel) :
        Base(kernel) {
    }
    
    /*!
     * Checks whether \c p lies on \c surface
     * \param p The point to test
     * \param surface The surface
     * \return (true) if the \c p lies on \c surface
     */
    result_type operator()(const Point_2& p, const Surface_3& surface) const {
        
        CERR("\nqkvais_on_3: p: " << p << "\n surface: " << surface << "\n");
        
        CGAL_precondition_msg(p.is_finite(), 
                              "Is_on_3: Point at inf not supported");
        CGAL_precondition_msg(!p.is_z_at_infinity(),
                              "Is_on_3: Point with |z|=oo not supported");
        
        result_type res = false;
        
        if (surface == p.surface()) {
            res = true;
        } else {
            
            typedef typename Curved_kernel_via_analysis_2l::Surface_pair_3
                Surface_pair_3;

            Surface_pair_3 pair = 
                Surface_pair_3::surface_pair_cache()(
                        std::make_pair(p.surface(), surface)
                );
            
            typedef typename 
                Surface_pair_3::Restricted_cad_3::Z_stack Z_stack;
            
            Z_stack z_stack = pair.z_stack_for(
                    p.projected_point()
            );

            int level = z_stack.z_level_of_sheet(p.surface(), p.sheet());
            
            int level_s = z_stack.level_of_surface_in_z_cell(surface, level);

            res = (level_s != -1);
        }
        
        CERR("result: " << res << "\n");
        
        return res;
    }
};

//!\brief Tests whether a point lies on a supporting curve
template < class CurvedKernelViaAnalysis_2l >
class Is_on_2: public Curved_kernel_via_analysis_2_Functors::
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2l > {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;

    //! the base type
    typedef 
    Curved_kernel_via_analysis_2_Functors::
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2l >
    Base;
    
    CGAL_CKvA_2l_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef bool result_type;
    
    //! standard constructor
    Is_on_2(Curved_kernel_via_analysis_2l *kernel) :
        Base(kernel) {
    }
    
    /*!
     * Checks whether \c p lies on \c c 
     * \param p The point to test
     * \param c The curve
     * \return (true) if the \c p lies on \c c
     */
    result_type operator()(const Point_2& p, const Curve_analysis_2& c) const {

        CERR("\nqkvais_on_2: p: " << p << "\n curve: " << c << "\n");
        // FUTURE TODO implement Is_on_2 with Curve_2 == Surface_3 (eriC)
        CGAL_error_msg("Is_on_2 not implemented for Surfaces");
        
        result_type res = false;
        
        CERR("result: " << res << "\n");
        
        return res;
    }
};

#undef CGAL_CKvA_2l_GRAB_BASE_FUNCTOR_TYPES 

} //  Curved_kernel_via_analysis_2l_Functors

} // namespace CGALi
 
CGAL_END_NAMESPACE

#endif // CGAL_CURVED_KERNEL_VIA_ANALYSIS_2L_FUNCTORS_H
// EOF
