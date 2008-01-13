// TODO: Add licence
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
//
// ============================================================================

#ifndef CGAL_SURFACE_ARC_2L_H
#define CGAL_SURFACE_ARC_2L_H

/*! \file Surface_arc_2l.h
 *  \brief defines class \c Surface_arc_2l
 *  
 *  arcs on surfaces, lifted from 2D.
 */

#include <CGAL/basic.h>

#include <CGAL/Curved_kernel_via_analysis_2l/Surface_point_2l.h>
#include <CGAL/Curved_kernel_via_analysis_2/Arc_2.h>
#include <CGAL/Curved_kernel_via_analysis_2.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

// pre-declaration
template < class CurvedKernelViaAnalysis_2l, class SurfacePair_3, class Rep_ >
class Surface_arc_2l;

// TODO documentation
template < class CurvedKernelViaAnalysis_2l, class SurfacePair_3 >
class Surface_arc_2l_rep : 
        public Arc_2_base_rep< CurvedKernelViaAnalysis_2l > {

public:

    //! this type's template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;

    //! this instance's second template parameter
    typedef SurfacePair_3 Surface_pair_3;

    //! the instance itself
    typedef 
    Surface_arc_2l_rep< Curved_kernel_via_analysis_2l, Surface_pair_3 >
    Self;
    
    //! the base type
    typedef Arc_2_base_rep< Curved_kernel_via_analysis_2l > Base;
    

    //! type of surface
    typedef typename Surface_pair_3::Surface_3 Surface_3;

    //!\name Constructors
    //!@{

    Surface_arc_2l_rep() :
        _m_sheet(-1), _m_sheet_min(-1), _m_sheet_max(-1), 
        _m_is_z_vertical(false) {
    }

    //!@}
    
    
protected:
    //! supporting surface
    mutable Surface_3 _m_surface;

    //! sheet number of arc
    mutable int _m_sheet;

    //! sheet number of arc at x-min point
    mutable int _m_sheet_min;
    
    //! sheet number of arcs at x-max point
    mutable int _m_sheet_max;
    
    //! indicates whether arc is vertical
    mutable int _m_is_z_vertical;
    
    // befriending the handle
    friend class 
    Surface_arc_2l< Curved_kernel_via_analysis_2l, Surface_pair_3, Self >;

};


//! represents an xy-monotone arc on a surface
template < 
  class CurvedKernelViaAnalysis_2l, 
  class SurfacePair_3,
  class Rep_ = 
    CGALi::Surface_arc_2l_rep< CurvedKernelViaAnalysis_2l, SurfacePair_3 >
>
class Surface_arc_2l :
    public CGALi::Arc_2_base< 
        CurvedKernelViaAnalysis_2l, 
        //Surface_arc_2l< CurvedKernelViaAnalysis_2l, SurfacePair_3 >, 
        typename CurvedKernelViaAnalysis_2l::Arc_2,
        Rep_ 
    > 
{

public:

    //!\name Public types
    //!@{
    
    //! this type's first template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;
    
    //! this type's second template parameter
    typedef SurfacePair_3 Surface_pair_3;

    //! this type's third template parameter
    typedef Rep_ Rep;

    //! the class itself
    typedef 
    Surface_arc_2l< Curved_kernel_via_analysis_2l, Surface_pair_3, Rep >
    Self;
    
    //! the base class
    typedef CGALi::Arc_2_base< Curved_kernel_via_analysis_2l, 
                               typename Curved_kernel_via_analysis_2l::Arc_2, 
                               Rep > 
    Base;
    
    //! type of planar point
    typedef typename Curved_kernel_via_analysis_2l::Point_2::Base 
    Projected_point_2;
    
    //! type of planar arc
    typedef Base Projected_arc_2;

    //! type of surface
    typedef typename Surface_pair_3::Surface_3 Surface_3;
    
    //! type of planar point
    typedef typename Curved_kernel_via_analysis_2l::Point_2 Surface_point_2l;

    //!@}
    
    //!\name Constructors
    //!@{

    /*!\brief
     * Default constructor
     */
    Surface_arc_2l() : 
        Base() {   
    }

protected:
    
    /*!\brief
     * constructs an arc from a given represenation
     */
    Surface_arc_2l(Rep rep) : 
        Base(rep) { 
    }
    
    /*!\brief
     * Standard constructor for an bounded arc on xy-monotone part
     * of the surface.
     * It represents the arc on \c surface covertical to \c arc which
     * lies on \c sheet of the xy-monotone subsurface.
     *
     * \pre arc.curve_end(MIN) = p
     * \pre arc.curve_end(MAX) = q
     */
    Surface_arc_2l(Curved_kernel_via_analysis_2l *kernel,
                   const Projected_arc_2& arc, 
                   const Surface_point_2l& p,
                   const Surface_point_2l& q,
                   const Surface_3& surface,
                   int sheet, int sheet_p, int sheet_q) :
        Base(arc) 
        // TODO rebind, i.e., replace p and q in "arc"
    {
        _set_ckva(kernel);

        // TODO remove location tests for contraction!
        CGAL_precondition(
                arc.location(CGAL::ARR_MIN_END) == CGAL::ARR_INTERIOR
        );
        CGAL_precondition(
                p.compare_xy(arc.curve_end(CGAL::ARR_MIN_END)) ==
                CGAL::EQUAL
        );
        CGAL_precondition(
                arc.location(CGAL::ARR_MAX_END) == CGAL::ARR_INTERIOR
        );
        CGAL_precondition(
                q.compare_xy(arc.curve_end(CGAL::ARR_MAX_END)) ==
                CGAL::EQUAL
        );
        
        //this->ptr()->_projected_segment = seg;
        this->ptr()->_m_surface = surface;
        CGAL_precondition(sheet >= 0);
        // TODO add precond CGAL_precondition(sheet < #sheets over arc);
        this->ptr()->_m_sheet = sheet;

        // TODO add precond CGAL_precondition(sheet < #sheets over min);
        // TODO add precond CGAL_precondition(sheet < #sheets over max);
        CGAL_precondition_code(
                // TODO add sanity checks for sheet_at_min/max wrt to sheet
                // use "adjacency information"
        );
        this->ptr()->_m_sheet_min = sheet_p;
        this->ptr()->_m_sheet_max = sheet_q;
    }

    /*!\brief
     * Standard constructor for a ray on xy-monotone part
     * of the surface.
     * It represents the arc on \c surface covertical to \c arc which
     * lies on \c sheet of the xy-monotone subsurface.
     *
     * \pre arc.curve_end(MIN) = p || arc.curve_end(MAX) == p
     */
    Surface_arc_2l(Curved_kernel_via_analysis_2l *kernel,
                   const Projected_arc_2& arc, 
                   const Surface_point_2l& p,
                   const Surface_3& surface,
                   int sheet, int sheet_p) :
        // TODO rebind, i.e., replace p
        Base(arc) 
    {
        _set_ckva(kernel);
        
        bool min_finite = 
            (arc.curve_end(CGAL::ARR_MIN_END) == CGAL::ARR_INTERIOR);
        CGAL_precondition_code(
                bool max_finite = 
                (arc.curve_end(CGAL::ARR_MAX_END) == CGAL::ARR_INTERIOR);
        );
        CGAL_precondition(min_finite || max_finite && 
                          !(!min_finite && !max_finite));
        
        CGAL_precondition_code(
                if (min_finite) {
                    CGAL_precondition(
                            p.compare_xy(arc.curve_end(CGAL::ARR_MIN_END)) ==
                            CGAL::EQUAL
                    );
                } else {
                    CGAL_precondition(
                            p.compare_xy(arc.curve_end(CGAL::ARR_MAX_END)) ==
                            CGAL::EQUAL
                    );
                }
        );
        
        //this->ptr()->_projected_segment = seg;
        this->ptr()->_m_surface = surface;
        CGAL_precondition(sheet >= 0);
        // TODO add precond CGAL_precondition(sheet < #sheets over arc);
        this->ptr()->_m_sheet = sheet;

        // TODO add precond CGAL_precondition(sheet < #sheets over min/ax);
        CGAL_precondition_code(
                // TODO add sanity checks for sheet_at_min/max wrt to sheet
                // use "adjacency information"
        );
        if (min_finite) {
            this->ptr()->_m_sheet_min = sheet_p;
            this->ptr()->_m_sheet_max = sheet;
        } else {
            this->ptr()->_m_sheet_min = sheet;
            this->ptr()->_m_sheet_max = sheet_p;
        }
    }
    

    /*!\brief
     * Standard constructor for a branch on xy-monotone part
     * of the surface.
     * It represents the arc on \c surface covertical to \c arc which
     * lies on \c sheet of the xy-monotone subsurface.
     *
     * \pre arc.curve_end(MIN) = p || arc.curve_end(MAX) == p
     */
    Surface_arc_2l(Curved_kernel_via_analysis_2l *kernel,
                   const Projected_arc_2& arc, 
                   const Surface_3& surface,
                   int sheet) :
        // TODO rebind?
        Base(arc) 
    {
        _set_ckva(kernel);

        bool min_finite = 
            (arc.curve_end(CGAL::ARR_MIN_END) == CGAL::ARR_INTERIOR);
        bool max_finite = 
            (arc.curve_end(CGAL::ARR_MAX_END) == CGAL::ARR_INTERIOR);
        CGAL_precondition(!min_finite && !max_finite);
        
        //this->ptr()->_projected_segment = seg;
        this->ptr()->_m_surface = surface;
        CGAL_precondition(sheet >= 0);
        // TODO add precond CGAL_precondition(sheet < #sheets over arc);
        this->ptr()->_m_sheet = sheet;
        this->ptr()->_m_sheet_min = sheet;
        this->ptr()->_m_sheet_max = sheet;
    }
    
    // constructors for vertical arcs
    
    //! represents a bounded vertical arc
    Surface_arc_2l(Curved_kernel_via_analysis_2l *kernel,
                   const Surface_point_2l& p,
                   const Surface_point_2l& q,
                   const Surface_3& surface) :
        Base() 
        // TODO base would be degenerate which is not allowed, so we have
        // to deal with it throughout the whole class, i.e., 
        // rep requires to store a single point p.base() 
        // whenever Base is default constructed. 
        // This point must be used, as the default constructed
        // Base is useless.
    {
        _set_ckva(kernel);
        
        CGAL_precondition(p.compare_xy(q) == CGAL::EQUAL);
        // TODO check that surface has a vertical line through p and q
        this->ptr()->_m_is_z_vertical = true;
        this->ptr()->_m_surface = surface;
    }

    //! represents a vertical ray
    Surface_arc_2l(Curved_kernel_via_analysis_2l *kernel,
                   const Surface_point_2l p,
                   CGAL::Arr_curve_end inf_end,
                   const Surface_3& surface) :
        Base() 
        // TODO see above
    {
        _set_ckva(kernel);

        // TODO chack that surface has a vertical line through p
        this->ptr()->_m_is_z_vertical = true;
        this->ptr()->_m_surface = surface;
        // TODO make use of inf_end using private constructors of 
        // Surface_point_2l
    }

    //! represents a vertical branch
    Surface_arc_2l(Curved_kernel_via_analysis_2l *kernel,
                   const Projected_point_2& p,
                   const Surface_3& surface) :
        Base()
        // TODO see above
    {
        _set_ckva(kernel);

        // TODO chack that surface has a vertical line through p
        this->ptr()->_m_is_z_vertical = true;
        this->ptr()->_m_surface = surface;
        // TODO set curve-ends to -oo and +oo using private constructors of 
        // Surface_point_2l
    }
    
    //!@}

public:
    // TODO missing constructors for arcs whose projection is a bounded arc
    // but whose ends (at least one or both) approach a z-vertical asympote

    
    // TODO 
    // access to curve_end(CGAL::Arr_curve_end);

    //!\name Access functions
    //!@{

    /*!\brief
     * returns projected arc
     */
    Projected_arc_2 projected_arc() const {
        return *dynamic_cast< const Projected_arc_2* >(this);
    }

    /*\brief
     * returns the supporting surfaces of 3d-segment
     */
    Surface_3 surface() const {
        return this->ptr()->_m_surface;
    }

    /*!\brief
     * returns the sheet of the 3d-segment
     *
     * \pre !is_z_vertical()
     */
    int sheet() const {
        CGAL_precondition(!is_z_vertical());
        return this->ptr()->_m_sheet;
    }

    /*!\brief
     * returns the sheet of the 3d-segment at min or max end
     *
     * \pre !is_z_vertical()
     */
    int sheet(CGAL::Arr_curve_end end) const {
        CGAL_precondition(!is_z_vertical());
        return (end == CGAL::ARR_MIN_END ? this->ptr()->_m_sheet_min :
                this->ptr()->_m_sheet_max);
    }
    
    /*!\brief
     * return whether arc is z-vertical
     */
    bool is_z_vertical() const {
        return this->ptr()->_m_is_z_vertical;
    }

    //!@}

};    




// FUTURE-TODO a surface can have a vertical plane here over the arc_2
// -> Surface_patch_2l{};


} // namespace CGALi

CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_ARC_2L_H
// EOF
