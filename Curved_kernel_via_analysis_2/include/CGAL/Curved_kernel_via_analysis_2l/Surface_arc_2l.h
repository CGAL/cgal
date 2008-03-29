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

#include <SoX/GAPS/Restricted_cad_3.h>
#include <SoX/GAPS/Restricted_cad_3_accessor.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

// pre-declaration
template < class CurvedKernelViaAnalysis_2l, class SurfacePair_3, class Rep_ >
class Surface_arc_2l;

//! representation class for arc on a surface
template < class CurvedKernelViaAnalysis_2l, class SurfacePair_3 >
class Surface_arc_2l_rep : 
        public Arc_2_rep< CurvedKernelViaAnalysis_2l > {

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
    typedef Arc_2_rep< Curved_kernel_via_analysis_2l > Base;
    
    //! typedef of Point_2 (Kernel_point_2)
    typedef typename Curved_kernel_via_analysis_2l::Point_2 Point_2;

    //! type of projected kernel
    typedef typename 
    Curved_kernel_via_analysis_2l::Curved_kernel_via_analysis_2
    Projected_kernel_2;
    
    //! type of projected point
    typedef typename Projected_kernel_2::Point_2 Projected_point_2;

    //! type of projected point
    typedef typename Projected_kernel_2::Arc_2 Projected_arc_2;
    
    //! type of surface
    typedef typename Surface_pair_3::Surface_3 Surface_3;

    //!\name Constructors
    //!@{

    //! standard constructor
    Surface_arc_2l_rep() :
        _m_sheet(-1), _m_sheet_min(-1), _m_sheet_max(-1), 
        _m_is_z_vertical(false) {
    }

    //! constructor for vertical arcs
    Surface_arc_2l_rep(const Point_2& p, const Point_2& q,
                       const Surface_3& surface) :
        _m_projected_point(p.projected_point()),
        _m_surface(surface),
        _m_sheet(-1), _m_sheet_min(-1), _m_sheet_max(-1),
        _m_is_z_vertical(true) {

        CGAL_precondition(p.compare_xy(q) == CGAL::EQUAL);
        CGAL::Comparison_result cmp = p.compare_xyz(p, true);
        CGAL_precondition(cmp != CGAL::EQUAL);
        if (cmp == CGAL::LARGER) {
            Base::_m_min = q;
            Base::_m_max = p;
        }
        // else
        Base::_m_min = p;
        Base::_m_max = q;
    }
    
    //!@}
    
    
protected:
    //! projected arc
    mutable boost::optional< Projected_arc_2 > _m_projected_arc;

    //! projected point (for z-vertical arcs)
    mutable boost::optional< Projected_point_2 > _m_projected_point;

    //! supporting surface
    mutable Surface_3 _m_surface;

    //! sheet number of arc
    mutable int _m_sheet;

    //! sheet number of arc at x-min point
    mutable int _m_sheet_min;
    
    //! sheet number of arcs at x-max point
    mutable int _m_sheet_max;
    
    //! indicates whether arc is vertical
    mutable bool _m_is_z_vertical;
    
    // befriending the handle
    friend class 
    Surface_arc_2l< Curved_kernel_via_analysis_2l, Surface_pair_3, Self >;

    // for setting sheets + projected_arcs
    friend class Curved_kernel_via_analysis_2l::Trim_2;
    friend class Curved_kernel_via_analysis_2l::Split_2;
    friend class Curved_kernel_via_analysis_2l::Merge_2;
};


//! represents an xy-monotone arc on a surface
template < 
  class CurvedKernelViaAnalysis_2l, 
  class SurfacePair_3,
  class Rep_ = 
    CGALi::Surface_arc_2l_rep< CurvedKernelViaAnalysis_2l, SurfacePair_3 >
>
class Surface_arc_2l : public 
   CurvedKernelViaAnalysis_2l::Curved_kernel_via_analysis_2::Arc_2::
   template rebind< CurvedKernelViaAnalysis_2l, Rep_ >::Other {
    
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
    
    //! type of projected kernel
    typedef typename 
    Curved_kernel_via_analysis_2l::Curved_kernel_via_analysis_2
    Projected_kernel_2;
    
    //! type of projected point
    typedef typename Projected_kernel_2::Point_2 Projected_point_2;

    //! type of projected point
    typedef typename Projected_kernel_2::Arc_2 Projected_arc_2;

    //! type of surface
    typedef typename Surface_pair_3::Surface_3 Surface_3;
    
    //! type of planar point
    typedef typename Curved_kernel_via_analysis_2l::Point_2 Surface_point_2l;

    //! typedef of Kernel_arc_2
    typedef typename Curved_kernel_via_analysis_2l::Arc_2 Kernel_arc_2;
    
    //! type of rebinding
    typedef typename Projected_arc_2::
    template rebind < Curved_kernel_via_analysis_2l, Rep > Rebind;
    
    //! the base class
    typedef typename Rebind::Other Base;

    //!@}

public:    
    //!\name Standard constructors
    //!@{

    /*!\brief
     * Default constructor
     */
    Surface_arc_2l() : 
        Base() {   
    }
    
    //!@}

public:
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
    Surface_arc_2l(const Projected_arc_2& arc, 
                   const Surface_point_2l& p,
                   const Surface_point_2l& q,
                   const Surface_3& surface,
                   int sheet, int sheet_p, int sheet_q) :
        Base(Rebind()(arc, p, q)) {
        
        this->copy_on_write();
        
        CGAL_precondition(arc.is_finite(CGAL::ARR_MIN_END));
        CGAL_precondition(arc.is_finite(CGAL::ARR_MAX_END));

        CGAL_precondition(
                p.projected_point().
                compare_xy(arc.curve_end(CGAL::ARR_MIN_END)) ==
                CGAL::EQUAL
        );
        CGAL_precondition(
                q.projected_point().
                compare_xy(arc.curve_end(CGAL::ARR_MAX_END)) ==
                CGAL::EQUAL
        );
        
        this->ptr()->_m_projected_arc = arc;
        
        this->ptr()->_m_surface = surface;

        CGAL_precondition_code(
                int number_of_sheets = -1;
                typedef typename Surface_pair_3::Restricted_cad_3
                Restricted_cad_3;
                typedef SoX::Restricted_cad_3_accessor< Restricted_cad_3 > 
                Accessor;
                Restricted_cad_3 cad =
                Restricted_cad_3::cad_cache()(surface);
        );
        
        CGAL_precondition(sheet >= 0);
        CGAL_precondition_code(
                Projected_point_2 pt = Accessor::point_in_interior(arc);
                number_of_sheets = cad.z_stack_at(pt).number_of_z_cells();
        );
        CGAL_precondition(sheet < number_of_sheets);

        this->ptr()->_m_sheet = sheet;

        CGAL_precondition(sheet_p >= 0);
        CGAL_precondition_code(
                number_of_sheets = 
                cad.z_stack_at(p.projected_point()).number_of_z_cells();
        );
        CGAL_precondition(sheet_p < number_of_sheets);
        CGAL_precondition_code(
                // TODO add sanity checks for sheet_at_min/max wrt sheet (eriC
                // use "adjacency information" 
        );

        this->ptr()->_m_sheet_min = sheet_p;
        
        CGAL_precondition(sheet_q >= 0);
        CGAL_precondition_code(
                number_of_sheets = 
                cad.z_stack_at(q.projected_point()).number_of_z_cells();
        );
        CGAL_precondition(sheet_q < number_of_sheets);
        CGAL_precondition_code(
                // TODO add sanity checks for sheet_at_min/max wrt sheet (eriC
                // use "adjacency information" 
        );
        
        this->ptr()->_m_sheet_max = sheet_q;
        
        CGAL_postcondition(this->_check_surface_arc_interior());
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
    Surface_arc_2l(const Projected_arc_2& arc, 
                   const Surface_point_2l& p,
                   CGAL::Arr_curve_end z_inf_end_other,
                   const Surface_3& surface,
                   int sheet, int sheet_p) :
        Base(Rebind()(arc,
                      (arc.curve_end(CGAL::ARR_MIN_END).compare_xy(
                              p.projected_point()
                      ) == CGAL::EQUAL ?
                       p : 
                       Surface_point_2l(Rebind()(arc, CGAL::ARR_MIN_END), 
                                        surface, z_inf_end_other)
                      ),
                      (arc.curve_end(CGAL::ARR_MAX_END).compare_xy(
                              p.projected_point()
                      ) == CGAL::EQUAL ?
                       p : 
                       Surface_point_2l(Rebind()(arc, CGAL::ARR_MAX_END), 
                                        surface, z_inf_end_other)
                      )
             )
        ) {
        
        this->copy_on_write();

        CGAL_precondition(arc.is_finite(CGAL::ARR_MIN_END));
        CGAL_precondition(arc.is_finite(CGAL::ARR_MAX_END));
        
        this->ptr()->_m_projected_arc = arc;

        this->ptr()->_m_surface = surface;

        bool p_at_min = 
            p.projected_point().
            compare_xy(arc.curve_end(CGAL::ARR_MIN_END)) ==
            CGAL::EQUAL;

        CGAL_precondition_code(
                int number_of_sheets = -1;
                typedef typename Surface_pair_3::Restricted_cad_3
                Restricted_cad_3;
                typedef SoX::Restricted_cad_3_accessor< Restricted_cad_3 > 
                Accessor;
                Restricted_cad_3 cad =
                Restricted_cad_3::cad_cache()(surface);
        );
        
        CGAL_precondition(sheet >= 0);
        CGAL_precondition_code(
                Projected_point_2 pt = Accessor::point_in_interior(arc);
                number_of_sheets = cad.z_stack_at(pt).number_of_z_cells();
        );
        CGAL_precondition(sheet < number_of_sheets);

        this->ptr()->_m_sheet = sheet;

        CGAL_precondition(sheet_p >= 0);
        CGAL_precondition_code(
                number_of_sheets = 
                cad.z_stack_at(p.projected_point()).number_of_z_cells();
        );
        CGAL_precondition(sheet_p < number_of_sheets);
        
        CGAL_precondition_code(
                // TODO add sanity checks for sheet_at_min/max wrt sheet (eriC
                // use "adjacency information" 
        );
        
        if (p_at_min) {
            this->ptr()->_m_sheet_min = sheet_p;
            this->ptr()->_m_sheet_max = sheet;
        } else {
            this->ptr()->_m_sheet_min = sheet;
            this->ptr()->_m_sheet_max = sheet_p;
        }

        CGAL_postcondition(this->_check_surface_arc_interior());
    }
    
    /*!\brief
     * Standard constructor for an unbounded arc on xy-monotone part
     * of the surface with z-asympotic behavior at both ends.
     * It represents the arc on \c surface covertical to \c arc which
     * lies on \c sheet of the xy-monotone subsurface.
     */
    Surface_arc_2l(const Projected_arc_2& arc, 
                   CGAL::Arr_curve_end z_inf_end_p,
                   CGAL::Arr_curve_end z_inf_end_q,
                   const Surface_3& surface,
                   int sheet) :
        Base(Rebind()(arc, 
                      Surface_point_2l(arc.curve_end(CGAL::ARR_MIN_END), 
                                       surface, z_inf_end_p),
                      Surface_point_2l(arc.curve_end(CGAL::ARR_MAX_END), 
                                       surface, z_inf_end_q))
        ) {
        
        this->copy_on_write();
        
        CGAL_precondition(arc.is_finite(CGAL::ARR_MIN_END));
        CGAL_precondition(arc.is_finite(CGAL::ARR_MAX_END));
        
        this->ptr()->_m_projected_arc = arc;

        this->ptr()->_m_surface = surface;

        CGAL_precondition_code(
                int number_of_sheets = -1;
                typedef typename Surface_pair_3::Restricted_cad_3
                Restricted_cad_3;
                typedef SoX::Restricted_cad_3_accessor< Restricted_cad_3 > 
                Accessor;
                Restricted_cad_3 cad =
                Restricted_cad_3::cad_cache()(surface);
        );
        
        CGAL_precondition(sheet >= 0);
        CGAL_precondition_code(
                Projected_point_2 pt = Accessor::point_in_interior(arc);
                number_of_sheets = cad.z_stack_at(pt).number_of_z_cells();
        );
        CGAL_precondition(sheet < number_of_sheets);

        CGAL_precondition_code(
                // TODO add sanity checks for sheet_at_min/max wrt sheet (eriC
                // use "adjacency information" 
        );

        this->ptr()->_m_sheet = sheet;
        this->ptr()->_m_sheet_min = sheet;
        this->ptr()->_m_sheet_max = sheet;

        CGAL_postcondition(this->_check_surface_arc_interior());
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
    Surface_arc_2l(const Projected_arc_2& arc, 
                   const Surface_point_2l& p,
                   const Surface_3& surface,
                   int sheet, int sheet_p) :
        Base(Rebind()(arc, 
                      (arc.is_finite(CGAL::ARR_MIN_END) ?
                       p :
                       Surface_point_2l(Rebind()(arc, CGAL::ARR_MIN_END), 
                                        surface, sheet)),
                      (arc.is_finite(CGAL::ARR_MAX_END) ? 
                       p :
                       Surface_point_2l(Rebind()(arc, CGAL::ARR_MAX_END),
                                        surface, sheet)))
        ) {
        
        this->copy_on_write();
        
        this->ptr()->_m_projected_arc = arc;

        this->ptr()->_m_surface = surface;
        
        bool p_at_min = arc.is_finite(CGAL::ARR_MIN_END);
        CGAL_precondition_code(
                bool p_at_max = (arc.is_finite(CGAL::ARR_MAX_END));
        );
        CGAL_precondition(p_at_min || p_at_max && !(p_at_min && p_at_max));
        
        CGAL_precondition_code(
                int number_of_sheets = -1;
                typedef typename Surface_pair_3::Restricted_cad_3
                Restricted_cad_3;
                typedef SoX::Restricted_cad_3_accessor< Restricted_cad_3 > 
                Accessor;
                Restricted_cad_3 cad =
                Restricted_cad_3::cad_cache()(surface);
        );
        
        CGAL_precondition(sheet >= 0);
        CGAL_precondition_code(
                Projected_point_2 pt = Accessor::point_in_interior(arc);
                number_of_sheets = 
                cad.z_stack_at(pt).number_of_z_cells();
        );
        CGAL_precondition(sheet < number_of_sheets);
        
        this->ptr()->_m_sheet = sheet;
        
        CGAL_precondition(sheet_p >= 0);
        CGAL_precondition_code(
                number_of_sheets = 
                cad.z_stack_at(p.projected_point()).number_of_z_cells();
        );
        CGAL_precondition(sheet_p < number_of_sheets);
        
        CGAL_precondition_code(
                // TODO add sanity checks for sheet_at_min/max wrt sheet (eriC)
                // use "adjacency information"
        );
        
        if (p_at_min) {
            this->ptr()->_m_sheet_min = sheet_p;
            this->ptr()->_m_sheet_max = sheet;
        } else {
            this->ptr()->_m_sheet_min = sheet;
            this->ptr()->_m_sheet_max = sheet_p;
        }

        CGAL_postcondition(this->_check_surface_arc_interior());
    }

    /*!\brief
     * Standard constructor for a branch on a xy-monotone part
     * of the surface with a z-vertical asymptotic behaviour at the projected
     * finite end.
     * It represents the arc on \c surface covertical to \c arc which
     * lies on \c sheet of the xy-monotone subsurface
     */
    Surface_arc_2l(const Projected_arc_2& arc, 
                   CGAL::Arr_curve_end z_inf_end,
                   const Surface_3& surface,
                   int sheet) :
        Base(Rebind()(arc, 
                      (arc.is_finite(CGAL::ARR_MIN_END) ? 
                       Surface_point_2l(arc.curve_end(CGAL::ARR_MIN_END), 
                                        surface, z_inf_end) :
                       Surface_point_2l(Rebind()(arc, CGAL::ARR_MIN_END),
                                        surface, sheet)),
                      (arc.is_finite(CGAL::ARR_MAX_END) ? 
                       Surface_point_2l(arc.curve_end(CGAL::ARR_MAX_END), 
                                        surface, z_inf_end) :
                       Surface_point_2l(Rebind()(arc, CGAL::ARR_MAX_END), 
                                        surface, sheet)))
        ) {
        
        this->copy_on_write();
        
        bool min_finite = (arc.is_finite(CGAL::ARR_MIN_END));
        CGAL_precondition_code(
                bool max_finite = (arc.is_finite(CGAL::ARR_MAX_END));
        );
        CGAL_precondition(min_finite || max_finite && 
                          !(!min_finite && !max_finite));
        
        this->ptr()->_m_projected_arc = arc;
        
        this->ptr()->_m_surface = surface;

        CGAL_precondition_code(
                int number_of_sheets = -1;
                typedef typename Surface_pair_3::Restricted_cad_3
                Restricted_cad_3;
                typedef SoX::Restricted_cad_3_accessor< Restricted_cad_3 > 
                Accessor;
                Restricted_cad_3 cad =
                Restricted_cad_3::cad_cache()(surface);
        );
        
        CGAL_precondition(sheet >= 0);
        CGAL_precondition_code(
                Projected_point_2 pt = Accessor::point_in_interior(arc);
                number_of_sheets = 
                cad.z_stack_at(pt).number_of_z_cells();
        );
        CGAL_precondition(sheet < number_of_sheets);

        this->ptr()->_m_sheet = sheet;

        CGAL_precondition_code(
                // TODO add sanity checks for sheet_at_min/max wrt sheet (eriC)
                // use "adjacency information"
        );
        
        this->ptr()->_m_sheet_min = sheet;
        this->ptr()->_m_sheet_max = sheet;

        CGAL_postcondition(this->_check_surface_arc_interior());
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
    Surface_arc_2l(const Projected_arc_2& arc, 
                   const Surface_3& surface,
                   int sheet) :
        Base(Rebind()(arc,
                      Surface_point_2l(Rebind()(arc, CGAL::ARR_MIN_END), 
                                       surface, sheet),
                      Surface_point_2l(Rebind()(arc, CGAL::ARR_MAX_END),
                                       surface, sheet))
        ) {
        
        this->copy_on_write();
        
        CGAL_precondition(!arc.is_finite(CGAL::ARR_MIN_END));
        CGAL_precondition(!arc.is_finite(CGAL::ARR_MAX_END));
        
        this->ptr()->_m_projected_arc = arc;
        
        this->ptr()->_m_surface = surface;

        CGAL_precondition_code(
                int number_of_sheets = -1;
                typedef typename Surface_pair_3::Restricted_cad_3
                Restricted_cad_3;
                typedef SoX::Restricted_cad_3_accessor< Restricted_cad_3 > 
                Accessor;
                Restricted_cad_3 cad =
                Restricted_cad_3::cad_cache()(surface);
        );
        CGAL_precondition(sheet >= 0);
        CGAL_precondition_code(
                Projected_point_2 pt = Accessor::point_in_interior(arc);
                number_of_sheets = 
                cad.z_stack_at(pt).number_of_z_cells();
        );
        
        this->ptr()->_m_sheet = sheet;
        
        this->ptr()->_m_sheet_min = sheet;
        this->ptr()->_m_sheet_max = sheet;

        CGAL_postcondition(this->_check_surface_arc_interior());
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
    Surface_arc_2l(const Surface_point_2l& p,
                   const Surface_point_2l& q,
                   const Surface_3& surface) :
        Base(Rep(p, q, surface)) {
        
        CGAL_precondition(p.projected_point().
                          compare_xy(q.projected_point()) == CGAL::EQUAL);
        CGAL_precondition_code(
                typedef typename Surface_pair_3::Restricted_cad_3
                Restricted_cad_3;
                Restricted_cad_3 cad =
                Restricted_cad_3::cad_cache()(surface);
        );
        CGAL_precondition(
                cad.supports_vertical_line_at(p.projected_point(), surface)
        );

        CGAL_postcondition(this->_check_surface_arc_interior());
    }

    //! represents a vertical ray
    Surface_arc_2l(const Surface_point_2l& p,
                   CGAL::Arr_curve_end z_inf_end,
                   const Surface_3& surface) :
        Base(Rep(p, 
                 Surface_point_2l(p.projected_point(), surface, z_inf_end), 
                 surface
             )
        ) {
        
        CGAL_precondition_code(
                typedef typename Surface_pair_3::Restricted_cad_3
                Restricted_cad_3;
                Restricted_cad_3 cad =
                Restricted_cad_3::cad_cache()(surface);
        );
        CGAL_precondition(
                cad.supports_vertical_line_at(p.projected_point(), surface)
        );

        CGAL_postcondition(this->_check_surface_arc_interior());
    }

    //! represents a vertical branch
    Surface_arc_2l(const Projected_point_2& p,
                   const Surface_3& surface) :
        Base(Rep(Surface_point_2l(p, surface, CGAL::ARR_MIN_END), 
                 Surface_point_2l(p, surface, CGAL::ARR_MAX_END), 
                 surface
             )
        ) {

        CGAL_precondition_code(
                typedef typename Surface_pair_3::Restricted_cad_3
                Restricted_cad_3;
                Restricted_cad_3 cad =
                Restricted_cad_3::cad_cache()(surface);
        );
        CGAL_precondition(cad.supports_vertical_line_at(p, surface));

        CGAL_postcondition(this->_check_surface_arc_interior());
    }
    
    //!@}

protected:
    //!\name Constructors for rebind/replace_endpoints
    //!@{
    
    /*!\brief
     * constructs an arc from a given represenation
     */
    Surface_arc_2l(Rep rep) : 
        Base(rep) { 
    }

    //!@}

protected:
    //!\name Sanity checks
    //!@{

    //! checks whether arc is constructed properly
    bool _check_surface_arc_interior() {
        // TODO implement _check_surface_arc_interior
        return true;
    }

    //!@}

public:

    //!\name Access functions
    //!@{

    /*!\brief
     * returns projected arc
     */
    Projected_arc_2 projected_arc() const {
        CGAL_precondition(!this->is_z_vertical());
        if (!this->ptr()->_m_projected_arc) {
            CGAL_precondition(dynamic_cast< const Kernel_arc_2* >(this));
            this->ptr()->_m_projected_arc = 
                typename Kernel_arc_2::Rebind()(
                        *dynamic_cast< const Kernel_arc_2* >(this)
                );
        }
        return *(this->ptr()->_m_projected_arc);
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
     * returns the sheet of the 3d-segment at a given projected point \c pt
     *
     * \pre !is_z_vertical()
     */
    int sheet(const Projected_point_2& pt) const {
        CGAL_precondition(!is_z_vertical());
        CGAL_precondition(
                this->projected_arc().compare_y_at_x(pt) == CGAL::EQUAL
        );
        if (this->sheet() != this->sheet(CGAL::ARR_MIN_END) && 
            this->projected_arc().is_finite(CGAL::ARR_MIN_END) && 
            this->projected_arc().curve_end(CGAL::ARR_MIN_END) == pt) {
            return this->ptr()->_m_sheet_min;
        } else if (this->sheet() != this->sheet(CGAL::ARR_MAX_END) && 
                   this->projected_arc().is_finite(CGAL::ARR_MAX_END) && 
                   this->projected_arc().curve_end(CGAL::ARR_MAX_END) == pt) {
            return this->ptr()->_m_sheet_max;
        }
        return this->ptr()->_m_sheet;
    }
    
    /*!\brief
     * return whether arc is z-vertical
     */
    bool is_z_vertical() const {
        return this->ptr()->_m_is_z_vertical;
    }

    /*!\brief
     * returns projected point of vertical arc
     *
     * \pre is_z_vertical()
     */
    Projected_point_2 projected_point() const {
        CGAL_precondition(this->is_z_vertical());
        CGAL_precondition(this->ptr()->_m_projected_point);
        return *(this->ptr()->_m_projected_point);
    }
    
    //!@}


protected:
    //!\name Proteced members
    //!@{

    /*!\brief 
     * replaces this arc's end-points by \c src and \c tgt with arcnos
     * \c arcno_min and \c arcno_max.
     * 
     * new curve ends are sorted lexicographical in case of need; 
     * all preconditions must be checked by the caller
     */
    std::pair< Kernel_arc_2, CGAL::Comparison_result > 
    _replace_endpoints(
            const Surface_point_2l& p1, const Surface_point_2l& p2,
            int arcno1, int arcno2,
            int sheet1, int sheet2) const {
        
        CERR("\n_sa_2l_replace_endpoints\n");    

        if (is_z_vertical()) {
            Rep rep(*this->ptr());
            CGAL::Comparison_result cmp = p1.compare_xyz(p2);
            if (cmp == CGAL::LARGER) {
                rep._m_min = p2;
                rep._m_max = p1;
            } else {
                rep._m_min = p1;
                rep._m_max = p2;
            }
            return std::make_pair(Kernel_arc_2(rep), cmp);
        } 
        
        // else
        std::pair< Kernel_arc_2, CGAL::Comparison_result >
            replaced = Base::_replace_endpoints(
                    p1, p2, arcno1, arcno2
            );

        if (replaced.second == CGAL::LARGER) {
            std::swap(sheet1, sheet2);
        }
        
        if (sheet1 >= 0) {
            replaced.first.ptr()->_m_sheet_min = sheet1;
        }
        if (sheet2 >= 0) {
            replaced.first.ptr()->_m_sheet_max = sheet2;
        }
        
        replaced.first.ptr()->_m_projected_arc = boost::none;
        
        return replaced;
    }

    //!@}

public:
    //!\name IO
    //!@{
    
    //! write represenation to \c os
    void write(std::ostream& os) const { 
        os << "Arc_2l(";
        if (this->is_z_vertical()) {
            os << "Point_2(" << this->projected_point() << "), ";
            os << "MinPoint(" << this->ptr()->_m_min << "), ";
            os << "MaxPoint(" << this->ptr()->_m_max << "), ";
            os << "Surface(" << this->surface() << ", Z-VERT)";
        } else {
            os << "Arc_2(" << this->projected_arc() << "), ";
            os << "MinPoint(" << this->ptr()->_m_min << "), ";
            os << "MaxPoint(" << this->ptr()->_m_max << "), ";
            os << "Surface(" << this->surface() << ", " 
               << this->sheet() << ", " << this->sheet(CGAL::ARR_MIN_END)
               << ", " << this->sheet(CGAL::ARR_MAX_END)
               << ")";
        }
        os << std::flush;
    }

    //!@}

    //!\name Friends
    //!@{

    //! for _replace_points
    friend class Curved_kernel_via_analysis_2l::Trim_2;
    friend class Curved_kernel_via_analysis_2l::Split_2;
    friend class Curved_kernel_via_analysis_2l::Merge_2;
   
    //! for replace endpoints
    friend class Self::Rebind;
    
    //! for rebind
    friend class Self::Rebind::Other;
    
    //!@}
};    

/*!\relates Surface_arc_2l
 * \brief 
 * output operator
 */
template < class CurvedKernelViaAnalysis_2l, class SurfacePair_3, class Rep_ >
std::ostream& operator<< (
        std::ostream& os,
        const 
        Surface_arc_2l< CurvedKernelViaAnalysis_2l, SurfacePair_3, Rep_ >& 
        arc) {
    
    arc.write(os);
    
    return os;
}

} // namespace CGALi

CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_ARC_2L_H
// EOF
