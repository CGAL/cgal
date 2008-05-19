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
//
// ============================================================================

#ifndef CGAL_SURFACE_POINT_2L_H
#define CGAL_SURFACE_POINT_2L_H

/*! \file Surface_point_2l.h
 *  \brief defines class \c Surface_point_2l
 *  
 *  Kernel for generic points and arcs on surfaces, lifted from 2D.
 */

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>

#include <iostream>
#include <boost/optional.hpp>

#include <CGAL/Curved_kernel_via_analysis_2.h>
#include <CGAL/Curved_kernel_via_analysis_2/Point_2.h>

#include <SoX/GAPS/Restricted_cad_3.h>
#include <SoX/GAPS/Restricted_cad_3_accessor.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

// pre-declaration
template < class CurvedKernelViaAnalysis_2l, class SurfacePair_3, class Rep_ >
class Surface_point_2l;

//! representation type of Surface_point_2l
template < class CurvedKernelViaAnalysis_2l, class SurfacePair_3 >
class Surface_point_2l_rep : 
        public Point_2_rep< CurvedKernelViaAnalysis_2l > {

public:

    //! this instance's first template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;

    //! this instance's second template parameter
    typedef SurfacePair_3 Surface_pair_3;

    //! the instance itself
    typedef 
    Surface_point_2l_rep< Curved_kernel_via_analysis_2l, Surface_pair_3 > Self;

    //! type of base class
    typedef Point_2_rep< Curved_kernel_via_analysis_2l > Base;
    
    //! type of projected kernel
    typedef typename 
    Curved_kernel_via_analysis_2l::Curved_kernel_via_analysis_2
    Projected_kernel_2;
    
    //! type of projected point
    typedef typename Projected_kernel_2::Point_2 Projected_point_2;
    
    //! type of surface
    typedef typename Surface_pair_3::Surface_3 Surface_3;

    //! type of cgal's inexact kernel
    typedef CGAL::Cartesian< double > Kernel;

    //! type of double approximation
    typedef CGAL::Point_3< Kernel > Approximation_3;
    
    //!\name Constructors
    //!@{

    //! default constructor
    Surface_point_2l_rep() :
        Base(), _m_sheet(-1) {
    }
    
    //!@}
public:
    //! projected point
    mutable boost::optional< Projected_point_2 > _m_projected_point;
    
    //! supporting surface
    mutable Surface_3 _m_surface;
    
    //! sheet number of point
    mutable int _m_sheet;

    //! if not z-finite, this members stores whether at -oo or +oo
    mutable boost::optional< CGAL::Arr_curve_end > _m_z_inf_end;
    
    //! approximation
    mutable boost::optional< Approximation_3 > _m_approximation;

    // befriending the handle
    friend class 
    Surface_point_2l< Curved_kernel_via_analysis_2l, Surface_pair_3, Self >;
};


//! represents a point on a surface
template < 
  class CurvedKernelViaAnalysis_2l, 
  class SurfacePair_3,
  class Rep_ = 
    CGALi::Surface_point_2l_rep< CurvedKernelViaAnalysis_2l, SurfacePair_3 >
 >
class Surface_point_2l : public 
    CurvedKernelViaAnalysis_2l::Curved_kernel_via_analysis_2::Point_2::
    template rebind< CurvedKernelViaAnalysis_2l, Rep_  >::Other
{
public:

    //!\name Public types
    //!@{

    //! this instance's first template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;
    
    //! this instance's second template parameter
    typedef SurfacePair_3 Surface_pair_3;

    //! this instance's third template parameter
    typedef Rep_ Rep;

    //! the instance itself
    typedef 
    Surface_point_2l< Curved_kernel_via_analysis_2l, Surface_pair_3, Rep > 
    Self;
    
    //! type of projected kernel
    typedef typename 
    Curved_kernel_via_analysis_2l::Curved_kernel_via_analysis_2
    Projected_kernel_2;
    
    //! type of projected point
    typedef typename Projected_kernel_2::Point_2 Projected_point_2;

    //! type of surface
    typedef typename Surface_pair_3::Surface_3 Surface_3;

    //! the rebinding
    typedef typename Projected_point_2::template 
    rebind< Curved_kernel_via_analysis_2l, Rep  > Rebind;

    //! the base type
    typedef typename Rebind::Other Base;

    //! type of Curve_analysis
    typedef typename Curved_kernel_via_analysis_2l::Curve_analysis_2
    Curve_analysis_2;

    //! type of kernel point
    typedef typename Curved_kernel_via_analysis_2l::Point_2 Kernel_point_2;
    
    //! type of Approximation
    typedef typename Rep::Approximation_3 Approximation_3;

    //!@}

public:

    //!\name Simple Constructors
    //!@{

    /*!\brief
     * Default constructor
     */
    Surface_point_2l() : 
        Base() {   
    }

    //!@}

    //!\name Usual constructors
    //!@{

    //!\brief Constructs point on \c sheet of \c surface above \c point
    //!\pre sheet >= 0
    Surface_point_2l(const Projected_point_2& pt, 
                     const Surface_3& surface, 
                     int sheet) :
        Base(Rebind()(pt)) {
        
        this->copy_on_write();
        
        this->ptr()->_m_projected_point = pt;
        
        this->ptr()->_m_surface = surface;

        CGAL_precondition(sheet >= 0);
        CGAL_precondition_code((
            {
                if (pt.is_finite()) {
                    typedef typename Surface_pair_3::Restricted_cad_3
                        Restricted_cad_3;
                    Restricted_cad_3 cad =
                        Restricted_cad_3::cad_cache()(surface);
                    
                    typedef typename 
                        Surface_pair_3::Restricted_cad_3::Z_stack Z_stack;
                    int number_of_sheets = 
                        cad.z_stack_for(pt).number_of_z_cells();
                    CGAL_precondition(sheet < number_of_sheets);
                } else {
                    // TODO add test for number of sheets in designated face
                }
            })
        );

        this->ptr()->_m_sheet = sheet;
    }

    //!@}

protected:

    //!\name Constructors for special cases
    //!@{

    //!\brief Constructs point at z=+oo/-oo depending on  \c inf_end 
    //! of \c surface above \c point
    //!\pre sheet >= 0
    Surface_point_2l(const Projected_point_2& pt, 
                     const Surface_3& surface,
                     CGAL::Arr_curve_end inf_end) :
        Base(Rebind()(pt)) {

        this->copy_on_write();
        
        // TODO add preconditions?
        // surface has a vertical line, or surface has asymptotes
        
        this->ptr()->_m_projected_point = pt;

        this->ptr()->_m_surface = surface;
        
        this->ptr()->_m_z_inf_end = inf_end;
    }
    
    //!@}
    
protected:
    //!\name Constructor for rebind
    //!@{
    
    /*!\brief
     * constructs from a given represenation
     */
    Surface_point_2l(Rep rep) :
        Base(rep) {
    }

    //!@}
    
public:
    //!\name Access functions
    //!@{

    /*!\brief
     * returns projected point
     */
    inline
    Projected_point_2 projected_point() const {
        if (!this->ptr()->_m_projected_point) {
            CGAL_precondition(dynamic_cast< const Kernel_point_2* >(this));
            this->ptr()->_m_projected_point = 
                typename Kernel_point_2::Rebind()(
                        *dynamic_cast< const Kernel_point_2* >(this)
                );
        }
        return *this->ptr()->_m_projected_point;
    }

    /*\brief
     * returns the supporting surfaces of 3d-point
     */
    inline
    Surface_3 surface() const {
        return this->ptr()->_m_surface;
    }

    /*!\brief
     * returns the sheet number of the 3d-point
     */
    inline
    int sheet() const {
        return this->ptr()->_m_sheet;
    }
    
    /*!\brief
     * returns whether point's z-coordinate is infinite
     */
    inline bool is_z_at_infinity() const {
        return this->ptr()->_m_z_inf_end;
    }

    /*!\brief
     * returns whether point's z-coordinate is infinite
     */
    inline CGAL::Arr_curve_end z_infinity() const {
        return *this->ptr()->_m_z_inf_end;
    }
    
    //!@}
    
#define CGAL_CKvA_2l_GRAB_CK_FUNCTOR_FOR_POINT(X, Y, Z) \
    typename Curved_kernel_via_analysis_2l::X Y = \
         Curved_kernel_via_analysis_2l::instance().Z(); \


public:
    //!\name Comparisons
    //!@{
    
    //!\brief compares two points xyz-lexicographically
    //!
    //!\pre compared points have finite x/y-coordinates
    inline
    CGAL::Comparison_result compare_xyz(const Kernel_point_2& q, 
                                        bool equal_xy = false) const {
        CGAL_precondition(this->is_finite());
        CGAL_precondition(q.is_finite());
        
        CGAL_CKvA_2l_GRAB_CK_FUNCTOR_FOR_POINT(Compare_xyz_3, 
                                               compare_xyz_3,
                                               compare_xyz_3_object);
        CGAL_precondition(dynamic_cast< const Kernel_point_2* >(this));
        return compare_xyz_3(
                *dynamic_cast< const Kernel_point_2* >(this), q, equal_xy
        );
    }

    //!@}

#undef CGAL_CKvA_2l_GRAB_CK_FUNCTOR_FOR_POINT
    
    //!\name Approximation
    //!@{
    
    // returns an non-robust approximation of the point
    Approximation_3 to_double() const {
        
        if (!this->ptr()->_m_approximation) {

            typedef typename Curve_analysis_2::Algebraic_curve_kernel_2::
                Boundary Rational;
            
            double x, y, z;
            
            // X
            {
                x = CGAL::to_double(this->x()); // TODO replace by CGAL::
            }

            Rational bound(2e-10); // TODO correct threshold?
            // Y
            {
                typename Curve_analysis_2::Algebraic_curve_kernel_2::
                    Xy_coordinate_2 xy = 
                    this->curve().status_line_at_exact_x(this->x()).
                    algebraic_real_2(this->arcno());
                
                typename Curve_analysis_2::Algebraic_curve_kernel_2::
                    Lower_boundary_y_2 lower_boundary_y;
                typename Curve_analysis_2::Algebraic_curve_kernel_2::
                    Upper_boundary_y_2 upper_boundary_y;
                typename Curve_analysis_2::Algebraic_curve_kernel_2::
                    Refine_y_2 refine_y;
                
                while (upper_boundary_y(xy) - lower_boundary_y(xy) > bound) {
                    refine_y(xy);
                }
                
                y = CGAL::to_double(lower_boundary_y(xy));
            }
            // TODO replace xy by planar approximation (renderer, Pavel)
            // Z
            {
                typedef typename Surface_pair_3::Restricted_cad_3
                    Restricted_cad_3;
                typedef typename Surface_pair_3::Z_at_xy_isolator
                    Z_at_xy_isolator;
                Restricted_cad_3 cad =
                    Restricted_cad_3::cad_cache()(this->surface());
                boost::optional< Z_at_xy_isolator > isolator =
                    cad.isolator_at(this->projected_point(),
                                    this->surface());
                CGAL_assertion(isolator);
                
                while (isolator->length(this->sheet()) > bound) {
                    isolator->refine_interval(this->sheet());
                }
                
                z = CGAL::to_double(isolator->left_boundary(this->sheet()));
            }
            
            this->ptr()->_m_approximation = Approximation_3(x,y,z);
        }
        
        CGAL_postcondition(this->ptr()->_m_approximation);
        return *this->ptr()->_m_approximation;
    }

    //!@}

public:
    //!\name IO
    //!@{
    
    //! write represenation to \c os
    void write(std::ostream& os) const {
        os << "Point_2l(";
        os << "Point2(" << this->projected_point() << ", loc=" 
           << this->location() << "), ";
        os << "Surface(" << this->surface() << ", ";
        if (this->is_z_at_infinity()) {
            if (this->z_infinity() == CGAL::ARR_MIN_END) {
                os << "@-oo";
            } else {
                os << "@+oo";
            }
        } else {
            os << this->sheet();
        }
        os << ")";
        os << ")" << std::flush;
    }

    //!@}

    //!\name Friends
    //!@{

    //! for special arc constructors
    friend class Curved_kernel_via_analysis_2l::Arc_2;

    //! for rebind
    friend class Self::Rebind;

    //!@}

}; // Surface_point_2l


/*!\relates Surface_point_2l
 * \brief 
 * output operator
 */
template < class CurvedKernelViaAnalysis_2l, class SurfacePair_3 >
std::ostream& operator<< (
        std::ostream& os,
        const 
        Surface_point_2l< CurvedKernelViaAnalysis_2l, SurfacePair_3 >& 
        pt) {
    
    pt.write(os);
    
    return os;
}

} // namespace CGALi

CGAL_END_NAMESPACE


#endif // CGAL_SURFACE_POINT_2L_H
// EOF
