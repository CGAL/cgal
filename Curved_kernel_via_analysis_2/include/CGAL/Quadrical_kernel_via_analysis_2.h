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

#ifndef CGAL_QUADRICAL_KERNEL_VIA_ANALYSIS_2_H
#define CGAL_QUADRICAL_KERNEL_VIA_ANALYSIS_2_H

/*! \file Quadrical_kernel_via_analysis_2.h
 *  \brief defines class \c Quadrical_kernel_via_analysis_2
 *  
 *  Kernel for lifted generic points and arcs on embedded on a quadric
 */

#include <CGAL/basic.h>

#include <CGAL/Curved_kernel_via_analysis_2.h>
#include <CGAL/Curved_kernel_via_analysis_2/Curved_kernel_via_analysis_2_functors.h>

#include <CGAL/Curved_kernel_via_analysis_2l/Surface_point_2l.h>
#include <CGAL/Curved_kernel_via_analysis_2l/Surface_arc_2l.h>
#include <CGAL/Curved_kernel_via_analysis_2l/Curved_kernel_via_analysis_2l_functors.h>

#include <QdX/Quadric_pair_3.h>
#include <QdX/Quadric_3_z_at_xy_isolator_traits.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

// pre-declaration
template < class QuadricalKernelViaAnalysis_2, class SurfacePair_3 >
class Quadric_point_2;


template < class QuadricalKernelViaAnalysis_2, class SurfacePair_3 >
class Quadric_point_2_rep : 
    public Surface_point_2l_rep< QuadricalKernelViaAnalysis_2, SurfacePair_3 > 
{
public:
    //! this instance's first template parameter
    typedef QuadricalKernelViaAnalysis_2 Quadrical_kernel_via_analysis_2;

    //! this instance's second template parameter
    typedef SurfacePair_3 Surface_pair_3;

    //! the instance itself
    typedef 
    Surface_point_2l_rep< Quadrical_kernel_via_analysis_2, Surface_pair_3 > 
    Self;
    
    //! base class
    typedef 
    Surface_point_2l_rep< Quadrical_kernel_via_analysis_2, Surface_pair_3 >
    Base;
    
    //! type of curve
    typedef typename Base::Xy_coordinate_2 Xy_coordinate_2;

    //!\name Constructors
    //!@{
    
    //! default constructor
    Quadric_point_2_rep() :
        Base() {
    }

    //! standard constructor 
    Quadric_point_2_rep(const Xy_coordinate_2& xy) :
        Base(xy) {
    }
    
    //!@}

protected:

    // befriending the handle
    friend class 
    Quadric_point_2< Quadrical_kernel_via_analysis_2, Surface_pair_3 >;
};


//! represent point on a quadric
template < class QuadricalKernelViaAnalysis_2, class SurfacePair_3 >
class Quadric_point_2 : 
    public CGALi::Surface_point_2l< 
        QuadricalKernelViaAnalysis_2, 
        SurfacePair_3,
        CGALi::Quadric_point_2_rep< 
            QuadricalKernelViaAnalysis_2, SurfacePair_3 
        > > {
public:
    //! this instance's first template parameter
    typedef QuadricalKernelViaAnalysis_2 Quadrical_kernel_via_analysis_2;

    //! this instance's second template parameter
    typedef SurfacePair_3 Surface_pair_3;

    //! the instance itself
    typedef 
    Quadric_point_2< Quadrical_kernel_via_analysis_2, Surface_pair_3 > Self;
    
    //! the type of the representation
    typedef 
    CGALi::Quadric_point_2_rep< 
    Quadrical_kernel_via_analysis_2, Surface_pair_3 > Rep;

    //! the base type
    typedef CGALi::Surface_point_2l< 
    Quadrical_kernel_via_analysis_2, Surface_pair_3 , Rep > 
    Base;

    //! type of surface
    typedef typename Surface_pair_3::Surface_3 Surface_3;
    
    //! type of planar point
    typedef typename Base::Projected_point_2 Projected_point_2;
    
public:
    //!\name Standard constructors
    //!@{
    
    /*!\brief
     * Default constructor
     */
    Quadric_point_2() : 
        Base() {   
    }
    
    //!@}

public:
    //!\name Usual constructors 
    //!@{
    
    //!\brief Constructs point on \c sheet of \c surface above \c point
    //!\pre sheet >= 0
    Quadric_point_2(const Projected_point_2& pt, 
                    const Surface_3& surface, 
                    int sheet) :
        Base(pt, surface, sheet) {
        CGAL_precondition(sheet < 2);
    }
    
    //!@}
    
protected:
    //!\name special constructor for point with z=+-oo
    //!@{
    
    //!\brief Constructs point with z = -oo or z = +oo depending on 
    //! \c z_inf_end of \c surface above \c point
    //!\pre sheet >= 0
    Quadric_point_2(const Projected_point_2& pt, 
                    const Surface_3& surface,
                    CGAL::Arr_curve_end z_inf_end) :
        Base(pt, surface, z_inf_end) {
    }
    
    //!@}

protected: 
    //!\name Constructors for rebind
    //!@{
    
    /*!\brief
     * constructs from a given represenation
     */
    Quadric_point_2(Rep rep) :
        Base(rep) {
    }
    //!@}
    
    //! for constructint points
    friend class Quadrical_kernel_via_analysis_2::Construct_point_2;

    //! for rebind
    friend class Base::Rebind;

    //! for arc
    friend class Quadrical_kernel_via_analysis_2::Arc_2;

    //! for point constructors in arc constructors
    friend class CGALi::Surface_arc_2l< Quadrical_kernel_via_analysis_2, 
        Surface_pair_3,
        CGALi::Surface_arc_2l_rep< Quadrical_kernel_via_analysis_2, 
        Surface_pair_3 > >;

};

/*!\relates Quadric_point_2
 * \brief 
 * output operator
 */
template < class QuadricalKernelViaAnalysis_2, class SurfacePair_3 >
std::ostream& operator<< (
        std::ostream& os,
        const 
        Quadric_point_2< QuadricalKernelViaAnalysis_2, SurfacePair_3 >& 
        pt) {
    
    pt.write(os);
    
    return os;
}


// pre-declaration
template < class QuadricalKernelViaAnalysis_2, class SurfacePair_3 >
class Quadric_arc_2;

//! representation class for arcs on a quadric
template < class QuadricalKernelViaAnalysis_2, class SurfacePair_3 >
class Quadric_arc_2_rep : 
      public Surface_arc_2l_rep< QuadricalKernelViaAnalysis_2, SurfacePair_3 >
{

protected:
    
    //! this type's first template parameter
    typedef QuadricalKernelViaAnalysis_2 Quadrical_kernel_via_analysis_2;

    //! this instance's second template parameter
    typedef SurfacePair_3 Surface_pair_3;

    //! the instance itself
    typedef 
    Surface_arc_2l_rep< Quadrical_kernel_via_analysis_2, Surface_pair_3 >
    Self;
    
    // the base type
    typedef 
    Surface_arc_2l_rep< Quadrical_kernel_via_analysis_2, Surface_pair_3 > 
    Base;
    
protected:

    // befriending the handle
    friend class 
    Quadric_arc_2< Quadrical_kernel_via_analysis_2, Surface_pair_3 >;
};

//! represents xy-monotone arc on a quadric
template < class QuadricalKernelViaAnalysis_2, class SurfacePair_3 >
class Quadric_arc_2 :
    public CGALi::Surface_arc_2l< 
        QuadricalKernelViaAnalysis_2, 
        SurfacePair_3,
        CGALi::Surface_arc_2l_rep< QuadricalKernelViaAnalysis_2, 
        SurfacePair_3 > > {

public:

    //! this type's first template parameter
    typedef QuadricalKernelViaAnalysis_2 Quadrical_kernel_via_analysis_2;
    
    //! this instance's second template parameter
    typedef SurfacePair_3 Surface_pair_3;

    //! the class itself
    typedef Quadric_arc_2< Quadrical_kernel_via_analysis_2, Surface_pair_3 > 
    Self;
    
    //! the representation
    typedef
    CGALi::Surface_arc_2l_rep< Quadrical_kernel_via_analysis_2, 
                               Surface_pair_3 > 
    Rep;
    
    //! the base class
    typedef 
    CGALi::Surface_arc_2l< Quadrical_kernel_via_analysis_2, Surface_pair_3 > 
    Base;
    
    //! type of surface
    typedef typename Surface_pair_3::Surface_3 Surface_3;
    
    //! type of planar point
    typedef typename Base::Projected_point_2 Projected_point_2;

    //! type of planar arc
    typedef typename Base::Projected_arc_2 Projected_arc_2;
    
    //! type of surface point
    typedef typename Quadrical_kernel_via_analysis_2::Point_2 Quadric_point_2;

public:
    //!\name Simple constructors
    //!@{

    /*!\brief
     * Default constructor
     */
    Quadric_arc_2() : 
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
    Quadric_arc_2(const Projected_arc_2& arc, 
                  const Quadric_point_2& p,
                  const Quadric_point_2& q,
                  const Surface_3& surface,
                  int sheet, int sheet_p, int sheet_q) :
        Base(arc, p, q, surface, sheet, sheet_p, sheet_q) {
        CGAL_precondition(sheet < 2);
        CGAL_precondition(sheet_p < 2);
        CGAL_precondition(sheet_q < 2);
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
    Quadric_arc_2(const Projected_arc_2& arc, 
                  const Quadric_point_2& p,
                  CGAL::Arr_curve_end z_inf_end_other,
                  const Surface_3& surface,
                  int sheet, int sheet_p) :
        Base(arc, p, z_inf_end_other, surface, sheet, sheet_p) {
        CGAL_precondition(sheet < 2);
        CGAL_precondition(sheet_p < 2);
    }
    
    /*!\brief
     * Standard constructor for an unbounded arc on xy-monotone part
     * of the surface with z-asympotic behavior at both ends.
     * It represents the arc on \c surface covertical to \c arc which
     * lies on \c sheet of the xy-monotone subsurface.
     */
    Quadric_arc_2(const Projected_arc_2& arc, 
                  CGAL::Arr_curve_end z_inf_end_p,
                  CGAL::Arr_curve_end z_inf_end_q,
                  const Surface_3& surface,
                  int sheet) :
        Base(arc, z_inf_end_p, z_inf_end_q, surface, sheet) {
        CGAL_precondition(sheet < 2);
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
    Quadric_arc_2(const Projected_arc_2& arc, 
                  const Quadric_point_2& p,
                  const Surface_3& surface,
                  int sheet, int sheet_p) :
        Base(arc, p, surface, sheet, sheet_p) {
        CGAL_precondition(sheet < 2);
        CGAL_precondition(sheet_p < 2);
    }

    /*!\brief
     * Standard constructor for a branch on a xy-monotone part
     * of the surface with a z-vertical asymptotic behaviour at the projected
     * finite end.
     * It represents the arc on \c surface covertical to \c arc which
     * lies on \c sheet of the xy-monotone subsurface
     */
    Quadric_arc_2(const Projected_arc_2& arc, 
                  CGAL::Arr_curve_end z_inf_end,
                  const Surface_3& surface,
                  int sheet) :
        Base(arc, z_inf_end, surface, sheet) {
        CGAL_precondition(sheet < 2);
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
    Quadric_arc_2(const Projected_arc_2& arc, 
                  const Surface_3& surface,
                  int sheet) :
        Base(arc, surface, sheet) {
        CGAL_precondition(sheet < 2);
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
    Quadric_arc_2(const Quadric_point_2& p,
                  const Quadric_point_2& q,
                  const Surface_3& surface) :
        Base(p, q, surface) {
    }
    
    //! represents a vertical ray
    Quadric_arc_2(const Quadric_point_2& p,
                  CGAL::Arr_curve_end z_inf_end,
                  const Surface_3& surface) :
        Base(p, z_inf_end, surface) {
    }

    //! represents a vertical branch
    Quadric_arc_2(const Projected_point_2& p,
                  const Surface_3& surface) :
        Base(p, surface) {
    }
    
    //!@}

protected:
    //!\name Constructors for rebind/replace_endpoints
    //!@{
    
    /*!\brief
     * constructs an arc from a given represenation
     */
    Quadric_arc_2(Rep rep) : 
        Base(rep) { 
    }

    //!@}
    
public:
    //!\name Static members
    //!@{
    inline
    static bool can_intersect_only_at_curve_ends(const Quadric_arc_2& cv1,
                                                 const Quadric_arc_2& cv2) {
        int min_sheet_at_min = std::min(
                (cv1.is_finite(CGAL::ARR_MIN_END) ? 
                 cv1.sheet(CGAL::ARR_MIN_END) : cv1.sheet()),
                (cv2.is_finite(CGAL::ARR_MIN_END) ? 
                 cv2.sheet(CGAL::ARR_MIN_END) : cv2.sheet())
        );
        int min_sheet_at_max = std::min(
                (cv1.is_finite(CGAL::ARR_MAX_END) ? 
                 cv1.sheet(CGAL::ARR_MAX_END) : cv1.sheet()),
                (cv2.is_finite(CGAL::ARR_MAX_END) ? 
                 cv2.sheet(CGAL::ARR_MAX_END) : cv2.sheet())
        );
        
        return (cv1.sheet() != cv2.sheet() && 
                min_sheet_at_min == min_sheet_at_max && min_sheet_at_min == 0
        );
    }

    //!@}

    // friends
    //! for constructors
    friend class Quadrical_kernel_via_analysis_2::Construct_arc_2;

    //! for replace endpoints
    friend class Base::Rebind::Other;

    //! for rebind
    friend class Base::Rebind;

    friend class Quadrical_kernel_via_analysis_2::Arc_2::Base;
};   

#define CGAL_CKvA_2l_GRAB_BASE_FUNCTOR_TYPES \
    typedef typename Curved_kernel_via_analysis_2l::Curve_2 Curve_2; \
    typedef typename Curved_kernel_via_analysis_2l::Point_2 Point_2; \
    typedef typename Curved_kernel_via_analysis_2l::Arc_2 Arc_2; \


namespace Quadrical_kernel_via_analysis_2_Functors {

template < class CurvedKernelViaAnalysis_2l >
class Compare_x_on_identification_2 : public 
// Remark: Does not exist in CKvA_2
Curved_kernel_via_analysis_2_Functors::
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2l > {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;

    //! the base type
    typedef Curved_kernel_via_analysis_2_Functors::
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2l >
    Base;

    CGAL_CKvA_2l_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<2>            Arity;
    
    //! standard constructor
    Compare_x_on_identification_2(Curved_kernel_via_analysis_2l *kernel) :
        Base(kernel) {
    }
        
    /*!
     * Compare the x-coordinates of two points on the identification
     * \param p1 The first point.
     * \param p2 The second point.
     * \return LARGER if x(p1) > x(p2);
     *         SMALLER if x(p1) \< x(p2);
     *         EQUAL if x(p1) = x(p2).
     */
    result_type operator()(const Point_2 &p1, const Point_2 &p2) const {
        return Curved_kernel_via_analysis_2l::instance().
            kernel().compare_x_2_object()
            (p1.x(), p2.x());
    }
};

template < class CurvedKernelViaAnalysis_2l >
class Compare_xy_2 : 
        public CurvedKernelViaAnalysis_2l::Base::Compare_xy_2 {
    
public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;

    //! the base type
    typedef typename Curved_kernel_via_analysis_2l::Base::Compare_xy_2 Base;

    CGAL_CKvA_2l_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<2>            Arity;
    
    //! standard constructor
    Compare_xy_2(Curved_kernel_via_analysis_2l *kernel) :
        Base(kernel) {
    }
        
    /*!
     * Compare the coordinates of two points lexicographically
     * \param p1 The first point.
     * \param p2 The second point.
     * \return LARGER if x(p1) >lex x(p2);
     *         SMALLER if x(p1) \<lex x(p2);
     *         EQUAL if x(p1) = x(p2).
     */
    result_type operator()(const Point_2& p1, const Point_2& p2,
                           bool equal_x = false) const {
        
        CGAL::Comparison_result res = 
            (equal_x ? CGAL::EQUAL : 
             Curved_kernel_via_analysis_2l::instance().
             compare_x_2_object()(p1, p2)
            );
        
        if (res != CGAL::EQUAL) {
            // do nothing
        } else if (p1.sheet() != p2.sheet()) {
            res = CGAL::compare(p1.sheet(), p2.sheet());
        } else {
            Base base_compare_xy(this->_ckva());
            res = base_compare_xy(p1, p2, true);
            if (p1.sheet() == 1 && p2.sheet() == 1) {
                res = -res;
            }
        }
        return res;
    }
}; // Compare_xy_2

template < class CurvedKernelViaAnalysis_2l >
class Compare_y_near_boundary_2 : 
        public CurvedKernelViaAnalysis_2l::Base::Compare_y_near_boundary_2 {
public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;

    //! the base type
    typedef typename 
    Curved_kernel_via_analysis_2l::Base::Compare_y_near_boundary_2
    Base;

    CGAL_CKvA_2l_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<3>            Arity;
    

    //! standard constructor
    Compare_y_near_boundary_2(Curved_kernel_via_analysis_2l *kernel) :
        Base(kernel) {
    }

    /*! Compare the y-coordinates of 2 lines at their ends near the boundary
     * of the parameter space at x = +/- oo.
     * \param cv1 the first arc.
     * \param cv2 the second arc.
     * \param ce the line end indicator.
     * \return the second comparison result.
     * \pre the ce ends of the lines xcv1 and xcv2 lie either on the left
     * boundary or on the right boundary of the parameter space.
     */
    result_type operator()(const Arc_2& cv1, const Arc_2& cv2,
                           CGAL::Arr_curve_end ce) const {
        
        CERR("\nqkva_compare_y_near_boundary;\ncv1: " << cv1 << ";\ncv2: " <<
             cv2 << ";\nce: " << ce << "\n");
        
        CGAL::Comparison_result res = CGAL::EQUAL;
        
        CGAL_precondition(
                cv1.location(ce) == CGAL::ARR_LEFT_BOUNDARY ||
                cv1.location(ce) == CGAL::ARR_RIGHT_BOUNDARY
        );
        CGAL_precondition(
                cv2.location(ce) == CGAL::ARR_LEFT_BOUNDARY ||
                cv2.location(ce) == CGAL::ARR_RIGHT_BOUNDARY
        );

        int sheet1 = cv1.sheet();
        int sheet2 = cv2.sheet();

        Base base_compare_y_near_boundary(this->_ckva());
        typename Curved_kernel_via_analysis_2l::Base::Compare_y_at_x_right_2
            base_compare_y_at_x_right(this->_ckva());
        typename Curved_kernel_via_analysis_2l::Base::Compare_y_at_x_left_2
            base_compare_y_at_x_left(this->_ckva());
        
        if (sheet1 != sheet2) {
            res = CGAL::compare(sheet1, sheet2);
        } else {
            if (!cv1.is_finite(ce)) {
                res = base_compare_y_near_boundary(cv1, cv2, ce);
                if (sheet1 == 1) {
                    CGAL_assertion(sheet2 == 1);
                    res = -res;
                }
            } else {
                if (ce == CGAL::ARR_MIN_END) {
                    res = base_compare_y_at_x_right(
                            cv1, cv2, cv1.curve_end(CGAL::ARR_MIN_END)
                    );
                } else {
                    res = base_compare_y_at_x_left(
                            cv1, cv2, cv1.curve_end(CGAL::ARR_MAX_END)
                    );
                }
                // already reversed the case sheet1 == sheet2 == 1
            }
        }
        
        CERR("result: " << res << "\n");
        return res;
    }
}; // Compare_y_near_boundary_2

template < class CurvedKernelViaAnalysis_2l >
class Compare_y_at_x_2 : 
        public CurvedKernelViaAnalysis_2l::Base::Compare_y_at_x_2 {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;

    //! the base type
    typedef typename 
    Curved_kernel_via_analysis_2l::Base::Compare_y_at_x_2 Base;

    CGAL_CKvA_2l_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<2>            Arity;
    
    //! standard constructor
    Compare_y_at_x_2(Curved_kernel_via_analysis_2l *kernel) :
        Base(kernel) {
    }

    /*!
     * Return the location of the given point with respect to the input curve.
     * \param cv The curve.
     * \param p The point.
     * \pre p is in the x-range of cv.
     * \return SMALLER if y(p) \< cv(x(p)), i.e. the point is below the curve;
     *         LARGER if y(p) > cv(x(p)), i.e. the point is above the curve;
     *         EQUAL if p lies on the curve.
     */
    result_type operator()(const Point_2& p, const Arc_2& cv) const {
     
        CERR("\nqkva_compare_y_at_x;\np: " << p << ";\ncv: " << cv << "\n"); 
        CGAL::Comparison_result res = CGAL::EQUAL;
        
        // FUTURE TODO p can lie on boundary

        int sp = p.sheet();
        int sa = cv.sheet();
        if (cv.is_finite(CGAL::ARR_MIN_END) && 
            p == cv.curve_end(CGAL::ARR_MIN_END)) {
            sa = cv.sheet(CGAL::ARR_MIN_END);
        } else if (cv.is_finite(CGAL::ARR_MAX_END) && 
                   p == cv.curve_end(CGAL::ARR_MAX_END)) {
            sa = cv.sheet(CGAL::ARR_MAX_END);
        }  
        if (sa != sp) {
            res = CGAL::compare(sp, sa);
        } else {
            Base base_compare_y_at_x(this->_ckva());
            res = base_compare_y_at_x(p, cv);
            if (sa == 1) {
                CGAL_assertion(sp == 1);
                res = -res;
            }
        }
        
        CERR("result: " << res << "\n");
        return res;
    }
}; // Compare_y_at_x_2


template < class CurvedKernelViaAnalysis_2l >
class Compare_y_at_x_left_2 :
        public CurvedKernelViaAnalysis_2l::Base::Compare_y_at_x_left_2 {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;

    //! the base type
    typedef typename Curved_kernel_via_analysis_2l::Base::Compare_y_at_x_left_2
    Base;

    CGAL_CKvA_2l_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<3>            Arity;
    
    //! standard constructor
    Compare_y_at_x_left_2(Curved_kernel_via_analysis_2l *kernel) :
        Base(kernel) {
    }

    /*!
     * Compares the y value of two x-monotone curves immediately to the left
     * of their intersection point. If one of the curves is vertical
     * (emanating downward from p), it's always considered to be below the
     * other curve.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param p The intersection point.
     * \pre The point p lies on both curves, and both of them must be also be
     *      defined (lexicographically) to its left.
     * \return The relative position of cv1 with respect to cv2 immdiately to
     *         the left of p: SMALLER, LARGER or EQUAL.
     */
    result_type operator() (const Arc_2& cv1, const Arc_2& cv2,
                            const Point_2& p) const {

        CERR("\nqkva_compare_y_at_x_left(cv2);\ncv1: " << cv1 << ";\ncv2: " <<
             cv2 << ";\np: " << p << "\n");
        
        CGAL::Comparison_result res = CGAL::EQUAL;

        int sheet1 = cv1.sheet();
        int sheet2 = cv2.sheet();
        
        if (sheet1 != sheet2) {
            res = CGAL::compare(sheet1, sheet2);
        } else {
            Base base_compare_y_at_x_left(this->_ckva());
            res = base_compare_y_at_x_left(cv1, cv2, p);
            if (sheet1 == 1) {
                CGAL_assertion(sheet2 == 1);
                res = -res;
            }
        }
        
        CERR("result: " << res << "\n");
        return res;
    }
}; // Compare_y_at_x_2_left


template < class CurvedKernelViaAnalysis_2l >
class Compare_y_at_x_right_2 :
        public CurvedKernelViaAnalysis_2l::Base::Compare_y_at_x_right_2 {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;

    //! the base type
    typedef typename 
    Curved_kernel_via_analysis_2l::Base::Compare_y_at_x_right_2
    Base;

    CGAL_CKvA_2l_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<3>            Arity;
    
    //! standard constructor
    Compare_y_at_x_right_2(Curved_kernel_via_analysis_2l *kernel) :
        Base(kernel) {
    }

    /*!
     * Compares the y value of two x-monotone curves immediately to the right
     * of their intersection point. If one of the curves is vertical
     * (emanating downward from p), it's always considered to be below the
     * other curve.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param p The intersection point.
     * \pre The point p lies on both curves, and both of them must be also be
     *      defined (lexicographically) to its right.
     * \return The relative position of cv1 with respect to cv2 immdiately to
     *         the right of p: SMALLER, LARGER or EQUAL.
     */
    result_type operator() (const Arc_2& cv1, const Arc_2& cv2,
                            const Point_2& p) const {

        CERR("\nqkva_compare_y_at_x_right(cv2);\ncv1: " << cv1 
             << ";\ncv2: " << cv2 << ";\np: " << p << "\n");

        CGAL::Comparison_result res = CGAL::EQUAL;
        
        int sheet1 = cv1.sheet();
        int sheet2 = cv2.sheet();
        
        if (sheet1 != sheet2) {
            res = CGAL::compare(sheet1, sheet2);
        } else {
            Base base_compare_y_at_x_right(this->_ckva());
            res = base_compare_y_at_x_right(cv1, cv2, p);
            if (sheet1 == 1) {
                CGAL_assertion(sheet2 == 1);
                res = -res;
            }
        }
        
        CERR("result: " << res << "\n");
        return res;
    }
}; // Compare_y_at_x_2_right


template < class CurvedKernelViaAnalysis_2l >
class Do_overlap_2 :
        public CurvedKernelViaAnalysis_2l::Base::Do_overlap_2 {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;

    //! the base type
    typedef typename Curved_kernel_via_analysis_2l::Base::Do_overlap_2 Base;

    CGAL_CKvA_2l_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type 
    typedef bool result_type;
    typedef Arity_tag<2> Arity;
    
    //! standard constructor
    Do_overlap_2(Curved_kernel_via_analysis_2l *kernel) :
        Base(kernel) {
    }
    
    /*!\brief
     * Check whether two given curves overlap, i.e., they have infinitely
     * many intersection points
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the curves overlap; (false) otherwise.
     */
    bool operator()(const Arc_2& cv1, const Arc_2& cv2) const {
    
        CERR("\nqkva_do_overlap:\ncv1: " << cv1 << ";\ncv2: " << cv2 << "\n");
        
        int sheet1 = cv1.sheet();
        int sheet2 = cv2.sheet();

        bool res = (sheet1 == sheet2);
        
        if (res) {
            Base base_do_overlap(this->_ckva());
            res = base_do_overlap(cv1, cv2);
        }

        CERR("result: " << res << "\n");
        return res;
    }
}; // Do_overlap_2


//!\brief Tests two objects, whether they are equal
template < class CurvedKernelViaAnalysis_2l >
class Equal_2 :
        public CurvedKernelViaAnalysis_2l::Base::Equal_2 {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;

    //! the base type
    typedef typename Curved_kernel_via_analysis_2l::Base::Equal_2 Base;

    CGAL_CKvA_2l_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef bool result_type;
    typedef Arity_tag<2> Arity;
    
    //! standard constructor
    Equal_2(Curved_kernel_via_analysis_2l *kernel) :
        Base(kernel) {
    }
    
    /*!
     * Check if the two points are the same.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return (true) if the two point are the same; (false) otherwise.
     */
    result_type operator()(const Point_2& p1, const Point_2& p2) const {

        CERR("\nqkva_equal;\np1: " << p1 << ";\np2: " << p2 << "\n");

        bool res = (Curved_kernel_via_analysis_2l::instance().
                    compare_xy_2_object()(p1, p2) == 
                    CGAL::EQUAL);

        CERR("result: " << res << "\n");
        return res;
    }
     
    /*!
     * Check if the two x-monotone curves are the same (have the same graph).
     * \param cv1 The first 
     *        curve(_m_curved_kernel->kernel().compare_xy_2_object()
     (p1.xy(), p2.xy()));.
     * \param cv2 The second curve.
     * \return (true) if the two curves are the same; (false) otherwise.
     */
    result_type operator()(const Arc_2& cv1, const Arc_2& cv2) const {
        
        CERR("\nqkva_equal;\ncv1: " << cv1 << ";\ncv2: " << cv2 << "\n");
        
        int sheet1 = cv1.sheet();
        int sheet2 = cv2.sheet();
        
        bool res = (sheet1 == sheet2);
        
        if (res) {
            Base base_equal(this->_ckva());
            res = base_equal(cv1, cv2);
        }
        
        CERR("result: " << res << "\n");
        return res;
    }
}; // Equal_2


//! checks wether and how two arcs are intersection - with first filtering
template < class CurvedKernelViaAnalysis_2l >
class Intersect_2 :
        public CurvedKernelViaAnalysis_2l::Base::Intersect_2 {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;

    //! the base type
    typedef typename Curved_kernel_via_analysis_2l::Base::Intersect_2 Base;

    CGAL_CKvA_2l_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef std::iterator< output_iterator_tag, CGAL::Object > result_type;
    typedef Arity_tag<3> Arity;    
    
    //! standard constructor
    Intersect_2(Curved_kernel_via_analysis_2l *kernel) :
        Base(kernel) {
    }
    
    /*!
     * Find all intersections of the two given curves and insert them to the 
     * output iterator. If two arcs intersect only once, only a single will be
     * placed to the iterator. Type of output iterator is \c CGAL::Object 
     * containing either an \c Arc_2 object (overlap) or a \c Point_2 object
     * with multiplicity (point-wise intersections)
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param oi The output iterator.
     * \return The past-the-end iterator.
     */
    template < class OutputIterator >
    OutputIterator operator()(const Arc_2& cv1, const Arc_2& cv2,
                              OutputIterator oi) const {

        CERR("\nqkva_intersect;\ncv1: " << cv1 << ";\ncv2: " << cv2);

        CGAL_precondition(cv1.surface() == cv2.surface());
        int sheet1 = cv1.sheet();
        int sheet2 = cv2.sheet();

        // handle special case of two segments on same curve and at endpoints
        if ((sheet1 == sheet2 && cv1.curve().id() == cv2.curve().id() &&
             !cv1.do_overlap(cv2)) || 
            Arc_2::can_intersect_only_at_curve_ends(cv1, cv2)) {
            
            std::list< std::pair< Point_2, unsigned int > > ipoints;
            Arc_2::_intersect_at_endpoints(
                    cv1, cv2, std::back_inserter(ipoints)
            );
            
            for (typename std::list< std::pair< Point_2, unsigned int > >::
                     const_iterator it = ipoints.begin();
                 it != ipoints.end(); it++) {
                
                *oi++ = CGAL::make_object(*it);
            }

        } else if (sheet1 == sheet2) {
            
            // call projected intersection
            std::list< CGAL::Object > tmp;
            
            Base base_intersect(this->_ckva());
            base_intersect(cv1, cv2, std::back_inserter(tmp));
            
            for (std::list< CGAL::Object >::const_iterator it = tmp.begin();
                 it != tmp.end(); it++) {
                
                typedef typename Curved_kernel_via_analysis_2l::
                    Curved_kernel_via_analysis_2 Projected_kernel_2;
                
                typedef typename Projected_kernel_2::Point_2 P_point_2;
                typedef typename Projected_kernel_2::Arc_2 P_arc_2;

                typedef typename Base::Point_2 Base_point_2;
                typedef typename Base::Arc_2 Base_arc_2;

                Base_arc_2 base_arc;

                if (CGAL::assign(base_arc, *it)) {
                    
                    P_arc_2 p_arc = typename Base::Arc_2::Rebind()(base_arc);
                    
                    // lift overlapping arcs
                    int sheet_min = sheet1;
                    int sheet_max = sheet1;
                    Point_2 pt_min;
                    Point_2 pt_max;
                    
                    bool min_finite = p_arc.is_finite(CGAL::ARR_MIN_END);
                    bool max_finite = p_arc.is_finite(CGAL::ARR_MAX_END);
                    
                    if (min_finite) {
                        if (p_arc.curve_end(CGAL::ARR_MIN_END) ==
                            cv1.projected_arc().curve_end(CGAL::ARR_MIN_END)){
                            pt_min = cv1.curve_end(CGAL::ARR_MIN_END);
                            sheet_min = cv1.sheet(CGAL::ARR_MIN_END);
                        } else if (p_arc.curve_end(CGAL::ARR_MIN_END) ==
                                   cv2.projected_arc().
                                   curve_end(CGAL::ARR_MIN_END)){
                            pt_min = cv2.curve_end(CGAL::ARR_MIN_END);
                            sheet_min = cv2.sheet(CGAL::ARR_MIN_END);
                        }
                    }
                    if (max_finite) {
                        if (p_arc.curve_end(CGAL::ARR_MAX_END) ==
                            cv1.projected_arc().
                            curve_end(CGAL::ARR_MAX_END)){
                            pt_max = cv1.curve_end(CGAL::ARR_MAX_END);
                            sheet_max = cv1.sheet(CGAL::ARR_MAX_END);
                        } else if (p_arc.curve_end(CGAL::ARR_MAX_END) ==
                                   cv2.projected_arc().
                                   curve_end(CGAL::ARR_MAX_END)){
                            pt_max = cv2.curve_end(CGAL::ARR_MAX_END);
                            sheet_max = cv2.sheet(CGAL::ARR_MAX_END);
                        }
                    }
                    
                    // Remark: Current ellipsoidal implemenation does not
                    //         allow z-asymptotes. So no need to handle them
                    //         currently. Only the "usual" cases for 
                    //         surface_arcs remain.
                    Arc_2 arc;

                    if (min_finite && max_finite) {
                        arc = Arc_2(
                                p_arc, 
                                pt_min, pt_max,
                                cv1.surface(),
                                sheet1,
                                sheet_min, sheet_max
                        );
                    } else if (!min_finite && !max_finite) {
                        arc = Arc_2(p_arc, cv1.surface(), sheet1);
                    } else {
                        arc = Arc_2(p_arc, 
                                    (min_finite ? pt_min : pt_max),
                                    cv1.surface(), 
                                    sheet1,
                                    (min_finite ? sheet_min : sheet_max));
                    }
                    
                } else {
                    
                    std::pair< Base_point_2, unsigned int > b_pt;
                    CGAL_assertion_code(bool check =)
                        CGAL::assign(b_pt, *it);
                    CGAL_assertion(check);
                    
                    P_point_2 p_pt = 
                        typename Base::Point_2::Rebind()(b_pt.first);
                    
                    typename 
                        Curved_kernel_via_analysis_2l::Construct_point_on_arc_2
                        construct_point_on_arc = 
                        Curved_kernel_via_analysis_2l::instance().
                        construct_point_on_arc_2_object();
                    
                    Point_2 pt = 
                        construct_point_on_arc(
                                p_pt.x(), 
                                p_pt.curve(), 
                                p_pt.arcno(),
                                (p_pt.curve().id() == cv1.curve().id() ?
                                 cv1 : cv2)
                        );
                    
                    *oi++ = CGAL::make_object(std::make_pair(pt, b_pt.second));
                }
            }
        }
        return oi;
    }

}; // Intersect_2;


template < class CurvedKernelViaAnalysis_2l >
class Trim_2 :
        public CurvedKernelViaAnalysis_2l::Base::Trim_2 {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;

    //! the base type
    typedef typename Curved_kernel_via_analysis_2l::Base::Trim_2 Base;

    CGAL_CKvA_2l_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef Arc_2 result_type;
    typedef Arity_tag<3> Arity;
    
    //! standard constructor
    Trim_2(Curved_kernel_via_analysis_2l *kernel) :
        Base(kernel) {
    }
    
    /*!\brief
     * Returns a 
     * \param cv curve ??
     * \param p ??
     * \param q ??
     * \return (true) if the curves overlap; (false) otherwise.
     */
    Arc_2 operator()(const Arc_2& cv, const Point_2& p, const Point_2& q) {
        
        CERR("\nqkva_trim:\ncv: " << cv << ";\np: " << p 
             << ";\nq: "<< q << "\n");
        CGAL_precondition(
                !Curved_kernel_via_analysis_2l::instance().
                equal_2_object()(p, q)
        );
        CGAL_precondition(cv.compare_y_at_x(p) == CGAL::EQUAL);
        CGAL_precondition(cv.compare_y_at_x(q) == CGAL::EQUAL);  
        CGAL_precondition(p.is_finite());
        CGAL_precondition(q.is_finite());
        Arc_2 arc = cv._replace_endpoints(
                p, q, 
                (cv.is_vertical() ? -1 : cv.arcno(p.x())),
                (cv.is_vertical() ? -1 : cv.arcno(q.x())),
                (cv.is_z_vertical() ? -1 : cv.sheet(p.projected_point())),
                (cv.is_z_vertical() ? -1 : cv.sheet(q.projected_point()))
        ).first;

        CERR("result: " << arc << "\n");
        return arc;
    }
};


template < class CurvedKernelViaAnalysis_2l >
class Split_2 :
        public CurvedKernelViaAnalysis_2l::Base::Split_2 {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;

    //! the base type
    typedef typename Curved_kernel_via_analysis_2l::Base::Split_2 Base;

    CGAL_CKvA_2l_GRAB_BASE_FUNCTOR_TYPES;

    //! the result type
    typedef void result_type;
    typedef Arity_tag<4> Arity;    
    
    //! standard constructor
    Split_2(Curved_kernel_via_analysis_2l *kernel) :
        Base(kernel) {
    }
    
    /*!
     * Split a given x-monotone curve at a given point into two sub-curves.
     * \param cv The curve to split
     * \param p The split point.
     * \param c1 Output: The left resulting subcurve (p is its right endpoint)
     * \param c2 Output: The right resulting subcurve (p is its left endpoint)
     * \pre p lies on cv but is not one of its end-points.
     */
    void operator()(const Arc_2& cv, const Point_2 & p,
                    Arc_2& c1, Arc_2& c2) const {
        
        CERR("\nqkva_split:\ncv: " << cv << ";\np: " << p << "\n");
        
        CGAL_precondition(cv.compare_y_at_x(p) == CGAL::EQUAL);
        CGAL_precondition(p.is_finite());
        // check that p is not an end-point of the arc
        CGAL_precondition_code(
                cv._same_arc_compare_xy(cv._minpoint(), p) != CGAL::EQUAL &&
                cv._same_arc_compare_xy(cv._maxpoint(), p) != CGAL::EQUAL);
        
        c1 = cv._replace_endpoints(
                cv._minpoint(), p, 
                -1, (cv.is_vertical() ? -1 : cv.arcno()),
                -1, (cv.is_z_vertical() ? -1 : cv.sheet())
        ).first;
        c2 = cv._replace_endpoints(
                p, cv._maxpoint(), 
                (cv.is_vertical() ? -1 : cv.arcno()), -1,
                (cv.is_z_vertical() ? -1 : cv.sheet()), -1
        ).first;
        
        CERR("result:\nc1: " << c1 << ";\nc2: " << c2 << "\n");
    }
};

template < class CurvedKernelViaAnalysis_2l >
class Are_mergeable_2 :
        public CurvedKernelViaAnalysis_2l::Base::Are_mergeable_2 {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;

    //! the base type
    typedef typename Curved_kernel_via_analysis_2l::Base::Are_mergeable_2
    Base;

    CGAL_CKvA_2l_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef bool result_type;
    typedef Arity_tag<2> Arity;    
    
    //! standard constructor
    Are_mergeable_2(Curved_kernel_via_analysis_2l *kernel) :
        Base(kernel) {
    }
    
    /*!\brief
     * Check whether two given curves (arcs) are mergeable
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the two arcs are mergeable, i.e., they are supported
     * by the same curve and share a common endpoint; (false) otherwise.
     */
    bool operator()(const Arc_2& cv1, const Arc_2& cv2) const {
    
        CERR("\nqkva_are_mergeable:\ncv1: " << cv1 
             << ";\ncv2: " << cv2 << "\n");
        
        int sheet1 = cv1.sheet();
        int sheet2 = cv2.sheet();

        bool res = true;
        
        if (sheet1 != sheet2 && cv1.curve().id() == cv2.curve().id()) {
            res = false;
        } else if (Arc_2::can_intersect_only_at_curve_ends(cv1, cv2)) {
            res = false;
        }
        
        if (res) {
            Base base_are_mergeable(this->_ckva());
            res = base_are_mergeable(cv1, cv2);
        }
        
        CERR("result: " << res << "\n");
        return res;
    }
};


template < class CurvedKernelViaAnalysis_2l >
class Merge_2 :
        public CurvedKernelViaAnalysis_2l::Base::Merge_2 {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;

    //! the base type
    typedef typename Curved_kernel_via_analysis_2l::Base::Merge_2 Base;

    CGAL_CKvA_2l_GRAB_BASE_FUNCTOR_TYPES;

    //! the result type
    typedef void result_type;
    typedef Arity_tag<2> Arity;    
    
    //! standard constructor
    Merge_2(Curved_kernel_via_analysis_2l *kernel) :
        Base(kernel) {
    }

    /*!\brief
     * Merge two given x-monotone curves into a single one
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param c Output: The resulting curve.
     * \pre The two curves are mergeable, that is they are supported by the
     *      same curve and share a common endpoint.
     */  
    void operator()(const Arc_2& cv1, const Arc_2& cv2, Arc_2& c) const {
    
        CERR("\nqkva_merge:\ncv1: " << cv1 << ";\ncv2: " << cv2 << "\n");
        CGAL_precondition(cv1.are_mergeable(cv2));
        Arc_2::simplify(cv1, cv2);
        
        Point_2 src, tgt;
        int arcno_s = -1, arcno_t = -1;
        int sheet_s = -1, sheet_t = -1;
        bool replace_src; // true if cv2 < *this otherwise *this arc < cv2 arc
        // arcs are mergeable => they have one common finite end-point
        replace_src = (cv1._minpoint() == cv2._maxpoint());
        src = (replace_src ? cv2._minpoint() : cv1._minpoint());
        tgt = (replace_src ? cv1._maxpoint() : cv2._maxpoint());
              
        if (!cv1.is_vertical()) {
            arcno_s = (replace_src ? cv2.arcno(CGAL::ARR_MIN_END) :
                       cv1.arcno(CGAL::ARR_MIN_END));
            arcno_t = (replace_src ? cv1.arcno(CGAL::ARR_MAX_END) :
                       cv2.arcno(CGAL::ARR_MAX_END));
        }
        if (!cv1.is_z_vertical()) {
            sheet_s = (replace_src ? cv2.sheet(CGAL::ARR_MIN_END) :
                       cv1.sheet(CGAL::ARR_MIN_END));
            sheet_t = (replace_src ? cv1.sheet(CGAL::ARR_MAX_END) :
                       cv2.sheet(CGAL::ARR_MAX_END));
        }

        Arc_2 arc = cv1._replace_endpoints(src, tgt, 
                                           arcno_s, arcno_t,
                                           sheet_s, sheet_t).first;
        c = arc;

        CERR("result: " << c << "\n");
    }
};


template < class CurvedKernelViaAnalysis_2l >
class Make_x_monotone_2 :
        public CurvedKernelViaAnalysis_2l::Base::Make_x_monotone_2 {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;

    //! the base type
    typedef typename Curved_kernel_via_analysis_2l::Base::Make_x_monotone_2
    Base;

    CGAL_CKvA_2l_GRAB_BASE_FUNCTOR_TYPES;
    
    typedef typename 
    Curved_kernel_via_analysis_2l::Curved_kernel_via_analysis_2
    Curved_kernel_via_analysis_2;
    

    typedef QdX::Quadric_3_z_at_xy_isolator_traits< 
        Curved_kernel_via_analysis_2, Curve_2 // Remark: Curve_2 == Surface_3
    > Z_at_xy_traits;
    typedef QdX::Quadric_pair_3< Z_at_xy_traits > Surface_pair_3;
    
    //! the result type 
    typedef std::iterator< output_iterator_tag, CGAL::Object > result_type;
    typedef Arity_tag<2> Arity;   
    
    //! standard constructor
    Make_x_monotone_2(Curved_kernel_via_analysis_2l *kernel) :
        Base(kernel) {
    } 

    /*!
     * decompose a given arc into list of x-monotone pieces 
     * (subcurves) and insert them to the output iterator. Since \c Arc_2 
     * is by definition x-monotone, an input arc is passed to the 
     * output iterator directly. 
     * \param cv The curve.
     * \param oi The output iterator, whose value-type is Object. 
     * The returned objects are all wrappers X_monotone_curve_2 objects.
     * \return The past-the-end iterator.
     */
    template < class OutputIterator >
    OutputIterator operator()(const Arc_2& cv, OutputIterator oi) const {
    
        *oi++ = CGAL::make_object(cv);
        return oi;
    }
    
    /*!
     * decompose a given curve into list of x-monotone pieces 
     * (subcurves) and insert them to the output iterator. 
     * \param cv The curve.
     * \param oi The output iterator, whose value-type is Object. 
     * The returned objects are all wrappers X_monotone_curve_2 objects.
     * \return The past-the-end iterator.
     */
    template < class OutputIterator >
    OutputIterator operator()(const Curve_2& cv, OutputIterator oi) const {
        
        CERR("\nqkva_make_x_monotone:\ncv: " << cv << "\n");
        
        // construct surface pair 
        Surface_pair_3 pair = Surface_pair_3::surface_pair_cache()(
                std::make_pair(
                        Curved_kernel_via_analysis_2l::instance().reference(),
                        cv
                )
        );
        
        // compute lifted points and arcs
        std::list< Point_2 > points;
        std::list< Arc_2 > arcs;
        pair.template spatial_intersections< Curved_kernel_via_analysis_2l >(
                this->_ckva()->reference(),
                std::back_inserter(arcs),
                std::back_inserter(points)
        );

        // output as CGAL::Objects
        for (typename std::list< Point_2 >::const_iterator pit = 
                 points.begin();
             pit != points.end(); pit++) {
            *oi++ = CGAL::make_object(*pit);
        }

        for (typename std::list< Arc_2 >::const_iterator ait = arcs.begin();
             ait != arcs.end(); ait++) {
            *oi++ = CGAL::make_object(*ait);
        }

        return oi;
    }
}; // Make_x_monotone_2

#undef CGAL_CKvA_2l_GRAB_BASE_FUNCTOR_TYPES

} // Quadrical_kernel_via_analysis_2_functors

} // namespace CGALi

//! basic kernel to maintain points and arcs on a quadric
template < class CurvedKernelViaAnalysis_2, class SurfacePair_3 >
class Quadrical_kernel_via_analysis_2 :
  public CurvedKernelViaAnalysis_2::
  template rebind<
    Quadrical_kernel_via_analysis_2< 
      CurvedKernelViaAnalysis_2, SurfacePair_3 
    >,
    CGALi::Quadric_point_2< 
      Quadrical_kernel_via_analysis_2<
        CurvedKernelViaAnalysis_2, SurfacePair_3
      >,
      SurfacePair_3
    >,
    CGALi::Quadric_arc_2< 
      Quadrical_kernel_via_analysis_2<
        CurvedKernelViaAnalysis_2, SurfacePair_3
      >,
      SurfacePair_3
    >
  >::Other
{
public:
    //! \name public typedefs
    //!@{
    
    //! this instance's first template argument
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;
    
    //! this instance's second template parameter
    typedef SurfacePair_3 Surface_pair_3;
     
    //! myself
    typedef Quadrical_kernel_via_analysis_2< 
        Curved_kernel_via_analysis_2, Surface_pair_3 
    > 
    Self;

    //! type of curve kernel
    typedef typename 
    Curved_kernel_via_analysis_2::Curve_kernel_2 Curve_kernel_2;

    //! type of curve analysis
    typedef typename Curve_kernel_2::Curve_analysis_2 Curve_analysis_2;

    //!@}

    //!\name embedded types  for \c Arrangement_2 package
    //!@{

    //! type of surface
    typedef typename Surface_pair_3::Surface_3 Surface_3;

    //! type of curve_2
    typedef Surface_3 Curve_2;
        
    //! type of a point on generic curve
    typedef CGALi::Quadric_point_2< Self, Surface_pair_3 > Point_2; 

    //! type of an arc on generic curve
    typedef CGALi::Quadric_arc_2< Self, Surface_pair_3 > Arc_2; 

    //! type of weakly x-monotone arc for \c ArrangementTraits_2
    typedef Arc_2 X_monotone_curve_2;

    //! tag specifies which boundary functors are implemented
    typedef CGAL::Arr_all_boundary_tag Boundary_category;

    //!@}
    
public:
    //!\name base type
    //!@{

    //! the base type
    typedef typename Curved_kernel_via_analysis_2::
    template rebind< Self, Point_2, Arc_2 >::Other Base;
    
    //!@}

public:
    //! \name Constructors
    //!@{

    //! default constructor
    Quadrical_kernel_via_analysis_2() :
        Base() {
    }
    
    //! standard constructor
    Quadrical_kernel_via_analysis_2(const Surface_3& reference) :
        Base(),
        _m_reference(reference) {
    }
    
    //! construct using specific \c Curve_kernel_2 instance (for controlling)
    Quadrical_kernel_via_analysis_2(const Curve_kernel_2& kernel,
                                    const Surface_3& reference) :
        Base(kernel),
        _m_reference(reference) {
    }
    
    //!@}

    //!\name Access members
    //!@{

    //! returns the reference surface
    inline
    const Surface_3& reference() {
        return _m_reference;
    }
    //!@}

public:
    //!\name embedded constructions and predicates 
    //!@{
    
    //! type of Construct_point_2 functor
    typedef 
    CGALi::Curved_kernel_via_analysis_2l_Functors::Construct_point_2l<Self> 
    Construct_point_2;
    //! returns an instance of Construct_point_2 functor
    Construct_point_2 construct_point_2_object() const { 
        return Construct_point_2(&Self::instance());
    }
    
    //! type of Construct_projected_point_2 functor
    typedef 
    typename Curved_kernel_via_analysis_2::Construct_point_2 
    Construct_projected_point_2;
    
    //! returns an instance of Construct_projected_point_2 functor
    Construct_projected_point_2 construct_projected_point_2_object() const { 
        return 
            Curved_kernel_via_analysis_2(this->kernel()).
            construct_point_2_object();
    }

    //! type of Construct_arc_2 functor
    typedef 
    CGALi::Curved_kernel_via_analysis_2l_Functors::Construct_arc_2l<Self> 
    Construct_arc_2;
    //! returns an instance of Construct_arc_2 functor
    Construct_arc_2 construct_arc_2_object() const { 
        return Construct_arc_2(&Self::instance());
    }

    //! type of Construct_projected_arc_2 functor
    typedef 
    typename Curved_kernel_via_analysis_2::Construct_arc_2 
    Construct_projected_arc_2;
    //! returns an instance of Construct_projected_arc_2 functor
    Construct_projected_arc_2 construct_projected_arc_2_object() const { 
        return
            Curved_kernel_via_analysis_2(this->kernel()).
            construct_arc_2_object();
    }
    
// declares curved kernel functors, for each functor defines a member function
// returning an instance of this functor
#define CGAL_CKvA_2l_functor_pred(Y, Z) \
    typedef CGALi::Curved_kernel_via_analysis_2l_Functors::Y< Self > Y; \
    Y Z() const { return Y(&Self::instance()); }
    
#define CGAL_CKvA_2l_functor_cons(Y, Z) CGAL_CKvA_2l_functor_pred(Y, Z)

public:
    
    CGAL_CKvA_2l_functor_cons(Construct_point_on_arc_2,
                              construct_point_on_arc_2_object);
    
    CGAL_CKvA_2l_functor_pred(Compare_xyz_3, compare_xyz_3_object);
    
    
#undef CGAL_CKvA_2l_functor_pred
#undef CGAL_CKvA_2l_functor_cons
    
// declares curved kernel functors, for each functor defines a member function
// returning an instance of this functor
#define CGAL_QKvA_2_functor_pred(Y, Z) \
    typedef CGALi::Quadrical_kernel_via_analysis_2_Functors::Y< Self > Y; \
    Y Z() const { return Y(&Self::instance()); }

#define CGAL_QKvA_2_functor_cons(Y, Z) CGAL_QKvA_2_functor_pred(Y, Z)

    CGAL_QKvA_2_functor_pred(Compare_x_on_identification_2, 
                             compare_x_on_identification_2_object);
    
    CGAL_QKvA_2_functor_pred(Compare_xy_2, compare_xy_2_object);
    
    CGAL_QKvA_2_functor_pred(Compare_y_near_boundary_2, 
                             compare_y_near_boundary_2_object);
    CGAL_QKvA_2_functor_pred(Compare_y_at_x_2, compare_y_at_x_2_object);
    CGAL_QKvA_2_functor_pred(Compare_y_at_x_left_2, 
                             compare_y_at_x_left_2_object);
    CGAL_QKvA_2_functor_pred(Compare_y_at_x_right_2, 
                             compare_y_at_x_right_2_object);
    
    CGAL_QKvA_2_functor_pred(Equal_2, equal_2_object);

    CGAL_QKvA_2_functor_pred(Do_overlap_2, do_overlap_2_object);

    CGAL_QKvA_2_functor_cons(Intersect_2, intersect_2_object);

    CGAL_QKvA_2_functor_pred(Trim_2, trim_2_object);
    CGAL_QKvA_2_functor_pred(Split_2, split_2_object);

    CGAL_QKvA_2_functor_pred(Are_mergeable_2, are_mergeable_2_object);
    CGAL_QKvA_2_functor_pred(Merge_2, merge_2_object);

    CGAL_QKvA_2_functor_cons(Make_x_monotone_2, make_x_monotone_2_object);

    // TODO implement Is_on_2 (Pavel)
    
    //!@}
    
#undef CGAL_QKvA_2_functor_pred
#undef CGAL_QKvA_2_functor_cons

protected:
    //!\name Data members
    
    //! the reference surface
    Surface_3 _m_reference;

}; // class Quadrical_kernel_via_analysis_2

CGAL_END_NAMESPACE

#endif // CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_H
// EOF
