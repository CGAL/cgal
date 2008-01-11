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

#include <CGAL/Curved_kernel_via_analysis_2l.h>

#include <QdX/gfx_utils.h>

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
    // TODO add data
    //! double approxximation
    boost::optional< QdX::Gfx_point_3 > _m_gfx_point;
    
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
        > 
>
{
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
    
    //!\name Constructors
    //!@{
    
    /*!\brief
     * Default constructor
     */
    Quadric_point_2() : 
        Base() {   
    }

    /*!\brief
     * constructs from a given represenation
     */
    Quadric_point_2(Rep rep) :
        Base(rep) {
    }
    
protected:
    //!\brief Constructs point on \c sheet of \c surface above \c point
    //!\pre sheet >= 0
    Quadric_point_2(Quadrical_kernel_via_analysis_2 *kernel,
                     const Projected_point_2& pt, 
                     const Surface_3& surface, 
                     int sheet) :
        Base(kernel, pt, surface, sheet) {
        CGAL_precondition(sheet < 2);
    }
    
    //!@}
    
public:
    //!\name IO
    //!@{
    
    //! write represenation to \c os
    void write(std::ostream& os) const {
        os << Base(*this) << " " 
           << "Surface(" << this->surface() << ", " 
           << this->sheet() 
           << ")";
    }
    
    //!@}

    friend class Quadrical_kernel_via_analysis_2::Construct_point_2;
};

/*!\relates Quadric_point_2
 * \brief 
 * output operator
 */
template < class QuadricalKernelViaAnalysis_2, class SurfacePair_3 >
std::ostream& operator<< (
        std::ostream& os,
        const 
        Quadric_point_2<QuadricalKernelViaAnalysis_2, SurfacePair_3 >& 
        pt) {
    
    pt.write(os);
    
    return os;
}



// pre-declaration
template < class QuadricalKernelViaAnalysis_2, class SurfacePair_3 >
class Quadric_arc_2;

// TODO documentation
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
    // TODO add data
    //! gfx approx
    // TODO

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
        SurfacePair_3 > 
    > 
{

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

    //!\name Constructors
    //!@{

    /*!\brief
     * Default constructor
     */
    Quadric_arc_2() : 
        Base() {   
    }

    // TODO check what to do with this c'tor
    /*!\brief
     * constructs an arc from a given represenation
     */
    Quadric_arc_2(Rep rep) : 
        Base(rep) { 
    }

protected:
    /*!\brief
     * constructs an arc on \c sheet of surface \c surface, 
     * whose projection is \c arc with given \c source and \c target.
     *
     * \pre levels must be valid
     */
    Quadric_arc_2(Quadrical_kernel_via_analysis_2 *kernel,
                  const Projected_arc_2& arc, 
                  const Quadric_point_2& p,
                  const Quadric_point_2& q,
                  const Surface_3& surface,
                  int sheet, int sheet_p, int sheet_q) :
        Base(kernel, arc, p. q. surface, sheet, sheet_p, sheet_q) {
#if 0 // TODO check what todo with this
        if (seg.is_vertical() && level > 0) {
            this->ptr()->is_reversed_ = !this->ptr()->is_reversed_;
        }
#endif
        CGAL_precondition(sheet < 2);
        CGAL_precondition(sheet_p < 2);
        CGAL_precondition(sheet_q < 2);
    }    

    /*!\brief
     * Standard constructor for a ray on xy-monotone part
     * of the surface.
     * It represents the arc on \c surface covertical to \c arc which
     * lies on \c sheet of the xy-monotone subsurface.
     *
     * \pre arc.curve_end(MIN) = p || arc.curve_end(MAX) == p
     */
    Quadric_arc_2(Quadrical_kernel_via_analysis_2 *kernel,
                  const Projected_arc_2& arc, 
                  const Quadric_point_2& p,
                  const Surface_3& surface,
                  int sheet, int sheet_p) :
        Base(kernel, arc, p, surface, sheet, sheet_p) {
        
        CGAL_precondition(sheet < 2);
        CGAL_precondition(sheet_p < 2);
    }
    

    /*!\brief
     * Standard constructor for a branch on xy-monotone part
     * of the surface.
     * It represents the arc on \c surface covertical to \c arc which
     * lies on \c sheet of the xy-monotone subsurface.
     *
     * \pre arc.curve_end(MIN) = p || arc.curve_end(MAX) == p
     */
    Quadric_arc_2(Quadrical_kernel_via_analysis_2 *kernel,
                  const Projected_arc_2& arc, 
                  const Surface_3& surface,
                  int sheet) :
        Base(kernel, arc, surface, sheet) {
        CGAL_precondition(sheet < 2);
    }
    
    // constructors for vertical arcs
    
    //! represents a bounded vertical arc
    Quadric_arc_2(Quadrical_kernel_via_analysis_2 *kernel,
                  const Quadric_point_2& p,
                  const Quadric_point_2& q,
                  const Surface_3& surface) :
        Base(kernel, p, q, surface) {

    }

    //! represents a vertical ray
    Quadric_arc_2(Quadrical_kernel_via_analysis_2 *kernel,
                  const Quadric_point_2 p,
                  CGAL::Arr_curve_end inf_end,
                  const Surface_3& surface) :
        Base(kernel, p. inf_end, surface) {
        
    }

    //! represents a vertical branch
    Quadric_arc_2(Quadrical_kernel_via_analysis_2 *kernel,
                  const Projected_point_2& p,
                  const Surface_3& surface) :
        Base(kernel, p, surface) {
    }
    
    //!@}

};   


namespace Quadrical_kernel_via_analysis_2_Functors {

template <class CurvedKernel_2>
class Compare_x_on_identification_2
{
    typedef typename CurvedKernel_2::Point_2 Point_2;
    typedef typename CurvedKernel_2::Arc_2 Arc_2;
    
public:
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<2>            Arity;
    
    //! standard constructor
    Compare_x_on_identification_2(CurvedKernel_2 *kernel) :
        _m_kernel(kernel) {
        CGAL_assertion(kernel != NULL);
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
        return _m_kernel->kernel().compare_x_2_object()
            (p1.x(), p2.x());
    }
    
private:
    //! pointer to \c CurvedKernel_2 ?
    CurvedKernel_2 *_m_kernel;
};


template <class CurvedKernel_2>
class Compare_xy_2
{
    typedef typename CurvedKernel_2::Point_2 Point_2;
    typedef typename CurvedKernel_2::Arc_2 Arc_2;
    
public:
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<2>            Arity;
    
    //! standard constructor
    Compare_xy_2(CurvedKernel_2 *kernel) :
        _m_kernel(kernel) {
        CGAL_assertion(kernel != NULL);
    }
        
    /*!
     * Compare the coordinates of two points lexicographically
     * \param p1 The first point.
     * \param p2 The second point.
     * \return LARGER if x(p1) >lex x(p2);
     *         SMALLER if x(p1) \<lex x(p2);
     *         EQUAL if x(p1) = x(p2).
     */
    template < class Point_2_ >
    result_type operator()(const Point_2_& p1_, const Point_2_& p2_,
                           bool equal_x = false) const {
        
        if (dynamic_cast<const Point_2*>((&p1_))) {
            const Point_2& p1 = *dynamic_cast<const Point_2*>((&p1_));
            CGAL_precondition(dynamic_cast<const Point_2*>((&p2_)));
            const Point_2& p2 = *dynamic_cast<const Point_2*>((&p2_));
            
            CGAL::Comparison_result res = 
                (equal_x ? CGAL::EQUAL : 
                 _m_kernel->kernel().compare_x_2_object()(p1.x(), p2.x())
                );
            
            if (res != CGAL::EQUAL) {
                // do nothing
            } else if (p1.sheet() != p2.sheet()) {
                res = CGAL::compare(p1.sheet(), p2.sheet());
            } else {
                res = _m_kernel->kernel().compare_xy_2_object()(
                        p1.xy(), p2.xy(), true
                );
                if (p1.sheet() == 1 && p2.sheet() == 1) {
                    res = -res;
                }
            }
            return res;
        } else {
            CGAL_precondition(!dynamic_cast<const Point_2*>((&p1_)));
            CGAL_precondition(
                    dynamic_cast<const typename Point_2::Projected_point_2*>
                    ((&p1_))
            );
            const typename Point_2::Projected_point_2& p1 = 
                *dynamic_cast<const typename Point_2::Projected_point_2*>
                ((&p1_));
            CGAL_precondition(!dynamic_cast<const Point_2*>((&p2_)));
            CGAL_precondition(
                    dynamic_cast<const typename Point_2::Projected_point_2*>
                    ((&p2_))
            );
            const typename Point_2::Projected_point_2& p2 = 
                *dynamic_cast<const typename Point_2::Projected_point_2*>
                ((&p2_));
            
            CGAL::Comparison_result res = 
                (equal_x ? CGAL::EQUAL : 
                 _m_kernel->kernel().compare_x_2_object()(p1.x(), p2.x())
                );
            return res;
        }
    }
    
private:
    //! pointer to \c CurvedKernel_2 ?
    CurvedKernel_2 *_m_kernel;
}; // Compare_xy_2


template < class CurvedKernel_2 >
class Compare_y_near_boundary_2 :  
    public Curved_kernel_via_analysis_2_Functors::
           Compare_y_near_boundary_2< CurvedKernel_2 >
{
   
    typedef typename CurvedKernel_2::Point_2 Point_2;
    typedef typename CurvedKernel_2::Arc_2 Arc_2;

public:
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<3>            Arity;
    
    typedef Curved_kernel_via_analysis_2_Functors::
    Compare_y_near_boundary_2< CurvedKernel_2 > Base;

    //! standard constructor
    Compare_y_near_boundary_2(CurvedKernel_2 *) {
    }

    /*! Compare the y-coordinates of 2 lines at their ends near the boundary
     * of the parameter space at x = +/- oo.
     * \param xcv1 the first arc.
     * \param xcv2 the second arc.
     * \param ce the line end indicator.
     * \return the second comparison result.
     * \pre the ce ends of the lines xcv1 and xcv2 lie either on the left
     * boundary or on the right boundary of the parameter space.
     */
    template < class Arc_2_ >
    result_type operator()(const Arc_2_& cv1, const Arc_2_& cv2,
                           CGAL::Arr_curve_end ce) const {
        
        CERR("\nquadric_compare_y_near_boundary; cv1: " << cv1 << "; cv2: " <<
             cv2 << "; end: " << ce << "\n");
        
        CGAL::Comparison_result res = CGAL::EQUAL;
        
        CGAL_precondition(dynamic_cast<const Arc_2*>((&cv1)));
        const Arc_2& arc1 = *dynamic_cast<const Arc_2*>((&cv1));
        CGAL_precondition(dynamic_cast<const Arc_2*>((&cv2)));
        const Arc_2& arc2 = *dynamic_cast<const Arc_2*>((&cv2));
        
        CGAL_precondition(
                arc1.location(ce) == CGAL::ARR_LEFT_BOUNDARY ||
                arc1.location(ce) == CGAL::ARR_RIGHT_BOUNDARY
        );
        CGAL_precondition(
                arc2.location(ce) == CGAL::ARR_LEFT_BOUNDARY ||
                arc2.location(ce) == CGAL::ARR_RIGHT_BOUNDARY
        );

        int s1 = cv1.sheet();
        int s2 = cv2.sheet();
        
        if (s1 != s2) {
            res = CGAL::compare(s1, s2);
        } else {
            bool unbounded_end = false;
            if (ce == CGAL::ARR_MIN_END) {
                unbounded_end = (arc1._minpoint().arc_rep() != NULL);
            } else {
                unbounded_end = (arc1._maxpoint().arc_rep() != NULL);
            }
            if (unbounded_end) {
                res = Base::operator()(arc1, arc2, ce);
                if (s1 == 1) {
                    CGAL_assertion(s2 == 1);
                    res = -res;
                }
            } else {
                if (ce == CGAL::ARR_MIN_END) {
                    res = 
                        Base::_m_curved_kernel->compare_y_at_x_right_2_object()
                        (arc1, arc2, arc1._minpoint());
                } else {
                    res = 
                        Base::_m_curved_kernel->compare_y_at_x_left_2_object()
                        (arc1, arc2, arc1._maxpoint());
                }
                // already reversed the case s1 == s2 == 1
            }
        }

        CERR("result: " << res << "\n");
        return res;
    }
}; // Compare_y_near_boundary_2


template < class CurvedKernel_2 >
class Compare_y_at_x_2 :  
    public Curved_kernel_via_analysis_2_Functors::
           Compare_y_at_x_2< CurvedKernel_2 >
{
   
    typedef typename CurvedKernel_2::Point_2 Point_2;
    typedef typename CurvedKernel_2::Arc_2 Arc_2;

public:
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<2>            Arity;
    
    typedef Curved_kernel_via_analysis_2_Functors::
    Compare_y_at_x_2< CurvedKernel_2 > Base;

    //! standard constructor
    Compare_y_at_x_2(CurvedKernel_2 *) {
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
    template < class Point_2_, class Arc_2_ >
    result_type operator()(const Point_2_& p, const Arc_2_& cv) const {
     
        CERR("\ncompare_y_at_x; p: " << p << ";\n cv:" << cv << "\n"); 
        CGAL::Comparison_result res = CGAL::EQUAL;
        
        CGAL_precondition(dynamic_cast<const Point_2*>((&p)));
        const Point_2& pt = *dynamic_cast<const Point_2*>((&p));
        CGAL_precondition(dynamic_cast<const Arc_2*>((&cv)));
        const Arc_2& arc = *dynamic_cast<const Arc_2*>((&cv));
        
        // FUTURE TODO p can lie on boundary

        int sp = pt.sheet();
        int sa = arc.sheet();
        
        if (sa != sp) {
            res = CGAL::compare(sp, sa);
        } else {
            res = Base::operator()(p, cv);
            if (sa == 1) {
                CGAL_assertion(sp == 1);
                res = -res;
            }
        }
        
        CERR("result: " << res << "\n");
        return res;
    }
}; // Compare_y_at_x_2

// TODO missing functors
// knows_y/reverse_y
// - Compare_y_at_x_left_2
// - Compare_y_at_x_right_2

// may_have_common_part/intersect_only_at_ends
// - Do_overlap_2
// - Equal_2
// - Are_mergeable_2
// - (Merge_2 - rewrite sheets?)

// ???
// - (Is_bounded_2)
// - (Is_on_2)

//! checks wether and how two arcs are intersection - with first filtering
template < class CurvedKernel_2 >
class Intersect_2 : 
        public Curved_kernel_via_analysis_2_Functors::
            Intersect_2< CurvedKernel_2 > {

    typedef typename CurvedKernel_2::Point_2 Point_2;

public:
    typedef typename 
    Curved_kernel_via_analysis_2_Functors::Intersect_2< CurvedKernel_2 > Base;

    typedef std::iterator<output_iterator_tag, CGAL::Object> result_type;
    typedef Arity_tag<3> Arity;    
    
    //! standard constructor
    Intersect_2(CurvedKernel_2 *kernel) :
        Base(kernel) {
        CGAL_assertion(kernel != NULL);
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
    template < class Arc_2_, class OutputIterator >
    OutputIterator operator()(const Arc_2_& cv1, const Arc_2_& cv2,
                              OutputIterator oi) const {

        CERR("\nquadric_2_intersect; cv1: " << cv1 
             << ";\n cv2:" << cv2 << "");
        
        // TODO check whether possible to intersect

        // call projected intersection
        std::list< CGAL::Object > tmp;
        Base::operator()(cv1, cv2, std::back_inserter(tmp));
        for (std::list< CGAL::Object >::const_iterator it = tmp.begin();
             it != tmp.end(); it++) {
            // TODO lift objects!
            // *oi++ = *it;
        }
        return oi;
    }
}; // Intersect_2;


template < class CurvedKernel_2>
class Make_x_monotone_2 
{
    typedef typename CurvedKernel_2::Curve_2 Curve_2;
    typedef typename CurvedKernel_2::Point_2 Point_2;
    typedef typename CurvedKernel_2::Arc_2 Arc_2;
   
public:
    typedef std::iterator<output_iterator_tag, CGAL::Object> result_type;
    typedef Arity_tag<2> Arity;   
    
    //! standard constructor
    Make_x_monotone_2(CurvedKernel_2 *kernel) :
        _m_curved_kernel(kernel) {
        CGAL_assertion(kernel != NULL);
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
    template<class OutputIterator>
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
    // TODO: move this to separate file Arr_kernel_traits_2.h ?
    template<class OutputIterator>
    OutputIterator operator()(const Curve_2& cv, OutputIterator oi) const {
    
        // TODO compute surface pair 
        
        // lift segments

        // TODO maybe adapt it to use QdX::P_curve_2
        return oi;
    }
    
protected:
    //! pointer to \c CurvedKernel_2 ?
    CurvedKernel_2 *_m_curved_kernel;
}; // Make_x_monotone_2


} // Quadrical_kernel_via_analysis_2_functors

} // namespace CGALi

//! basic kernel to maintain points and arcs on a quadric
template < class CurveKernel_2, class SurfacePair_3 >
class Quadrical_kernel_via_analysis_2 :
  public CGALi::Curved_kernel_via_analysis_2_base < CurveKernel_2 >,
  public CGALi::Curved_kernel_via_analysis_2_functors < 
    Quadrical_kernel_via_analysis_2< CurveKernel_2, SurfacePair_3 >,
     typename CurveKernel_2::Curve_2,
    CGALi::Quadric_point_2< 
      Quadrical_kernel_via_analysis_2< CurveKernel_2, SurfacePair_3 >,
      SurfacePair_3
    >,
    CGALi::Quadric_arc_2< 
      Quadrical_kernel_via_analysis_2< CurveKernel_2, SurfacePair_3 >,
      SurfacePair_3
    >
  > 
{
public:
    //! \name public typedefs
    //!@{
    
    //! this instance's first template argument
    typedef CurveKernel_2 Curve_kernel_2;

   //! this instance's second template parameter
    typedef SurfacePair_3 Surface_pair_3;
     

    //! myself
    typedef Quadrical_kernel_via_analysis_2< Curve_kernel_2, Surface_pair_3 > 
    Self;
    
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
    //!\name embedded constructions and predicates 
    //!@{
    
    typedef 
    CGALi::Curved_kernel_via_analysis_2l_Functors::Construct_point_2l<Self> 
    Construct_point_2;

    Construct_point_2 construct_point_2_object() const { 
        return Construct_point_2(
                (Quadrical_kernel_via_analysis_2 *)this
        ); 
    }

    typedef 
    CGALi::Curved_kernel_via_analysis_2_Functors::Construct_point_2<Self,
       typename Point_2::Projected_point_2 > 
    Construct_projected_point_2;
    
    Construct_projected_point_2 construct_projected_point_2_object() const { 
        return Construct_projected_point_2(
                (Quadrical_kernel_via_analysis_2 *)this
        ); 
    }

    typedef 
    CGALi::Curved_kernel_via_analysis_2l_Functors::Construct_arc_2l<Self> 
    Construct_arc_2;

    Construct_arc_2 construct_arc_2_object() const { 
        return Construct_arc_2(
                (Quadrical_kernel_via_analysis_2 *)this
        ); 
    }

// declares curved kernel functors, for each functor defines a member function
// returning an instance of this functor
#define CGAL_CKvA_2l_functor_pred(Y, Z) \
    typedef CGALi::Curved_kernel_via_analysis_2l_Functors::Y<Self> Y; \
    Y Z() const { return Y((Quadrical_kernel_via_analysis_2 *)this); }
    
#define CGAL_CKvA_2l_functor_cons(Y, Z) CGAL_CKvA_2l_functor_pred(Y, Z)

public:
    
    CGAL_CKvA_2l_functor_cons(Construct_point_on_arc_2,
                              construct_point_on_arc_2_object);
    
#undef CGAL_CKvA_2l_functor_pred
#undef CGAL_CKvA_2l_functor_cons
    
// declares curved kernel functors, for each functor defines a member function
// returning an instance of this functor
#define CGAL_QKvA_2_functor_pred(Y, Z) \
    typedef CGALi::Quadrical_kernel_via_analysis_2_Functors::Y<Self> Y; \
    Y Z() const { return Y((Quadrical_kernel_via_analysis_2 *)this); }

#define CGAL_QKvA_2_functor_cons(Y, Z) CGAL_QKvA_2_functor_pred(Y, Z)

    CGAL_QKvA_2_functor_pred(Compare_x_on_identification_2, 
                             compare_x_on_identification_2_object);
    
    CGAL_QKvA_2_functor_pred(Compare_xy_2, compare_xy_2_object);

    CGAL_QKvA_2_functor_cons(Intersect_2, intersect_2_object);

    CGAL_QKvA_2_functor_cons(Make_x_monotone_2, make_x_monotone_2_object);
    
    //!@}
    
#undef CGAL_QKvA_2_functor_pred
#undef CGAL_QKvA_2_functor_cons
   
protected:
    //!\name Protected internal types

    //!@{
    //! class collecting basic types
    typedef CGALi::Curved_kernel_via_analysis_2_base < CurveKernel_2 >
    Base_kernel;

    //! class collecting basic types
    typedef CGALi::Curved_kernel_via_analysis_2_functors < 
            Self, Curve_2, Point_2, Arc_2
    >  
    Base_functors;
    
    //!@}

public:
    //! \name Constructors
    //!@{

    //! default constructor
    Quadrical_kernel_via_analysis_2() :
        Base_kernel() {
    }
    
    //! standard constructor
    Quadrical_kernel_via_analysis_2(const Surface_3& reference) :
        Base_kernel(),
        _m_reference(reference) {
    }
    
    //! construct using specific \c Curve_kernel_2 instance (for controlling)
    Quadrical_kernel_via_analysis_2(const Curve_kernel_2& kernel,
                                    const Surface_3& reference) :
        Base_kernel(kernel),
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
 
protected:
    //!\name Data members
    Surface_3 _m_reference;

}; // class Quadrical_kernel_via_analysis_2

CGAL_END_NAMESPACE

#endif // CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_H
// EOF
