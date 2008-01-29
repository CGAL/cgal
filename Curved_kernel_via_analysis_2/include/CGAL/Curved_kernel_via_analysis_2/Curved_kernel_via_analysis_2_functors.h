
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//                 Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#ifndef CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_FUNCTORS_H
#define CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_FUNCTORS_H

/*! \file Curved_kernel_via_analysis_2/curved_kernel_via_analysis_2_functors.h
 *  \brief defines Curved_kernel_via_analysis_2 function objects + class
 */

#include <CGAL/Curved_kernel_via_analysis_2/Make_x_monotone_2.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

#ifndef CERR
//#define CKvA_DEBUG_PRINT_CERR
#ifdef CKvA_DEBUG_PRINT_CERR
#define CERR(x) std::cout << x
#else
#define CERR(x) static_cast<void>(0)
#endif
#endif

namespace Curved_kernel_via_analysis_2_Functors {

//! base functor class
template < class CurvedKernelViaAnalysis_2 > 
class Curved_kernel_via_analysis_2_functor_base {

public:
    //!\name Public types

    //!@{
    //! this instance's template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;
    
    //! the curve type
    typedef typename Curved_kernel_via_analysis_2::Curve_2 Curve_2;

    //! the point type
    typedef typename Curved_kernel_via_analysis_2::Point_2 Point_2;

    //! the arc type
    typedef typename Curved_kernel_via_analysis_2::Arc_2 Arc_2;
    
    //! the x-coordinate type
    typedef typename Point_2::X_coordinate_1 X_coordinate_1;

    //!@}

    //!\name Constructors
    //!@{

    Curved_kernel_via_analysis_2_functor_base(
            Curved_kernel_via_analysis_2 *kernel) :
        _m_curved_kernel(kernel) {
        CGAL_precondition(kernel != NULL);
    }

    //!@}

protected:
    //!\name Access members
    //!@{
    
    //! return pointer to curved kernel
    Curved_kernel_via_analysis_2* _ckva() const {
        return _m_curved_kernel;
    }

    //!@}
    
protected:
    //!\name Data members
    //!@{
    
    //! pointer to \c Curved_kernel_via_analysis_2 ?
    Curved_kernel_via_analysis_2 *_m_curved_kernel;
    
    //!@}
};


#define CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES \
    typedef typename Base::Curve_2 Curve_2; \
    typedef typename Base::Point_2 Point_2; \
    typedef typename Base::Arc_2 Arc_2; \
    typedef typename Base::X_coordinate_1 X_coordinate_1; \

// end define



template < class CurvedKernelViaAnalysis_2 >
class Construct_point_2 : public 
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef 
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;
    
    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef Point_2 result_type;

    //! standard constructor
    Construct_point_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }

    //!\brief constructs a finite point with x-coordinate
    //! \c x on curve \c c with arc number \c arcno
    //!
    //! implies no boundary conditions in x/y
    Point_2 operator()(const X_coordinate_1& x, const Curve_2& c, int arcno) {
        Point_2 pt(x, c, arcno);
        return pt;
    }
};


//!\brief Functor to construct point on an arc
//! \c x on curve \c c with arc number \c arcno
//!
//! implies no boundary conditions in x/y
template < class CurvedKernelViaAnalysis_2 >
class Construct_point_on_arc_2 : public 
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef 
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;
    
    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef Point_2 result_type;
    
    //! standard constructor
    Construct_point_on_arc_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }
    
    //! constructs points at x 
    Point_2 operator()(
            const X_coordinate_1& x,
            const Curve_2& c, int arcno,
            const Arc_2& arc) {
        CGAL_assertion(c.id() == arc.curve().id());
        CGAL_assertion(arcno == arc.arcno(x));

        typename Curved_kernel_via_analysis_2::Construct_point_2 
            construct_point = Curved_kernel_via_analysis_2::instance().
            construct_point_2_object();
        
        Point_2 pt = construct_point(x, c, arcno);
        
        // here we can modify the point, if we want to
        return pt;
    }
};

 
template < class CurvedKernelViaAnalysis_2 >
class Construct_arc_2 : public 
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef 
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;
    
    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef Arc_2 result_type;

    //! standard constructor
    Construct_arc_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }

    //! constructing non-vertical arcs
    //!@{
    
    //! \brief 
    //! constructs an arc with two finite end-points, supported by curve \c c
    //! with \c arcno (segment)  
    //! 
    //! \c arcno_p and \c arcno_q define arcnos of \c p and \c q w.r.t. 
    //! the curve \c c
    //!
    //! \pre p.x() != q.x()
    Arc_2 operator()(const Point_2& p, const Point_2& q, const Curve_2& c,
                     int arcno, int arcno_p, int arcno_q) {
        Arc_2 arc(p, q, c, arcno, arcno_p, arcno_q);
        return arc;
    }
    
    /*!\brief
     * constructs an arc with one finite end-point \c origin and one
     * x-infinite end, supported by curve \c c with \c arcno (ray I)
     *
     * \c inf_end defines whether the ray emanates from +/- x-infinity, 
     * \c arcno_o defines an arcno of point \c origin w.r.t. curve \c c
     */
    Arc_2 operator()(const Point_2& origin, CGAL::Arr_curve_end inf_end, 
                     const Curve_2& c, int arcno, int arcno_o) {
        Arc_2 arc(origin, inf_end, c, arcno, arcno_o);
        return arc;
    }

    /*!\brief
     * constructs an arc with one finite end-point \c origin and one asymtpotic
     * (y-infinite) end given by x-coordinate \c asympt_x (ray II)
     *
     * \c inf_end specifies +/-oo an asymptotic end is approaching, \c arcno_o
     * defines an arcno of point \c origin (arcno of asymptotic end is the
     * same as \c arcno )
     * \pre origin.x() != asympt_x
     */
    Arc_2 operator()(const Point_2& origin, const X_coordinate_1& asympt_x, 
                     CGAL::Arr_curve_end inf_end, const Curve_2& c, int arcno, 
                     int arcno_o) {
        Arc_2 arc(origin, asympt_x, inf_end, c, arcno, arcno_o);
        return arc;
    }

    /*!\brief
     * constructs an arc with two x-infinite ends supported by curve \c c
     * with \c arcno (branch I)
     */
    Arc_2 operator()(const Curve_2& c, int arcno) {
        Arc_2 arc(c, arcno);
        return arc;
    }
    
    /*!\brief
     * constructs an arc with two asymptotic ends defined by \c asympt_x1 and
     * \c asympt_x2 respectively, supported by curve \c c with \c arcno
     * (branch II)
     *
     * \c inf_end1/2 define +/-oo the repspective asymptotic end is approaching
     * \pre asympt_x1 != asympt_x2
     */
    Arc_2 operator()(const X_coordinate_1& asympt_x1, 
                     const X_coordinate_1& asympt_x2, 
                     CGAL::Arr_curve_end inf_end1, 
                     CGAL::Arr_curve_end inf_end2,
                     const Curve_2& c, int arcno) {
        Arc_2 arc(asympt_x1, asympt_x2, inf_end1, inf_end2, c, arcno);
        return arc;
    }
    
    /*!\brief
     * constructs an arc with one x-infinite end and one asymptotic end 
     * defined by x-coordinate \c asympt_x supported by curve \c c with 
     * \c arcno (branch III)
     *
     * \c inf_endx specifies whether the branch goes to +/- x-infinity,
     * \c inf_endy specifies +/-oo the asymptotic end approaches
     */
    Arc_2 operator()(CGAL::Arr_curve_end inf_endx, 
                     const X_coordinate_1& asympt_x,
                     CGAL::Arr_curve_end inf_endy, 
                     const Curve_2& c, int arcno) {
        Arc_2 arc(inf_endx, asympt_x, inf_endy, c, arcno);
        return arc;
    }
    
    //!@}
    
    //!\name constructing vertical arcs
    //!@{
    
    //! \brief 
    //! constructs a vertcial arc with two finite end-points \c p and \c q ,
    //! supported by curve \c c (vertical segment)
    //! 
    //! \pre p != q && p.x() == q.x()
    //! \pre c must have a vertical component at this x
    Arc_2 operator()(const Point_2& p, const Point_2& q, const Curve_2& c) {
        Arc_2 arc(p,q,c);
        return arc;
    }
    
    /*!\brief
     * constructs a vertical arc with one finite end-point \c origin and one
     * y-infinite end, supported by curve \c c (vertical ray)
     *
     * \c inf_end defines whether the ray emanates from +/- y-infninty, 
     * \pre c must have a vertical line component at this x
     */
    Arc_2 operator()(const Point_2& origin, CGAL::Arr_curve_end inf_end,
                     const Curve_2& c) {
        
        Arc_2 arc(origin, inf_end, c);
        return arc;
    }
    
    /*!\brief
     * constructs a vertical arc with two y-infinite ends, at x-coordinate 
     * \c x , supported by curve \c c (vertical branch)
     * 
     * \pre c must have a vertical line component at this x
     */
    Arc_2 operator()(const X_coordinate_1& x, const Curve_2& c) {
        Arc_2 arc(x, c);
        return arc;
    }
};


template < class CurvedKernelViaAnalysis_2 >
class Is_vertical_2 : public 
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef 
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;
    
    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;

    //! the result type
    typedef bool result_type;
    typedef Arity_tag<1> Arity;
    
    //! standard constructor
    Is_vertical_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }

    /*!
     * Check whether the given x-monotone curve is a vertical segment.
     * \param cv The curve.
     * \return (true) if the curve is a vertical segment; 
     * (false) otherwise.
     */
    result_type operator()(const Arc_2& cv) const {
        return cv.is_vertical();
    }
};


template < class CurvedKernelViaAnalysis_2 >
class Is_bounded_2 : public 
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef 
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;
    
    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef bool result_type;
    typedef Arity_tag<2> Arity;

    Is_bounded_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }
    
    /*! Is the end of an x-monotone curve bounded?
     * \param xcv The x-monotone curve.
     * \param ce The end of xcv identifier.
     * \return true is the curve end is bounded, and false otherwise
     */
    result_type operator()(const Arc_2& cv, Arr_curve_end ce) const {
        return (cv.is_finite(ce));
    }
};


template < class CurvedKernelViaAnalysis_2 >
class Parameter_space_in_x_2 : public 
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef 
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;
    
    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef CGAL::Arr_parameter_space result_type;
    typedef Arity_tag<2> Arity;
    
    //! standard constructor
    Parameter_space_in_x_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }
    
    /*! Obtains the parameter space at the end of a line along the x-axis.
     * \param xcv the line
     * \param ce the line end indicator:
     *     ARR_MIN_END - the minimal end of xc or
     *     ARR_MAX_END - the maximal end of xc
     * \return the parameter space at the ce end of the line xcv.
     *   ARR_LEFT_BOUNDARY  - the line approaches the identification arc from
     *                        the right at the line left end.
     *   ARR_INTERIOR       - the line does not approache the 
     *                        identification arc.
     *   ARR_RIGHT_BOUNDARY - the line approaches the identification arc from
     *                        the left at the line right end.
     */
    result_type operator()(const Arc_2& cv, CGAL::Arr_curve_end end) const {

        CGAL::Arr_parameter_space loc = cv.location(end);
        if(loc == CGAL::ARR_LEFT_BOUNDARY || loc == CGAL::ARR_RIGHT_BOUNDARY)
            return loc;
        return CGAL::ARR_INTERIOR;
    }
};


template < class CurvedKernelViaAnalysis_2 >
class Parameter_space_in_y_2 : public 
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {
   
public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef 
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;
    
    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef CGAL::Arr_parameter_space result_type;
    typedef Arity_tag<2> Arity;
    
    //! standard constructor
    Parameter_space_in_y_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }

    /*! Obtains the parameter space at the end of a line along the y-axis .
     * Note that if the line end coincides with a pole, then unless the line
     * coincides with the identification arc, the line end is considered to
     * be approaching the boundary, but not on the boundary.
     * If the line coincides with the identification arc, it is assumed to
     * be smaller than any other object.
     * \param xcv the line
     * \param ce the line end indicator:
     *     ARR_MIN_END - the minimal end of xc or
     *     ARR_MAX_END - the maximal end of xc
     * \return the parameter space at the ce end of the line xcv.
     *   ARR_BOTTOM_BOUNDARY  - the line approaches the south pole at the line
     *                          left end.
     *   ARR_INTERIOR         - the line does not approach boundary
     *   ARR_TOP_BOUNDARY     - the line approaches the north pole at the line
     *                          right end.
     */
    result_type operator()(const Arc_2& cv, CGAL::Arr_curve_end end) const {
            
        CGAL::Arr_parameter_space loc = cv.location(end);
        if(loc == CGAL::ARR_BOTTOM_BOUNDARY || loc == CGAL::ARR_TOP_BOUNDARY)
            return loc;
        return CGAL::ARR_INTERIOR;
    }
};


template < class CurvedKernelViaAnalysis_2 >
class Construct_min_vertex_2 : public 
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {
   
public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef 
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;
    
    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef Point_2 result_type;
    typedef Arity_tag<1> Arity;
    
    //! standard constructor
    Construct_min_vertex_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }

    /*!
     * Get the left end-point of the x-monotone curve
     * \param cv The curve.
     * \pre corresponding end-point must not lie at infinity
     * \return The left endpoint.
     */
    result_type operator()(const Arc_2& cv) const {
    
        return cv.curve_end(CGAL::ARR_MIN_END);
    }
};


template < class CurvedKernelViaAnalysis_2 >
class Construct_max_vertex_2 : public 
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {
   
public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef 
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;
    
    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef Point_2 result_type;
    typedef Arity_tag<1> Arity;
    
    //! standard constructor
    Construct_max_vertex_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }

    /*!
     * Get the right endpoint of the x-monotone curve (segment).
     * \param cv The curve.
     * \pre corresponding end-point must not lie at infinity
     * \return The right endpoint.
     */
    result_type operator()(const Arc_2& cv) const {
        
        return cv.curve_end(CGAL::ARR_MAX_END);
    }
};


template <class CurvedKernelViaAnalysis_2>
class Compare_x_2 : public 
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {
    
public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef 
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<2>            Arity;
    
    //! standard constructor
    Compare_x_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }
    
    /*!
     * Compare the x-coordinates of two points.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return LARGER if x(p1) > x(p2);
     *         SMALLER if x(p1) \< x(p2);
     *         EQUAL if x(p1) = x(p2).
     */
    result_type operator()(const Point_2 &p1, const Point_2 &p2) const {
        return Curved_kernel_via_analysis_2::instance().
            kernel().compare_x_2_object()
            (p1.x(), p2.x());
    }
};


template < class CurvedKernelViaAnalysis_2 >
class Compare_xy_2 : public 
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef 
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<2>            Arity;
    
    //! standard constructor
    Compare_xy_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }

    /*!
     * Compares two points lexigoraphically: by x, then by y.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return LARGER if x(p1) > x(p2), or if x(p1) = x(p2) and y(p1) > y(p2);
     *         SMALLER if x(p1) \< x(p2), or if x(p1) = x(p2) and 
     *                   y(p1) \< y(p2);
     *         EQUAL if the two points are equal.
     */
    result_type operator()(const Point_2& p1, const Point_2& p2,
                           bool equal_x = false) const {

        result_type res =
            (Curved_kernel_via_analysis_2::instance().
             kernel().compare_xy_2_object()
             (p1.xy(), p2.xy(), equal_x));
        return res;
    }
};


template <class CurvedKernelViaAnalysis_2>
class Compare_x_near_boundary_2 : public 
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef 
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<4>            Arity;

    Compare_x_near_boundary_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }

    /*! Compare the x-coordinate of a point with the x-coordinate of
     * a line end near the boundary at y = +/- oo.
     * \param p the point direction.
     * \param xcv the line, the endpoint of which is compared.
     * \param ce the line-end indicator -
     *            ARR_MIN_END - the minimal end of xc or
     *            ARR_MAX_END - the maximal end of xc.
     * \return the comparison result:
     *         SMALLER - x(p) \< x(xc, ce);
     *         EQUAL   - x(p) = x(xc, ce);
     *         LARGER  - x(p) > x(xc, ce).     
     * \pre p lies in the interior of the parameter space.
     * \pre the ce end of the line xcv lies on a boundary.
     */
    result_type operator()(const Point_2& p, const Arc_2& cv,
                           CGAL::Arr_curve_end ce) const {

        CERR("\ncompare_x_near_boundary: p: " << p << "\n cv: " <<
             cv << "; curve_end: " << ce << "\n");
        
        // this curve end has boundary only in y
        CGAL_precondition(cv.is_on_bottom_top(cv.location(ce)));
        if (cv.is_singular()) // the curve end goes to singularity => x-order
            return CGAL::EQUAL; // doesn't matter   

        CGAL::Comparison_result res = 
            Curved_kernel_via_analysis_2::instance().
            kernel().compare_x_2_object()(
                    p.x(), cv.curve_end_x(ce)
            );
        // for vertical arcs equality of x-coordinates means overlapping
        // in case of discontinuity => x-comparison is enough
        if (res != CGAL::EQUAL || cv.is_vertical() || cv.is_on_disc()) {
            CERR("result: " << res << "\n");
            return res;
        }
        // look at the side from which the vertical asymptote is approached 
        res = (ce == CGAL::ARR_MIN_END ? CGAL::SMALLER : CGAL::LARGER);
        CERR("result: " << res << "\n");
        return res;
    }

    /*! Compare the x-coordinates of 2 arcs ends near the boundary of the
     * parameter space at y = +/- oo.
     * \param xcv1 the first arc.
     * \param ce1 the first arc end indicator -
     *            ARR_MIN_END - the minimal end of xcv1 or
     *            ARR_MAX_END - the maximal end of xcv1.
     * \param xcv2 the second arc.
     * \param ce2 the second arc end indicator -
     *            ARR_MIN_END - the minimal end of xcv2 or
     *            ARR_MAX_END - the maximal end of xcv2.
     * \return the second comparison result:
     *         SMALLER - x(xcv1, ce1) \< x(xcv2, ce2);
     *         EQUAL   - x(xcv1, ce1) = x(xcv2, ce2);
     *         LARGER  - x(xcv1, ce1) > x(xcv2, ce2).
     * \pre the ce1 end of the line xcv1 lies on a boundary.
     * \pre the ce2 end of the line xcv2 lies on a boundary.
     */
    result_type operator()(const Arc_2& cv1, CGAL::Arr_curve_end ce1,
                           const Arc_2& cv2, CGAL::Arr_curve_end ce2) const {

        CERR("\ncompare_x_near_boundary: cv1: " << cv1 << "\n cv2: " <<
            cv2 << "; end1: " << ce1 << "; end2: " << ce2 << "\n");
        /*CGAL::Arr_boundary_type bnd1 = boundary(end1), 
            bnd2 = cv2.boundary(ce2);*/
        CGAL::Arr_parameter_space loc1 = cv1.location(ce1), 
            loc2 = cv2.location(ce2);
        CGAL_precondition(cv1.is_on_bottom_top(loc1) && 
                          cv1.is_on_bottom_top(loc2));
        
        if (cv1.is_singular() != cv1.is_singular()) {
            // only one curve end lies at singularity (another at +/-oo)
            CGAL_error_msg("SINGULARITY + INF comparison is not yet \
                implemented");
        }
        
        
        CGAL::Comparison_result res;
        if (cv1.is_singular() && cv1.is_singular()) {
            if (loc1 < loc2) {
                return CGAL::SMALLER;
            }
            if (loc1 > loc2) {
                return CGAL::LARGER;
            }
            // both ends lie at the same singularity => need special handling
            // but x-order doesn't matter
        } else  { // establish x-order
            res = Curved_kernel_via_analysis_2::instance().
                kernel().compare_x_2_object()(
                    cv1.curve_end_x(ce1),
                    cv2.curve_end_x(ce2)
            );
            // x-coordinate comparison is enough for these cases
            // we assume that either both curve ends lie on disc or neither of
            // them
            if (res != CGAL::EQUAL || 
                (cv1.is_on_disc() && cv1.is_on_disc())) {
                CERR("result: " << res << "\n");
                return res;
            }
        }    
        // now we either +/-oo case: ARR_MIN_END > vertical > ARR_MAX_END
        // or both ends lie at the same singularity: these cases can be 
        // handled simultaneously  
        if (cv1.is_vertical()) {
            if (!cv2.is_vertical()) {
                return (ce2 == CGAL::ARR_MIN_END ? 
                        CGAL::SMALLER : CGAL::LARGER);
            }
            // both are vertical
            if (loc1 == loc2) { // both ends converge to the same infinity
                return CGAL::EQUAL;
            }
            return (loc1 == CGAL::ARR_BOTTOM_BOUNDARY ? 
                    CGAL::SMALLER : CGAL::LARGER);
        } 
        
        if (cv2.is_vertical()) {
            return (ce1 == CGAL::ARR_MIN_END ? CGAL::LARGER : CGAL::SMALLER);
        }
        
        // otherwise: both ends have asymptotic behaviour or singularity
        if (ce1 == ce2) { // both ends approach asymptote from one side
            
            if(loc1 == loc2) { // need special y-comparison
                X_coordinate_1 x0(cv1.curve_end_x(ce1));
                res = cv1._compare_arc_numbers(
                        cv2, CGAL::ARR_INTERIOR, x0, 
                        (ce1 == CGAL::ARR_MIN_END ? 
                         CGAL::POSITIVE : CGAL::NEGATIVE)
                );
                if ((ce1 == CGAL::ARR_MAX_END &&
                     loc1 == CGAL::ARR_TOP_BOUNDARY) ||
                    (ce1 == CGAL::ARR_MIN_END &&
                     loc1 == CGAL::ARR_BOTTOM_BOUNDARY)) {
                    res = -res;
                }
                CERR("result: " << res << "\n");
                return res;
            }
            // else: order can be determined without y-comparison
            //(loc1 == CGAL::ARR_BOTTOM_BOUNDARY ? CGAL::SMALLER :
            res = CGAL::EQUAL;
                   //CGAL::LARGER);
            CERR("result: " << res << "\n");
            return res;
        }
        // curve ends approach vertical asymptote (or singularity) from
        // different sides => no comparisons required
        res =  (ce1 == CGAL::ARR_MIN_END ? CGAL::LARGER : CGAL::SMALLER);
        CERR("result: " << res << "\n");
        return res;
    }
};


template < class CurvedKernelViaAnalysis_2 >
class Compare_y_near_boundary_2 : public 
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef 
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<3>            Arity;
    
    //! standard constructor
    Compare_y_near_boundary_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
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
    result_type operator()(const Arc_2& cv1, const Arc_2& cv2,
                           CGAL::Arr_curve_end ce) const {\

        CERR("\ncompare_y_near_boundary; cv1: " << cv1 << "; cv2: " <<
             cv2 << "; end: " << ce << "\n");

        CGAL::Arr_parameter_space loc1 = cv1.location(ce);
        CGAL_precondition(cv1.is_on_left_right(loc1) &&
                          loc1 == cv2.location(ce));
        // comparing ids is the same as calling is_identical() ??
        if (cv1.id() == cv2.id()) {
            return CGAL::EQUAL;
        } 
        
        // in this setting same handling as for +/-oo ?
        CGAL::Comparison_result res = cv1._compare_arc_numbers(cv2, loc1);
        CERR("result: " << res << "\n");
        return res;
    }
};


template < class CurvedKernelViaAnalysis_2 >
class Compare_y_at_x_2 : public 
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {
   
public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef 
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<2>            Arity;
    
    //! standard constructor
    Compare_y_at_x_2(Curved_kernel_via_analysis_2 *kernel) :
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
     
#if 0

        CERR("\ncompare_y_at_x; p: " << p << ";\n cv:" << cv << "\n");
        CGAL::Arr_parameter_space loc1 = cv.location(CGAL::ARR_MIN_END),
            loc2 = cv.location(CGAL::ARR_MAX_END);/*, locp = p.location();*/
        
        /*CGAL::Boundary_type bndp_x = p.boundary_in_x(), 
            bndp_y = p.boundary_in_y(), 
            bnd1_x = boundary_in_x(CGAL::ARR_MIN_END), 
            bnd2_x = boundary_in_x(CGAL::ARR_MAX_END), 
            bnd1_y = boundary_in_y(CGAL::ARR_MIN_END),
            bnd2_y = boundary_in_y(CGAL::ARR_MAX_END);*/
        
        //CGAL_precondition(!(is_infinite(bndp_x) || is_infinite(bndp_y)));
        // handle special case when a curve end coincides with p at singularity
        /*if((bndp_x == CGAL::AFTER_SINGULARITY && bnd1_x == bndp_x) ||
            (bndp_x == CGAL::BEFORE_SINGULARITY && bnd2_x == bndp_x) ||
           (bndp_y == CGAL::AFTER_SINGULARITY && bnd1_y == bndp_y) ||
            (bndp_y == CGAL::BEFORE_SINGULARITY && bnd2_y == bndp_y))
            return CGAL::EQUAL;
        CGAL_precondition_msg(!is_singular(bndp_x), "Target point is not "
            "within the arc's x-range"); 
                
        if(is_singular(bndp_y)) {// singularity in y is always in x-range 
             if(bndp_y < CGAL::NO_BOUNDARY)
                 return CGAL::SMALLER;
             return CGAL::LARGER; // bndp_y > 0
        }*/
        bool eq_min = false, eq_max = false, in_x_range = true;
        /*if(is_on_disc(bndp_x)) {
            eq_min = (bndp_x < CGAL::NO_BOUNDARY && bnd1_x == bndp_x);
            eq_max = (bndp_x > CGAL::NO_BOUNDARY && bnd2_x == bndp_x);
            // report x-range assert violation if the point lies on disc in
            // x but neither of arc's ends do
            if(!(eq_min || eq_max))
                CGAL_error_msg("Target point is not within the arc's x-range");
        } else // we should be able to access x-coord when point is on disc */
            in_x_range = cv.is_in_x_range(p.x(), &eq_min, &eq_max);

        CGAL_precondition(in_x_range); // check x-range
        /*if(is_on_disc(bndp_y)) {
            if((eq_min && bndp_y < CGAL::NO_BOUNDARY && bnd1_y == bndp_y) ||
               (eq_max && bndp_y > CGAL::NO_BOUNDARY && bnd2_y == bndp_y))
               return CGAL::EQUAL;
             // otherwise handle by the boundary type
             if(bndp_y < CGAL::NO_BOUNDARY) 
                 return CGAL::SMALLER;
             return CGAL::LARGER; // bndp_y > 0
        }*/
        
        if (cv.is_vertical()) {
            if (cv.is_interior(loc1)) {
                // for vertical arcs we can ask for .xy() member
                if (Curved_kernel_via_analysis_2::instance().
                    compare_xy_2_object()(
                            p, cv._minpoint(), true
                    ) == CGAL::SMALLER) {
                    return CGAL::SMALLER;
                }
            }
            if (cv.is_interior(loc2)) {
                if (Curved_kernel_via_analysis_2::instance().
                    compare_xy_2_object()(
                            p, cv._maxpoint(), true
                    ) == CGAL::LARGER) {
                    return CGAL::LARGER;
                }
            }
            return CGAL::EQUAL; // p lies on a vertical arc
        }
        if (eq_min && loc1 != CGAL::ARR_INTERIOR) {
            return (loc1 == CGAL::ARR_BOTTOM_BOUNDARY ? 
                    CGAL::LARGER : CGAL::SMALLER);
        }
        if (eq_max && loc2 != CGAL::ARR_INTERIOR) {
            return (loc2 == CGAL::ARR_BOTTOM_BOUNDARY ? 
                    CGAL::LARGER : CGAL::SMALLER);
        }
        // what remains to be handled ?    
        /*if(is_on_disc(bndp_x)) { 
            // the point and a respective curve end lie on disc in x => need
            // comparison at x-infinity; 
            // since we compare point agains the arc: reverse the result
            return (- _compare_arc_numbers(p.xy(), bnd1_x));
        }*/
        // otherwise return reversed y-order of this arc and point p
        CGAL::Comparison_result res;
        if (eq_min) {
            res = Curved_kernel_via_analysis_2::instance().
                compare_xy_2_object()(
                    p, cv._minpoint(), true
            );
        } else if (eq_max) {
            res = Curved_kernel_via_analysis_2::instance().
                compare_xy_2_object()(
                    p, cv._maxpoint(), true
            );
        } else {
            res = -cv._compare_arc_numbers(p.xy(), CGAL::ARR_INTERIOR, p.x());
        }
        CERR("cmp result: " << res << "\n");
        return res;

#else

        CERR("\ncompare_y_at_x; p: " << p << ";\n cv:" << cv << "\n");
        CGAL::Arr_parameter_space loc1 = cv.location(CGAL::ARR_MIN_END),
            loc2 = cv.location(CGAL::ARR_MAX_END);/*, locp = p.location();*/
        bool eq_min, eq_max;

        CGAL_assertion_code (
           bool in_x_range =
        );
        cv.is_in_x_range(p.x(), &eq_min, &eq_max);
        CGAL_assertion(in_x_range);

        if (cv.is_vertical()) {
            if (cv.is_interior(loc1)) {
                // for vertical arcs we can ask for .xy() member
                if (Curved_kernel_via_analysis_2::instance().
                    compare_xy_2_object()(
                            p, cv._minpoint(), true
                    ) == CGAL::SMALLER) {
                    return CGAL::SMALLER;
                }
            }
            if (cv.is_interior(loc2)) {
                if (Curved_kernel_via_analysis_2::instance().
                    compare_xy_2_object()(
                            p, cv._maxpoint(), true
                    ) == CGAL::LARGER) {
                    return CGAL::LARGER;
                }
            }
            return CGAL::EQUAL; // p lies on a vertical arc
        }
        CGAL::Comparison_result res;
        if(eq_min) {
            res = Curved_kernel_via_analysis_2::instance().
                compare_xy_2_object()(
                        p, cv._minpoint(), true
                );
        } else if(eq_max) {
            res = Curved_kernel_via_analysis_2::instance().
                compare_xy_2_object()(
                        p, cv._maxpoint(), true
                );
            
        } else {
            Point_2 point_on_s
                = Curved_kernel_via_analysis_2::instance().
                construct_point_on_arc_2_object()
                ( p.x(), 
                  cv.curve(), 
                  cv.arcno(),
                  cv );
            res = Curved_kernel_via_analysis_2::instance().
                compare_xy_2_object()(p, point_on_s, true);
        }
        CERR("cmp result: " << res << "\n");
        return res;

#endif

    }
};


template < class CurvedKernelViaAnalysis_2 >
class Compare_y_at_x_left_2 : public 
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {
    
public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef 
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<3>            Arity;
    
    //! standard constructor
    Compare_y_at_x_left_2(Curved_kernel_via_analysis_2 *kernel) :
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

        CERR("\ncompare_y_at_x_left(cv2); cv1: " << cv1 << "; cv2: " <<
            cv2 << "; p: " << p << "\n");

        CGAL_precondition_code(
        CGAL::Arr_parameter_space locp = p.location();
        // ensure that p lies on both arcs and doesn't lie on the negative 
        // boundary
        CGAL_precondition(locp != CGAL::ARR_LEFT_BOUNDARY && 
                          cv1.compare_y_at_x(p) == CGAL::EQUAL && 
                          cv2.compare_y_at_x(p) == CGAL::EQUAL);
        );
        // check whether both arcs indeed lie to the left of p
        CGAL_precondition(
                (cv1.is_vertical() && 
                 cv1.location(CGAL::ARR_MIN_END) == 
                 CGAL::ARR_BOTTOM_BOUNDARY) ||
                cv1._same_arc_compare_xy(cv1._minpoint(), p) == CGAL::SMALLER
        );
        CGAL_precondition(
                (cv2.is_vertical() &&
                 cv2.location(CGAL::ARR_MIN_END) == CGAL::ARR_BOTTOM_BOUNDARY) 
                ||
                cv2._same_arc_compare_xy(cv2._minpoint(), p) == CGAL::SMALLER
        );
        if (cv1.is_vertical()) {
            // if both are vertical (they overlap), we return EQUAL
            if(cv2.is_vertical()) {
                return CGAL::EQUAL;
            }
            // a vertical arc is always smaller than the arc extending to the
            // left
            return CGAL::SMALLER;
        } 
        // a vertical arc is always smaller than the arc extending to the left;
        // due to the order, we have to return the opposite
        if (cv2.is_vertical()) {
            return CGAL::LARGER;
        }
            
        if (cv1.is_singular()) {// singularity in y
            CGAL_error_msg("Handling singularity in y is not yet implemented");
        }
        
        // vertical line immediately to the left of p: if p lies on boundary
        // get the vertical line over the last interval; otherwise
        // obtain the interval w.r.t. point's x-coordinate (this also valid
        // for discontinuity in y)
        /*if(bndp_x == CGAL::BEFORE_SINGULARITY || 
           bndp_x == CGAL::BEFORE_DISCONTINUITY)
           return _compare_arc_numbers(cv2, bndp_x);
        else*/

        CGAL::Comparison_result res = 
            cv1._compare_arc_numbers(cv2,
                                     CGAL::ARR_INTERIOR, 
                                     p.x(), CGAL::NEGATIVE);
        CERR("result: " << res << "\n");
        return res;
    }
};


template < class CurvedKernelViaAnalysis_2 >
class Compare_y_at_x_right_2 : public 
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {
   
public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef 
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<3>            Arity;
    
    //! standard constructor
    Compare_y_at_x_right_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }

     /*!
     * Compares the y value of two x-monotone curves immediately 
     * to the right of their intersection point. If one of the curves is
     * vertical (emanating upward from p), it's always considered to be above
     * the other curve.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param p The intersection point.
     * \pre The point p lies on both curves, and both of them must be 
     * also be defined (lexicographically) to its right.
     * \return The relative position of cv1 with respect to 
     * cv2 immdiately to the right of p: SMALLER, LARGER or EQUAL.
     */
    result_type operator()(const Arc_2& cv1, const Arc_2& cv2,
                           const Point_2& p) const {

        CERR("\ncompare_y_at_x_right; cv1: " << cv1 << "; cv2: " <<
            cv2 << "; p: " << p << "\n");
        
        CGAL_precondition_code(
                CGAL::Arr_parameter_space locp = p.location();
        );
        // ensure that p lies on both arcs and doesn't lie on the positive
        // boundary
        CGAL_precondition(locp != CGAL::ARR_RIGHT_BOUNDARY);
        CGAL_precondition(cv1.compare_y_at_x(p) == CGAL::EQUAL);
        CGAL_precondition(cv2.compare_y_at_x(p) == CGAL::EQUAL);
        
        // check whether both arcs indeed lie to the left of p    
        CGAL_precondition(
                (cv1.is_vertical() && 
                 cv1.location(CGAL::ARR_MAX_END) == CGAL::ARR_TOP_BOUNDARY) ||
                cv1._same_arc_compare_xy(p, cv1._maxpoint()) == CGAL::SMALLER
        );
        CGAL_precondition(
                (cv2.is_vertical() &&
                 cv2.location(CGAL::ARR_MAX_END) == CGAL::ARR_TOP_BOUNDARY) ||
                cv2._same_arc_compare_xy(p, cv2._maxpoint()) == CGAL::SMALLER
        );
        
        if (cv1.is_vertical()) {
            // if both are vertical (they overlap), we return EQUAL
            if (cv2.is_vertical()) {
                return CGAL::EQUAL;
            }
            // a vertical arc is always LARGER than arc extending to the
            // right
            return CGAL::LARGER;
        } 
        // a vertical arc is always LARGER than arc extending to the right; 
        // due to the order, we have to return the opposite
        if (cv2.is_vertical()) {
            return CGAL::SMALLER;
        }
        
        if (cv1.is_singular()) {// singularity in y
            CGAL_error_msg("Handling singularity in y is not yet \
                implemented");
        }
            
        // vertical line immediately to the right of p: if p lies on boundary
        // get the vertical line over the first interval; otherwise
        // obtain the interval w.r.t. point's x-coordinate (this also valid
        // for discontinuity in y)        
        /*if(bndp_x == CGAL::AFTER_SINGULARITY || 
                bndp_x == CGAL::AFTER_DISCONTINUITY)
           return _compare_arc_numbers(cv2, bndp_x);
        else*/
        CGAL::Comparison_result res = 
            cv1._compare_arc_numbers(cv2,
                                     CGAL::ARR_INTERIOR, 
                                     p.x(), 
                                     CGAL::POSITIVE);
        CERR("result: " << res << "\n");
        return res;
    }
};


template < class CurvedKernelViaAnalysis_2 >
class Is_in_x_range_2 : public 
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {
   
public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef 
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef bool result_type;
    typedef Arity_tag<2>            Arity;
    
    //! standard constructor
    Is_in_x_range_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }
    
    /*!\brief
     * Check whether a given point lies within the curve's x-range
     * \param cv The curve.
     * \param p the point
     * \return (true) if p lies in arc's x-range; (false) otherwise.
     */
    bool operator()(const Arc_2& cv, const Point_2& p) const {
        return cv.is_in_x_range(p);
    }
};


//!\brief Tests two objects, whether they are equal
template < class CurvedKernelViaAnalysis_2 >
class Equal_2 : public 
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef 
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef bool result_type;
    typedef Arity_tag<2> Arity;
    
    //! standard constructor
    Equal_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }
    
    /*!
     * Check if the two points are the same.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return (true) if the two point are the same; (false) otherwise.
     */
    result_type operator()(const Point_2& p1, const Point_2& p2) const {
        return (Curved_kernel_via_analysis_2::instance().
                compare_xy_2_object()(p1, p2) == 
                CGAL::EQUAL);
    }
     
    /*!
     * Check if the two x-monotone curves are the same (have the same graph).
     * \param cv1 The first 
     *        curve(this->_ckva()->kernel().compare_xy_2_object()
             (p1.xy(), p2.xy()));.
     * \param cv2 The second curve.
     * \return (true) if the two curves are the same; (false) otherwise.
     */
    result_type operator()(const Arc_2& cv1, const Arc_2& cv2) const {

        if (cv1.is_identical(cv2)) {
            return true;
        }
        // only one of the arcs is vertical => not equal
        if (cv1.is_vertical() != cv2.is_vertical()) {
            return false;
        }
        
        // distinct supporting curves implies inequality, provided the
        // coprimality condition is satisfied
        if (!cv1.curve().is_identical(cv2.curve())) {
            return false;
        }
        
        // here either both or none of the arcs are vertical, check for arcnos
        // equality
        if (!cv1.is_vertical() && cv1.arcno() != cv2.arcno()) {
            return false;
        }
        // otherwise compare respective curve ends: supporting curves and 
        // arcnos are equal => the curve ends belong to the same arc
        return ((cv1._same_arc_compare_xy(cv1._minpoint(), cv2._minpoint()) ==
                 CGAL::EQUAL &&
                 cv1._same_arc_compare_xy(cv1._maxpoint(), cv2._maxpoint()) ==
                 CGAL::EQUAL));
    }
    
};


template < class CurvedKernelViaAnalysis_2 >
class Do_overlap_2 : public 
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {
    
public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef 
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef bool result_type;
    typedef Arity_tag<2> Arity;
    
    //! standard constructor
    Do_overlap_2(Curved_kernel_via_analysis_2 *kernel) :
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
    
        CERR("\ndo_overlap\n");
        if (cv1.is_identical(cv2)) {
            return true;
        }
            
        Arc_2::simplify(cv1, cv2);
        // arcs with coprime support do not overlap    
        if (!cv1.curve().is_identical(cv2.curve())) {
            return false;
        }
        
        if (cv1.is_vertical() != cv2.is_vertical()) {
            return false; // only one arc is vertical => can't overlap
        }
        if (cv1.is_vertical()) { // both arcs are vertical
            // check for x-coordinates equality
            if (Curved_kernel_via_analysis_2::instance().
                kernel().compare_x_2_object()(
                        cv1._minpoint().x(),
                        cv2._minpoint().x()) != CGAL::EQUAL) {
                return false;
            }
            // compare y-coordinates of min curve ends
            switch(cv1._same_arc_compare_xy(
                           cv1._minpoint(), cv2._minpoint(), true)
            ) {
            case CGAL::EQUAL: // this->source == cv2->source => overlap !
                return true;            
            case CGAL::SMALLER: // this->source < cv2->source
                // check whether this->target > cv2->source
                return (cv1._same_arc_compare_xy(
                                cv1._maxpoint(), cv2._minpoint(), 
                                true) == CGAL::LARGER);
            case CGAL::LARGER: // this->source > cv2->source
                // check whether this->source < cv2->target
                return (cv1._same_arc_compare_xy(
                                cv1._minpoint(), cv2._maxpoint(), 
                                true) == CGAL::SMALLER);
            }
        }
        // processing non-vertical arcs
        if (cv1.arcno() != cv2.arcno()) {
            return false;
        }
        /* At this point, we have two non-vertical arcs supported by the same
         * curve with equal arc numbers in their interior. They do overlap if
         * their x-ranges overlap. Compare only x-coordinates */
        switch (cv1._same_arc_compare_xy(
                        cv1._minpoint(), cv2._minpoint(), false, true)
        ) {
        case CGAL::EQUAL: // this->source == cv2->source => overlap !
            return true;            
        case CGAL::SMALLER: // this->source < cv2->source
            // check whether this->target > cv2->source
            return (cv1._same_arc_compare_xy(
                            cv1._maxpoint(), cv2._minpoint(), false, 
                            true) == CGAL::LARGER);
        case CGAL::LARGER: // this->source > cv2->source
            // check whether this->source < cv2->target
            return (cv1._same_arc_compare_xy(
                            cv1._minpoint(), cv2._maxpoint(), false, 
                            true) == CGAL::SMALLER);
        }
        CGAL_error_msg("bogus comparison result");
        return false;
    }
};


template < class CurvedKernelViaAnalysis_2 >
class Intersect_2 : public 
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef 
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef std::iterator< output_iterator_tag, CGAL::Object > result_type;
    typedef Arity_tag<3> Arity;    
    
    //! standard constructor
    Intersect_2(Curved_kernel_via_analysis_2 *kernel) :
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
        
        CERR("\nintersect; cv1: " << cv1 
             << ";\n cv2:" << cv2 << "");
        
        // if arcs overlap, just store their common part, otherwise compute
        // point-wise intersections
        std::vector< Arc_2 > common_arcs;
        if (cv1._trim_if_overlapped(cv2, std::back_inserter(common_arcs))) {
            typename std::vector< Arc_2 >::const_iterator it;
            for(it = common_arcs.begin(); it < common_arcs.end(); it++) {
                *oi++ = CGAL::make_object(*it);
            }
            return oi; 
        }
        // process non-ov erlapping case        
        typedef std::pair< Point_2, unsigned int > Point_and_mult;
        typedef std::vector< Point_and_mult > Point_vector;
        Point_vector vec;
        typename Point_vector::const_iterator it;
        Arc_2::_intersection_points(cv1, cv2, std::back_inserter(vec));

        //std::cout << "results\n";
        for (it = vec.begin(); it != vec.end(); it++) {
            *oi++ = CGAL::make_object(*it);
        }
        return oi;
    }

};


template < class CurvedKernelViaAnalysis_2 >
class Trim_2 : public 
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {
   
public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef 
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef Arc_2 result_type;
    typedef Arity_tag<3> Arity;
    
    //! standard constructor
    Trim_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }
    
    /*!\brief
     * Returns a 
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the curves overlap; (false) otherwise.
     */
    Arc_2 operator()(const Arc_2& cv, const Point_2& p, const Point_2& q) {
    
        CERR("trim\n");

        CGAL_precondition(p.location()==CGAL::ARR_INTERIOR);
        CGAL_precondition(q.location()==CGAL::ARR_INTERIOR);
        
        CGAL_precondition(
                Curved_kernel_via_analysis_2::instance().
                compare_xy_2_object()(p, q) != CGAL::EQUAL
        );
        CGAL_precondition(cv.compare_y_at_x(p) == CGAL::EQUAL);
        CGAL_precondition(cv.compare_y_at_x(q) == CGAL::EQUAL);  

        return cv._trim(p, q);
    }
};

template < class CurvedKernelViaAnalysis_2 >
class Split_2 : public 
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {
   
public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef 
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef void result_type;
    typedef Arity_tag<4> Arity;    
    
    //! standard constructor
    Split_2(Curved_kernel_via_analysis_2 *kernel) :
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
        
        CGAL_precondition(cv.compare_y_at_x(p) == CGAL::EQUAL);
        // check that p is not an end-point of the arc
        CGAL_precondition_code(
                cv._same_arc_compare_xy(cv._minpoint(), p) != CGAL::EQUAL &&
                cv._same_arc_compare_xy(cv._maxpoint(), p) != CGAL::EQUAL);
        
        CERR("\nsplit\n");
        c1 = cv._replace_endpoints(
                cv._minpoint(), p, -1, (cv.is_vertical() ? -1 : cv.arcno())
        );
        c2 = cv._replace_endpoints(
                p, cv._maxpoint(), (cv.is_vertical() ? -1 : cv.arcno()), -1
        );
    }
};


template < class CurvedKernelViaAnalysis_2 >
class Are_mergeable_2 : public 
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {
   
public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef 
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef bool result_type;
    typedef Arity_tag<2> Arity;    
    
    //! standard constructor
    Are_mergeable_2(Curved_kernel_via_analysis_2 *kernel) :
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
    
        CERR("\nare_mergeable\n");
        
        if (cv1.do_overlap(cv2)) {// if arcs overlap they are not mergeable
            return false;   // this also call simplify
        }
        
        // touch in at most one point now and supporting curves are simplified
        // both arcs must be either vertical or not
        if (!cv1.curve().is_identical(cv2.curve()) || 
            cv1.is_vertical() != cv2.is_vertical()) {
            return false;
        }
        // if both arcs are non-vertical => they must have equal arcnos
        // to be mergeable            
        if (!cv1.is_vertical() && cv1.arcno() != cv2.arcno()) {
            return false;
        }
        // for non-vertical arcs arc numbers are equal => can use same_arc_cmp
        bool max_min = 
            (cv1._same_arc_compare_xy(cv1._maxpoint(), cv2._minpoint()) == 
             CGAL::EQUAL),
            min_max = false;
        if (!max_min) { // both cases cannot happen simultaneously
            min_max = 
                (cv1._same_arc_compare_xy(cv1._minpoint(), cv2._maxpoint()) == 
                 CGAL::EQUAL);
            if (!min_max) { // arcs have no common end-point => not mergeable
                return false;
            }
        }
        // check that the common point is not an event point
        if (cv1.is_vertical()) { // both arcs are vertical 
            Point_2 common = (max_min ? cv1._maxpoint() : cv1._minpoint());
            // a common end must be a finite point
            CGAL_precondition(cv1.is_interior(common.location()));
            // check that there are no other non-vertical branches coming 
            // through this point
            typedef typename 
                Curved_kernel_via_analysis_2::Curve_kernel_2::Curve_analysis_2 
                Curve_analysis_2;
            Curve_analysis_2 ca_2(cv1.curve());
            typename Curve_analysis_2::Status_line_1 cv_line = 
                ca_2.status_line_for_x(common.x());
            CGAL_assertion(cv_line.is_event()); // ??
            // we are not allowed to use number_of_incident_branches()
            // since the common point might be supported by different curve, 
            // and therefore its arcno might be not valid for *this arc
            for (int k = 0; k < cv_line.number_of_events(); k++) {
                // use a temporary object for comparison predicate
                typename Point_2::Xy_coordinate_2
                    tmp(common.x(), cv1.curve(), k);
                if (Curved_kernel_via_analysis_2::instance().
                    kernel().compare_xy_2_object()(
                            common.xy(), tmp) == 
                    CGAL::EQUAL) {
                    return false;
                }
            }
        } else if (cv1.interval_id() != cv2.interval_id()) {
            return false; // non-vertical case
        }
        
        return true;
    }
};

template < class CurvedKernelViaAnalysis_2 >
class Merge_2 : public 
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {
   
public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef 
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef void result_type;
    typedef Arity_tag<2> Arity;    
    
    //! standard constructor
    Merge_2(Curved_kernel_via_analysis_2 *kernel) :
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
    
        CERR("merge\n");
        CGAL_precondition(cv1.are_mergeable(cv2));
        Arc_2::simplify(cv1, cv2);
        
        Point_2 src, tgt;
        int arcno_s = -1, arcno_t = -1;
        bool replace_src; // true if cv2 < *this otherwise *this arc < cv2 arc
        // arcs are mergeable => they have one common finite end-point
        replace_src = 
            (cv1._same_arc_compare_xy(cv1._minpoint(), cv2._maxpoint()) == 
             CGAL::EQUAL);
        src = (replace_src ? cv2._minpoint() : cv1._minpoint());
        tgt = (replace_src ? cv1._maxpoint() : cv2._maxpoint());
              
        if (!cv1.is_vertical()) {
            arcno_s = (replace_src ? cv2.arcno(CGAL::ARR_MIN_END) :
                       cv1.arcno(CGAL::ARR_MIN_END));
            arcno_t = (replace_src ? cv1.arcno(CGAL::ARR_MAX_END) :
                       cv2.arcno(CGAL::ARR_MAX_END));
        }
        Arc_2 arc = cv1._replace_endpoints(src, tgt, arcno_s, arcno_t);
        // arc.set_boundaries_after_merge(*this, s); - no need to, since
        // boundaries are stored in Point_2 type and will be copied implicitly
        
        c = arc;
    }
};


//!\brief Tests whether a point lies on a supporting curve
template < class CurvedKernelViaAnalysis_2 >
class Is_on_2 : public 
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef 
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef bool result_type;
    typedef Arity_tag<2> Arity;
    
    //! standard constructor
    Is_on_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }
    
    /*!
     * Checks whether \c p lies on \c c 
     * \param p1 The point to test
     * \param p2 The curve
     * \return (true) if the \c p lies on \c c
     */
    result_type operator()(const Point_2& p, const Curve_2& c) const {
        result_type res = 
            (Curved_kernel_via_analysis_2::instance().
             kernel().sign_at_2_object()(c, p.xy())
             == CGAL::ZERO);
        return res;
    }
};



template < class CurvedKernelViaAnalysis_2>
class Make_x_monotone_2 : public 
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef 
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef std::iterator<output_iterator_tag, CGAL::Object> result_type;
    typedef Arity_tag<2> Arity;   
    
    //! standard constructor
    Make_x_monotone_2(Curved_kernel_via_analysis_2 *kernel) :
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
    // FUTURE TODO: move this to separate file Arr_kernel_traits_2.h ?
    template < class OutputIterator >
    OutputIterator operator()(const Curve_2& cv, OutputIterator oi) const {
    
        CGAL::CGALi::Make_x_monotone_2<Curved_kernel_via_analysis_2>
            make_x_monotone(&Curved_kernel_via_analysis_2::instance());
        return make_x_monotone(cv, oi);
    }
};

#undef CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES 

} // namespace Curved_kernel_via_analysis_2_Functors

//! collects main set of functors of a curved kernel
template < 
    class CurvedKernelViaAnalysis_2, 
    class Curve_2_, 
    class Point_2_, 
    class Arc_2_
>
class Curved_kernel_via_analysis_2_functors {
public:
    //! this instance's first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! this instance's second template parameter
    typedef Curve_2_ Curve_2;
    
    //! this instance's third template parameter
    typedef Point_2_ Point_2;
    
    //! this instance's fourth template parameter
    typedef Arc_2_ Arc_2;
    
    // declares curved kernel functors, 
    // for each functor defines a member function
// returning an instance of this functor
#define CGAL_CKvA_2_functor_pred(Y, Z) \
    typedef \
    Curved_kernel_via_analysis_2_Functors::Y< Curved_kernel_via_analysis_2 > \
    Y; \
    Y Z() const { return Y(&Curved_kernel_via_analysis_2::instance()); }

#define CGAL_CKvA_2_functor_cons(Y, Z) CGAL_CKvA_2_functor_pred(Y, Z)

    CGAL_CKvA_2_functor_pred(Compare_x_2, compare_x_2_object);  
    CGAL_CKvA_2_functor_pred(Compare_xy_2, compare_xy_2_object);
    CGAL_CKvA_2_functor_pred(Is_vertical_2, is_vertical_2_object); 
    CGAL_CKvA_2_functor_pred(Is_bounded_2, is_bounded_2_object);
    CGAL_CKvA_2_functor_pred(Parameter_space_in_x_2, 
                             parameter_space_in_x_2_object);
    CGAL_CKvA_2_functor_pred(Parameter_space_in_y_2, 
                             parameter_space_in_y_2_object);
    CGAL_CKvA_2_functor_cons(Construct_min_vertex_2,
                             construct_min_vertex_2_object);
    CGAL_CKvA_2_functor_cons(Construct_max_vertex_2,
                             construct_max_vertex_2_object);
    CGAL_CKvA_2_functor_pred(Compare_x_near_boundary_2,
                             compare_x_near_boundary_2_object);
    CGAL_CKvA_2_functor_pred(Compare_y_near_boundary_2,
                             compare_y_near_boundary_2_object);
    CGAL_CKvA_2_functor_pred(Compare_y_at_x_2, compare_y_at_x_2_object);   
    CGAL_CKvA_2_functor_pred(Compare_y_at_x_left_2,
                             compare_y_at_x_left_2_object);
    CGAL_CKvA_2_functor_pred(Compare_y_at_x_right_2,
                             compare_y_at_x_right_2_object);
    CGAL_CKvA_2_functor_pred(Equal_2, equal_2_object); 
    CGAL_CKvA_2_functor_pred(Is_in_x_range_2, is_in_x_range_2_object);
    CGAL_CKvA_2_functor_pred(Do_overlap_2, do_overlap_2_object);
    CGAL_CKvA_2_functor_cons(Intersect_2, intersect_2_object);
    CGAL_CKvA_2_functor_cons(Trim_2, trim_2_object);
    CGAL_CKvA_2_functor_cons(Split_2, split_2_object);  
    CGAL_CKvA_2_functor_pred(Are_mergeable_2, are_mergeable_2_object); 
    CGAL_CKvA_2_functor_cons(Merge_2, merge_2_object); 


    CGAL_CKvA_2_functor_pred(Is_on_2, is_on_2_object); 
    CGAL_CKvA_2_functor_cons(Make_x_monotone_2, make_x_monotone_2_object);


#undef CGAL_CKvA_2_functor_pred
#undef CGAL_CKvA_2_functor_cons

}; // Curved_kernel_via_analysis_2_functors

} // namespace CGALi


CGAL_END_NAMESPACE

#endif // CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_FUNCTORS_H
