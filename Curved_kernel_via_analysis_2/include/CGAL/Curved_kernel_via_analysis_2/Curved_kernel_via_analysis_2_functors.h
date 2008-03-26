// TODO add licence
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

/*! \file Curved_kernel_via_analysis_2/Curved_kernel_via_analysis_2_functors.h
 *  \brief defines Curved_kernel_via_analysis_2 function objects + class
 */

#include <CGAL/basic.h>
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

/*!\brief
 * Collects main functors for Curved_kernel_via_analysis_2
 */
namespace Curved_kernel_via_analysis_2_Functors {

/*!\brief
 * Base functor class for inheritance
 */
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

    //! type of curve kernel
    typedef typename Curved_kernel_via_analysis_2::Curve_kernel_2 
    Curve_kernel_2;

    //! type of curve analaysis
    typedef typename Curve_kernel_2::Curve_analysis_2 Curve_analysis_2;
    
    //! the x-coordinate type
    typedef typename Point_2::X_coordinate_1 X_coordinate_1;

    //!@}

public:
    //!\name Constructors
    //!@{
    
    /*!\brief
     * Constructs base functor
     *
     * \param kernel The kernel instance
     */
    Curved_kernel_via_analysis_2_functor_base(
            Curved_kernel_via_analysis_2 *kernel) :
        _m_curved_kernel(kernel) {
        CGAL_precondition(kernel != NULL);
    }

    //!@}

protected:
    //!\name Access members
    //!@{
    
    /*!\brief 
     * Return pointer to curved kernel
     * \return Pointer to stored kernel
     */
    Curved_kernel_via_analysis_2* _ckva() const {
        return _m_curved_kernel;
    }

    //!@}
    
protected:
    //!\name Data members
    //!@{
    
    //! stores pointer to \c Curved_kernel_via_analysis_2
    Curved_kernel_via_analysis_2 *_m_curved_kernel;
    
    //!@}
};


#define CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES \
    typedef typename Base::Curve_2 Curve_2; \
    typedef typename Base::Point_2 Point_2; \
    typedef typename Base::Arc_2 Arc_2; \
    typedef typename Base::Curve_analysis_2 Curve_analysis_2; \
    typedef typename Base::X_coordinate_1 X_coordinate_1; \

// end define


/*!\brief 
 * Functor to construct a point on a curve
 */
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

    /*!\brief 
     * Standard constructor
     *
     * \param kernel The kernel
     */
    Construct_point_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }

    /*!\brief
     * Constructs a interior point 
     * 
     * \param x x-coordinate
     * \param c The supporting curve
     * \param arcno Arcnumber on curve
     * \return The constructed point
     */
    Point_2 operator()(const X_coordinate_1& x, 
                       const Curve_analysis_2& c, 
                       int arcno) {
        Point_2 pt(x, c, arcno);
        return pt;
    }
};


/*!\brief
 * Functor to construct point on an arc
 */
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
    
    /*!\brief 
     * Standard constructor
     *
     * \param kernel The kernel
     */
    Construct_point_on_arc_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }
    
    /*!\brief
     * Constructs point on an arc
     * 
     * \param x x-coordinate of point
     * \param c The supporting curve
     * \param arcno Arcnumber on curve
     * \param arc Can be used to query meta data
     * \return The constructed point
     */
    Point_2 operator()(
            const X_coordinate_1& x,
            const Curve_analysis_2& c, int arcno,
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


/*!\brief
 * Functor to construct an x-monotone arc 
 */
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

    /*!\brief 
     * Standard constructor
     *
     * \param kernel The kernel
     */
    Construct_arc_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }

    //!\name Constructing non-vertical arcs
    //!@{
    
    /*!\brief 
     * Constructs a non-vertical arc with two interior end-points (segment).
     * 
     * \param p first endpoint
     * \param q second endpoint
     * \param c The supporting curve
     * \param arcno The arcnumber wrt \c c in the interior of the arc
     * \param arcno_p The arcnumber wrt \c c of the arc at \c p
     * \param arcno_q The arcnumber wrt \c c of the arc at \c q
     * \returns The constructed segment
     * 
     * \pre p.x() != q.x()
     *
     */
    Arc_2 operator()(const Point_2& p, const Point_2& q, 
                     const Curve_analysis_2& c,
                     int arcno, int arcno_p, int arcno_q) {
        Arc_2 arc(p, q, c, arcno, arcno_p, arcno_q);
        return arc;
    }
    
    /*!\brief
     * Constructs a non-vertical arc with one interior end-point and whose
     * other end approaches the left or right boundary of the parameter space 
     * (ray I)
     *
     * \param origin The interior end-point of the ray
     * \param inf_end Defining whether the arc emanates from the left or right
     *        boundary
     * \param c The supporting curve
     * \param arcno The arcnumber wrt \c c in the interior of the arc
     * \param arcno_o The arcnumber wrt \c c of the arc at \c origin 
     * \return The constructed ray
     */
    Arc_2 operator()(const Point_2& origin, CGAL::Arr_curve_end inf_end, 
                     const Curve_analysis_2& c, int arcno, int arcno_o) {
        Arc_2 arc(origin, inf_end, c, arcno, arcno_o);
        return arc;
    }

    /*!\brief
     * Constructs a non-vertical arc with one interior end-point and whose 
     * other end approaches a vertical asymptote (ray II)
     *
     * \param origin The interior end-point
     * \param asympt_x The x-coordinate of the vertical asymptote
     * \param inf_end Arc is approaching the bottom or top boundary
     * \param c The supporting curve
     * \param arcno The arcnumber wrt \c c in the interior of the arc
     * \param arcno_o The arcnumber wrt \c c of the arc at \c origin 
     * \return The constructed ray
     *
     * \pre origin.x() != asympt_x
     */
    Arc_2 operator()(const Point_2& origin, const X_coordinate_1& asympt_x, 
                     CGAL::Arr_curve_end inf_end, 
                     const Curve_analysis_2& c, int arcno, 
                     int arcno_o) {
        Arc_2 arc(origin, asympt_x, inf_end, c, arcno, arcno_o);
        return arc;
    }

    /*!\brief
     * Constructs a non-vertical arc with two non-interior ends at 
     * the left and right boundary (branch I)
     *
     * \param c The supporting curve
     * \param arcno The arcnumber wrt to \c c in the interior of the arc
     * \return The constructed branch
     */
    Arc_2 operator()(const Curve_analysis_2& c, int arcno) {
        Arc_2 arc(c, arcno);
        return arc;
    }
    
    /*!\brief
     * Constructs a non-vertical arc with two ends approaching vertical
     * asymptotes (branch II).
     *
     * \param asympt_x1 The x-coordinate of the first asymptote
     * \param inf_end1 Arc is approaching the bottom or top boundary at 
     *                 \c asympt_x1
     * \param asympt_x2 The x-coordinate of the second asymptote
     * \param inf_end2 Arc is approaching the bottom or top boundary at 
     *                 \c asympt_x2
     * \return The constructed branch
     *
     * \pre asympt_x1 != asympt_x2
     */
    Arc_2 operator()(const X_coordinate_1& asympt_x1, 
                     CGAL::Arr_curve_end inf_end1,
                     const X_coordinate_1& asympt_x2, 
                     CGAL::Arr_curve_end inf_end2,
                     const Curve_analysis_2& c, int arcno) {
        Arc_2 arc(asympt_x1, inf_end1, asympt_x2, inf_end2, c, arcno);
        return arc;
    }
    
    /*!\brief
     * Construct a non-vertical arc with one left- or right-boundary end 
     * and one end that approaches a vertical asymptote (branch III)
     *
     * \param inf_endx Defining whether the arc emanates from the left or right
     *        boundary
     * \param asympt_x The x-coordinate of the asymptote
     * \param inf_endy Arc is approaching the bottom or top boundary at
     *                 asympt_x
     * \return The constructed branch
     */
    Arc_2 operator()(CGAL::Arr_curve_end inf_endx, 
                     const X_coordinate_1& asympt_x,
                     CGAL::Arr_curve_end inf_endy, 
                     const Curve_analysis_2& c, int arcno) {
        Arc_2 arc(inf_endx, asympt_x, inf_endy, c, arcno);
        return arc;
    }
    
    //!@}
    //!\name Constructing vertical arcs
    //!@{
    
    /*!\brief 
     * Constructs a vertical arc with two interior end-points 
     * (vertical segment)
     *
     * \param p The first end-point
     * \param q The second end-point
     * \param c The supporting curve
     * \return The constructed arc
     * 
     * \pre p != q && p.x() == q.x()
     * \pre c must have a vertical component at this x
     */
    Arc_2 operator()(const Point_2& p, const Point_2& q, 
                     const Curve_analysis_2& c) {
        Arc_2 arc(p,q,c);
        return arc;
    }
    
    /*!\brief
     * Constructs a vertical arc with one interior end-point and 
     * one that reaches the bottom or top boundary (vertical ray)
     *
     * \param origin The interior end-point
     * \param inf_end Ray emanates from bottom or top boundary
     * \return The constructed ray
     *
     * \pre c must have a vertical line component at this x
     */
    Arc_2 operator()(const Point_2& origin, CGAL::Arr_curve_end inf_end,
                     const Curve_analysis_2& c) {
        
        Arc_2 arc(origin, inf_end, c);
        return arc;
    }
    
    /*!\brief
     * Constructs a vertical arc that connects bottom with top boundary
     * (vertical branch)
     *
     * \param x The x-coordinate of the branch
     * \return The constructed branch
     * 
     * \pre c must have a vertical line component at this x
     */
    Arc_2 operator()(const X_coordinate_1& x, const Curve_analysis_2& c) {
        Arc_2 arc(x, c);
        return arc;
    }
};

/*!\brief
 * Functor that checks whether a given arc is vertical
 */
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
    
    //! the arity of the functor
    typedef Arity_tag<1> Arity;
    
    /*!\brief 
     * Standard constructor
     *
     * \param kernel The kernel
     */
    Is_vertical_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }

    /*!\brief
     * Check whether the given x-monotone arc \c cv is a vertical segment.
     *
     * \param cv The curve.
     * \return \c true if the curve is a vertical segment, \c false otherwise.
     */
    result_type operator()(const Arc_2& cv) const {
        return cv.is_vertical();
    }
};

/*!\brief
 * Functor that checks whether a given arc end is interior
 */
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
    
    //! the arity of the functor
    typedef Arity_tag<2> Arity;

    Is_bounded_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }
    
    /*! Is the end of an arc interior.
     *
     * \param cv The arc
     * \param ce The end of cv identifier.
     * \return \c true is the arc end is interior, and \c false otherwise
     */
    result_type operator()(const Arc_2& cv, Arr_curve_end ce) const {
        return (cv.is_finite(ce));
    }
};


/*!\brief
 * Functor computing parameter space in x for arc
 */
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

    //! the arity of the functor
    typedef Arity_tag<2> Arity;
    
    /*!\brief 
     * Standard constructor
     *
     * \param kernel The kernel
     */
    Parameter_space_in_x_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }
    
    /*! Obtains the parameter space in x at the end of an arc.
     *
     * \param cv The arc
     * \param ce the arc end indicator:
     *     ARR_MIN_END - the minimal end of cv
     *     ARR_MAX_END - the maximal end of cv
     * \return the parameter space at the \c ce end of \c cv
     *   ARR_LEFT_BOUNDARY  - the arc approaches the left boundary of the 
     *                        parameter space
     *   ARR_INTERIOR       - the arc does not approach a boundary of the
     *                        parameter space
     *   ARR_RIGHT_BOUNDARY - the arc approaches the right boundary of the
     *                        parameter space
     */
    result_type operator()(const Arc_2& cv, CGAL::Arr_curve_end ce) const {

        CGAL::Arr_parameter_space loc = cv.location(ce);
        if(loc == CGAL::ARR_LEFT_BOUNDARY || loc == CGAL::ARR_RIGHT_BOUNDARY)
            return loc;
        return CGAL::ARR_INTERIOR;
    }
};

// Functor computing parameter space in y for arc
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

    //! the arity of the functor
    typedef Arity_tag<2> Arity;
    
    /*!\brief 
     * Standard constructor
     *
     * \param kernel The kernel
     */
    Parameter_space_in_y_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }

    /*! Obtains the parameter space in y at the end of an arc.
     * 
     * \param cv The arc
     * \param ce the arc end indicator:
     *     ARR_MIN_END - the minimal end of cv
     *     ARR_MAX_END - the maximal end of cv
     * \return the parameter space at the \c ce end of \c cv
     *   ARR_BOTTOM_BOUNDARY- the arc approaches the bottom boundary of the 
     *                        parameter space
     *   ARR_INTERIOR       - the arc does not approach a boundary of the
     *                        parameter space
     *   ARR_TOP_BOUNDARY   - the arc approaches the top boundary of the
     *                        parameter space
     */
    result_type operator()(const Arc_2& cv, CGAL::Arr_curve_end ce) const {
            
        CGAL::Arr_parameter_space loc = cv.location(ce);
        if(loc == CGAL::ARR_BOTTOM_BOUNDARY || loc == CGAL::ARR_TOP_BOUNDARY)
            return loc;
        return CGAL::ARR_INTERIOR;
    }
};

/*!\brief
 * Functor constructing minimum point of an arc (if interior)
 */
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

    //! the arity of the functor
    typedef Arity_tag<1> Arity;
    
    /*!\brief 
     * Standard constructor
     *
     * \param kernel The kernel
     */
    Construct_min_vertex_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }

    /*!\brief
     * Get the minimum end-point of the arc
     *
     * \param cv The arc.
     * \return The minimum end-point.
     *
     * \pre minimum end-point is interior
     */
    result_type operator()(const Arc_2& cv) const {
    
        return cv.curve_end(CGAL::ARR_MIN_END);
    }
};

/*!\brief
 * Functor constructing maximum point of an arc (if interior)
 */
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

    //! the arity of the functor
    typedef Arity_tag<1> Arity;
    
    /*!\brief 
     * Standard constructor
     *
     * \param kernel The kernel
     */
    Construct_max_vertex_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }

    /*!\brief
     * Get the maximum end-point of the arc.
     *
     * \param cv The arc.
     * \return The maximum end-point.
     *
     * \pre maximum end-point is interior
     */
    result_type operator()(const Arc_2& cv) const {
        
        return cv.curve_end(CGAL::ARR_MAX_END);
    }
};

/*!\brief
 * Functor that compares x-coordinates of two interior points
 */
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

    //! the arity of the functor
    typedef Arity_tag<2>            Arity;
    
    /*!\brief 
     * Standard constructor
     *
     * \param kernel The kernel
     */
    Compare_x_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }
    
    /*!\brief
     * Compare the x-coordinates of two points.
     *
     * \param p1 The first point.
     * \param p2 The second point.
     * \return CGAL::LARGER if x(p1) > x(p2);
     *         CGAL::SMALLER if x(p1) \< x(p2);
     *         CGAL::EQUAL if x(p1) = x(p2).
     */
    result_type operator()(const Point_2 &p1, const Point_2 &p2) const {
        return Curved_kernel_via_analysis_2::instance().
            kernel().compare_x_2_object()
            (p1.x(), p2.x());
    }
};

/*!\brief
 * Functor that compares coordinates of two interior points lexicographically
 */
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

    //! the arity of the functor
    typedef Arity_tag<2>            Arity;
    
    /*!\brief 
     * Standard constructor
     *
     * \param kernel The kernel
     */
    Compare_xy_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }

    /*!\brief 
     * Compares two points lexigoraphically: by x, then by y.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return CGAL::LARGER if x(p1) > x(p2), 
     *                      or if x(p1) = x(p2) and y(p1) > y(p2);
     *         CGAL::SMALLER if x(p1) \< x(p2), 
     *                       or if x(p1) = x(p2) and y(p1) \< y(p2);
     *         CGAL::EQUAL if the two points are equal.
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


/*!\brief
 * Functor that compares x-coordinates near the top or bottom boundary
 */
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

    //! the arity of the functor
    typedef Arity_tag<4>            Arity;

    Compare_x_near_boundary_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }

    /*!\brief Compare the x-coordinate of a point with the x-coordinate of
     * an arcend near the boundary at bottom or top boundary
     * 
     * \param p the point direction.
     * \param cv the arc, the endpoint of which is compared.
     * \param ce the arc-end indicator -
     *            ARR_MIN_END - the minimal end of cv or
     *            ARR_MAX_END - the maximal end of cv.
     * \return the comparison result:
     *         SMALLER - x(p) \< x(cv, ce);
     *         EQUAL   - x(p) = x(cv, ce);
     *         LARGER  - x(p) > x(cv, ce).
     *
     * \pre p lies in the interior of the parameter space.
     * \pre the ce end of the line cv lies on a boundary.
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

    /*! Compare the x-coordinates of 2 arcs ends near the top or bottom 
     * boundary of the parameter space
     * \param cv1 the first arc.
     * \param ce1 the first arc end indicator -
     *            ARR_MIN_END - the minimal end of cv1 or
     *            ARR_MAX_END - the maximal end of cv1.
     * \param cv2 the second arc.
     * \param ce2 the second arc end indicator -
     *            ARR_MIN_END - the minimal end of cv2 or
     *            ARR_MAX_END - the maximal end of cv2.
     * \return the second comparison result:
     *         SMALLER - x(cv1, ce1) \< x(cv2, ce2);
     *         EQUAL   - x(cv1, ce1) = x(cv2, ce2);
     *         LARGER  - x(cv1, ce1) > x(cv2, ce2).
     *
     * \pre the ce1 end of the arc cv1 lies on a boundary.
     * \pre the ce2 end of the arc cv2 lies on a boundary.
     */
    result_type operator()(const Arc_2& cv1, CGAL::Arr_curve_end ce1,
                           const Arc_2& cv2, CGAL::Arr_curve_end ce2) const {

        CERR("\ncompare_x_near_boundary: cv1: " << cv1 << "\n cv2: " <<
            cv2 << "; end1: " << ce1 << "; end2: " << ce2 << "\n");
        /*CGAL::Arr_boundary_type bnd1 = boundary(end1), 
            bnd2 = cv2.boundary(ce2);*/
        CGAL::Arr_parameter_space 
            loc1 = cv1.location(ce1), 
            loc2 = cv2.location(ce2);
        CGAL_precondition(cv1.is_on_bottom_top(loc1));
        CGAL_precondition(cv1.is_on_bottom_top(loc2));
        
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


/*!\brief
 * Functor that compares ends of arcs near left or right boundary
 */
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

    //! the arity of the functor
    typedef Arity_tag<3>            Arity;
    
    /*!\brief 
     * Standard constructor
     *
     * \param kernel The kernel
     */
    Compare_y_near_boundary_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }

    /*!\brief
     * Compare the y-coordinates of 2 arcs at their ends near the left
     * or right boundary of the parameter space
     * 
     * \param cv1 the first arc.
     * \param cv2 the second arc.
     * \param ce the arc end indicator.
     * \return CGAL::SMALLER is y(cv1,ce) < y(cv2,ce2) near left/right boundary
     *         CGAL::EQUAL is y(cv1,ce) = y(cv2,ce2) near left/right boundary,
     *         i.e., they overlap
     * \return CGAL::LARGER is y(cv1,ce) > y(cv2,ce2) near left/right boundary
     *
     * \pre the ce ends of the arcs cv1 and cv2 lie either on the left
     * boundary or on the right boundary of the parameter space.
     */
    result_type operator()(const Arc_2& cv1, const Arc_2& cv2,
                           CGAL::Arr_curve_end ce) const {
        
        CERR("\ncompare_y_near_boundary; cv1: " << cv1 << "; cv2: " <<
             cv2 << "; end: " << ce << "\n");
        
        CGAL::Arr_parameter_space loc1 = cv1.location(ce);
        CGAL_precondition(cv1.is_on_left_right(loc1));
        CGAL_precondition(loc1 == cv2.location(ce));
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


/*!\brief
 * Functor that computes relative vertical alignment of an interior point and 
 * an arc
 */
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

    //! the arity of the functor
    typedef Arity_tag<2>            Arity;
    
    /*!\brief 
     * Standard constructor
     *
     * \param kernel The kernel
     */
    Compare_y_at_x_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }

    /*!\brief
     * Return the relative vertical alignment of a point with an arc
     * 
     * \param p The point
     * \param cv The arc
     * \return CGAL::SMALLER if y(p) \< cv(x(p)), i.e., 
     *                       the point is below the arc;
     *         CGAL::LARGER if y(p) > cv(x(p)), i.e.,
     *                      the point is above the arc;
     *         CGAL::EQUAL if p lies on the arc
     *
     * \pre p is in the x-range of cv.
     *
     */
    result_type operator()(const Point_2& p, const Arc_2& cv) const {
        
        CERR("\ncompare_y_at_x; p: " << p << ";\n cv:" << cv << "\n");

        bool eq_min, eq_max;
        CGAL_assertion_code (
           bool in_x_range =
        )
        cv.is_in_x_range(p.x(), &eq_min, &eq_max);
        CGAL_assertion(in_x_range);

        if (cv.is_vertical()) {
            if (cv.is_finite(CGAL::ARR_MIN_END)) {
                // for vertical arcs we can ask for .xy() member
                if (Curved_kernel_via_analysis_2::instance().kernel().
                    compare_xy_2_object()(
                            p.xy(), cv._minpoint().xy(), true
                    ) == CGAL::SMALLER) {
                    return CGAL::SMALLER;
                }
            }
            if (cv.is_finite(CGAL::ARR_MAX_END)) {
                if (Curved_kernel_via_analysis_2::instance().kernel().
                    compare_xy_2_object()(
                            p.xy(), cv._maxpoint().xy(), true
                    ) == CGAL::LARGER) {
                    return CGAL::LARGER;
                }
            }
            return CGAL::EQUAL; // p lies on a vertical arc
        }
        CGAL::Comparison_result res;
        if (eq_min) {
            res = Curved_kernel_via_analysis_2::instance().kernel().
                compare_xy_2_object()(
                        p.xy(), cv._minpoint().xy(), true
                );
        } else if (eq_max) {
            res = Curved_kernel_via_analysis_2::instance().kernel().
                compare_xy_2_object()(
                        p.xy(), cv._maxpoint().xy(), true
                );
        } else {
            Point_2 point_on_s
                = Curved_kernel_via_analysis_2::instance().
                construct_point_on_arc_2_object()
                ( p.x(), 
                  cv.curve(), 
                  cv.arcno(),
                  cv );
            res = Curved_kernel_via_analysis_2::instance().kernel().
                compare_xy_2_object()(p.xy(), point_on_s.xy(), true);
        }
        CERR("cmp result: " << res << "\n");
        return res;
    }
};

/*!\brief
 * Functor that computes the relative vertical aligment of two arcs left
 * of a point
 */
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

    //! the arity of the functor
    typedef Arity_tag<3>            Arity;
    
    /*!\brief 
     * Standard constructor
     *
     * \param kernel The kernel
     */
    Compare_y_at_x_left_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }
    
    /*!\brief
     * Compares the relative vertical alignment of two arcs 
     * immediately to the left of one of their intersection points. 
     * 
     * If one of the arcs is vertical (emanating downward from p), it is 
     * always considered to be below the other curve.
     *
     * \param cv1 The first arc
     * \param cv2 The second arc
     * \param p The intersection point.
     * \return The relative vertical alignment of cv1 with respect to cv2 
     *         immediately to the left of p: CGAL::SMALLER, CGAL::LARGER or 
               CGAL::EQUAL.
     *
     * \pre The point p lies on both curves, and both of them must be also be
     *      defined (lexicographically) to its left.
     */
    result_type operator() (const Arc_2& cv1, const Arc_2& cv2,
                            const Point_2& p) const {

        CERR("\ncompare_y_at_x_left(cv2); cv1: " << cv1 << "; cv2: " <<
            cv2 << "; p: " << p << "\n");

        // ensure that p lies on both arcs 
        CGAL_precondition(cv1.compare_y_at_x(p) == CGAL::EQUAL);
        CGAL_precondition(cv2.compare_y_at_x(p) == CGAL::EQUAL);
        
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


/*!\brief
 * Functor that computes the relative vertical aligment of two arcs right
 * of a point
 */
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

    //! the arity of the functor
    typedef Arity_tag<3>            Arity;
    
    /*!\brief 
     * Standard constructor
     *
     * \param kernel The kernel
     */
    Compare_y_at_x_right_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }

    /*!\brief
     * Compares the relative vertical alignment of two arcs 
     * immediately to the right of one of their intersection points. 
     * 
     * If one of the arcs is vertical (emanating downward from p), it is 
     * always considered to be below the other curve.
     *
     * \param cv1 The first arc
     * \param cv2 The second arc
     * \param p The intersection point.
     * \return The relative vertical alignment of cv1 with respect to cv2 
     *         immediately to the right of p: CGAL::SMALLER, CGAL::LARGER or 
               CGAL::EQUAL.
     *
     * \pre The point p lies on both curves, and both of them must be also be
     *      defined (lexicographically) to its right.
     */
    result_type operator()(const Arc_2& cv1, const Arc_2& cv2,
                           const Point_2& p) const {

        CERR("\ncompare_y_at_x_right; cv1: " << cv1 << "; cv2: " <<
            cv2 << "; p: " << p << "\n");
        
        // ensure that p lies on both arcs and doesn't lie on the positive
        // boundary
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

/*!\brief
 * Functor that checks whether a point is in the x-range of an arc
 */
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

    //! the arity of the functor
    typedef Arity_tag<2>            Arity;
    
    /*!\brief 
     * Standard constructor
     *
     * \param kernel The kernel
     */
    Is_in_x_range_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }
    
    /*!\brief
     * Check whether a given point lies within the curve's x-range
     * 
     * \param cv The arc
     * \param p the point
     * \return \c true if p lies in arc's x-range; \c false otherwise.
     */
    bool operator()(const Arc_2& cv, const Point_2& p) const {
        return cv.is_in_x_range(p);
    }
};

/*!\brief
 * Tests two objects, whether they are equal
 */
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

    //! the arity of the functor
    typedef Arity_tag<2> Arity;
    
    /*!\brief 
     * Standard constructor
     *
     * \param kernel The kernel
     */
    Equal_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }
    
    /*!\brief
     * Check if the two points are the same
     *
     * \param p1 The first point.
     * \param p2 The second point.
     * \return \c true if the two point are the same; \c false otherwise.
     */
    result_type operator()(const Point_2& p1, const Point_2& p2) const {
        return (Curved_kernel_via_analysis_2::instance().
                compare_xy_2_object()(p1, p2) == 
                CGAL::EQUAL);
    }
     
    /*!\brief 
     * Check if the two arcs are the same 
     * 
     * \param cv1 The first arc
     * \param cv2 The second arc
     * \return \c true if the two curves are the same; \c false otherwise.
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


/*!\brief
 * Functor that checks whether two arcs overlap
 */
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

    //! the arity of the functor
    typedef Arity_tag<2> Arity;
    
    /*!\brief 
     * Standard constructor
     *
     * \param kernel The kernel
     */
    Do_overlap_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }
    
    /*!\brief
     * Check whether two given arcs overlap, i.e., they have infinitely
     * many intersection points
     *
     * \param cv1 The first arc
     * \param cv2 The second arc
     * \return \c true if the curves overlap; \c false otherwise
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


/*!\brief
 * Functor that computes the intersections of two arcs
 */
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

    //! the arity of the functor
    typedef Arity_tag<3> Arity;    
    
    /*!\brief 
     * Standard constructor
     *
     * \param kernel The kernel
     */
    Intersect_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }
    
    /*!\brief
     * Find all intersections of the two given arcs and insert them to the 
     * output iterator. 
     *
     * Type of output iterator is \c CGAL::Object 
     * containing either an \c Arc_2 object (overlap) or a \c 
     * std::pair< Point_2, unsigned int >, where the unsigned int denotes
     * the multiplicity of the zero-dimensional intersection (0 if unknown)
     *
     * \param cv1 The first arc
     * \param cv2 The second arc
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

        for (it = vec.begin(); it != vec.end(); it++) {
            *oi++ = CGAL::make_object(*it);
        }
        return oi;
    }

};


/*!\brief
 * Functors that trims an arc
 */
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

    //! the arity of the functor
    typedef Arity_tag<3> Arity;
    
    /*!\brief 
     * Standard constructor
     *
     * \param kernel The kernel
     */
    Trim_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }
    
    /*!\brief
     * Returns a trimmed version of an arc
     * 
     * \param cv The arc
     * \param p the new first endpoint
     * \param q the new second endpoint
     * \return The trimmed arc
     *
     * \pre p != q
     * \pre both points must be interior and must lie on \c cv
     */
    Arc_2 operator()(const Arc_2& cv, const Point_2& p, const Point_2& q) {
    
        CERR("trim\n");

        CGAL_precondition(p.location() == CGAL::ARR_INTERIOR);
        CGAL_precondition(q.location() == CGAL::ARR_INTERIOR);
        
        CGAL_precondition(
                !Curved_kernel_via_analysis_2::instance().
                equal_2_object()(p, q)
        );
        CGAL_precondition(cv.compare_y_at_x(p) == CGAL::EQUAL);
        CGAL_precondition(cv.compare_y_at_x(q) == CGAL::EQUAL);  

        return cv._trim(p, q);
    }
};


/*!\brief
 * Functor that splits a arc at an interior point
 */
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

    //! the arity of the functor
    typedef Arity_tag<4> Arity;    
    
    /*!\brief 
     * Standard constructor
     *
     * \param kernel The kernel
     */
    Split_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }
    
    /*!\brief
     * Split a given arc at a given point into two sub-arcs
     * 
     * \param cv The arc to split
     * \param p The split point
     * \param c1 Output: The left resulting subcurve (p is its right endpoint)
     * \param c2 Output: The right resulting subcurve (p is its left endpoint)
     * 
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
        ).first;
        c2 = cv._replace_endpoints(
                p, cv._maxpoint(), (cv.is_vertical() ? -1 : cv.arcno()), -1
        ).first;
    }
};

/*!\brief
 * Functor that computes whether two arcs are mergeable
 */
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

    //! the arity of the functor
    typedef Arity_tag<2> Arity;    
    
    /*!\brief 
     * Standard constructor
     *
     * \param kernel The kernel
     */
    Are_mergeable_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }
    
    /*!\brief
     * Check whether two given arcs are mergeable
     *
     * \param cv1 The first arc
     * \param cv2 The second arc
     * \return \c true if the two arcs are mergeable, i.e., they are supported
     * by the same curve and share a common endpoint; \c false otherwise.
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

/*!\brief
 * Functor that merges two arcs
 */
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

    //! the arity of the functor
    typedef Arity_tag<2> Arity;    
    
    /*!\brief 
     * Standard constructor
     *
     * \param kernel The kernel
     */
    Merge_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }

    /*!\brief
     * Merges two given arcs into a single one
     *
     * \param cv1 The first arc
     * \param cv2 The second arc
     * \param c Output: The resulting arc
     * 
     * \pre The two arcs are mergeable, that is they are supported by the
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
        Arc_2 arc = cv1._replace_endpoints(src, tgt, arcno_s, arcno_t).first;
        // arc.set_boundaries_after_merge(*this, s); - no need to, since
        // boundaries are stored in Point_2 type and will be copied implicitly
        
        c = arc;
    }
};


/*!\brief
 * Functor that computes whether a point lies on a supporting curve
 */
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

    //! the arity of the functor
    typedef Arity_tag<2> Arity;
    
    /*!\brief 
     * Standard constructor
     *
     * \param kernel The kernel
     */
    Is_on_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }
    
    /*!\brief
     * Checks whether \c p lies on \c c 
     * 
     * \param p The point to test
     * \param c The curve
     * \return \c true if the \c p lies on \c c, \c false otherwise
     */
    result_type operator()(const Point_2& p, const Curve_analysis_2& c) const {
        result_type res = 
            (Curved_kernel_via_analysis_2::instance().
             kernel().sign_at_2_object()(c, p.xy())
             == CGAL::ZERO);
        return res;
    }
};


/*!\brief 
 * Functor that decomposes curve into x-monotone arcs and isolated points
 */
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

    //! the arity of the functor
    typedef Arity_tag<2> Arity;   
    
    /*!\brief 
     * Standard constructor
     *
     * \param kernel The kernel
     */
    Make_x_monotone_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    } 

    /*!\brief
     * Decomposes a given arc into list of x-monotone arcs 
     * (subcurves) and insert them to the output iterator. Since \c Arc_2 
     * is by definition x-monotone, an input arc is passed to the 
     * output iterator directly. 
     *
     * \param cv The arc
     * \param oi The output iterator, whose value-type is Object
     * The returned objects are all wrappers Arc_2 objects
     * \return The past-the-end iterator
     */
    template < class OutputIterator >
    OutputIterator operator()(const Arc_2& cv, OutputIterator oi) const {
    
        *oi++ = CGAL::make_object(cv);
        return oi;
    }
    
    /*!\brief 
     * Decomposes a given curve into list of x-monotone arcs
     * (subcurves) and isolated points and insert them to the output iterator. 
     * 
     * \param cv The curve
     * \param oi The output iterator, whose value-type is Object.
     * The returned objects either wrapper Arc_2 or Point_2 objects
     * \return The past-the-end iterator
     */
    template < class OutputIterator >
    OutputIterator operator()(const Curve_2& cv, OutputIterator oi) const {
    
        CGAL::CGALi::Make_x_monotone_2< Curved_kernel_via_analysis_2 >
            make_x_monotone(&Curved_kernel_via_analysis_2::instance());
        return make_x_monotone(cv, oi);
    }
};

#undef CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES 

} // namespace Curved_kernel_via_analysis_2_Functors

/*!\brief
 * Collects main set of functors of a curved kernel
 */
template < 
    class CurvedKernelViaAnalysis_2, 
    class Curve_2_, 
    class Point_2_, 
    class Arc_2_
>
class Curved_kernel_via_analysis_2_functors {

public:
    //!\name Public types
    //!@{
    
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
    /*!\brief functor */ \
    typedef \
    Curved_kernel_via_analysis_2_Functors::Y< Curved_kernel_via_analysis_2 > \
    Y; \
    /*! returns instance of functor */ \
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
// EOF
