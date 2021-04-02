// Copyright (c) 2007,2008,2009,2010,2011 Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//                 Pavel Emeliyanenko <asm@mpi-sb.mpg.de>

#ifndef CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_FUNCTORS_H
#define CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_FUNCTORS_H

/*!\file include/CGAL/Curved_kernel_via_analysis_2/Curved_kernel_via_analysis_2_functors.h
 * \brief defines Curved_kernel_via_analysis_2 function objects + class
 */

#include <CGAL/config.h>
#include <CGAL/Curved_kernel_via_analysis_2/Make_x_monotone_2.h>
#include <CGAL/iterator.h>
namespace CGAL {

namespace internal {

#ifndef CERR
//#define CKvA_DEBUG_PRINT_CERR
#ifdef CKvA_DEBUG_PRINT_CERR
#define CERR(x) std::cerr << x
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
    typedef typename Point_2::Coordinate_1 Coordinate_1;
    //typedef typename Point_2::Y_coordinate_1 Y_coordinate_1;


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
        CGAL_precondition(kernel != nullptr);
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
    typedef typename Base::Point_2 Point_2; \
    typedef typename Base::Arc_2 Arc_2; \
    typedef typename Base::Curve_analysis_2 Curve_analysis_2; \
    typedef typename Base::Coordinate_1 Coordinate_1; \
    //typedef typename Base::Y_coordinate_1 Y_coordinate_1;

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

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES

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
     * Constructs an interior point
     *
     * \param x x-coordinate
     * \param c The supporting curve
     * \param arcno Arcnumber on curve
     * \return The constructed point
     */
    Point_2 operator()(const Coordinate_1& x,
                       const Curve_analysis_2& c,
                       int arcno) {
        Point_2 pt(x, c, arcno);
        return pt;
    }

    Point_2 operator()(const Coordinate_1& x,
                       const Arc_2& a) {
        CGAL_precondition(a.is_in_x_range(x));
        Point_2 pt(x, a.curve(), a.arcno(x));
        return pt;
    }


    template <typename T>
    Point_2 operator() (const T& x,
                        const T& y) {
        Point_2 pt(x,y);
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

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES

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
            const Coordinate_1& x,
            const Curve_analysis_2& c, int arcno, const Arc_2& arc
    ) {

        // avoid compiler warning
        (void)arc;

        //CGAL::set_pretty_mode(std::cerr);
        CERR("Construct_pt_on_arc: " << CGAL::to_double(x) << ", " << arcno <<
             ", " << c.id() <<  "\narc = " << arc << "\n");

        //CGAL_assertion(c.id() == arc.curve().id());
        //CGAL_assertion(arcno == arc.arcno(x));

        Point_2 pt = Base::_ckva()->construct_point_2_object()(x, c, arcno);

        // here we can modify the point wrt "data stored in arc",
        // if we want to
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

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES

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
    Arc_2 operator()(const Point_2& origin, const Coordinate_1& asympt_x,
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
    Arc_2 operator()(const Coordinate_1& asympt_x1,
                     CGAL::Arr_curve_end inf_end1,
                     const Coordinate_1& asympt_x2,
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
                     const Coordinate_1& asympt_x,
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
    Arc_2 operator()(const Coordinate_1& x, const Curve_analysis_2& c) {
        Arc_2 arc(x, c);
        return arc;
    }

    //!@}
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

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES

    //! the result type
    typedef bool result_type;

    //! the arity of the functor

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

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES

    //! the result type
    typedef Point_2 result_type;

    //! the arity of the functor

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

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES

    //! the result type
    typedef Point_2 result_type;

    //! the arity of the functor

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
 * Functor constructing an interior point of on an arc.
 */
template < class CurvedKernelViaAnalysis_2 >
class Construct_interior_vertex_2 : public
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES

    //! the result type
    typedef Point_2 result_type;

    /*!\brief
     * Standard constructor
     *
     * \param kernel The kernel
     */
    Construct_interior_vertex_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }

    /*!\brief
     * Get an interior point on the arc.
     *
     * \param arc The arc.
     * \return A point on the interior of the arc.
     */
  result_type operator()(const Arc_2& arc) const {
    Point_2 p = compute_interior_vertex(arc);
    CGAL_postcondition(this->_ckva()->is_on_2_object()(p,arc));
    return p;
  }


      private:

    result_type compute_interior_vertex(const Arc_2& arc) const {

      typedef typename Curved_kernel_via_analysis_2::Curve_2 Curve_analysis_2;

      typedef typename Curved_kernel_via_analysis_2::Curve_kernel_2
                                                      Algebraic_curve_kernel_2;
      typedef typename Algebraic_curve_kernel_2::Bound Bound;
      typedef typename Algebraic_curve_kernel_2::Coordinate_1 Coordinate_1;

      typedef CGAL::Polynomial< Bound > Poly_rat_1;
      typedef CGAL::Polynomial< Poly_rat_1 > Poly_rat_2;

      typedef CGAL::Polynomial_traits_d< Poly_rat_1 > PT_rat_1;
      typedef CGAL::Polynomial_traits_d< Poly_rat_2 > PT_rat_2;

      typedef CGAL::Fraction_traits< Poly_rat_2 > FT_2;

      if (!arc.is_vertical())
      {
        Bound x_coord = arc.boundary_in_x_range_interior();
        int arcno = arc.arcno();
        const Curve_analysis_2& ca = arc.curve();

        Point_2 p = Curved_kernel_via_analysis_2::instance().
          construct_point_on_arc_2_object()
          (Coordinate_1(x_coord), ca, arcno, arc);
        return p;
      }

      Bound y_coord = 0;
      if (arc.is_finite(ARR_MIN_END))
      {
        if (arc.is_finite(ARR_MAX_END))
        {
          // We need torefine the interval because there is a chance that
          // the low of the upper point is below the high of the lower point.
          y_coord = Curved_kernel_via_analysis_2::instance().kernel().
            bound_between_y_2_object() (arc.curve_end(ARR_MIN_END).xy(),
                                         arc.curve_end(ARR_MAX_END).xy());
        }
        else
        {
          std::pair<Bound,Bound> approx_pair =
            Curved_kernel_via_analysis_2::instance().kernel().
            approximate_relative_y_2_object()
            (arc.curve_end(ARR_MIN_END).xy(),4);
          y_coord = approx_pair.second + Bound(1);

        }
      }
      else
      {
        if (arc.is_finite(ARR_MAX_END))
        {
          std::pair<Bound,Bound> approx_pair =
            Curved_kernel_via_analysis_2::instance().kernel().
            approximate_relative_y_2_object()
            (arc.curve_end(ARR_MAX_END).xy(),4);
          y_coord = approx_pair.first-Bound(1);
        }
      }

      /*! \todo Try to remove this polynomial stuff */
      typename PT_rat_1::Construct_polynomial cp1;
      Poly_rat_2 poly2 = typename PT_rat_2::Construct_polynomial()
        (cp1(-y_coord), cp1(Bound(1)));

      typename FT_2::Denominator_type dummy;
      typename FT_2::Numerator_type curve_poly;
      typename FT_2::Decompose() (poly2, curve_poly, dummy);

      Curve_analysis_2 curve = Curved_kernel_via_analysis_2::instance().
        kernel().construct_curve_2_object()(curve_poly);
      Point_2 p =  Curved_kernel_via_analysis_2::instance().
        construct_point_on_arc_2_object()(arc.x(), curve, 0, arc);
      return p;
    }
};


/*!\brief
 * Functor that compares x-coordinates of two interior points
 */
template < class CurvedKernelViaAnalysis_2 >
class Compare_x_2 : public
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES

    //! the result type
    typedef CGAL::Comparison_result result_type;

    //! the arity of the functor

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
            kernel().compare_1_object()
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

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES

    //! the result type
    typedef CGAL::Comparison_result result_type;

    //! the arity of the functor

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
        // TODO add CGAL_precondition(p1.location() == CGAL::ARR_INTERIOR);
        // TODO add CGAL_precondition(p2.location() == CGAL::ARR_INTERIOR);

        CERR("\ncompare_xy; p1: " << p1
             << ";\n p2:" << p2 << "");

        if (p1.id() == p2.id()) {
            result_type res = CGAL::EQUAL;
            CERR("result: " << res << "\n");
            return res;
        }

        result_type res =
            (Curved_kernel_via_analysis_2::instance().
             kernel().compare_xy_2_object()
             (p1.xy(), p2.xy(), equal_x));

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

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES

    //! the result type
    typedef CGAL::Comparison_result result_type;

    //! the arity of the functor

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

        CGAL_precondition(p.is_finite());

        CERR("\ncompare_y_at_x; p: " << p << ";\n cv:" << cv << "\n");

        bool eq_min, eq_max;
        CGAL_assertion_code (
           bool in_x_range =
        )
        cv.is_in_x_range(p.x(), &eq_min, &eq_max);
        CGAL_assertion(in_x_range);

        // TODO replace with this->compare_xy_2_object()?
        typename Base::Curve_kernel_2::Compare_xy_2 cmp_xy(
            Base::_ckva()->kernel().compare_xy_2_object());

        if (cv.is_vertical()) {

            if (cv.is_finite(CGAL::ARR_MIN_END)) {

                // for vertical arcs we can ask for .xy() member
                if (cmp_xy(p.xy(), cv._minpoint().xy(), true) ==
                         CGAL::SMALLER) {
                    CERR("cmp result: " << CGAL::SMALLER << "\n");
                    return CGAL::SMALLER;
                }
            }

            if (cv.is_finite(CGAL::ARR_MAX_END)) {
                if (cmp_xy(p.xy(), cv._maxpoint().xy(), true) ==
                    CGAL::LARGER) {
                    CERR("cmp result: " << CGAL::LARGER << "\n");
                    return CGAL::LARGER;
                }
            }
            CERR("cmp result: " << CGAL::EQUAL << "\n");
            return CGAL::EQUAL; // p lies on a vertical arc
        }
        CGAL::Comparison_result res;
        if (eq_min) {
          if(cv._minpoint().is_finite()){
            res = cmp_xy(p.xy(), cv._minpoint().xy(), true);
          }else{
            res = (cv.location(CGAL::ARR_MIN_END)==ARR_TOP_BOUNDARY)?
              CGAL::SMALLER : CGAL::LARGER;
          }
        } else if (eq_max) {
          CGAL_precondition(p.is_finite());
          if(cv._maxpoint().is_finite()){
            res = cmp_xy(p.xy(), cv._maxpoint().xy(), true);
          }else{
            res = (cv.location(CGAL::ARR_MAX_END)==ARR_TOP_BOUNDARY)?
              CGAL::SMALLER : CGAL::LARGER;
          }
        } else {
          Point_2 point_on_s =
            Base::_ckva()->construct_point_on_arc_2_object()
            (p.x(), cv.curve(), cv.arcno(), cv );

          CGAL_precondition(p.is_finite());
          CGAL_precondition(point_on_s.is_finite());
            res = cmp_xy(p.xy(), point_on_s.xy(), true);
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

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES

    //! the result type
    typedef CGAL::Comparison_result result_type;

    //! the arity of the functor

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

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES

    //! the result type
    typedef CGAL::Comparison_result result_type;

    //! the arity of the functor

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

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES

    //! the result type
    typedef bool result_type;

    //! the arity of the functor

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

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES

    //! the result type
    typedef bool result_type;

    //! the arity of the functor

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

        Arc_2::simplify(cv1,cv2);
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

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES

    //! the result type
    typedef bool result_type;

    //! the arity of the functor

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
                kernel().compare_1_object()(
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
class Intersect_2 :
    public Curved_kernel_via_analysis_2_functor_base<CurvedKernelViaAnalysis_2>
{
public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
      Base;

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES

    typedef unsigned int                              Multiplicity;
    typedef std::pair<Point_2, Multiplicity>          Intersection_point;
    typedef boost::variant<Intersection_point, Arc_2> Intersection_result;

    //! the result type
    typedef CGAL::cpp98::iterator<std::output_iterator_tag, Intersection_result>
      result_type;

    //! the arity of the functor

    /*!\brief
     * Standard constructor
     *
     * \param kernel The kernel
     */
    Intersect_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }

    // TODO add operators for non-x-monotone arcs + curves (i.e., all combis)

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
        std::vector<Arc_2> arcs;
        if (cv1._trim_if_overlapped(cv2, std::back_inserter(arcs))) {
            for (const auto& item : arcs) *oi++ = Intersection_result(item);
            return oi;
        }
        // process non-ov erlapping case
        std::vector<Intersection_point> vec;
        Arc_2::_intersection_points(cv1, cv2, std::back_inserter(vec));
        for (const auto& item : vec) *oi++ = Intersection_result(item);
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

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES

    //! the result type
    typedef Arc_2 result_type;

    //! the arity of the functor

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

        CGAL_precondition(!Base::_ckva()->equal_2_object()(p, q));
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

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES

    //! the result type
    typedef void result_type;

    //! the arity of the functor

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
        CGAL_precondition(cv._same_arc_compare_xy(cv._minpoint(), p) != CGAL::EQUAL);
        CGAL_precondition(cv._same_arc_compare_xy(cv._maxpoint(), p) != CGAL::EQUAL);

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

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES

    //! the result type
    typedef bool result_type;

    //! the arity of the functor

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

            Curve_analysis_2 ca_2(cv1.curve());
            typename Curve_analysis_2::Status_line_1 cv_line =
                ca_2.status_line_for_x(common.x());
            CGAL_assertion(cv_line.is_event()); // ??

            // we are not allowed to use number_of_incident_branches()
            // since the common point might be supported by different curve,
            // and therefore its arcno might be not valid for *this arc
            for (int k = 0; k < cv_line.number_of_events(); k++) {
                // use a temporary object for comparison predicate
                typename Point_2::Coordinate_2 tmp(common.x(), cv1.curve(), k);
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

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES

    //! the result type
    typedef void result_type;

    //! the arity of the functor

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

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES

    //! the result type
    typedef bool result_type;

    //! the arity of the functor

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

    /*!\brief
     * Checks whether \c p lies on \c cv
     *
     * \param p The point to test
     * \param c The curve arc
     * \return \c true if the \c p lies on \c cv, \c false otherwise
     */
    result_type operator()(const Point_2& p, const Arc_2& cv) const {
      // fix by Michael Hemmer:
      if(cv.is_in_x_range(p.x())){
        return cv.compare_y_at_x(p) == CGAL::EQUAL;
      }
      return false;

      // This is the old code that seems to interpret the arc as open
      // In particular, it returns false for vertical arcs.
      // (Michael Hemmer)
//       bool is_left, is_right;
//       result_type res =
//         (cv.is_in_x_range(p.x(),&is_left,&is_right)) &&
//         !(is_left) &&
//         !(is_right) &&
//         (cv.compare_y_at_x(p) == CGAL::EQUAL);
//       return res;
    }

};


/*!\brief
 * Functor that decomposes curve into x-monotone arcs and isolated points
 */
template <class CurvedKernelViaAnalysis_2>
class Make_x_monotone_2 : public
Curved_kernel_via_analysis_2_functor_base<CurvedKernelViaAnalysis_2> {

public:
  //! this instance' first template parameter
  typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

  //! the base type
  typedef
  Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2>
    Base;

  CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES

  //! the result type
  typedef CGAL::cpp98::iterator< std::output_iterator_tag, CGAL::Object>
    result_type;

  //! the arity of the functor

  /*!\brief
   * Standard constructor
   *
   * \param kernel The kernel
   */
  Make_x_monotone_2(Curved_kernel_via_analysis_2* kernel) :
    Base(kernel)
  {}

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
  template <class OutputIterator>
  OutputIterator operator()(const Arc_2& cv, OutputIterator oi) const
  {
    typedef typename Curved_kernel_via_analysis_2::Point_2      Point_2;
    typedef typename Curved_kernel_via_analysis_2::Arc_2        Arc_2;
    typedef boost::variant<Point_2, Arc_2>      Make_x_monotone_result;
    *oi++ = Make_x_monotone_result(cv);
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
  template <class OutputIterator>
  OutputIterator operator()(const Curve_analysis_2& cv, OutputIterator oi)
    const
  {
    CGAL::internal::Make_x_monotone_2< Curved_kernel_via_analysis_2 >
      make_x_monotone(Base::_ckva());

    return make_x_monotone(cv, oi);
  }

  /*!\brief
   * Splits an input object \c obj into x-monotone arcs and isolated points
   *
   * \param obj the polymorph input object: can represet \c Point_2,
   * \c Arc_2, \c Non_x_monotone_arc_2 or \c Curve_analysis_2
   * \param oi Output iterator that stores CGAL::Object, which either
   *           encapsulates \c Point_2 or \c Arc_2
   * \return Past-the-end iterator of \c oi
   */
  template <class OutputIterator>
  OutputIterator operator()(CGAL::Object obj, OutputIterator oi)
  {
    typedef typename Curved_kernel_via_analysis_2::Point_2      Point_2;
    typedef typename Curved_kernel_via_analysis_2::Arc_2        Arc_2;
    typedef typename Curved_kernel_via_analysis_2::Non_x_monotone_arc_2
      Non_x_monotone_arc_2;
    typedef boost::variant<Point_2, Arc_2>
      Make_x_monotone_result;

    Curve_analysis_2 curve;
    Point_2 p;
    Arc_2 xcv;
    Non_x_monotone_arc_2 nxarc;

    if (CGAL::assign(curve, obj)) oi = (*this)(curve, oi);
    else if (CGAL::assign(nxarc, obj))
      std::cerr << "AU BACKE" << std::endl;
    //oi = std::transform(nxarc.begin(), nxarc.end(), oi,
    //      std::ptr_fun(CGAL::make_object<Arc_2>));
    else if (CGAL::assign(xcv, obj)) *oi++ = Make_x_monotone_result(xcv);
    else if (CGAL::assign(p, obj)) *oi++ = Make_x_monotone_result(p);
    else CGAL_error();
    return oi;
  }
};


// left-right

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

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES

    //! the result type
    typedef CGAL::Arr_parameter_space result_type;

    //! the arity of the functor

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

    result_type operator()(const Point_2& pt) const {
      return pt.location();
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

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES

    //! the result type
    typedef CGAL::Comparison_result result_type;

    //! the arity of the functor

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


      CGAL::Comparison_result res = CGAL::EQUAL;

      CGAL::Arr_parameter_space loc1 = cv1.location(ce);
      CGAL_precondition(cv1.is_on_left_right(loc1));
      CGAL_precondition(loc1 == cv2.location(ce));
      // comparing ids is the same as calling is_identical() ??
      if (cv1.id() == cv2.id()) {
        CERR("result: " << res << "\n"); // EQUAL
        return res;
      }

      // both points lie on the left-right-identification, i.e., equalx
      // TODO this interface is NOT official
      CGAL::Object obj1 =
        cv1.curve().asymptotic_value_of_arc(loc1, cv1.arcno());
      CGAL::Object obj2 =
        cv2.curve().asymptotic_value_of_arc(loc1, cv2.arcno());

      typename Point_2::Curved_kernel_via_analysis_2::Curve_kernel_2::
        Algebraic_real_1 y1, y2;
      CGAL::Arr_parameter_space ps1, ps2;

      if (CGAL::assign(ps1, obj1)) {
        if (CGAL::assign(ps2, obj2)) {
          res = CGAL::EQUAL;
        } else {
          CGAL_assertion(CGAL::assign(y2, obj2));
          res = (ps1 == CGAL::ARR_BOTTOM_BOUNDARY ?
                 CGAL::SMALLER : CGAL::LARGER);
        }
      } else {
        CGAL_assertion_code(bool check = )
          CGAL::assign(y1, obj1);
        CGAL_assertion(check);
        if (CGAL::assign(ps2, obj2)) {
          res = (ps2 == CGAL::ARR_TOP_BOUNDARY ?
                 CGAL::SMALLER : CGAL::LARGER);
        } else {
          CGAL_assertion_code(bool check = )
            CGAL::assign(y2, obj2);
          CGAL_assertion(check);

          // Remark: Is filtered
          res = Base::_ckva()->kernel().compare_1_object()(y1, y2);
        }
      }

      if (res != EQUAL) {
        CERR("result: " << res << "\n");
        return res;
      }

      CGAL_precondition(cv1.is_on_left_right(loc1));
      CGAL_precondition(loc1 == cv2.location(ce));
      // comparing ids is the same as calling is_identical() ??
      if (cv1.id() == cv2.id()) {
        CGAL::Comparison_result res = CGAL::EQUAL;
        CERR("result: " << res << "\n");
        return res;
      }

      // in this setting same handling as for +/-oo ?
      res = cv1._compare_arc_numbers(cv2, loc1);
      CERR("result: " << res << "\n");
      return res;
    }
};

// bottom-top

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

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES

    //! the result type
    typedef CGAL::Arr_parameter_space result_type;

    //! the arity of the functor

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

    result_type operator()(const Point_2& pt) const {
      return pt.location();
    }

};


/*!\brief
 * Functor that compares x-limits at the top or bottom boundary
 */
template < class CurvedKernelViaAnalysis_2 >
class Compare_x_at_limit_2 : public
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES

    //! the result type
    typedef CGAL::Comparison_result result_type;

    //! the arity of the functor

    Compare_x_at_limit_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }

    /*!\brief Compare the x-limit of a vertical with the x-limit of
     * an arc end near the boundary at bottom or top boundary
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


        CERR("\ncompare_x_at_limit: p: " << p << "\n cv: " <<
             cv << "; curve_end: " << ce << "\n");

        // this curve end has boundary only in y
        CGAL_precondition(cv.is_on_bottom_top(cv.location(ce)));
        if (cv.is_singular()) // the curve end goes to contraction => x-order
            return CGAL::EQUAL; // doesn't matter

        CGAL::Comparison_result res =
            Curved_kernel_via_analysis_2::instance().
            kernel().compare_1_object()(
                    p.x(), cv.curve_end_x(ce)
            );
        CERR("result: " << res << "\n");
        return res;
    }

    /*! Compare the x-limits of 2 arcs ends near the top or bottom
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

        CERR("\ncompare_x_at_limit: cv1: " << cv1 << "\n cv2: " <<
            cv2 << "; end1: " << ce1 << "; end2: " << ce2 << "\n");
        /*CGAL::Arr_boundary_type bnd1 = boundary(end1),
            bnd2 = cv2.boundary(ce2);*/
        CGAL_precondition_code(CGAL::Arr_parameter_space loc1 = cv1.location(ce1));
        CGAL_precondition_code(CGAL::Arr_parameter_space loc2 = cv2.location(ce2));
        CGAL_precondition(cv1.is_on_bottom_top(loc1));
        CGAL_precondition(cv1.is_on_bottom_top(loc2));

        if (cv1.is_singular() != cv1.is_singular()) {
            // only one curve end lies at singularity (another at +/-oo)
            CGAL_error_msg("SINGULARITY + INF comparison is not yet \
                implemented");
        }

        CGAL::Comparison_result res = Curved_kernel_via_analysis_2::instance().
          kernel().compare_1_object()(
                                      cv1.curve_end_x(ce1),
                                      cv2.curve_end_x(ce2)
                                      );
        CERR("result: " << res << "\n");
        return res;
    }
};


/*!\brief
 * Functor that compares x-coordinates near the top or bottom boundary
 */
template < class CurvedKernelViaAnalysis_2 >
class Compare_x_near_limit_2 : public
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES

    //! the result type
    typedef CGAL::Comparison_result result_type;

    //! the arity of the functor

    Compare_x_near_limit_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }

    /*! Compare the x-coordinates of 2 arcs ends near the top or bottom
     * boundary of the parameter space
     * \param cv1 the first arc.
     * \param cv2 the second arc.
     * \param ce the arc end indicator -
     *            ARR_MIN_END - the minimal end of curves or
     *            ARR_MAX_END - the maximal end of curves.
     * \return the second comparison result:
     *         SMALLER - x(cv1, ce) \< x(cv2, ce);
     *         EQUAL   - x(cv1, ce) = x(cv2, ce);
     *         LARGER  - x(cv1, ce) > x(cv2, ce).
     *
     * \pre the ce1 end of the arc cv1 lies on a boundary.
     * \pre the ce2 end of the arc cv2 lies on a boundary.
     * \pre both curve ends are on the same boundary
     */
    result_type operator()(const Arc_2& cv1, const Arc_2& cv2,
                           CGAL::Arr_curve_end ce) const {

        CERR("\ncompare_x_near_limit: cv1: " << cv1 << "\n cv2: " <<
            cv2 << "; ce: " << ce << "\n");

        CGAL::Arr_parameter_space loc1 = cv1.location(ce);
        CGAL_precondition_code(CGAL::Arr_parameter_space loc2 =
                               cv2.location(ce));
        CGAL_precondition(cv1.is_on_bottom_top(loc1));
        CGAL_precondition(cv1.is_on_bottom_top(loc2));
        CGAL_precondition(cv1.compare_x_at_limit(ce, cv2, ce) == CGAL::EQUAL);

        CGAL_precondition(loc1 == loc2);

        Coordinate_1 x0(cv1.curve_end_x(ce));
        CGAL::Comparison_result res = cv1._compare_arc_numbers(
                                       cv2, CGAL::ARR_INTERIOR, x0,
                                       (ce == CGAL::ARR_MIN_END ?
                                        CGAL::POSITIVE : CGAL::NEGATIVE)
                                       );
        if ((ce == CGAL::ARR_MAX_END &&
             loc1 == CGAL::ARR_TOP_BOUNDARY) ||
            (ce == CGAL::ARR_MIN_END &&
             loc1 == CGAL::ARR_BOTTOM_BOUNDARY)) {
          res = opposite(res);
        }
        CERR("result: " << res << "\n");
        return res;
    }
};


/*!\brief
 * Functor to construct a point on a curve
 */
template < class CurvedKernelViaAnalysis_2 >
class Compare_endpoints_xy_2 : public
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES

    //! the result type
    typedef CGAL::Comparison_result result_type;

    /*!\brief
     * Standard constructor
     *
     * \param kernel The kernel
     */
    Compare_endpoints_xy_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }

    /*!\brief
     * Compares endpoints of arc lexicographically
     *
     * \param arc the arc
     * \return The order of endpoints
     */
    CGAL::Comparison_result operator()(const Arc_2& xcv) {
      if (xcv.is_left_to_right()) {
        return (CGAL::SMALLER);
      } else {
        return (CGAL::LARGER);
      }
      return CGAL::EQUAL;
    }
};


/*!\brief
 * Functor to construct a point on a curve
 */
template < class CurvedKernelViaAnalysis_2 >
class Construct_opposite_2 : public
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES

    //! the result type
    typedef Arc_2 result_type;

    /*!\brief
     * Standard constructor
     *
     * \param kernel The kernel
     */
    Construct_opposite_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }

    /*!\brief
     * Construct an arc of opposite direction
     *
     * \param arc the arc to be reversed
     * \return The reversed arc
     */
    Arc_2 operator()(const Arc_2& xcv) {
      return xcv.flip();
    }
};



/*!\brief
 * Functor that computes the x-extreme points of a curve
 */
template < class CurvedKernelViaAnalysis_2 >
class X_extreme_points_2 : public
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES

    typedef typename Curve_analysis_2::Coordinate_2 Coordinate_2;

    //! the result type
    typedef CGAL::cpp98::iterator< std::output_iterator_tag, Coordinate_2 >
         result_type;

    //! the arity of the functor

    /*!\brief
     * Standard constructor
     *
     * \param kernel The kernel
     */
    X_extreme_points_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }

    /*!\brief
     */
    template < class OutputIterator >
    OutputIterator operator()(const Curve_analysis_2& ca,
                              OutputIterator oi) const {

        typedef typename Curve_analysis_2::Status_line_1 Status_line_1;

        int events = ca.number_of_status_lines_with_event();

        for ( int i = 0; i < events; i++ ) {

            Status_line_1 status_line = ca.status_line_at_event(i);

            int lifts = status_line.number_of_events();

            for( int j = 0; j < lifts; j++ ) {

                std::pair<int, int> arcs =
                     status_line.number_of_incident_branches(j);

                if (arcs.first == 0 || arcs.second == 0)
                    *oi++ = status_line.algebraic_real_2(j);
            }
        }
        return oi;
    }

};

/*!\brief
 * Functor that computes the y-extreme points of a curve
 */
template < class CurvedKernelViaAnalysis_2 >
class Y_extreme_points_2 : public
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;

    CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES

    typedef typename Curve_analysis_2::Coordinate_2 Coordinate_2;

    //! the result type
    typedef CGAL::cpp98::iterator< std::output_iterator_tag, Coordinate_2 >
             result_type;

    //! the arity of the functor

    /*!\brief
     * Standard constructor
     *
     * \param kernel The kernel
     */
    Y_extreme_points_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }

    /*!\brief
     */
    template < class OutputIterator >
    OutputIterator operator()(const Curve_analysis_2& ca,
                              OutputIterator oi) const {

        typedef typename Curve_analysis_2::Status_line_1 Status_line_1;

        typename Base::Curve_kernel_2 curve_kernel = Base::_ckva()->kernel();

        Curve_analysis_2 ca_yx
            = curve_kernel.swap_x_and_y_2_object() (ca);

        std::vector<Coordinate_2> y_critical_points;

        Base::_ckva()->x_extreme_points_2_object()(ca_yx,
             std::back_inserter(y_critical_points));

        for( typename std::vector<Coordinate_2>::iterator it
                 = y_critical_points.begin();
             it != y_critical_points.end();
             it++ ) {

            Coordinate_1 curr_x = curve_kernel.get_y_2_object()( *it );

            Status_line_1 status_line = ca.status_line_at_exact_x(curr_x);

            int lifts = status_line.number_of_events();

            for( int i = 0; i < lifts; i++ ) {
                Coordinate_2 lift_xy = status_line.algebraic_real_2(i);

                bool y_coordinate_found;

                while(true) {

                    if( curve_kernel.upper_boundary_y_2_object() (lift_xy) <
                        curve_kernel.lower_boundary_x_2_object() (*it) ) {
                        y_coordinate_found = false;
                        break;
                    }
                    if( curve_kernel.upper_boundary_y_2_object() (lift_xy) >=
                        curve_kernel.upper_boundary_x_2_object() (*it) ) {
                        y_coordinate_found = true;
                        break;
                    }

                    curve_kernel.refine_x_2_object() (*it);
                }

                if(y_coordinate_found) {
                    *oi++ = Coordinate_2(curr_x, ca, i);
                    break;
                }
            }
        }
        return oi;
    }

};


#undef CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES

} // namespace Curved_kernel_via_analysis_2_Functors

/*!\brief
 * Collects the main set of functors of a curved kernel
 */
template < class CurvedKernelViaAnalysis_2, class Dummy = void>
class Curved_kernel_via_analysis_2_functors {

public:
    //!\name Public types
    //!@{

    //! this instance's first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;
    //!@}
//     typedef Curved_kernel_via_analysis_2_functors<
//         CurvedKernelViaAnalysis_2 > Functor_base;

// declares curved kernel functors, for each functor defines a member function
// returning an instance of this functor
#define CGAL_CKvA_2_functor_pred(Y, Z) \
    typedef \
    Curved_kernel_via_analysis_2_Functors::Y< Curved_kernel_via_analysis_2 \
        > Y; \
    Y Z() const { return Y(&Curved_kernel_via_analysis_2::instance()); }

#define CGAL_CKvA_2_functor_cons(Y, Z) CGAL_CKvA_2_functor_pred(Y, Z)

    CGAL_CKvA_2_functor_pred(Compare_x_2, compare_x_2_object);
    CGAL_CKvA_2_functor_pred(Compare_xy_2, compare_xy_2_object);
    CGAL_CKvA_2_functor_pred(Is_vertical_2, is_vertical_2_object);

    CGAL_CKvA_2_functor_cons(Construct_min_vertex_2,
                             construct_min_vertex_2_object);
    CGAL_CKvA_2_functor_cons(Construct_max_vertex_2,
                             construct_max_vertex_2_object);
    CGAL_CKvA_2_functor_cons(Construct_interior_vertex_2,
                             construct_interior_vertex_2_object);

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

    // left-right
    CGAL_CKvA_2_functor_pred(Parameter_space_in_x_2,
                             parameter_space_in_x_2_object);
    CGAL_CKvA_2_functor_pred(Compare_y_near_boundary_2,
                             compare_y_near_boundary_2_object);

    // bottom-top
    CGAL_CKvA_2_functor_pred(Parameter_space_in_y_2,
                             parameter_space_in_y_2_object);
    CGAL_CKvA_2_functor_pred(Compare_x_at_limit_2,
                             compare_x_at_limit_2_object);
    CGAL_CKvA_2_functor_pred(Compare_x_near_limit_2,
                             compare_x_near_limit_2_object);

    CGAL_CKvA_2_functor_cons(X_extreme_points_2, x_extreme_points_2_object);
    CGAL_CKvA_2_functor_cons(Y_extreme_points_2, y_extreme_points_2_object);

    CGAL_CKvA_2_functor_cons(Construct_point_2,
                             construct_point_2_object);

    CGAL_CKvA_2_functor_cons(Construct_point_on_arc_2,
                             construct_point_on_arc_2_object);

    CGAL_CKvA_2_functor_cons(Construct_arc_2,
                             construct_arc_2_object);

    CGAL_CKvA_2_functor_pred(Compare_endpoints_xy_2,
                             compare_endpoints_xy_2_object);

    CGAL_CKvA_2_functor_cons(Construct_opposite_2,
                             construct_opposite_2_object);

#undef CGAL_CKvA_2_functor_pred
#undef CGAL_CKvA_2_functor_cons

}; // Curved_kernel_via_analysis_2_functors

} // namespace internal

} //namespace CGAL

#endif // CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_FUNCTORS_H
// EOF
