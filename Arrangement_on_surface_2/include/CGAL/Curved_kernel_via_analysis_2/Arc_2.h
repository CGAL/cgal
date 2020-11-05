// Copyright (c) 2007,2008,2009,2010,2011 Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//                 Pavel Emeliyanenko <asm@mpi-sb.mpg.de>

#ifndef CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_ARC_2_H
#define CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_ARC_2_H

/*!\file include/CGAL/Curved_kernel_via_analysis_2/Arc_2.h
 *\brief Defines class \c Arc_2 that represents an arc on a curve that
 * can be analyzed.
 */

#include <CGAL/config.h>
#include <CGAL/Handle_with_policy.h>

#include <iostream>
#include <boost/optional.hpp>
#include <boost/none.hpp>
#include <boost/type_traits/is_same.hpp>

#include <CGAL/Bbox_2.h>
#include <CGAL/Arr_enums.h>

#ifndef CGAL_CKvA_USE_CACHES
#define CGAL_CKvA_USE_CACHES 1
#endif

#include <CGAL/Curved_kernel_via_analysis_2/Point_2.h>

#include <CGAL/Curved_kernel_via_analysis_2/Sweep_curves_adapter_2.h>

namespace CGAL {

namespace internal {

#ifndef CKvA_CERR
//#define CKvA_DEBUG_PRINT_CERR
#ifdef CKvA_DEBUG_PRINT_CERR
#define CKvA_CERR(x) std::cerr << x
#else
#define CKvA_CERR(x) static_cast<void>(0)
#endif
#endif

/*!\brief
 * Default representation class for Arc_2
 */
template < class CurvedKernelViaAnalysis_2 >
class Arc_2_rep {

public:
    //!\name Public types
    //!@{

    //! this instance's first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! myself
    typedef Arc_2_rep< Curved_kernel_via_analysis_2 > Self;

    //! type of curve kernel
    typedef typename Curved_kernel_via_analysis_2::Curve_kernel_2
    Curve_kernel_2;

    //! type of curve that can be analzed
    typedef typename Curve_kernel_2::Curve_analysis_2 Curve_analysis_2;

    //! type of a point on a point that can be analyzed
    typedef typename Curved_kernel_via_analysis_2::Point_2 Point_2;

    //! type of boundary value in x-range of an arc
    typedef typename Curve_kernel_2::Bound Bound;

    //!@}

public:
    //!\name Constructors
    //!@{

    //! default constructor
    Arc_2_rep() :
        _m_arcno(-1), _m_arcno_min(-1), _m_arcno_max(-1),
        _m_is_vertical(false),
        _m_left_to_right(1) {
    }

    //! copy constructor
    Arc_2_rep(const Self& s):
        _m_min(s._m_min), _m_max(s._m_max),
        _m_support(s._m_support),
        _m_arcno(s._m_arcno),
        _m_arcno_min(s._m_arcno_min),
        _m_arcno_max(s._m_arcno_max),
        _m_is_vertical(s._m_is_vertical),
        _m_left_to_right(s._m_left_to_right) {
    }

    //! standard constructor
    Arc_2_rep(const Point_2& p, const Point_2& q, const Curve_analysis_2& c,
              int arcno = -1, int arcno_p = -1, int arcno_q = -1,
              bool is_vertical = false) :
        _m_min(p), _m_max(q),
        _m_support(c),
        _m_arcno(arcno), _m_arcno_min(arcno_p), _m_arcno_max(arcno_q),
        _m_is_vertical(is_vertical),
        _m_left_to_right(1) {

        // set end-point arcnos from segment's interior
        if (_m_arcno_min == -1) {
            _m_arcno_min = _m_arcno;
        }
        if (_m_arcno_max == -1) {
            _m_arcno_max = _m_arcno;
        }
    }

    //!@}

public:
    //!\name Data members
    //!@{

    //! minimal end-points of an arc
    mutable Point_2 _m_min;

    //! maximal end-points of an arc
    mutable Point_2 _m_max;

    //! supporting curve
    mutable Curve_analysis_2 _m_support;

    //! interior arcno
    mutable int _m_arcno;

    //! arcno at min
    mutable int _m_arcno_min;

    //! arcno at max
    mutable int _m_arcno_max;

    //! indicates whether arc is vertical
    mutable bool _m_is_vertical;

    //! stores a direction (left to right, right to left)
    mutable bool _m_left_to_right;

    //! stores the index of an interval this arc belongs to
    mutable boost::optional<int> _m_interval_id;

    //! stores boundary value in x-range of non-vertical interval
    mutable boost::optional< Bound > _m_boundary_in_interval;

    //! stores a bbox for an arc
    mutable boost::optional< CGAL::Bbox_2 > _m_bbox;

    //!@}
};


/*!\brief
 * Class defines an arc on a curve that can be analyzed
 *
 * An arc is either non-vertical or vertical. If it is non-vertical,
 * we can assign a constant arc number to its interior and we may can assign
 * (non-identical) arc numbers to its end-points. It depends on whether
 * the arc touches the boundary of the parameter space or not.
 * If it is vertical, no arc number is available.
 *
 * We distinguish between interior arcs, rays, and branches. An interior arc
 * lies completely in the interior of the parameter space,
 * while a ray has one end that lies on the boundary of the parameter space,
 * and a branch has two end that lie on the boundary.
 *
 */
template < class CurvedKernelViaAnalysis_2, class Rep_ >
class Arc_2 :
        public CGAL::Handle_with_policy< Rep_ > {

public:
    //!\name Public types
    //!@{

    //! this instance's first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! this instance's second template parameter
    typedef Rep_ Rep;

    //! this instance itself
    typedef Arc_2< Curved_kernel_via_analysis_2, Rep > Self;

    //! type of curve kernel
    typedef typename Curved_kernel_via_analysis_2::Curve_kernel_2
    Curve_kernel_2;

    //! type of an x-coordinate
    typedef typename Curve_kernel_2::Coordinate_1 Coordinate_1;

    //! type of an xy-coordinate
    typedef typename Curve_kernel_2::Coordinate_2 Coordinate_2;

    //! type of "rational" value in x-range
    typedef typename Curve_kernel_2::Bound Bound;

    //! type of analysis of a pair of curves
    typedef typename Curve_kernel_2::Curve_analysis_2 Curve_analysis_2;

    //! type of analysis of a pair of curves
    typedef typename Curve_kernel_2::Curve_pair_analysis_2
    Curve_pair_analysis_2;

    //! type of a kernel point
    typedef typename Curved_kernel_via_analysis_2::Point_2 Point_2;
    // Remark: Point_2 is already Kernel_point_2 -> no need to introduce it

    //! type of kernel arc
    typedef typename Curved_kernel_via_analysis_2::Arc_2 Kernel_arc_2;

    //! the handle superclass
    typedef ::CGAL::Handle_with_policy< Rep > Base;

    //!@}

    #if !defined(CGAL_NO_ASSERTIONS)
    static const bool Kernel_arc_2_equals_Arc_2 = boost::is_same<Arc_2, Kernel_arc_2>::value;
    #endif

public:
    //!\name Rebind
    //!{

    /*!\brief
     * An auxiliary structure for rebinding the arc with a NewCKvA_2 and a
     * NewRep
     */
    template < typename NewCKvA_2, typename NewRep >
    class rebind
    {
    public:
        //! this instance's first template parameter
        typedef NewCKvA_2 New_curved_kernel_via_analysis_2;

        //! this instance's second template parameter
        typedef NewRep New_rep;

        //! the rebound type
        typedef Arc_2< New_curved_kernel_via_analysis_2, NewRep > Other;

        //! surface point type
        typedef typename Other::Point_2 New_point_2;

        //! type of rebound arc
        typedef typename New_curved_kernel_via_analysis_2::Arc_2 Rebound_arc_2;

        /*!\brief
         * Constructs supporting arc of type \c Rebound_arc_2 from the
         * (possible unbounded) \c arc
         * of type \c Self and replaces \c min and \c max point by the
         * given instances.
         *
         * All known items of the base class rep will be copied.
         *
         * \param arc Input arc
         * \param min New endpoint at min
         * \param max New endpoint at max
         * \return An arc of type \c Rebound_arc_2
         */
        Rebound_arc_2 operator()(const Self& arc,
                                 const New_point_2& min,
                                 const New_point_2& max) {
            New_rep newrep;
            newrep._m_min = min;
            newrep._m_max = max;

            copy_members(arc, newrep);

            return Rebound_arc_2(newrep);
        }

        /*!\brief
         * Gives direct access to minimal or maximal endpoint of an arc,
         * even if it is not interior!
         *
         * \param arc Input arc
         * \param ce Specifies which end is queried
         * \return min or max point (might be not interior!)
         */
        Point_2 operator()(const Self& arc, CGAL::Arr_curve_end ce) {
            if (ce == CGAL::ARR_MIN_END) {
                return arc._minpoint();
            }
            // else
            return arc._maxpoint();
        }

        // TODO move to SfA_2l
        /*!\brief
         * Reverse rebind, that is extracts original arc type from a
         * rebound instance
         *
         * \param arc A rebound arc
         * \return An arc of type \c Self extracted from \c arc
         */
        Self operator()(const Rebound_arc_2& arc) {
            Rep rep;

            rep._m_min = typename New_point_2::Rebind()(arc._minpoint());
            rep._m_max = typename New_point_2::Rebind()(arc._maxpoint());

            rep._m_support = arc.ptr()->_m_support;

            rep._m_arcno = arc.ptr()->_m_arcno;
            rep._m_arcno_min = arc.ptr()->_m_arcno_min;
            rep._m_arcno_max = arc.ptr()->_m_arcno_max;

            rep._m_is_vertical = arc.ptr()->_m_is_vertical;
            rep._m_left_to_right = arc.ptr()->_m_left_to_right;

            rep._m_interval_id = arc.ptr()->_m_interval_id;

            rep._m_boundary_in_interval =
                arc.ptr()->_m_boundary_in_interval;

            return Self(rep);
        }

    protected:
        /*!\brief
         * copies main members to a rep
         *
         * \param arc Source arc
         * \param newrep Destination representation
         */
        void copy_members(const Self& arc, New_rep& newrep) {

            newrep._m_support = arc.ptr()->_m_support;

            newrep._m_arcno = arc.ptr()->_m_arcno;
            newrep._m_arcno_min = arc.ptr()->_m_arcno_min;
            newrep._m_arcno_max = arc.ptr()->_m_arcno_max;

            newrep._m_is_vertical = arc.ptr()->_m_is_vertical;
            newrep._m_left_to_right = arc.ptr()->_m_left_to_right;

            newrep._m_interval_id = arc.ptr()->_m_interval_id;

            newrep._m_boundary_in_interval =
                arc.ptr()->_m_boundary_in_interval;
        }
    };

    //!}

public:
    //!\name basic constructors
    //!@{

    /*!\brief
     * Default constructor
     */
    Arc_2() :
        Base(Rep()) {
    }

    /*!\brief
     * copy constructor
     */
#ifdef DOXYGEN_RUNNING
    Arc_2(const Self& a) :
        Base(static_cast<const Base&>(a)) {
    }
#endif
    //!@}

public:
    //!\name Constructors for non-vertical arcs
    //!@{

    /*!\brief
     * Constructs an arc with two interior end-points (segment).
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
    Arc_2(const Point_2& p, const Point_2& q, const Curve_analysis_2& c,
          int arcno, int arcno_p, int arcno_q) :
        Base(Rep(p, q, c, arcno, arcno_p, arcno_q)) {

        CGAL_precondition(!p.is_identical(q));
        CGAL_precondition(Curved_kernel_via_analysis_2::instance().
                          compare_x_2_object()(p,q) != CGAL::EQUAL);
        // preconditions for arc-numbers and event points (should the common
        // parts be moved to a dedicated method ?)
        CGAL_precondition(arcno >= 0);
        CGAL_precondition(arcno_p >= 0);
        CGAL_precondition(arcno_q >= 0);
        // check end-points arcnos validity and coprimality condition
        // for supporting curves
        _check_pt_arcno_and_coprimality(p, arcno_p, c);
        _check_pt_arcno_and_coprimality(q, arcno_q, c);
        _fix_curve_ends_order(); // lexicographical order of curve ends
    }

   /*!\brief
     * Constructs an arc with one interior end-point and another end
     * at the left or right boundary of the parameter space (ray I).
     *
     * \param origin The interior end-point of the ray
     * \param inf_end Defining whether the arcs emanates from the left or right
     *        boundary
     * \param c The supporting curve
     * \param arcno The arcnumber wrt \c c in the interior of the arc
     * \param arcno_o The arcnumber wrt \c c of the arc at \c origin
     * \return The constructed ray
     */
    Arc_2(const Point_2& origin, CGAL::Arr_curve_end inf_end,
          const Curve_analysis_2& c, int arcno, int arcno_o) :
        Base(Rep(origin, Point_2(inf_end, c, arcno), c, arcno, arcno_o)) {

        CGAL_precondition(arcno >= 0);
        CGAL_precondition(arcno_o >= 0);
        // check end-points arcnos validity and coprimality condition
        // for supporting curves

        _check_pt_arcno_and_coprimality(origin, arcno_o, c);
        _fix_curve_ends_order(); // lexicographical order of curve ends
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
    Arc_2(const Point_2& origin, const Coordinate_1& asympt_x,
          CGAL::Arr_curve_end inf_end, const Curve_analysis_2& c, int arcno,
          int arcno_o) :
        Base(Rep(origin, Point_2(asympt_x, c, inf_end), c, arcno, arcno_o)) {

        CGAL_precondition(
                Curved_kernel_via_analysis_2::instance().
                kernel().compare_1_object()(origin.x(), asympt_x)
                != CGAL::EQUAL);
        CGAL_precondition(arcno >= 0);
        CGAL_precondition(arcno_o >= 0);
        _check_pt_arcno_and_coprimality(origin, arcno_o, c);

        _fix_curve_ends_order(); // lexicographical order of curve ends
    }

    /*!\brief
     * Constructs a non-vertical arc with two non-interior ends at the
     * left and right boundary (branch I)
     *
     * \param c The supporting curve
     * \param arcno The arcnumber wrt to \c c in the interior of the arc
     * \return The constructed branch
     */
    Arc_2(const Curve_analysis_2& c, int arcno) :
        Base(Rep(Point_2(CGAL::ARR_MIN_END, c, arcno),
                 Point_2(CGAL::ARR_MAX_END, c, arcno), c, arcno)) {

        CGAL_precondition(arcno >= 0);
        _fix_curve_ends_order();
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
    Arc_2(const Coordinate_1& asympt_x1, CGAL::Arr_curve_end inf_end1,
          const Coordinate_1& asympt_x2, CGAL::Arr_curve_end inf_end2,
          const Curve_analysis_2& c, int arcno) :
        Base(Rep(Point_2(asympt_x1, c, inf_end1),
                 Point_2(asympt_x2, c, inf_end2),
                 c, arcno)) {

        CGAL_precondition(
                Curved_kernel_via_analysis_2::instance().
                kernel().compare_1_object()(asympt_x1, asympt_x2)
                != CGAL::EQUAL);
        CGAL_precondition(arcno >= 0);
        _fix_curve_ends_order();
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
    Arc_2(CGAL::Arr_curve_end inf_endx, const Coordinate_1& asympt_x,
          CGAL::Arr_curve_end inf_endy, const Curve_analysis_2& c, int arcno) :
        Base(Rep(Point_2(inf_endx, c, arcno),
                 Point_2(asympt_x, c, inf_endy), c, arcno)) {

        CGAL_precondition(arcno >= 0);
        _fix_curve_ends_order();
    }

    //!@}

public:
    //!\name Constructors for vertical arcs
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
    Arc_2(const Point_2& p, const Point_2& q, const Curve_analysis_2& c) :
        Base(Rep(p, q, c, -1, -1, -1, true)) {

        CGAL_precondition(!p.is_identical(q));
        CGAL_precondition(Curved_kernel_via_analysis_2::instance().
                          compare_x_2_object()(p,q) == CGAL::EQUAL);
        CGAL_precondition(Curved_kernel_via_analysis_2::instance().
                          compare_xy_2_object()(p, q, true) != CGAL::EQUAL);
        // check coprimality condition for supporting curves
        _check_pt_arcno_and_coprimality(p, -1, c);
        _check_pt_arcno_and_coprimality(p, -1, c);
        _fix_curve_ends_order();
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
    Arc_2(const Point_2& origin, CGAL::Arr_curve_end inf_end,
          const Curve_analysis_2& c) :
        Base(Rep(origin, Point_2(origin.x(), c, inf_end),
                 c, -1, -1, -1, true)) {

        // check coprimality condition for supporting curves
        _check_pt_arcno_and_coprimality(origin, -1, c);
        _fix_curve_ends_order();
    }

    /*!\brief
     * Constructs a vertical arc that connects bottom with top boundary
     * (vertical branch)
     *
     * \param x The x-coordinate of the arc
     * \return The constructed branch
     *
     * \pre c must have a vertical line component at this x
     */
    Arc_2(const Coordinate_1& x, const Curve_analysis_2& c) :
        Base(Rep(Point_2(x, c, CGAL::ARR_MIN_END),
                 Point_2(x, c, CGAL::ARR_MAX_END), c, -1, -1, -1, true)) {

        _fix_curve_ends_order();
    }

    //!@}

protected:
    //!\name Constructor for replace endpoints + rebind
    //!@{

    /*!\brief
     * Constructs an arc from a given representation used in rebind
     *
     * \param rep Input representation
     * \return The constructed arc
     */
    Arc_2(Rep rep) :
        Base(rep) {
    }

    //!@}

public:
    //!\name Destructors (can be many ?)
    //!@{

    /*!\brief
     * Standard destructor
     */
    virtual ~Arc_2() {
    }

    //!@}

#define CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_ARC(X, Y) \
    typename Curved_kernel_via_analysis_2::X Y = \
         Curved_kernel_via_analysis_2::instance().Y##_object(); \


public:
    //!\name Parameter space
    //!@{

    /*!\brief
     * location of arc's end
     *
     * \param ce The intended end
     * \return The location of arc's \c ce in parameterspace
     */
    CGAL::Arr_parameter_space location(CGAL::Arr_curve_end ce) const {
        if (ce == CGAL::ARR_MIN_END) {
            return _minpoint().location();
        }
        return _maxpoint().location();
    }

    /*!\brief
     *  Sets boundary type for an end of an arc
     *
     * It is supposed that the user thoroughly understands malicious
     * consequences that may result from the misuse of the location
     *
     * \param ce The intended end
     * \param loc The location to store
     */
    void set_location(CGAL::Arr_curve_end ce,
                      CGAL::Arr_parameter_space loc) const {
        (ce == CGAL::ARR_MIN_END ?
         _minpoint().set_location(loc) : _maxpoint().set_location(loc));
    }

    //!@}

    //!\name Access functions
    //!@{

    /*!\brief
     * Is a curve-end finite?
     *
     * \param ce The intended end
     * \return \c true, if finite, \c false, otherwise
     */
    bool is_finite(CGAL::Arr_curve_end ce) const {
        const Point_2& pt =
            (ce == CGAL::ARR_MIN_END ? _minpoint() : _maxpoint());
        return pt.is_finite();
    }

    /*!\brief
     * returns arc's interior curve end
     *
     * \param ce The intended end
     * \return The minimal point of the arc, or the maximal point of the arc
     *
     *  \pre accessed curve end has finite coordinates
     */
    Point_2 curve_end(CGAL::Arr_curve_end ce) const {
        const Point_2& pt =
            (ce == CGAL::ARR_MIN_END ? _minpoint() : _maxpoint());
#if !CGAL_ARRANGEMENT_ON_DUPIN_CYCLIDE
        CGAL_precondition(pt.location() == CGAL::ARR_INTERIOR ||
                          pt.is_finite());
#endif
        return pt;
    }

    /*!\brief
     * returns x-coordinate of arc's curve end
     *
     * \param ce The intended end
     * \return x-coordinate of arc's end at \c ce
     *
     * \pre accessed curve end has finite x-coordinate
     */
    inline
    Coordinate_1 curve_end_x(CGAL::Arr_curve_end ce) const {
        CGAL_precondition(
                !(ce == CGAL::ARR_MIN_END ? _minpoint().is_on_left_right() :
                  _maxpoint().is_on_left_right()));
        return (ce == CGAL::ARR_MIN_END ? _minpoint().x() : _maxpoint().x());
    }

    /*!\brief
     * supporting curve of the arc
     *
     * \return supporting curve of the arc
     */
    inline
    const Curve_analysis_2& curve() const {
        return this->ptr()->_m_support;
    }

    /*!\brief arc number in interior
     *
     * \return arc number
     *
     * \pre !is_vertical()
     */
    inline
    int arcno() const {
        CGAL_precondition(!is_vertical());
        return this->ptr()->_m_arcno;
    }

    /*!\brief
     * arc number of end of arc, which may be different from arc number in its
     * interior
     *
     * \param ce The intended end
     * \return Arc number of intended end
     * \pre !is_vertical()
     */
    inline
    int arcno(CGAL::Arr_curve_end ce) const {
        CGAL_precondition(!is_vertical());
        return (ce == CGAL::ARR_MIN_END ? this->ptr()->_m_arcno_min :
                this->ptr()->_m_arcno_max);
    }

    /*!\brief
     * arc number at given x-coordinate
     *
     * If \c x0 is equal to source's or target's x-coordinate,
     * then the arc number of that point is returned.
     * Otherwise the arc number of the interior is returned.
     *
     * \param x0 queried x-coordinate
     * \returns arcnumber at \c x0
     * \pre !is_vertical()
     * \pre \c x0 must be within arcs's x-range.
     */
    inline
    int arcno(const Coordinate_1& x0) const {
        CGAL_precondition(!is_vertical());
        CGAL_precondition(is_in_x_range(x0));

        if (this->ptr()->_m_arcno_min != this->ptr()->_m_arcno &&
            is_finite(CGAL::ARR_MIN_END) &&
            Curved_kernel_via_analysis_2::instance().
            kernel().compare_1_object()(x0, _minpoint().x()) ==
            CGAL::EQUAL) {
            return this->ptr()->_m_arcno_min;
        }
        if (this->ptr()->_m_arcno_max != this->ptr()->_m_arcno &&
            is_finite(CGAL::ARR_MAX_END) &&
            Curved_kernel_via_analysis_2::instance().
            kernel().compare_1_object()(x0, _maxpoint().x()) ==
            CGAL::EQUAL) {
            return this->ptr()->_m_arcno_max;
        }
        return this->ptr()->_m_arcno;
    }

    /*!\brief
     * checks if the arc is vertical
     *
     * \return \c true, if vertical, \c false, otherwise
     */
    inline
    bool is_vertical() const {
        return this->ptr()->_m_is_vertical;
    }

    /*!\brief
     * returns x-coordinate of vertical arc
     *
     * \return x-coordinate of line that contains vertical arc
     * \pre is_vertical
     */
    inline
    const Coordinate_1& x() const {
        CGAL_precondition(is_vertical());
        return _minpoint().x();
    }

    //!@}

    //!\name Direction
    //!@{

    /*!\brief
     * checks if the arc is direction left-to-right (min-to-max)
     *
     * \return \c true, if directed left-to-right, \c false, otherwise
     */
    inline
    bool is_left_to_right() const {
      return this->ptr()->_m_left_to_right;
    }

    /*! Flip the arc . */
    Arc_2 flip () const  {
      Arc_2   opp(*this);
      opp.copy_on_write();
      opp.ptr()->_m_left_to_right = !this->ptr()->_m_left_to_right;

      return opp;
    }

    //!@}

    //!\name Intervals
    //!@{

    /*!\brief
     * returns the index of an open interval between two events of the
     * curve the arc belongs to
     *
     * \return interval id of supporting curve for this arc
     * \pre !is_vertical()
     */
    inline
    int interval_id() const {
        CGAL_precondition(!is_vertical());
        if(!this->ptr()->_m_interval_id)
            this->ptr()->_m_interval_id = _compute_interval_id();
        return *(this->ptr()->_m_interval_id);
    }


    /*!\brief
     * returns boundary value in interior of x-range of non-vertical
     * interval
     *
     * \return a rational x-coordinate in the interior of the arc's x-range
     * \pre !is_vertical()
     */
    Bound boundary_in_x_range_interior() const {
        CGAL_precondition(!is_vertical());
        if(!this->ptr()->_m_boundary_in_interval) {
            this->ptr()->_m_boundary_in_interval =
                _compute_boundary_in_interval();
            CGAL_postcondition_code(
                    typename Curve_analysis_2::Status_line_1 cv_line =
                    curve().status_line_at_exact_x(
                            Coordinate_1(
                                    *this->ptr()->_m_boundary_in_interval
                            )
                    );
            );
            CGAL_postcondition(cv_line.index() == interval_id());
        }
        return *(this->ptr()->_m_boundary_in_interval);
    }

    //!@}

public:
    //! \name Shortcuts for code readability
    //!@{

    //! tests whether this boundary type represents +/-oo
    inline static bool is_infinite(/*CGAL::Arr_boundary_type bnd*/) {
        return false; //(bnd == CGAL::ARR_UNBOUNDED);
    }

    //! tests whether this boundary type represents a singularity
    inline static bool is_singular(/*CGAL::Arr_boundary_type bnd*/) {
        return false; //(bnd == CGAL::ARR_CONTRACTION);
    }

    //! tests whether this boundary type represents lying on discontinuity
    inline static bool is_on_disc(/*CGAL::Arr_boundary_type bnd*/) {
        return false; //(bnd == CGAL::ARR_IDENTIFICATION);
    }

    //! returns true if a parameter encodes an entity in the interior
    inline static bool is_interior(CGAL::Arr_parameter_space loc) {
        return (loc == CGAL::ARR_INTERIOR);
    }

    //! returns true if a parameter encodes bottom or top boundary placement
    inline static bool is_on_bottom_top(CGAL::Arr_parameter_space loc) {
        return (loc == CGAL::ARR_BOTTOM_BOUNDARY ||
                loc == CGAL::ARR_TOP_BOUNDARY);
    }

    //! returns true if a parameter encodes left or right boundary placement
    inline static bool is_on_left_right(CGAL::Arr_parameter_space loc) {
        return (loc == CGAL::ARR_LEFT_BOUNDARY ||
                loc == CGAL::ARR_RIGHT_BOUNDARY);
    }

    //!@}

public:
    //! \name Predicates
    //!@{

      /*!
     * Compare the relative x-limits of a vertical line at an interior point
     * and the arc's end on a bottom or top boundary
     *
     * \param p A reference point; we refer to a vertical line incident to p.
     * \param ce ARR_MIN_END if we refer to the arc's minimal end,
     *            ARR_MAX_END if we refer to its maximal end.
     * \return CGAL::SMALLER if p lies to the left of the arc;
     *         CGAL::LARGER  if p lies to the right of the arc;
     *         CGAL::EQUAL   in case of an overlap.
     *
     * \pre the arc's relevant end is on bottom or top boundary
     */
    CGAL::Comparison_result compare_x_at_limit(
            CGAL::Arr_curve_end ce,
            const Point_2& p
    ) const {

        CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_ARC(Compare_x_at_limit_2,
                                            compare_x_at_limit_2)
        // compare with nullptr, in order to avoid a performance warning with VC++
        CGAL_precondition(Kernel_arc_2_equals_Arc_2 ||
                          dynamic_cast< const Kernel_arc_2* >(this) != nullptr);
        return compare_x_at_limit_2(
                p, *dynamic_cast< const Kernel_arc_2* >(this), ce
        );
    }


    /*!\brief
     * Compare the relative x-limits of the curve end of \c *this
     * and \c cv2
     * \param ce1 ARR_MIN_END if we refer to this' minimal end,
     *             ARR_MAX_END if we refer to this' maximal end.
     * \param cv2 The second curve.
     * \param ce2 ARR_MIN_END if we refer to its minimal end,
     *             ARR_MAX_END if we refer to its maximal end.
     * \return CGAL::SMALLER if \c this lies to the left of cv2;
     *         CGAL::LARGER  if \c this lies to the right of cv2;
     *         CGAL::EQUAL   in case of an overlap.
     *
     * \pre the curve ends lie on the bottom or top boundary
     */
    CGAL::Comparison_result compare_x_at_limit(
            CGAL::Arr_curve_end ce1,
            const Kernel_arc_2& cv2, CGAL::Arr_curve_end ce2) const {

        CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_ARC(Compare_x_at_limit_2,
                                            compare_x_at_limit_2)
        CGAL_precondition(Kernel_arc_2_equals_Arc_2 ||
                          dynamic_cast< const Kernel_arc_2* >(this) != nullptr);
        return compare_x_at_limit_2(
                *dynamic_cast< const Kernel_arc_2* >(this), ce1, cv2, ce2
        );
    }

    /*!
     * Compare the relative x-positions of an interior point
     * and the arc's end on a bottom or top boundary
     *
     * \param p A reference point; we refer to a vertical line incident to p.
     * \param ce ARR_MIN_END if we refer to the arc's minimal end,
     *            ARR_MAX_END if we refer to its maximal end.
     * \return CGAL::SMALLER if p lies to the left of the arc;
     *         CGAL::LARGER  if p lies to the right of the arc;
     *         CGAL::EQUAL   in case of an overlap.
     *
     * \pre the arc's relevant end is on bottom or top boundary
     */
    CGAL::Comparison_result compare_x_near_limit(
            CGAL::Arr_curve_end ce,
            const Point_2& p
    ) const {

        CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_ARC(Compare_x_near_limit_2,
                                            compare_x_near_limit_2)
        CGAL_precondition(Kernel_arc_2_equals_Arc_2 ||
                          dynamic_cast< const Kernel_arc_2* >(this) != nullptr);
        return compare_x_near_limit_2(
                p, *dynamic_cast< const Kernel_arc_2* >(this), ce
        );
    }

    /*!\brief
     * Compare the relative x-positions of the curve end of \c *this
     * and \c cv2
     * \param ce1 ARR_MIN_END if we refer to this' minimal end,
     *             ARR_MAX_END if we refer to this' maximal end.
     * \param cv2 The second curve.
     * \param ce2 ARR_MIN_END if we refer to its minimal end,
     *             ARR_MAX_END if we refer to its maximal end.
     * \return CGAL::SMALLER if \c this lies to the left of cv2;
     *         CGAL::LARGER  if \c this lies to the right of cv2;
     *         CGAL::EQUAL   in case of an overlap.
     *
     * \pre the curve ends lie on the bottom or top boundary
     */
    CGAL::Comparison_result compare_x_near_limit(const Kernel_arc_2& cv2,
                                                 CGAL::Arr_curve_end ce) const {

        CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_ARC(Compare_x_near_limit_2,
                                            compare_x_near_limit_2)
        CGAL_precondition(Kernel_arc_2_equals_Arc_2 ||
                          dynamic_cast< const Kernel_arc_2* >(this) != nullptr);
        return compare_x_near_limit_2(*dynamic_cast< const Kernel_arc_2* >(this), cv2, ce);
    }

    /*!\brief
     * Compare the relative y-positions of two arcs whose ends approach
     * the left or right boundary from the same side
     *
     * \param cv2 The second arc
     * \param ce ARR_MIN_END if we compare near left boundary
     *            ARR_MAX_END if we compare near right boundary
     * \return CGAL::SMALLER if this arc lies below cv2;
     *         CGAL::LARGER if this arc lies above cv2;
     *         CGAL::EQUAL in case of an overlap.
     *
     * \pre The ends are defined on left or right boundary
     */
    CGAL::Comparison_result compare_y_near_boundary(
            const Kernel_arc_2& cv2,
            CGAL::Arr_curve_end ce
    ) const {

        CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_ARC(Compare_y_near_boundary_2,
                                            compare_y_near_boundary_2)
        CGAL_precondition(Kernel_arc_2_equals_Arc_2 ||
                          dynamic_cast< const Kernel_arc_2* >(this) != nullptr);
        return compare_y_near_boundary_2(
                *dynamic_cast< const Kernel_arc_2* >(this), cv2, ce
        );
    }

    /*!\brief
     * Compares the relative vertical alignment of a point with this arc
     *
     * \param p The point.
     * \return
     * CGAL::SMALLER if y(p) \< arc(x(p)), i.e. the point is below the arc;
     * CGAL::LARGER if y(p) > arc(x(p)), i.e. the point is above the arc;
     * CGAL::EQUAL if p lies on the arc.
     *
     * \pre p is in the x-range of the arc.
     */
    CGAL::Comparison_result compare_y_at_x(const Point_2& p) const {

        CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_ARC(Compare_y_at_x_2,
                                            compare_y_at_x_2)
        CGAL_precondition(Kernel_arc_2_equals_Arc_2 ||
                          dynamic_cast< const Kernel_arc_2* >(this) != nullptr);
        return compare_y_at_x_2(p, *dynamic_cast< const Kernel_arc_2* >(this));
    }

    /*!\brief
     * Compares the relative vertical aligment of this arc with a second
     * immediately to the left of one of their intersection points.
     *
     * If one of the curves is vertical (emanating downward from p),
     * it is always considered to be below the other curve.
     *
     * \param cv2 The second arc
     * \param p The intersection point
     *
     * \return The relative vertical alignment this arc with respect to cv2
     *         immediately to the left of p: SMALLER, LARGER or EQUAL.
     *
     * \pre The point p lies on both curves, and both of them must be also be
     *      defined (lexicographical) to their left.
     */
    CGAL::Comparison_result compare_y_at_x_left(const Kernel_arc_2& cv2,
                                                const Point_2 &p) const {

        CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_ARC(Compare_y_at_x_left_2,
                                            compare_y_at_x_left_2)
        CGAL_precondition(Kernel_arc_2_equals_Arc_2 ||
                          dynamic_cast< const Kernel_arc_2* >(this) != nullptr);
        return compare_y_at_x_left_2(
                *dynamic_cast< const Kernel_arc_2* >(this), cv2, p
        );
    }

    /*!\brief
     * Compares the relative vertical aligment of this arc with a second
     * immediately to the right of one of their intersection points.
     *
     * If one of the curves is vertical (emanating downward from p),
     * it is always considered to be below the other curve.
     *
     * \param cv2 The second arc
     * \param p The intersection point
     *
     * \return The relative vertical alignment this arc with respect to cv2
     *         immediately to the right of p: SMALLER, LARGER or EQUAL.
     *
     * \pre The point p lies on both curves, and both of them must be also be
     *      defined (lexicographical) to their right.
     */
    CGAL::Comparison_result compare_y_at_x_right(const Kernel_arc_2& cv2,
                                                 const Point_2 &p) const {

        CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_ARC(Compare_y_at_x_right_2,
                                            compare_y_at_x_right_2)
        CGAL_precondition(Kernel_arc_2_equals_Arc_2 ||
                          dynamic_cast< const Kernel_arc_2* >(this) != nullptr);
        return compare_y_at_x_right_2(
                *dynamic_cast< const Kernel_arc_2* >(this), cv2, p
        );
    }

    /*!\brief
     * Check if the given x-value is in the x-range of the arc inclusive.
     *
     * \param x The x-value.
     * \param *eq_min Output: Is this value equal to the x-coordinate of the
     *                       ARR_MIN_END point.
     * \param *eq_max Output: Is this value equal to the x-coordinate of the
     *                       ARR_MAX_END point.
     * \return \c true, if p.x() is in x-range of arc, \c false otherwise
     */
    bool is_in_x_range(const Coordinate_1& x,
                       bool *eq_min = nullptr, bool *eq_max = nullptr) const {

        if (eq_min != nullptr && eq_max != nullptr) {
            *eq_min = *eq_max = false;
        }

        if (is_vertical()) {
            if (x == this->x()) {
                if (eq_min != nullptr) {
                    *eq_min = true;
                }
                if (eq_max != nullptr) {
                    *eq_max = true;
                }
                return true;
            }
            // else
            return false;
        }

        // precomputations:
        CGAL::Comparison_result resmin = CGAL::LARGER;
        CGAL::Arr_parameter_space min_loc = location(CGAL::ARR_MIN_END);
        bool min_has_x =
            (is_finite(CGAL::ARR_MIN_END) ||
             min_loc == CGAL::ARR_BOTTOM_BOUNDARY ||
             min_loc == CGAL::ARR_TOP_BOUNDARY);
        if (min_has_x) {
            resmin = Curved_kernel_via_analysis_2::instance().
                kernel().compare_1_object()(x, _minpoint().x());
            if (eq_min != nullptr) { // TODO asymptotic end in x-range?
                *eq_min = (resmin == CGAL::EQUAL);
            }
        }

        CGAL::Comparison_result resmax = CGAL::SMALLER;
        CGAL::Arr_parameter_space max_loc = location(CGAL::ARR_MAX_END);
        bool max_has_x =
            (is_finite(CGAL::ARR_MAX_END) ||
             max_loc == CGAL::ARR_BOTTOM_BOUNDARY ||
             max_loc == CGAL::ARR_TOP_BOUNDARY);

        if (max_has_x) {
            resmax = Curved_kernel_via_analysis_2::instance().
                kernel().compare_1_object()(x, _maxpoint().x());
            if (eq_max != nullptr) { // TODO asymptotic end in x-range?
                *eq_max = (resmax == CGAL::EQUAL);
            }
        }

        bool res =
            (resmin != CGAL::SMALLER && resmax != CGAL::LARGER);
        return res;
    }

    /*!\brief
     * Checks whether an x-coordinate lies in the interiors of this arc's
     * x-range
     *
     * \param x The query coordinate
     * \return \c true, if \c x lies in the interior of this arc's x-range,
     * \c false otherwise
     */
    // TODO do we need this special method ?
    bool is_in_x_range_interior(const Coordinate_1& x) const
    {
        bool eq_min, eq_max;
        if (!is_in_x_range(x, &eq_min, &eq_max) || eq_min || eq_max) {
            return false;
        }
        return true;
    }

    /*!\brief
     * Checks whether a given arc is equal to this one
     *
     * \param cv2 The query arc
     * \return \c true iff this arc is equal to \c cv, \c false otherwise
     */
    bool is_equal(const Kernel_arc_2& cv2) const {

        CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_ARC(Equal_2,
                                            equal_2)

        CGAL_precondition(Kernel_arc_2_equals_Arc_2 ||
                          dynamic_cast< const Kernel_arc_2* >(this) != nullptr);
        return equal_2(*dynamic_cast< const Kernel_arc_2* >(this), cv2);
    }

    /*!\brief
     * checks whether this arcs overlaps with another
     *
     * \param cv2 The query arc
     * \return \c true, if both arcs have infinitely many intersection points,
     *         \c false otherwise
     */
    bool do_overlap(const Kernel_arc_2& cv2) const {

        CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_ARC(Do_overlap_2,
                                            do_overlap_2)
        CGAL_precondition(Kernel_arc_2_equals_Arc_2 ||
                          dynamic_cast< const Kernel_arc_2* >(this) != nullptr);
        return do_overlap_2(*dynamic_cast< const Kernel_arc_2* >(this), cv2);
    }

    /*!\brief
     * multiplicity of intersection
     *
     * The intersection multiplicity of \c *this and \c cv2 at point \c p is
     * returned.
     *
     * \param cv2 The second arc
     * \param p The intersection point
     * \return The multiplicity of the intersection at \c p
     * \pre \c p must be an intersection point.
     */
    int multiplicity_of_intersection(
            const Kernel_arc_2& cv2, const Point_2& p) const {

        // intersection point must lie in the interior of both arcs
        CGAL_precondition_code( // because of macro stupidity one needs
            bool eq_min1;       // to omit commas in declaration
            bool eq_max1;
            bool eq_min2;
            bool eq_max2;
        );
        CGAL_precondition(is_in_x_range(p.x(), &eq_min1, &eq_max1));
        CGAL_precondition(cv2.is_in_x_range(p.x(), &eq_min2, &eq_max2));
        CGAL_precondition(is_vertical() || (!eq_min1 && !eq_max1));
        CGAL_precondition(cv2.is_vertical() || (!eq_min2 && !eq_max2));

        // there must be an intersection at this point (in_x_range is checked
        // internally by compare_y_at_x() ?
        CGAL_expensive_precondition(compare_y_at_x(p) == CGAL::EQUAL &&
            cv2.compare_y_at_x(p) == CGAL::EQUAL);

        Kernel_arc_2::simplify(*dynamic_cast< const Kernel_arc_2*>(this), cv2);
        CGAL_precondition(!curve().is_identical(cv2.curve()));
        if (is_vertical() || cv2.is_vertical()) {
            CGAL_assertion(!(is_vertical() && cv2.is_vertical()));
            return 1;
        }

        Curve_pair_analysis_2 cpa_2 =
            Curved_kernel_via_analysis_2::instance().
            kernel().construct_curve_pair_2_object()(curve(), cv2.curve());

        typename Curve_pair_analysis_2::Status_line_1 cpv_line =
                cpa_2.status_line_for_x(p.x());

        CGAL_precondition(cpv_line.is_intersection());
        int j = cpv_line.event_of_curve(arcno(p.x()), curve()),
            mult = cpv_line.multiplicity_of_intersection(j);

        CGAL_postcondition(mult > 0);
        return mult;
    }

    //!@}

    //!\name Constructing functions
    //!@{

    /*!\brief
     * Find all intersections of this arc with another one and
     * insert them to the output iterator.
     *
     * Type of output iterator is \c CGAL::Object. It either contains
     * an \c Arc_2 object (overlap) or a
     * <tt>std::pair\<Point_2, unsigned int></tt> (intersection point +
     * multiplicity). A past-the-end iterator is returned.
     *
     * \param cv2 The second arc
     * \param oi The outputiterator
     * \return A past-the-end iterator of \c oi
     */
    template < class OutputIterator >
    OutputIterator intersections(const Kernel_arc_2& cv2,
                                 OutputIterator oi) const {

        CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_ARC(Intersect_2,
                                            intersect_2)

        CGAL_precondition(Kernel_arc_2_equals_Arc_2 ||
                          dynamic_cast< const Kernel_arc_2* >(this) != nullptr);
        return intersect_2(
                *dynamic_cast< const Kernel_arc_2* >(this), cv2, oi
        );
    }

    /*!\brief
     * Computes the next intersection of \c *this and \c cv2 right of \c p
     * in lexicographical order and returns it through \c intersection
     * argument
     *
     * intersect_right_of_point is not called when using sweep_curves() with
     * intersection dictionary and without validation of internal structures
     * (as is standard). Hence we can be lazy here for the moment
     * without losing performance.
     *
     * \param cv2 The second arc
     * \param p The minimal bound point
     * \param intersection The next intersection
     * \return \c true, if there is a next intersection and
     *         \c intersection has been set properly, \c false otherwise
     * \pre The arcs are not allowed to overlap
     */
    bool intersect_right_of_point(const Kernel_arc_2& cv2, const Point_2& p,
                                  Point_2& intersection) const {

        CGAL_precondition(!this->do_overlap(cv2));

        // TODO rewrite intersect_right_of_point (Pavel)
        // use static member for Intersect, Left & Right
        // with parameters for direction and where to stop
        typedef std::vector<std::pair<Point_2, int> > Point_container;
        Point_container tmp;
        _intersection_points(
                *dynamic_cast< const Kernel_arc_2*>(this), cv2,
                back_inserter(tmp)
        );
        typename Point_container::const_iterator it;
        for (it = tmp.begin(); it != tmp.end(); it++) {
            if(it->first > p) {
                intersection = it->first;
                return true;
            }
        }
        return false;
    }

    /*!\brief
     * Computes the next intersection of \c *this and \c cv2 left of \c p
     * in lexicographical order and returns it through \c intersection
     * argument
     *
     * intersect_right_of_point is not called when using sweep_curves() with
     * intersection dictionary and without validation of internal structures
     * (as is standard). Hence we can be lazy here for the moment
     * without losing performance.
     *
     * \param cv2 The second arc
     * \param p The maximal bound point
     * \param intersection The next intersection
     * \return \c true, if there is a next intersection
     *         and \c intersection has been set properly, \c false otherwise
     * \pre The arcs are not allowed to overlap
     */
    bool intersect_left_of_point(const Kernel_arc_2& cv2, const Point_2& p,
                                 Point_2& intersection) const {

        CGAL_precondition(!this->do_overlap(cv2));

        // TODO rewrite intersect_left_of_point (Pavel)
        // use static member for Intersect, Left & Right
        // with parameters for direction and where to stop
        typedef std::vector<std::pair<Point_2, int> > Point_container;
        Point_container tmp;
        _intersection_points(
                *dynamic_cast< const Kernel_arc_2*>(this), cv2,
                back_inserter(tmp)
        );
        typename Point_container::const_reverse_iterator it;
        for(it = tmp.rbegin(); it != tmp.rend(); it++) {
            if(it->first < p) {
                intersection = it->first;
                return true;
            }
        }
        return false;
    }

    /*!\brief
     * Returns a trimmed version of an arc
     *
     * \param p the new first endpoint
     * \param q the new second endpoint
     * \return The trimmed arc
     *
     * \pre p != q
     * \pre both points must be interior and must lie on \c cv
     */
    // do we need this method separetely ??
    Kernel_arc_2 trim(const Point_2& p, const Point_2& q) const {

        CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_ARC(Trim_2, trim_2)

        CGAL_precondition(Kernel_arc_2_equals_Arc_2 ||
                          dynamic_cast< const Kernel_arc_2* >(this) != nullptr);
        return trim_2(*dynamic_cast< const Kernel_arc_2* >(this), p, q);
    }

    /*!\brief
     * Split an arc at a given point into two sub-arcs
     *
     * \param p The split point
     * \param s1 Output: The left resulting sub-arc (p is its right endpoint)
     * \param s2 Output: The right resulting sub-arc (p is its left endpoint)
     *
     * \pre p lies on cv but is not one of its end-points.
     */
    void split(const Point_2& p, Kernel_arc_2& s1, Kernel_arc_2& s2) const {

        CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_ARC(Split_2,
                                            split_2)
        CGAL_precondition(Kernel_arc_2_equals_Arc_2 ||
                          dynamic_cast< const Kernel_arc_2* >(this) != nullptr);
        split_2(*dynamic_cast< const Kernel_arc_2* >(this), p, s1, s2);
    }

    /*!\brief
     * Check whether this arc can be merged with a second
     *
     * \param cv2 The second arc
     * \return \c true if the two arcs are mergeable, i.e., they are supported
     * by the same curve and share a common endpoint; \c false otherwise.
     */
    bool are_mergeable(const Kernel_arc_2& cv2) const {

        CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_ARC(Are_mergeable_2,
                                            are_mergeable_2)
        CGAL_precondition(Kernel_arc_2_equals_Arc_2 ||
                          dynamic_cast< const Kernel_arc_2* >(this) != nullptr);
        return are_mergeable_2(
                *dynamic_cast< const Kernel_arc_2* >(this), cv2
        );
    }

  /*!\brief
     * Merges this arc with a second
     *
     * \param cv2 The second arc
     * \return The resulting arc
     *
     * \pre The two arcs are mergeable, that is they are supported by the
     *      same curve and share a common endpoint.
     */
    Kernel_arc_2 merge(const Kernel_arc_2& cv2) const {

        CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_ARC(Merge_2, merge_2)
        Kernel_arc_2 tmp;
        CGAL_precondition(Kernel_arc_2_equals_Arc_2 ||
                          dynamic_cast< const Kernel_arc_2* >(this) != nullptr);
        merge_2(*dynamic_cast< const Kernel_arc_2* >(this), cv2, tmp);
        return tmp;
    }

    //!@}
    //!\name Simplification
    //!@{

    /*! \brief
     *  simplifies representation of \c cv and/or \c p in case they have
     *  non-coprime supporting curves.
     *
     *  \return \c true if simplification took place, \c false otherwise
     */
    static bool simplify(const Kernel_arc_2& cv, const Coordinate_2& p) {

        if (cv.curve().is_identical(p.curve())) {
            return false;
        }

        std::vector<Curve_analysis_2> parts_of_f, parts_of_g, common;

        if (Curved_kernel_via_analysis_2::instance().
            kernel().decompose_2_object()(
                    cv.curve(), p.curve(),
                    std::back_inserter(parts_of_f),
                    std::back_inserter(parts_of_g),
                    std::back_inserter(common))) {

            CGAL_assertion((parts_of_f.size() == 1 ||
                            parts_of_g.size() == 1) && common.size() == 1);
            if (parts_of_f.size() == 1) {
                cv._simplify_by(
                    Curved_kernel_via_analysis_2::instance().
                        kernel().construct_curve_pair_2_object()
                            (parts_of_f[0], common[0]));
            }
            if (parts_of_g.size() == 1) {
                p.simplify_by(Curved_kernel_via_analysis_2::instance().
                        kernel().construct_curve_pair_2_object()
                            (parts_of_g[0], common[0]));
            }
            return true;
        }
        return false;
    }

    /*!\brief
     * simplifies representation of \c cv1 and/or \c cv2 in case they have
     * non-coprime supporting curves.
     *
     *  \return \c true if simplification took place, \c false otherwise
     */
    static bool simplify(const Kernel_arc_2& cv1, const Kernel_arc_2& cv2) {

        if (cv1.curve().is_identical(cv2.curve())) {
            return false;
        }

        std::vector<Curve_analysis_2> parts_of_f, parts_of_g, common;

        if (Curved_kernel_via_analysis_2::instance().
            kernel().decompose_2_object()(
                    cv1.curve(), cv2.curve(),
                    std::back_inserter(parts_of_f),
                    std::back_inserter(parts_of_g),
                    std::back_inserter(common))) {
            CGAL_assertion((parts_of_f.size() == 1 ||
                       parts_of_g.size() == 1) && common.size() == 1);
            if (parts_of_f.size() == 1) {
                cv1._simplify_by(Curved_kernel_via_analysis_2::instance().
                        kernel().construct_curve_pair_2_object()
                            (parts_of_f[0], common[0]));
            }
            if (parts_of_g.size() == 1) {
                cv2._simplify_by(Curved_kernel_via_analysis_2::instance().
                        kernel().construct_curve_pair_2_object()
                            (parts_of_g[0], common[0]));
            }
            return true;
        }
        return false;
    }

    //!@}

protected:
    //!\name Trimming
    //!@{

    /*!\brief
     * Returns a trimmed version of an arc (internal version that does not use
     * functor)
     *
     * \param p the new first endpoint
     * \param q the new second endpoint
     * \return The trimmed arc
     *
     * \pre p != q
     * \pre both points must be interior and must lie on \c cv
     */
    // TODO implement in functor?
    Kernel_arc_2 _trim(const Point_2& p, const Point_2& q) const {

        if (p.location() == CGAL::ARR_INTERIOR &&
            q.location() == CGAL::ARR_INTERIOR) {

            return _replace_endpoints(p, q,
                    (is_vertical() ? -1 : arcno(p.x())),
                    (is_vertical() ? -1 : arcno(q.x()))).first;
        }

        if (p.location() != CGAL::ARR_INTERIOR &&
            q.location() != CGAL::ARR_INTERIOR)
            return static_cast<const Kernel_arc_2&>(*this);

        Kernel_arc_2 left_arc, right_arc;
        CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_ARC(Split_2, split_2)

        if (p.location() != CGAL::ARR_INTERIOR &&
            q.location() == CGAL::ARR_INTERIOR) {

            split_2(static_cast<const Kernel_arc_2&>(*this), q, left_arc,
                    right_arc);
            return left_arc;
        }
        // if (p.location() == CGAL::ARR_INTERIOR &&
        //   q.location() != CGAL::ARR_INTERIOR)

        split_2(static_cast<const Kernel_arc_2&>(*this), p, left_arc,
                    right_arc);
        return right_arc;
    }

public:

    /*!\brief
     * Trims this arc and \c cv2 to the common x-range, if it is non-trivial
     *
     * \param cv2 the second arc
     * \param trimmed1 Output: trimmed version of \c *this to joint x-range of
     *                 \c *this and \c cv2
     * \param trimmed1 Output: trimmed version of \c cv2 to joint x-range of
     *                 \c *this and \c cv2
     * \return \c true, if \c *this and \c cv2 share a non-trivial
     *         common x-range, \c false otherwise
     */
    bool trim_by_arc(const Kernel_arc_2& cv2, Kernel_arc_2& trimmed1,
                     Kernel_arc_2& trimmed2) const {

        CGAL_precondition(Kernel_arc_2_equals_Arc_2 ||
                          dynamic_cast< const Kernel_arc_2* >(this) != nullptr);

        const Kernel_arc_2& cv1 = static_cast< const Kernel_arc_2& >(*this);

        Point_2 common_left, common_right;

        bool joint = cv1._joint_x_range(cv2, common_left, common_right);

        if (!joint) {
            return false;
        }

        typename Curve_kernel_2::Compare_1 compare_x(
                Curved_kernel_via_analysis_2::instance().
                    kernel().compare_1_object());

        CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_ARC(Construct_point_on_arc_2,
                                            construct_point_on_arc_2)

        Point_2 left1, left2;

        if (common_left.location() != CGAL::ARR_LEFT_BOUNDARY) {
            if ((cv1.location(CGAL::ARR_MIN_END) !=
                 CGAL::ARR_LEFT_BOUNDARY)  &&
                (compare_x(cv1.curve_end_x(CGAL::ARR_MIN_END),
                           common_left.x()) == CGAL::EQUAL) ) {
                left1 = cv1._minpoint();
            } else {
                left1 = construct_point_on_arc_2(common_left.x(),
                                               cv1.curve(),
                                               cv1.arcno(),
                                               cv1);
            }
            if ((cv2.location(CGAL::ARR_MIN_END) !=
                 CGAL::ARR_LEFT_BOUNDARY)  &&
                (compare_x(cv2.curve_end_x(CGAL::ARR_MIN_END),
                           common_left.x()) == CGAL::EQUAL) ) {
                left2 = cv2._minpoint();
            } else {
                left2 = construct_point_on_arc_2(common_left.x(),
                                               cv2.curve(),
                                               cv2.arcno(),
                                               cv2);
            }
        } else {
            left1 = cv1._minpoint();
            left2 = cv2._minpoint();
        }


        Point_2 right1, right2;

        if (common_right.location() != CGAL::ARR_RIGHT_BOUNDARY) {

            if ((cv1.location(CGAL::ARR_MAX_END) !=
                 CGAL::ARR_RIGHT_BOUNDARY)  &&
                (compare_x(cv1.curve_end_x(CGAL::ARR_MAX_END),
                           common_right.x()) == CGAL::EQUAL) ) {
                right1 = cv1._maxpoint();
            } else {
                right1 = construct_point_on_arc_2(common_right.x(),
                                                cv1.curve(),
                                                cv1.arcno(),
                                                cv1);
            }
            if ((cv2.location(CGAL::ARR_MAX_END) !=
                 CGAL::ARR_RIGHT_BOUNDARY)  &&
                (compare_x(cv2.curve_end_x(CGAL::ARR_MAX_END),
                           common_right.x()) == CGAL::EQUAL) ) {
                right2 = cv2._maxpoint();
            } else {
                right2 = construct_point_on_arc_2(common_right.x(),
                                                cv2.curve(),
                                                cv2.arcno(),
                                                cv2);
            }

        } else {
            right1 = cv1._maxpoint();
            right2 = cv2._maxpoint();
        }

        trimmed1 = cv1._trim(left1, right1);
        trimmed2 = cv2._trim(left2, right2);

        return joint;
    }

    //!@}

protected:
    //!\name Protected helper methods
    //!@{

    /*!\brief
     * function to ensure lexicographical order of the curve ends
     *
     * must be called once from constructor
     */
    void _fix_curve_ends_order() {
        CGAL::Comparison_result res =
            _same_arc_compare_xy(_minpoint(), _maxpoint());
        // curve ends cannot be identical
        CGAL_precondition(res != CGAL::EQUAL);
        if(res == CGAL::LARGER) { // swap curve ends and corresponding arcnos
            std::swap(this->ptr()->_m_min, this->ptr()->_m_max);
            std::swap(this->ptr()->_m_arcno_min, this->ptr()->_m_arcno_max);
        }
        // for non-vertical arcs check arcno constancy in the arc's interior
        // for vertical arcs check that there are no intersection points
        // between curve ends
        _check_arc_interior();
    }

    // p.curve() <-> p.arcno()
    // c <-> arcno_on_c
    /*!\brief
     * establishes preconditions that point \c pt lies on the curve
     * \c c with arc number \c arcno_on_c, also checks that point's supporting
     * curve and \c c are coprime
     *
     * \param pt Given point
     * \param arcno_on_c Arcno on curve
     * \param c Supporting curve
     */
  void _check_pt_arcno_and_coprimality(const Point_2&
                                         CGAL_precondition_code(pt),
                                       int CGAL_precondition_code(arcno_on_c),
                                       const Curve_analysis_2&
                                         CGAL_precondition_code(c)) const
  {
        CGAL_precondition_code(

        if (!c.is_identical(pt.curve())) {
            // -1 defines that no arcnos preconditions need to be established
            if (arcno_on_c != -1) {
                typename Curve_pair_analysis_2::Status_line_1
                    cpv_line;
                Curve_pair_analysis_2 cpa_2 =
                    Curved_kernel_via_analysis_2::instance().
                      kernel().construct_curve_pair_2_object()(pt.curve(), c);

                cpv_line = cpa_2.status_line_for_x(pt.x());
                CGAL_precondition(cpv_line.event_of_curve(pt.arcno(),
                                                          pt.curve())
                    == cpv_line.event_of_curve(arcno_on_c, c));
            }
            std::vector< Curve_analysis_2 > dummy[3];
            // ensure that curves are not decomposable
            CGAL_precondition(!Curved_kernel_via_analysis_2::instance().
                              kernel().decompose_2_object()(
                                      c, pt.curve(),
                                      std::back_inserter(dummy[0]),
                                      std::back_inserter(dummy[1]),
                                      std::back_inserter(dummy[2]))
            );
        } else if (arcno_on_c != -1) {
            CGAL_precondition(pt.arcno() == arcno_on_c);
        }
        );
    }

    /*!\brief
     * establishes preconditions to ensure that there are no event
     * points in the arc's interior (only at source and target) and its arc
     * number is constant
     *
     * \pre before calling this method source and target must be sorted
     * using \c _fix_curve_ends_order()
     */
    void _check_arc_interior() const {

#if !(defined(CGAL_KERNEL_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
        || defined(NDEBUG))

        if(is_vertical()) {
            Coordinate_1 x0 = _minpoint().x();
            typename Curve_analysis_2::Status_line_1 cv_line;
            cv_line = curve().status_line_for_x(x0);
            CGAL_precondition(cv_line.is_event());
            CGAL_precondition(cv_line.covers_line());

            // check that there are no intersections between min and max
            // curve ends
            bool inf_src =
                (_minpoint().location() == CGAL::ARR_BOTTOM_BOUNDARY);
            bool inf_tgt =
                (_maxpoint().location() == CGAL::ARR_TOP_BOUNDARY);
            // either no events over this line or the vertical line has at
            // least one finite end
            CGAL_precondition(cv_line.number_of_events() == 0 ||
                !(inf_src && inf_tgt));

            typename Curve_kernel_2::Compare_xy_2 cmp_xy(
                Curved_kernel_via_analysis_2::instance().
                    kernel().compare_xy_2_object());

            for(int k = 0; k < cv_line.number_of_events(); k++) {
            // TODO: replace by _compare_arc_numbers !! (Pavel)
          // no way since _compare_arc_numbers compares only against *this arc

              Coordinate_2 tmp(x0, curve(), k);
                bool res1 = true, res2 = true;
                if(!inf_src)
                    res1 = (cmp_xy(_minpoint().xy(), tmp, true) ==
                         CGAL::SMALLER);
                if(!inf_tgt)
                    res2 = (cmp_xy(tmp, _maxpoint().xy(), true) ==
                         CGAL::SMALLER);
                CGAL_precondition_msg(!(res1 && res2),
                  "Events are not allowed in the interior of a vertical arc!");
            }
            return;
        }

        typename Curve_analysis_2::Status_line_1 src_line, tgt_line,
            tmp;
        bool inf_src = (_minpoint().location() == CGAL::ARR_LEFT_BOUNDARY),
             inf_tgt = (_maxpoint().location() == CGAL::ARR_RIGHT_BOUNDARY);
        src_line = (inf_src ? curve().status_line_of_interval(0) :
            curve().status_line_for_x(_minpoint().x()));
        tgt_line = (inf_tgt ? curve().status_line_of_interval(
            curve().number_of_status_lines_with_event()) :
            curve().status_line_for_x(_maxpoint().x()));

        int src_idx = src_line.index(), tgt_idx = tgt_line.index(),
            diff = tgt_idx - src_idx;
        bool no_events_between = true;
        // it's supposed that arcs are not degenerate but lexicographic
        // order may not be established
        if(src_line.is_event())
            no_events_between = (tgt_line.is_event() ? (diff == 1) :
                (diff == 0)||(diff == 1));
        else
            no_events_between = (tgt_line.is_event() ? (diff == 0)||
                (diff == -1) : (diff == 0));

        if(!no_events_between) {
            // iterate through all events between source and target
            // to check that all events points lie above our arc
            int m_src_idx = src_idx + (src_line.is_event() ? 1 : 0),
                m_tgt_idx = tgt_idx - 1, low = m_src_idx, high = m_tgt_idx;
            int i, j;
            if(low > high) // do we need to check it ?
                std::swap(low, high);
            std::pair<int, int> ipair;
            for(i = low; i <= high; i++) {
                tmp = curve().status_line_at_event(i);
                for(j = 0; j < tmp.number_of_events(); j++) {
                    ipair = tmp.number_of_incident_branches(j);
                    if(ipair.first != 1||ipair.second != 1)
                        break;
                }
                // there must be at least one event and arcno() is not LARGER
                // than this event index
                CGAL_precondition(j < tmp.number_of_events() && arcno() <= j);
            }
        }
        // check validity of the curve-ends arcnos

        const typename Curved_kernel_via_analysis_2::
            Curve_interval_arcno_cache& map_interval_arcno =
            Curved_kernel_via_analysis_2::instance().interval_arcno_cache();

        if (src_line.is_event()) {
            CGAL_precondition(map_interval_arcno(src_line, 0,
                arcno()).first == this->ptr()->_m_arcno_min);
        } else {
            CGAL_precondition(arcno() == this->ptr()->_m_arcno_min);
        }
        if (tgt_line.is_event()) {
            CGAL_precondition(map_interval_arcno(tgt_line, 1,
                arcno()).first == this->ptr()->_m_arcno_max);
        } else {
            CGAL_precondition(arcno() == this->ptr()->_m_arcno_max);
        }
#endif
    }

    /*!\brief
     * compares y-coordinates of two arcs over an open (or closed)
     * interval or at exact x-coordinate
     *
     * \c where specifies whether to compare at negative/positive boundary or
     * at finite point. if \c where = ARR_INTERIOR \c perturb defines to
     * compare slightly to the left, on, or to the right of \c x0
     *
     * \param cv2 the second arc
     * \param where the location in parameter space
     * \param x0 The x-coordinate
     * \param perturb determines whether to pertub slightly to the left/right
     * \return the relative vertical alignment
     *
     * \pre !is_on_bottom_top(where)
     */
    CGAL::Comparison_result _compare_arc_numbers(
            const Kernel_arc_2& cv2,
            CGAL::Arr_parameter_space where,
            Coordinate_1 x0 = Coordinate_1(),
            CGAL::Sign perturb = CGAL::ZERO) const {

        CGAL_precondition(!is_on_bottom_top(where));
        CGAL_assertion(Kernel_arc_2_equals_Arc_2 ||
                       dynamic_cast< const Kernel_arc_2*>(this) != nullptr);
        Kernel_arc_2::simplify(*dynamic_cast< const Kernel_arc_2*>(this), cv2);
        if(curve().is_identical(cv2.curve()))
            return CGAL::sign(arcno() - cv2.arcno());
        return _compare_coprime(cv2, where, x0, perturb);
    }

    /*!\brief
     * computes vertical ordering of \c *this and \c cv2
     * having coprime supporting curves
     *
     * \param cv2 the second arc
     * \param where the location in parameter space
     * \param x0 The x-coordinate
     * \param perturb determines whether to pertub slightly to the left/right
     * \return the relative vertical alignment
     */
     CGAL::Comparison_result _compare_coprime(
            const Kernel_arc_2& cv2,
            CGAL::Arr_parameter_space where,
            Coordinate_1 x0,
            CGAL::Sign perturb) const {

#ifdef CKvA_DEBUG_PRINT_CERR
        CKvA_CERR("\n_compare_coprime; this: "
             << *dynamic_cast< const Kernel_arc_2*>(this)
             << "; g: " << cv2.curve().polynomial_2()
             << "; arcno_on_g: " << cv2.arcno() << "; where: " << where
        );
        if (where == CGAL::ARR_INTERIOR) {
            CKvA_CERR("; x = " << CGAL::to_double(x0));
        }
        CKvA_CERR("\n");
#endif

        typename Curve_pair_analysis_2::Status_line_1 cpv_line;
        Curve_pair_analysis_2 cpa_2 = Curved_kernel_via_analysis_2::instance().
                        kernel().construct_curve_pair_2_object()
                            (curve(), cv2.curve());

        if(where == CGAL::ARR_INTERIOR)
            cpv_line = cpa_2.status_line_for_x(x0, perturb);
        else
            cpv_line = cpa_2.status_line_of_interval(
                    // TODO don't mix up location (where) and finiteness!
                    where == CGAL::ARR_LEFT_BOUNDARY ? 0 :
                    cpa_2.number_of_status_lines_with_event());

        CGAL::Sign res =
            CGAL::sign(cpv_line.event_of_curve(arcno(), curve()) -
                       cpv_line.event_of_curve(cv2.arcno(), cv2.curve()));
        CKvA_CERR("result: " << res << "\n");
        return res;
    }

    /*\brief
     * internal comparison of two curve ends "lying" on the same arc
     *
     * since points are supposed to lie on the same arc, converging to the
     * boundary implies equality
     *
     * \param p first endpoint
     * \param q second endpint
     * \param equal_x \c true indicates to skip the comparison by x
     * \param only_x \c true indicates to report only the comparison by x
     * \returns the result of the queried comparison
     */
    CGAL::Comparison_result _same_arc_compare_xy(
            const Point_2& p,
            const Point_2& q,
            bool equal_x = false,
            bool only_x = false) const {

        CKvA_CERR("\n_same_arc_compare_xy; this: "
             << *dynamic_cast< const Kernel_arc_2*>(this)
             << "; p: " << p
             << "; q: " << q
             << "; equal_x: " << equal_x
             << "; only_x: " << only_x
             << "\n"
        );

        CGAL::Comparison_result res;

        if (p.is_identical(q)) {
            res = CGAL::EQUAL;
            CKvA_CERR("result1: " << res << "\n");
            return res;
        }

        CGAL::Arr_parameter_space locp = p.location(), locq = q.location();
        if (!equal_x || only_x) {

            if (!p.is_on_left_right() && !q.is_on_left_right()) {
                // both xs are finite: require x-comparisons
                res = Curved_kernel_via_analysis_2::instance().
                    compare_x_2_object()(p, q);
                if (res != CGAL::EQUAL) {
                    CKvA_CERR("result2: " << res << "\n");
                    return res;
                }
            } else if(locp != locq) {
                CGAL_assertion(p.is_on_left_right() || q.is_on_left_right());
                // at least one of the points lies at infty: suffice to cmp
                // boundaries
                if (locp == CGAL::ARR_LEFT_BOUNDARY) {
                    res = CGAL::SMALLER;
                    CKvA_CERR("result3: " << res << "\n");
                    return res;
                } else if (locp == CGAL::ARR_RIGHT_BOUNDARY) {
                    res = CGAL::LARGER;
                    CKvA_CERR("result4: " << res << "\n");
                    return res;
                } else if (locq == CGAL::ARR_LEFT_BOUNDARY) {
                    res = CGAL::LARGER;
                    CKvA_CERR("result5: " << res << "\n");
                    return res;
                } else if (locq == CGAL::ARR_RIGHT_BOUNDARY) {
                    res = CGAL::SMALLER;
                    CKvA_CERR("result6: " << res << "\n");
                    return res;
                }
            } // else: proceed to y-comparison
        }
        if (only_x) {
            res = CGAL::EQUAL;
            CKvA_CERR("result7: " << res << "\n");
            return res;
        }
        if (locp == locq) {
            if(locp != CGAL::ARR_INTERIOR) {
                res = CGAL::EQUAL; // both points are at the same inf in y
                CKvA_CERR("result8: " << res << "\n");
                return res;
            }
            // compare only y-values;
            res = Curved_kernel_via_analysis_2::instance().
                compare_xy_2_object()(p, q, true);
            CKvA_CERR("result9: " << res << "\n");
            return res;
        }
        // here: locp != locq && one of them is at inf y
        if (locp == CGAL::ARR_INTERIOR) {
            res = (locq == CGAL::ARR_BOTTOM_BOUNDARY ?
                   CGAL::LARGER : CGAL::SMALLER);
            CKvA_CERR("result10: " << res << "\n");
            return res;
        }
        // here: locp != locq && locp is at infty
        res = (locp == CGAL::ARR_BOTTOM_BOUNDARY ?
               CGAL::SMALLER : CGAL::LARGER);
        CKvA_CERR("result11: " << res << "\n");
        return res;
    }

    /*!\brief
     * min end-point of this arc (provided for code readability)
     *
     * \return min endpoint of arc (may lie on a boundary!)
     */
    inline
    const Point_2& _minpoint() const {
        return this->ptr()->_m_min;
    }

    /*!\brief
     * max end-point of this arc (provided for code readability)
     *
     * \return max endpoint of arc (may lie on a boundary!)
     */
    inline
    const Point_2& _maxpoint() const {
        return this->ptr()->_m_max;
    }

    /*!\brief
     * computes this arc's interval index
     *
     * \pre !is_vertical()
     */
    int _compute_interval_id() const {
        CGAL_precondition(!is_vertical());
        // we are interested in interval "to the right"
        CGAL::Arr_parameter_space min_loc = location(CGAL::ARR_MIN_END);
        bool min_has_x =
            (is_finite(CGAL::ARR_MIN_END) ||
             min_loc == CGAL::ARR_BOTTOM_BOUNDARY ||
             min_loc == CGAL::ARR_TOP_BOUNDARY);

        if (!min_has_x) {
            return 0;
        }
        // else
        typename Curve_analysis_2::Status_line_1 cv_line =
            curve().status_line_for_x(_minpoint().x(), CGAL::POSITIVE);
        return cv_line.index();
    }

    /*!\brief
     * computes this rational value in the interiors of the arc's x-range
     *
     * \pre !is_vertical()
     */
    Bound _compute_boundary_in_interval() const {
        CGAL_precondition(!is_vertical());
        // a curve end at negative boundary => 0th interval

        Bound res(0);

        typename Curve_kernel_2::Approximate_relative_1 approx_x;

        typename Curve_kernel_2::Bound_between_1 bound_between_x;

        CGAL::Arr_parameter_space min_loc = location(CGAL::ARR_MIN_END);
        bool min_has_x =
          (is_finite(CGAL::ARR_MIN_END) ||
           min_loc == CGAL::ARR_BOTTOM_BOUNDARY ||
           min_loc == CGAL::ARR_TOP_BOUNDARY);

        CGAL::Arr_parameter_space max_loc = location(CGAL::ARR_MAX_END);
        bool max_has_x =
          (is_finite(CGAL::ARR_MAX_END) ||
           max_loc == CGAL::ARR_BOTTOM_BOUNDARY ||
           max_loc == CGAL::ARR_TOP_BOUNDARY);

        if (min_has_x) {
          Coordinate_1 min_x = _minpoint().x();
          if (max_has_x) {
            Coordinate_1 max_x = _maxpoint().x();
            res = bound_between_x(min_x, max_x);
          } else {
            std::pair<Bound,Bound> min_pair=approx_x(min_x,4);
            res = min_pair.second + Bound(1);
          }
        } else {
          if (max_has_x) {
            Coordinate_1 max_x = _maxpoint().x();
            std::pair<Bound,Bound> max_pair=approx_x(max_x,4);
            res = max_pair.first - Bound(1);
          } else {
            // res stays 0
          }
        }
        CGAL_postcondition(is_in_x_range_interior(Coordinate_1(res)));
        return res;
    }

    /*!\brief
     * Replaces this arc's end-points by \c p1 and \c p2 with arcnos
     * \c arcno1 and \c arcno2.
     *
     * new curve ends are sorted lexicographical in case of need;
     * all preconditions must be checked by the caller
     *
     * \param p1 new first endpoint
     * \param p2 new second endpoint
     * \param arcno1 new first arcno (at \c p1)
     * \param arcno1 new second arcno (at \c p2)
     * \return pair whose first entry represent the refined arc, and whose
     *         second entry reports the lexicographic comparison of p1 and p2
     */
    std::pair< Kernel_arc_2, CGAL::Comparison_result >
    _replace_endpoints(
            const Point_2& p1, const Point_2& p2,
            int arcno1 = -1, int arcno2 = -1) const {

        CKvA_CERR("\n_replace_endpoints\n");

        Rep rep(*(this->ptr()));
        rep._m_min = p1;
        rep._m_max = p2;
        if (!is_vertical()) {
            if (arcno1 >= 0) {
                rep._m_arcno_min = arcno1;
            }
            if (arcno2 >= 0) {
                rep._m_arcno_max = arcno2;
            }
        }

        CGAL::Comparison_result cmp = _same_arc_compare_xy(p1,p2);
        if (cmp == CGAL::LARGER) {
            std::swap(rep._m_min, rep._m_max);
            std::swap(rep._m_arcno_min, rep._m_arcno_max);
        }
        /* no need to recompute location since they are set during
           construction of respective curve ends */
        rep._m_is_vertical = this->ptr()->_m_is_vertical;
        rep._m_left_to_right = this->ptr()->_m_left_to_right;

        rep._m_interval_id = boost::none;
        rep._m_boundary_in_interval = boost::none;

        return std::make_pair(Kernel_arc_2(rep), cmp);
    }

    /*!\brief
     * Simplifies representation of the arc !! DEPRECATED FUNCTION !!
     *
     * Given a decomposition of the arcs's supporting curve into a pair of two
     * curves \c cpa_2, we search for a curve this arc lies on and reset arc's
     * supporting curve and arcnos appropriately.
     *
     * \param cpa_2 analysis of curve pair that should be used
     *              in simplification
     * \pre \c cpa_2 must correspond to a decomposition of this arc's
     * supporting curve
     */
    void _simplify_by(const Curve_pair_analysis_2& cpa_2) const {

        typedef typename Curve_analysis_2::Polynomial_2 Polynomial_2;
        Polynomial_2 f = curve().polynomial_2();
        CGAL_precondition_code(
             Polynomial_2 mult = cpa_2.curve_analysis(0).polynomial_2() *
                    cpa_2.curve_analysis(1).polynomial_2();
             typename CGAL::Polynomial_traits_d<Polynomial_2>::Total_degree
                deg;
        );
        // common parts and full parts
        CGAL_precondition(CGAL::resultant(mult, f).degree() < 1);
        CGAL_precondition(mult.degree() == f.degree());
        CGAL_precondition(deg(mult) == deg(f));

        Coordinate_1 x0;
        if(is_vertical()) {
            // processing vertical arcs: search for supporting curve which has
            // vertical line at this x0 (must be exactly 1 curve)
            x0 = _minpoint().x();
            Curve_analysis_2 ca_2(cpa_2.curve_analysis(0));
            if(ca_2.status_line_for_x(x0).covers_line())
                this->ptr()->_m_support = ca_2;
            else {
                ca_2 = cpa_2.curve_analysis(1);
                CGAL_assertion(ca_2.status_line_for_x(x0).covers_line());
                this->ptr()->_m_support = ca_2;
            }
            return;
        }

        // processing non-vertical arcs
        typename Curve_pair_analysis_2::Status_line_1 cpv_line;
        std::pair<int, int> ipair;
        // preserve original supporting curve
        Curve_analysis_2 orig_curve(curve());

        // TODO do we mean location of is_finite?
        bool inf1_x = (_minpoint().location() == CGAL::ARR_LEFT_BOUNDARY);
        bool curve_idx;
        if(!inf1_x) {
            x0 = _minpoint().x();
            cpv_line = cpa_2.status_line_for_x(x0, CGAL::POSITIVE);
        } else
            cpv_line = cpa_2.status_line_of_interval(0);

        CGAL_precondition_code(
            typename Curve_analysis_2::Status_line_1
                cv_line = (inf1_x ? orig_curve.status_line_of_interval(0) :
                        orig_curve.status_line_for_x(x0, CGAL::POSITIVE));
        );
        CGAL_precondition(cpv_line.number_of_events() ==
            cv_line.number_of_events());

        { // search for new supporting curve and new arcno
            // since supporting curve was decomposed in two parts, arcno
            // represents y-position here
            ipair = cpv_line.curves_at_event(arcno());
            // this must be 1-curve event
            CGAL_assertion(!(ipair.first != -1&&ipair.second != -1));
            this->ptr()->_m_arcno = (ipair.first != -1 ? ipair.first :
                ipair.second);
            curve_idx = (ipair.first == -1);
            this->ptr()->_m_support = cpa_2.curve_analysis(curve_idx);
        }
        // search for source arcno
        /////////////// ATTENTION: this only holds for 2D plane topology !!
        ///////////////////////////////////////////////////////////////////
        // TODO do we mean location of is_finite?
        if(_minpoint().location() == CGAL::ARR_INTERIOR)  {

            cpv_line = cpa_2.status_line_for_x(x0);
            CGAL_precondition(cpv_line.number_of_events() ==
                    orig_curve.status_line_for_x(x0).number_of_events());
            ipair = cpv_line.curves_at_event(this->ptr()->_m_arcno_min);
            if(ipair.first != -1 && ipair.second != -1)
                // choose simpler supporting curve

              this->ptr()->_m_arcno_min = ((curve_idx) ?
                                           ipair.second : ipair.first);
            else {
                CGAL_assertion(ipair.first != -1||ipair.second != -1);
                this->ptr()->_m_arcno_min = (ipair.first != -1 ?
                    ipair.first : ipair.second);
            }
        } else // for infinite curve end arcno equals to interior arcno
            this->ptr()->_m_arcno_min = arcno();

        // search for new target arcno
        /////////////// ATTENTION: this only holds for 2D plane topology !!
        ///////////////////////////////////////////////////////////////////
        // TODO do we mean location of is_finite?
        if(_maxpoint().location() == CGAL::ARR_INTERIOR) {

            x0 = _maxpoint().x();
            cpv_line = cpa_2.status_line_for_x(x0);
            CGAL_precondition(cpv_line.number_of_events() ==
                    orig_curve.status_line_for_x(x0).number_of_events());

            ipair = cpv_line.curves_at_event(this->ptr()->_m_arcno_max);
            if(ipair.first != -1 && ipair.second != -1)
                // choose simpler supporting curve (the one which matches
                //  interior arcno)
                this->ptr()->_m_arcno_max = (curve_idx ?
                    ipair.second : ipair.first);
            else {
                CGAL_assertion(ipair.first != -1||ipair.second != -1);
                this->ptr()->_m_arcno_max = (ipair.first != -1 ?
                    ipair.first : ipair.second);
            }
        } else // for infinite curve end arcno equals to interior arcno
            this->ptr()->_m_arcno_max = arcno();

        // invalidate curve-specific data
        this->ptr()->_m_interval_id = boost::none;
        this->ptr()->_m_boundary_in_interval = boost::none;
    }
    //!@}

protected:
    //!\name Protected intersection methods
    //!@{

    /*!\brief
     * returns \c true if the two arcs \c *this and \c cv2 overlap,
     * overlapping part(s) are inserted to the output iterator \c oi
     * (of type \c Kernel_arc_2 ); if no overlapping parts found -
     * returns \c false
     *
     * \param cv2 The second arc
     * \param oi Report overlapping parts to this output iterator
     * \return \c true, if there was an overlap, \c false otherwise
     */
    template < class OutputIterator >
    bool _trim_if_overlapped(const Kernel_arc_2& cv2, OutputIterator oi) const
    {

        CKvA_CERR("\n_trim_if_overlapped: this: "
             << *dynamic_cast< const Kernel_arc_2*>(this) << "; and "
             << cv2 << "\n");
        // one arc is vertical and the other one is not, or x-ranges are not
        // overlapping => quit
        if (is_vertical() != cv2.is_vertical()) {
            return false;
        }

        if (is_vertical()) { // here process vertical case
            // check for x-coordinates equality
            if (Curved_kernel_via_analysis_2::instance().
                compare_x_2_object()(
                    _minpoint(),
                    cv2._minpoint()) != CGAL::EQUAL) {
                return false;
            }
            Kernel_arc_2::simplify(
                    *dynamic_cast< const Kernel_arc_2*>(this), cv2
            );
            // coprime support => no overlaps
            if(!curve().is_identical(cv2.curve()))
                return false;

            // LARGER source and smaller target
            Point_2 src = (_same_arc_compare_xy(_minpoint(), cv2._minpoint(),
                 true) == CGAL::LARGER ? _minpoint() : cv2._minpoint()),
                    tgt = (_same_arc_compare_xy(_maxpoint(), cv2._maxpoint(),
                 true)  == CGAL::SMALLER ? _maxpoint() : cv2._maxpoint());
            // vertical arcs do not overlap
            if(_same_arc_compare_xy(src, tgt, true) != CGAL::SMALLER)
                return false;
            // construct a common part
            *oi++ = (_replace_endpoints(src, tgt, -1, -1).first);
            return true;
        }
        // ask for joint x-range of two arcs
        // (LARGER source & smaller target curve ends)
        Point_2 src, tgt;
        if (!_joint_x_range(cv2, src, tgt)) {
            return false;
        }

        if (curve().is_identical(cv2.curve())) {
            if(arcno() != cv2.arcno()) // arcnos are not equal => no overlaps
                return false;
            int a_min = (src.is_on_left_right() ? -1 : arcno(src.x())),
                a_max = (tgt.is_on_left_right() ? -1 : arcno(tgt.x()));
            // construct a common  part
            *oi++ = _replace_endpoints(src, tgt, a_min, a_max).first;
            return true;
        }

        // we are left with two non-vertical arcs whose supporting curves
        // are different => look for overlapping parts of the curves
        typedef std::vector<std::pair<Curve_analysis_2, int> >
            Curve_arcno_container;
        typedef std::vector<Curve_analysis_2> Curve_container;
        Curve_container parts_f, parts_g, common;

        if (!Curved_kernel_via_analysis_2::instance().
            kernel().decompose_2_object()(
                    curve(), cv2.curve(),
                    std::back_inserter(parts_f),
                    std::back_inserter(parts_g),
                    std::back_inserter(common))) {
            return false; // supporting curves are coprime => quit
        }
        Coordinate_1 x0;
        bool yes = false, inf_x = src.is_on_left_right();
        if(!inf_x) // choose a target x-coordinate from the joint x-range
            x0 = src.x();
        std::pair<int, int> ipair;
        Curve_pair_analysis_2 cpa_2;
        Curve_arcno_container found, overlaps;

        CKvA_CERR("_trim_if_overlapped: non-coprime supporting curves\n");

        typename Curve_pair_analysis_2::Status_line_1 cpv_line;
        // iterate to find all overlapping parts
        typename Curve_container::const_iterator it_parts, it_com;
        for (it_com = common.begin(); it_com != common.end(); it_com++) {
            for(it_parts = parts_f.begin(); it_parts != parts_f.end();
                    it_parts++) {

                cpa_2 = Curved_kernel_via_analysis_2::instance().
                        kernel().construct_curve_pair_2_object()
                            (*it_com, *it_parts);
                cpv_line = (inf_x ? cpa_2.status_line_of_interval(0) :
                    cpa_2.status_line_for_x(x0, CGAL::POSITIVE));
                // no intersections at this curve pair => skip it
                if(arcno() >= cpv_line.number_of_events())
                    continue;
                ipair = cpv_line.curves_at_event(arcno(),*it_com,*it_parts);
                // this must be 1-curve event: is this true ???
                CGAL_assertion(!(ipair.first != -1&&ipair.second != -1));
                if(ipair.first != -1) // lies on a common part
                    found.push_back(std::make_pair(*it_com, ipair.first));
            }
        }

        // now iterate over all "suspicious" common parts to find real overlaps
        typename Curve_arcno_container::const_iterator it_found;
        for (it_found = found.begin(); it_found != found.end(); it_found++) {
            for (it_parts = parts_g.begin(); it_parts != parts_g.end();
                 it_parts++) {

                cpa_2 = Curved_kernel_via_analysis_2::instance().
                        kernel().construct_curve_pair_2_object()
                            (it_found->first, *it_parts);

                cpv_line = (inf_x ? cpa_2.status_line_of_interval(0) :
                    cpa_2.status_line_for_x(x0, CGAL::POSITIVE));
                // no intersections at this curve pair => skip it
                if(cv2.arcno() >= cpv_line.number_of_events())
                    continue;
                ipair = cpv_line.curves_at_event(cv2.arcno(),
                                                 it_found->first,
                                                 *it_parts);
                // this must be 1-curve event: is this true ???
                CGAL_assertion(!(ipair.first != -1&&ipair.second != -1));
                if(ipair.first == -1 || ipair.first == it_found->second)
                    continue;
                // lies on a common part and arcnos are the same: VUALA!!!
                // here we need to "clip" [src.x(), tgt.x()] w.r.t. the
                // defining x-range of a common part *it_found.. how ?
                yes = true; // we've got it!
                // now construct a common arc
                Rep rep(*(this->ptr()));
                rep._m_min = src;
                rep._m_max = tgt;
                rep._m_support = it_found->first;
                rep._m_arcno = it_found->second;
                rep._m_arcno_min = rep._m_arcno_max = rep._m_arcno;

                if(!inf_x) {
                    int a = arcno(src.x());
                    if(a != arcno()) {
                        cpv_line = cpa_2.status_line_for_x(src.x());
                        ipair = cpv_line.curves_at_event(a,
                                                         it_found->first,
                                                         *it_parts);
                        // should ultimately lie on the common curve ?
                        CGAL_assertion(ipair.first != -1);
                        rep._m_arcno_min = ipair.first;
                    }
                }
                if(!tgt.is_on_left_right()) {
                    int a = arcno(tgt.x());
                    if(a != arcno()) {
                        cpv_line = cpa_2.status_line_for_x(tgt.x());
                        ipair = cpv_line.curves_at_event(a,
                                                         it_found->first,
                                                         *it_parts);
                        // should ultimately lie on the common curve ?
                        CGAL_assertion(ipair.first != -1);
                        rep._m_arcno_max = ipair.first;
                    }
                }
                *oi++ = Kernel_arc_2(rep);
            }
        }
        return yes;
    }

    /*!\brief
     * computes zero-dimensional intersections of \c cv1 with \c cv2.
     *
     * Intersection points
     * are inserted to the output iterator \c oi as objects of type
     * <tt>std::pair<Point_2, unsigned int></tt> (intersection point +
     * multiplicity)
     *
     * \param cv1 the first arc
     * \param cv2 the second arc
     * \param oi reporting zero-dimensional intersections through this output
     *        iterator
     * \pre !cv1.do_overlap()
     */
    template < class OutputIterator >
    static OutputIterator _intersection_points(
            const Kernel_arc_2& cv1, const Kernel_arc_2& cv2,
            OutputIterator oi) {

        // handle a special case when two arcs are supported by the same
        // curve => only end-point intersections

        CKvA_CERR("\nintersection_points\n");
        Kernel_arc_2::simplify(cv1, cv2);
        if (cv1.curve().is_identical(cv2.curve())) {
            return _intersect_at_endpoints(cv1, cv2, oi);
        }

        // else general case: distinct supporting curves
        return _intersect_coprime_support(cv1, cv2, oi);
    }

    /*!\brief
     * computes intersection of two arcs meeting only at their curve ends.
     *
     * Intersection points are returned in the output interator \c oi as object
     * of type std::pair<Point_2, int> (intersection + multiplicity)
     *
     * \param cv1 the first arc
     * \param cv2 the second arc
     * \param oi reporting zero-dimensional intersections through this output
     *        iterator
     *
     */
    template < class OutputIterator >
    static OutputIterator _intersect_at_endpoints(const Kernel_arc_2& cv1,
                                                  const Kernel_arc_2& cv2,
                                                  OutputIterator oi) {

        CKvA_CERR("\n_intersect_at_endpoints\n");

        CGAL_precondition(!cv1.do_overlap(cv2));
        /* Since *this and cv2 do not overlap and cannot contain singularities
         * in the interior, the only remaining candidates for intersections are
         * their finite endpoints (if any), for vertical arcs as well.
         */
        /*CGAL::Bound_type bnd_x, bnd_y,
            bnd1_x = cv2.boundary_in_x(CGAL::ARR_MIN_END),
            bnd1_y = cv2.boundary_in_y(CGAL::ARR_MIN_END),
            bnd2_x = cv2.boundary_in_x(CGAL::ARR_MAX_END),
            bnd2_y = cv2.boundary_in_y(CGAL::ARR_MAX_END);*/

        // TODO do we mean location of is_finite?
        bool f2_min = (cv2._minpoint().location() == CGAL::ARR_INTERIOR),
             f2_max = (cv2._maxpoint().location() == CGAL::ARR_INTERIOR);
        if(!(f2_min || f2_max)) // neither of curve ends is finite =>
            return oi;          // no intersections

        Point_2 pt;

        CGAL::Arr_curve_end end = CGAL::ARR_MIN_END;

        while(1) {
            CGAL::Arr_parameter_space loc = cv1.location(end);
            //bnd_x = boundary_in_x(end), bnd_y = boundary_in_y(end);
            if(loc != CGAL::ARR_INTERIOR)
                goto Lendloop;
            pt = cv1.curve_end(end);
            // easy case: intersection at singularity doesn't require to
            // compare x/y-coordinates
            /*if(is_singular(bnd_x)) {
                if(bnd1_x == bnd_x || bnd2_x == bnd_x)
                    *oi++ = std::make_pair(pt, 0);

            } else if(is_singular(bnd_y)) {
                if(bnd1_y == bnd_y || bnd2_y == bnd_y)
                    *oi++ = std::make_pair(pt, 0);

            } else if(is_on_disc(bnd_x)) {

    // CONFUSION: if bndx != bnd1_x should we compare ys at -oo
    // or at +oo ? or is this true for discontinuity:
    // 0th interval == the last interval ? (i.e. intervals are mirrored ?)
    // what if both conditions are satisfied at a time ? duplicates ?
                if(bnd1_x == CGAL::AFTER_DISCONTINUITY &&
                    _compare_arc_numbers(cv2, bnd1_x) == CGAL::EQUAL)
                    *oi++ = std::make_pair(pt, 0);

                if(bnd2_x == CGAL::BEFORE_DISCONTINUITY &&
                    _compare_arc_numbers(cv2, bnd2_x) == CGAL::EQUAL)
                    *oi++ = std::make_pair(pt, 0);

            } else if(is_on_disc(bnd_y)) {
                  // disc in y: compare only x-coordinates !
    // what if both conditions are satisfied at a time ? duplicates ?

                if(bnd1_y == CGAL::AFTER_DISCONTINUITY &&
                    kernel_2.compare_1_object()(pt.x(), _minpoint().x()) ==
                        CGAL::EQUAL)
                    *oi++ = std::make_pair(pt, 0);

                if(bnd2_y == CGAL::BEFORE_DISCONTINUITY &&
                    kernel_2.compare_1_object()(pt.x(), _maxpoint().x()) ==
                        CGAL::EQUAL)
                    *oi++ = std::make_pair(pt, 0);
              // ordinar normal case:
              // selection is exclusive since arcs cannot intersect twice
              // at the same finite end-point
              } else*/ if((f2_min && pt == cv2._minpoint()) ||
                          (f2_max && pt == cv2._maxpoint())) {
                  *oi++ = std::make_pair(pt, 0);
              }
        Lendloop:
            if (end == CGAL::ARR_MAX_END) {
                break;
            }
            end = CGAL::ARR_MAX_END;
        }
        return oi;
    }

    /*!\brief
     * computes a joint x-range of two arcs and returns \c true
     * if arcs' x-ranges overlap; otherwise returns \c false
     *
     * \param cv2 The second arc
     * \param pt_low Output: Point indicating the lower bound of the the joint
     *        x-range
     * \param pt_high Output: Point indicating the upper bound of the the joint
     *        x-range
     * \return \c true, if arcs overlap, \c false otherwise
     *
     * \pre both arcs are not vertical
     */
    bool _joint_x_range(const Kernel_arc_2& cv2, Point_2& pt_low,
                        Point_2& pt_high) const {

        CKvA_CERR("\n_joint_x_range\n");

        CGAL_precondition(!is_vertical());
        CGAL_precondition(!cv2.is_vertical());

        Point_2 pt1 = _minpoint(), pt2 = cv2._minpoint();
        Point_2 low = pt2, high;
        // find intersection x-range: larger source & smaller target
        if (pt1.location() != CGAL::ARR_LEFT_BOUNDARY) {
            if (pt2.location() != CGAL::ARR_LEFT_BOUNDARY) {
                low = (Curved_kernel_via_analysis_2::instance().
                       compare_x_2_object()(pt1, pt2) ==
                       CGAL::LARGER ? pt1 : pt2);
            } else {
                low = pt1;
            }
        }
        pt1 = _maxpoint(), pt2 = cv2._maxpoint(), high = pt2;
        if (pt1.location() != CGAL::ARR_RIGHT_BOUNDARY) {
            if(pt2.location() != CGAL::ARR_RIGHT_BOUNDARY) {
                high = (Curved_kernel_via_analysis_2::instance().
                        compare_x_2_object()(pt1, pt2) ==
                        CGAL::SMALLER ? pt1 : pt2);
            } else {
                high = pt1;
            }
        }
        if (!low.is_on_left_right() && !high.is_on_left_right() &&
            Curved_kernel_via_analysis_2::instance().
            compare_x_2_object()(low, high) !=
            CGAL::SMALLER) {// disjoint x-ranges
            return false;
        }
        pt_low = low;
        pt_high = high;

        return true;
    }

    /*!\brief
     * computes zero-dimensional
     * intersections of two arcs having coprime supporting curves
     *
     * intersection points are inserted to the output iterator \c oi as objects
     * of type <tt>std::pair<Point_2, unsigned int></tt> (intersection point +
     * multiplicity)
     *
     * \param cv1 the first arc
     * \param cv2 the second arc
     * \param oi reporting zero-dimensional intersections through this output
     *        iterator
     */
    template <class OutputIterator>
    static OutputIterator _intersect_coprime_support(const Kernel_arc_2& cv1,
                                                     const Kernel_arc_2& cv2,
                                                     OutputIterator oi) {
        // vertical arcs: the interesting case is when only one of the arcs is
        // vertical - otherwise there is no intersection (different x-coords),
        // or they overlap (not allowed), or they touch at the end-points
        // (already tested)

        CKvA_CERR("\n_intersect_coprime_support: " << cv1 <<
            " and " << cv2 << "\n");

        if (cv1.is_vertical() || cv2.is_vertical()) {
            CGAL_assertion(cv1.is_vertical() != cv2.is_vertical());
            // due to coprimality condition, supporting curves are different =>
            // they have no common vertical line therefore there is no
            // intersection
            const Kernel_arc_2& vert = (cv1.is_vertical() ? cv1 : cv2),
                nonvert = (cv1.is_vertical() ? cv2 : cv1);
            Coordinate_1 x = vert._minpoint().x();
            // vertical arc does not lie within another arc's x-range => no
            // intersections
            if (!nonvert.is_in_x_range(x)) {
                return oi;
            }
            typename Curved_kernel_via_analysis_2:: Construct_point_on_arc_2
                construct_point_on_arc =
                Curved_kernel_via_analysis_2::instance().
                construct_point_on_arc_2_object();


            Point_2 xy = construct_point_on_arc(
                    x, nonvert.curve(), nonvert.arcno(x), nonvert
            );
            if (vert.compare_y_at_x(xy) == CGAL::EQUAL) {
                *oi++ = std::make_pair(xy, 1);
            }
            return oi;
        }

        Point_2 low_x, high_x;
        // x-ranges are disjoint => nothing to do
        if (!cv1._joint_x_range(cv2, low_x, high_x)) {
            return oi;
        }
        bool inf_low = low_x.is_on_left_right(),
            inf_high = high_x.is_on_left_right();
        Curve_analysis_2 f = cv1.curve(), g = cv2.curve();
        Curve_pair_analysis_2 cpa_2 =
            Curved_kernel_via_analysis_2::instance().
                kernel().construct_curve_pair_2_object()(f, g);
        int low_idx = 0,
            high_idx = cpa_2.number_of_status_lines_with_event()-1;

        bool index_at_event_min=false;
        bool index_at_event_max=false;

        typename Curve_pair_analysis_2::Status_line_1 line;
        if(!inf_low) {
            line = cpa_2.status_line_for_x(low_x.x());
            low_idx = line.index();
            index_at_event_min=line.is_event();
            if(index_at_event_min) {
                if((cv1._minpoint().is_on_bottom_top() &&
                    low_x.x() == cv1._minpoint().x()) ||
                   (cv2._minpoint().is_on_bottom_top() &&
                    low_x.x() == cv2._minpoint().x())) {
                 // hack: no intersection with asymptotic end
                  low_idx++;
                  index_at_event_min=false;
                }
            }
        }

        if(!inf_high) {
            line = cpa_2.status_line_for_x(high_x.x());
            high_idx = line.index();
            index_at_event_max=line.is_event();
            if(!index_at_event_max) {
              high_idx--;
            } else if((cv1._maxpoint().is_on_bottom_top() &&
                high_x.x() == cv1._maxpoint().x()) ||
                (cv2._maxpoint().is_on_bottom_top() &&
                 high_x.x() == cv2._maxpoint().x())) {
              // hack: no intersection with asymptotic end
              high_idx--;
              index_at_event_max=false;
            }
        }

        // run over all event points within the joint x-range of two arcs
        // looking whether a particular event is made of both curves, i.e.,
        // grabbing all 2-curve events
        int arcno1, arcno2, mult;

        typename CGAL::Polynomial_traits_d<
            typename Curve_kernel_2::Polynomial_2>::Total_degree deg;

        bool which_curve = (deg(f.polynomial_2()) < deg(g.polynomial_2()));
        for(int i = low_idx; i <= high_idx; i++) {
            typename Curve_pair_analysis_2::Status_line_1 tmp =
                cpa_2.status_line_at_event(i);
            if(!tmp.is_intersection())
                continue;

            Coordinate_1 x0 = tmp.x();
            if((i == low_idx && index_at_event_min) ||
               (i == high_idx && index_at_event_max)) {
                arcno1 = cv1.arcno(x0);
                arcno2 = cv2.arcno(x0);
                mult = 0; // intersection at end-point
            } else {
                arcno1 = cv1.arcno();
                arcno2 = cv2.arcno();
                mult = -1; // need to compute
            }

            int pos = tmp.event_of_curve(arcno1, f);
            if (pos != tmp.event_of_curve(arcno2, g)) {
                continue;
            }
            if (mult == -1) {
                mult = tmp.multiplicity_of_intersection(pos);
            }

            // pick up the curve with lower degree
            typename Curved_kernel_via_analysis_2::Construct_point_on_arc_2
                construct_point_on_arc =
                Curved_kernel_via_analysis_2::instance().
                construct_point_on_arc_2_object();

            if (which_curve) {
                Point_2 p = construct_point_on_arc(
                        x0, cv1.curve(), arcno1, cv1
                );
                *oi++ = std::make_pair(p, mult);
            } else {
                Point_2 p = construct_point_on_arc(
                        x0, cv2.curve(), arcno2, cv2
                );
                *oi++ = std::make_pair(p, mult);
            }
        }
        return oi;
    }

    #undef CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_ARC
    //!@}


    //!\name Approximation
    //!@{

private:

  std::pair<double,double> y_interval_for_curve_end(
      const Arc_2& arc,
      CGAL::Arr_curve_end end,
      long prec)
    const {

    double PINF = std::numeric_limits<double>::infinity();
    double MINF = -PINF;

    std::pair< Bound, Bound > y_approx;

    switch (this->location(end)) {

    case (CGAL::ARR_TOP_BOUNDARY): {
      return std::make_pair(PINF, PINF); // early exit
    }
    case (CGAL::ARR_BOTTOM_BOUNDARY): {
      return std::make_pair(MINF, MINF); // early exit
    }
    case(CGAL::ARR_LEFT_BOUNDARY):
    case(CGAL::ARR_RIGHT_BOUNDARY): {

      CGAL::Object obj = this->curve().asymptotic_value_of_arc(
          this->location(end), this->arcno()
      );

      CGAL::Arr_parameter_space ps;
      Coordinate_1 asym_info;

      if (CGAL::assign(ps, obj)) {
        if (ps == CGAL::ARR_BOTTOM_BOUNDARY) {
          return std::make_pair(MINF, MINF); // early exit
        } else {
          CGAL_assertion(ps == CGAL::ARR_TOP_BOUNDARY);
          return std::make_pair(PINF, PINF); // early exit
        }

      } else {

        CGAL_assertion_code(bool check =)
          CGAL::assign(asym_info, obj);
        CGAL_assertion(check);

        y_approx = Curved_kernel_via_analysis_2::instance().
          kernel().approximate_absolute_1_object()(asym_info, prec);
      }
      break;
    }
    case (CGAL::ARR_INTERIOR): {

      y_approx =
        Curved_kernel_via_analysis_2::instance().
        kernel().approximate_absolute_y_2_object()(
            arc.curve_end(end).xy(), prec
        );
      break;

    }
    default: {
       CGAL_error();
    }
    } // switch

    return std::make_pair(CGAL::to_double(y_approx.first),
                          CGAL::to_double(y_approx.second));
  }

public:

    /*!\brief
     * bbounding box for arc
     */
    CGAL::Bbox_2 bbox() const {
      if (!this->ptr()->_m_bbox) {

        double PINF = std::numeric_limits<double>::infinity();
        double MINF = -PINF;

        double xmin;
        double xmax;

        double ymin = PINF; // correct as we modify shrink below
        double ymax = MINF; // correct as we modify shrink below

        // TODO choose precision
        long prec = 53;

        std::pair< double, double > y_dapprox;

        // left end

        // xmin for left
        if (this->location(CGAL::ARR_MIN_END) == CGAL::ARR_INTERIOR ||
            this->location(CGAL::ARR_MIN_END) == CGAL::ARR_BOTTOM_BOUNDARY ||
            this->location(CGAL::ARR_MIN_END) == CGAL::ARR_TOP_BOUNDARY) {

          std::pair< Bound, Bound > x_approx =
            Curved_kernel_via_analysis_2::instance().
            kernel().approximate_absolute_1_object()(
                this->curve_end_x(CGAL::ARR_MIN_END), prec
            );

          xmin = CGAL::to_double(x_approx.first);

        } else {

          // left end can only lie on LEFT BOUNDARY
          xmin = MINF;

        }

        // ymin/ymax for left
        y_dapprox = y_interval_for_curve_end(*this, CGAL::ARR_MIN_END, prec);

        // adapt y-interval
        ymin = (CGAL::min)(ymin, y_dapprox.first);
        ymax = (CGAL::max)(ymax, y_dapprox.second);

        // right end

        // xmax for right
        if (this->location(CGAL::ARR_MAX_END) == CGAL::ARR_INTERIOR ||
            this->location(CGAL::ARR_MAX_END) == CGAL::ARR_BOTTOM_BOUNDARY ||
            this->location(CGAL::ARR_MAX_END) == CGAL::ARR_TOP_BOUNDARY) {

          std::pair< Bound, Bound > x_approx =
            Curved_kernel_via_analysis_2::instance().
            kernel().approximate_absolute_1_object()(
                this->curve_end_x(CGAL::ARR_MAX_END), prec
            );

          xmax = CGAL::to_double(x_approx.second);

        } else {

          // right end can only lie on RIGHT BOUNDARY
          xmax = PINF;

        }

        // ymin/ymax for right
        y_dapprox = y_interval_for_curve_end(*this, CGAL::ARR_MAX_END, prec);

        // adapt y-interval
        ymin = (CGAL::min)(ymin, y_dapprox.first);
        ymax = (CGAL::max)(ymax, y_dapprox.second);

        // search local extrema on a non-vertical arc

        if (!this->is_vertical()) {

          // TODO remove algebraic notation (Mult, Solve_2)

          typedef typename Curve_kernel_2::Multiplicity_type Multiplicity_type;

          std::vector< std::pair< Coordinate_2, Multiplicity_type > > pts;

          // TODO can this->curve()->polynomial be not squarefree?
          Curved_kernel_via_analysis_2::instance().
            kernel().solve_2_object()(
                this->curve(),
                Curved_kernel_via_analysis_2::instance().
                kernel().construct_curve_2_object()(
                    CGAL::differentiate(this->curve().polynomial_2(),0)
                    // ^^ f_x
                ),
                std::back_inserter(pts)
            );

          int n = static_cast<int>(pts.size());
          CKvA_CERR("check candidates for y-extremal points: #" << n );
          for (int i = 0; i < n; i++) {

            const Coordinate_2& curr_xy = pts[i].first;


#if 0
            // EBEB: Disabled this test as curr_xy's curve
            //       is not guaranteed to be wrt this->curve()

            CKvA_CERR("check if arcnos match: " <<
                 curr_xy << "; arc = " << *this << "\n\n");
            // this is the simpler test, thus we evaluate it first
            if (this->arcno() == curr_xy.arcno()) {
              CKvA_CERR("check if x-coordinate lies in interior: " <<
                   curr_xy << "; arc = " << *this << "\n\n");
              // this is the more sophisticated test, thus second
              if (this->is_in_x_range_interior(curr_xy.x())) {
                CKvA_CERR("update y coordinates");

                std::pair< Bound, Bound > xy_approx =
                  Curved_kernel_via_analysis_2::instance().
                  kernel().approximate_absolute_y_2_object()
                  (curr_xy, prec);

                // adapt y-interval
                ymin = (CGAL::min)(ymin,
                                 CGAL::to_double(xy_approx.first));
                ymax = (CGAL::max)(ymax,
                                 CGAL::to_double(xy_approx.second));
              }
            }
#else
            // this is the more sophisticated test, thus second
            if (this->is_in_x_range_interior(curr_xy.x())) {
              // TODO replace with is_on
              Point_2 curr_pt =
                Curved_kernel_via_analysis_2::instance().
                construct_point_2_object()(curr_xy.x(),
                                           curr_xy.curve(),
                                           curr_xy.arcno());
              if (this->compare_y_at_x(curr_pt) == CGAL::EQUAL) {

                CKvA_CERR("update y coordinates");

                std::pair< Bound, Bound > xy_approx =
                  Curved_kernel_via_analysis_2::instance().
                  kernel().approximate_absolute_y_2_object()
                  (curr_xy, prec);

                // adapt y-interval
                ymin = (CGAL::min)(ymin,
                                 CGAL::to_double(xy_approx.first));
                ymax = (CGAL::max)(ymax,
                                 CGAL::to_double(xy_approx.second));
              }
            }
#endif
          }
        }
        this->ptr()->_m_bbox = CGAL::Bbox_2(xmin, ymin, xmax, ymax);
      }

      return *(this->ptr()->_m_bbox);
    }

    //!@}


public:
    //!\name IO
    //!@{

    /*!\brief
     * output operator
     *
     * write arc to \c os
     */
    void write(std::ostream& os) const {

        switch (::CGAL::get_mode(os)) {
        case ::CGAL::IO::PRETTY:
            os << "arc@" << this->id() << "[(sup@" << this->curve().id();
            if (this->is_vertical()) {
              os << ", VERTICAL at x = "
                 << CGAL::to_double(curve_end_x(ARR_MAX_END)) ;
            } else {
                os << ", ARCNO=" << this->arcno(CGAL::ARR_MIN_END)
                   << "," << this->arcno()
                   << "," << this->arcno(CGAL::ARR_MAX_END);
            }
            os << "; l2r: " << this->is_left_to_right() << "); \n";
            os <<"min: " << this->_minpoint() << ";\n";
            os << "max: " << this->_maxpoint() << "]";

            break;

        case ::CGAL::IO::BINARY:
          std::cerr << "BINARY format not yet implemented" << std::endl;
        break;
        default:
          // ASCII
          os << "Arc_2(";
          os << this->ptr()->_m_min;
          os << ",";
          os << this->ptr()->_m_max;
          os << ",";
          os << this->ptr()->_m_support;
          os << ",";
          os << this->ptr()->_m_arcno;
          os << ",";
          os << this->ptr()->_m_arcno_min;
          os << ",";
          os << this->ptr()->_m_arcno_max;
          os << ",";
          os << this->ptr()->_m_is_vertical;
          os << ",";
          os << this->ptr()->_m_left_to_right;
          os << ")";
        }
    }


    /*!\brief
     * input operator
     *
     * read arc from \c is
     */
    void read(std::istream& is) {

      CGAL_precondition(CGAL::is_ascii(is));

      Rep rep;

      // read "Arc_2("
      swallow(is, 'A');
      swallow(is, 'r');
      swallow(is, 'c');
      swallow(is, '_');
      swallow(is, '2');
      swallow(is, '(');

      Point_2 min, max;

      // read values
      is >> rep._m_min;
      swallow(is, ',');
      is >> rep._m_max;
      swallow(is, ',');
      is >> rep._m_support;
      swallow(is, ',');
      is >> rep._m_arcno;
      swallow(is, ',');
      is >> rep._m_arcno_min;
      swallow(is, ',');
      is >> rep._m_arcno_max;
      swallow(is, ',');
      is >> rep._m_is_vertical;
      swallow(is, ',');
      is >> rep._m_left_to_right;

      // read the ')'
      swallow(is, ')');

      *this = Arc_2< Curved_kernel_via_analysis_2, Rep >(rep);
    }

    //!@}

    //! equality
    inline
    bool operator == (const Kernel_arc_2& arc2) const {
        return  is_equal(arc2);
    }

#if defined(_MSC_VER)
    // befriending the kernel point
    friend typename Curved_kernel_via_analysis_2::Point_2;

    // befriending the kernel arc
    friend typename Curved_kernel_via_analysis_2::Arc_2;

    // befriending the functors
#define CGAL_BEFRIEND_CKvA_2_FUNCTOR(Z) \
    friend typename Curved_kernel_via_analysis_2::Z; \
    friend typename Curved_kernel_via_analysis_2_Functors:: \
        Z<Curved_kernel_via_analysis_2>
#else // defined(_MSC_VER) || defined(__clang__) || defined(__INTEL_COMPILER)
    // befriending the kernel point
    //friend class Curved_kernel_via_analysis_2::Point_2;

    // befriending the kernel arc
    //friend class Curved_kernel_via_analysis_2::Arc_2;

    // befriending the functors
#define CGAL_BEFRIEND_CKvA_2_FUNCTOR(Z) \
    friend class Curved_kernel_via_analysis_2_Functors:: \
        Z<Curved_kernel_via_analysis_2>
#endif // defined(_MSC_VER) || defined(__clang__) || defined(__INTEL_COMPILER)


//Curved_kernel_via_analysis_2_functors<
  //              Curved_kernel_via_analysis_2> >;

    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Construct_arc_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Is_vertical_2);

    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Construct_min_vertex_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Construct_max_vertex_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Compare_y_at_x_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Compare_y_at_x_left_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Compare_y_at_x_right_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Is_in_x_range_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Equal_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Do_overlap_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Intersect_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Trim_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Split_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Are_mergeable_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Merge_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Is_on_2);

    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Parameter_space_in_x_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Compare_y_near_boundary_2);

    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Parameter_space_in_y_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Compare_x_at_limit_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Compare_x_near_limit_2);

#undef CGAL_BEFRIEND_CKvA_2_FUNCTOR

private:

    // type of CurveSweepTraits model
    typedef CGAL::Sweep_curves_adapter_2< Curved_kernel_via_analysis_2 > SCA_2;
    // befriend segment for Self::_intersection_points
    friend class internal::Generic_arc_2<SCA_2>;

    /*
    // befriend all functors
#define CGAL_BEFRIEND_SWEEP_CURVES_ADAPTER_2_FUNCTOR(Z) \
    friend class CGAL::Sweep_curves_functors::Z< SCA_2 >; \

    CGAL_BEFRIEND_SWEEP_CURVES_ADAPTER_2_FUNCTOR(Compare_xy_2)
    CGAL_BEFRIEND_SWEEP_CURVES_ADAPTER_2_FUNCTOR(Less_xy_2)
    CGAL_BEFRIEND_SWEEP_CURVES_ADAPTER_2_FUNCTOR(Compare_y_at_x_2)
    CGAL_BEFRIEND_SWEEP_CURVES_ADAPTER_2_FUNCTOR(Equal_y_at_x_2)
    CGAL_BEFRIEND_SWEEP_CURVES_ADAPTER_2_FUNCTOR(Multiplicity_of_intersection_2)
    CGAL_BEFRIEND_SWEEP_CURVES_ADAPTER_2_FUNCTOR(Compare_y_right_of_point_2)
    CGAL_BEFRIEND_SWEEP_CURVES_ADAPTER_2_FUNCTOR(Source_2)
    CGAL_BEFRIEND_SWEEP_CURVES_ADAPTER_2_FUNCTOR(Target_2)
    CGAL_BEFRIEND_SWEEP_CURVES_ADAPTER_2_FUNCTOR(Construct_segment_2)
    CGAL_BEFRIEND_SWEEP_CURVES_ADAPTER_2_FUNCTOR(Is_degenerate_2)
    CGAL_BEFRIEND_SWEEP_CURVES_ADAPTER_2_FUNCTOR(Do_overlap_2)
    CGAL_BEFRIEND_SWEEP_CURVES_ADAPTER_2_FUNCTOR(New_endpoints_2)
    CGAL_BEFRIEND_SWEEP_CURVES_ADAPTER_2_FUNCTOR(New_endpoints_opposite_2)
    CGAL_BEFRIEND_SWEEP_CURVES_ADAPTER_2_FUNCTOR(Intersect_2)
    CGAL_BEFRIEND_SWEEP_CURVES_ADAPTER_2_FUNCTOR(Intersect_right_of_point_2)
    CGAL_BEFRIEND_SWEEP_CURVES_ADAPTER_2_FUNCTOR(Make_x_monotone_2)

    #undef CGAL_BEFRIEND_SWEEP_CURVES_ADAPTER_2_FUNCTOR */

}; // class Arc_2

/*!\relates Arc_2
 * \brief
 * output operator
 *
 * writes \c arc to \c os
 */
template < class CurvedKernelViaAnalysis_2, class Rep_>
inline
std::ostream& operator<<(
    std::ostream& os,
    const Arc_2<CurvedKernelViaAnalysis_2, Rep_>& arc) {

  arc.write(os);
  return os;
}


//! \brief Reads the objects from stream.
template < class CurvedKernelViaAnalysis_2, class Rep_ >
std::istream& operator>> (
    std::istream& is,
    Arc_2< CurvedKernelViaAnalysis_2, Rep_ >& arc) {

  CGAL_precondition(CGAL::is_ascii(is));

  //typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;
  //typedef Rep_ Rep;

  arc.read(is);

  return is;
}



} // namespace internal

} //namespace CGAL

#endif // CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_ARC_2_H
// EOF
