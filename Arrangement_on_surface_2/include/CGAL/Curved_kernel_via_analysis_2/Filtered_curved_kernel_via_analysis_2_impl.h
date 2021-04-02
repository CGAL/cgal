// Copyright (c) 2008,2009,2010,2011 Max-Planck-Institute Saarbruecken (Germany),
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
//                 Pavel Emeliyanenko <asm@mpi-inf.mpg.de>
//                 Michael Kerber <mkerber@mpi-inf.mpg.de>

#ifndef CGAL_FILTERED_CURVED_KERNEL_VIA_ANALYSIS_2_IMPL_H
#define CGAL_FILTERED_CURVED_KERNEL_VIA_ANALYSIS_2_IMPL_H

/*!\file include/CGAL/Filtered_curved_kernel_via_analysis_2.h
 * \brief defines class \c Filtered_curved_kernel_via_analysis_2
 *
 * Defines points and arcs supported by curves that can be analyzed
 * and where some operations are filtered.
 */

#include <CGAL/config.h>
#include <CGAL/Curved_kernel_via_analysis_2/Curved_kernel_via_analysis_2_impl.h>
#include <CGAL/Curved_kernel_via_analysis_2/Curved_kernel_via_analysis_2_functors.h>

#include <CGAL/Bbox_2.h>

namespace CGAL {

#ifndef CKvA_CERR
//#define FCKvA_DEBUG_PRINT_CERR
#ifdef FCKvA_DEBUG_PRINT_CERR
#define CKvA_CERR(x) std::cout << x
#else
#define CKvA_CERR(x) static_cast<void>(0)
#endif
#endif

namespace internal {

namespace Filtered_curved_kernel_via_analysis_2_Functors {

#define CGAL_FILTERED_CKvA_2_GRAB_BASE_FUNCTOR_TYPES \
    typedef typename Curved_kernel_via_analysis_2::Curve_2 Curve_2; \
    typedef typename Curved_kernel_via_analysis_2::Point_2 Point_2; \
    typedef typename Curved_kernel_via_analysis_2::Arc_2 Arc_2; \


template < class CurvedKernelViaAnalysis_2, class FunctorBase >
class Compare_xy_2 :
        public FunctorBase::Compare_xy_2 {

public:
    //! this instance template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    CGAL_FILTERED_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;

#if DOXYGEN_RUNNING
    typedef typename Curved_kernel_via_analysis_2::Curve_2 Curve_2;

    //! type of point
    typedef typename Curved_kernel_via_analysis_2::Point_2 Point_2;

    //! type of arc
    typedef typename Curved_kernel_via_analysis_2::Arc_2 Arc_2;
#endif

    //! the base type
    typedef typename FunctorBase::Compare_xy_2 Base;

    //! the result type
    typedef CGAL::Comparison_result result_type;

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

        CKvA_CERR("\nfilteredcompare_xy_; p1: " << p1 << "; p2: " <<
             p2 << "\n");

        return Base::operator()(p1, p2, equal_x);
    }
};


// TODO implement Compare_y_limit_on_boundary_2


template < class CurvedKernelViaAnalysis_2, class FunctorBase >
class Compare_y_near_boundary_2 :
        public FunctorBase::Compare_y_near_boundary_2 {

public:
    //! this instance template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    CGAL_FILTERED_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;

#if DOXYGEN_RUNNING
    //! type of curve
    typedef typename Curved_kernel_via_analysis_2::Curve_2 Curve_2;

    //! type of point
    typedef typename Curved_kernel_via_analysis_2::Point_2 Point_2;

    //! type of arc
    typedef typename Curved_kernel_via_analysis_2::Arc_2 Arc_2;
#endif

    //! the bae type
    typedef typename FunctorBase::Compare_y_near_boundary_2
    Base;

    //! the result type
    typedef CGAL::Comparison_result result_type;

    //! standard constructor
    Compare_y_near_boundary_2(Curved_kernel_via_analysis_2 *kernel) :
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

        CKvA_CERR("\nfilteredcompare_y_near_boundary; cv1: " << cv1
                  << "; cv2: " << cv2 << "; end: " << ce << "\n");

        CGAL_assertion_code (
                CGAL::Arr_parameter_space loc1 = cv1.location(ce);
        )
        CGAL_precondition(Arc_2::is_on_left_right(loc1));
        CGAL_precondition(loc1 == cv2.location(ce));
        // comparing ids is the same as calling is_identical() ??
        if (cv1.id() == cv2.id()) {
            return CGAL::EQUAL;
        }

        typedef typename Arc_2::Curve_kernel_2 Curve_kernel_2;

        typedef typename Curve_kernel_2::Coordinate_1 Coordinate_1;

        CGAL::Object obj1, obj2;
        Coordinate_1 asym_info1, asym_info2;
        CGAL::Arr_parameter_space ps1, ps2;

        obj1 =
          cv1.curve().asymptotic_value_of_arc(cv1.location(ce), cv1.arcno());
        obj2 =
          cv2.curve().asymptotic_value_of_arc(cv2.location(ce), cv2.arcno());

        CGAL::Comparison_result filter_res = CGAL::EQUAL;

        if (CGAL::assign(ps1, obj1)) {
            if (CGAL::assign(ps2, obj2)) {
                if (ps1 == ps2) {
                    filter_res = CGAL::EQUAL;
                } else {
                    filter_res = (ps2 == CGAL::ARR_TOP_BOUNDARY ?
                                  CGAL::SMALLER : CGAL::LARGER);
                }
            } else {
                CGAL_assertion(CGAL::assign(asym_info2, obj2));
                filter_res = (ps1 == CGAL::ARR_TOP_BOUNDARY ?
                              CGAL::LARGER : CGAL::SMALLER);
            }
        } else {
            CGAL_assertion_code(bool check = )
                CGAL::assign(asym_info1, obj1);
            CGAL_assertion(check);
            if (CGAL::assign(ps2, obj2)) {
                filter_res = (ps2 == CGAL::ARR_TOP_BOUNDARY ?
                              CGAL::SMALLER : CGAL::LARGER);
            } else {
                CGAL_assertion_code(bool check = )
                    CGAL::assign(asym_info2, obj2);
                CGAL_assertion(check);
                filter_res = Base::_ckva()->kernel().compare_1_object()(
                        asym_info1, asym_info2
                );
            }
        }

        if (filter_res != CGAL::EQUAL) {
            CGAL_assertion_code(
            {
                Base base_compare_y_near_boundary(this->_ckva());

                CGAL::Comparison_result check_res =
                    base_compare_y_near_boundary(cv1, cv2, ce);
                CGAL_assertion(check_res == filter_res);
            }
            );
            return filter_res;
        }

        return Base::operator()(cv1, cv2, ce);
    }
};

template < class CurvedKernelViaAnalysis_2, class FunctorBase >
class Compare_y_at_x_2 :
        public FunctorBase::Compare_y_at_x_2 {

public:
    //! this instance template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    CGAL_FILTERED_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;

#if DOXYGEN_RUNNING
    //! type of curve
    typedef typename Curved_kernel_via_analysis_2::Curve_2 Curve_2;

    //! type of point
    typedef typename Curved_kernel_via_analysis_2::Point_2 Point_2;

    //! type of arc
    typedef typename Curved_kernel_via_analysis_2::Arc_2 Arc_2;
#endif

    //! the bae type
    typedef typename FunctorBase::Compare_y_at_x_2 Base;

    //! the result type
    typedef CGAL::Comparison_result result_type;

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

        return Base::operator()(p, cv);
    }
};

template < class CurvedKernelViaAnalysis_2, class FunctorBase >
class Compare_y_at_x_left_2 :
        public FunctorBase::Compare_y_at_x_left_2 {

public:
    //! this instance template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    CGAL_FILTERED_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;

#if DOXYGEN_RUNNING
    //! type of curve
    typedef typename Curved_kernel_via_analysis_2::Curve_2 Curve_2;

    //! type of point
    typedef typename Curved_kernel_via_analysis_2::Point_2 Point_2;

    //! type of arc
    typedef typename Curved_kernel_via_analysis_2::Arc_2 Arc_2;
#endif

    //! the bae type
    typedef typename FunctorBase::Compare_y_at_x_left_2 Base;

    //! the result type
    typedef CGAL::Comparison_result result_type;

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

        return Base::operator()(cv1, cv2, p);
    }
};

template < class CurvedKernelViaAnalysis_2, class FunctorBase >
class Compare_y_at_x_right_2 :
        public FunctorBase::Compare_y_at_x_right_2 {

public:
    //! this instance template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    CGAL_FILTERED_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;

#if DOXYGEN_RUNNING
    //! type of curve
    typedef typename Curved_kernel_via_analysis_2::Curve_2 Curve_2;

    //! type of point
    typedef typename Curved_kernel_via_analysis_2::Point_2 Point_2;

    //! type of arc
    typedef typename Curved_kernel_via_analysis_2::Arc_2 Arc_2;
#endif

    //! the bae type
    typedef typename FunctorBase::Compare_y_at_x_right_2 Base;

    //! the result type
    typedef CGAL::Comparison_result result_type;

    //! standard constructor
    Compare_y_at_x_right_2(Curved_kernel_via_analysis_2 *kernel) :
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

        CKvA_CERR("\ncompare_y_at_x_right(cv2); cv1: " << cv1 << "; cv2: " <<
            cv2 << "; p: " << p << "\n");

        return Base::operator()(cv1, cv2, p);
    }
};

template < class CurvedKernelViaAnalysis_2, class FunctorBase >
class May_have_intersection_2 :
        public Curved_kernel_via_analysis_2_Functors::
    Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {

public:
    //! this instance template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    CGAL_FILTERED_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;

    typedef Curved_kernel_via_analysis_2_Functors::
        Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 >
    Base;

#if DOXYGEN_RUNNING
    //! type of curve
    typedef typename Curved_kernel_via_analysis_2::Curve_2 Curve_2;

    //! type of point
    typedef typename Curved_kernel_via_analysis_2::Point_2 Point_2;

    //! type of arc
    typedef typename Curved_kernel_via_analysis_2::Arc_2 Arc_2;
#endif

private:

    typedef typename Curved_kernel_via_analysis_2::Curve_kernel_2
    Curve_kernel_2;

    typedef typename Curve_kernel_2::Curve_analysis_2 Curve_analysis_2;
    typedef typename Curve_kernel_2::Coordinate_1 Coordinate_1;
    typedef typename Curve_kernel_2::Coordinate_2 Coordinate_2;
    typedef typename Curve_kernel_2::X_real_traits_1 X_real_traits_1;
    typedef typename Curve_kernel_2::Y_real_traits_1 Y_real_traits_1;

    typename X_real_traits_1::Lower_boundary x_low;
    typename X_real_traits_1::Upper_boundary x_high;
    typename X_real_traits_1::Refine x_refine;

    typename Y_real_traits_1::Lower_boundary y_low;
    typename Y_real_traits_1::Upper_boundary y_high;
    typename Y_real_traits_1::Refine y_refine;

    typedef typename Coordinate_1::Rational Boundary;

public:
    typedef bool result_type;

    //! standard constructor
    May_have_intersection_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }

    /*!\brief
     * Checks whether \c cv1 and \c cv2 can have an intersection. If
     * not it certainly returns false, if possible, it return true.
     */
    bool operator()(const Arc_2& cv1, const Arc_2& cv2) const {

        Arc_2 trimmed_cv1, trimmed_cv2;

        if(! cv1.is_vertical() && ! cv2.is_vertical() ) {

            if(! cv1.trim_by_arc(cv2,trimmed_cv1,trimmed_cv2)) {
                return false;
            }

        } else {
            trimmed_cv1 = cv1;
            trimmed_cv2 = cv2;
        }

        std::list< CGAL::Bbox_2 > boxes1, boxes2;

        construct_covering_approximation(trimmed_cv1,
                                         std::back_inserter(boxes1));

        construct_covering_approximation(trimmed_cv2,
                                         std::back_inserter(boxes2));

        if (!boxes1.empty() && !boxes2.empty()) {
            // TODO better strategy than quadratic pair of for-loops (MK)
            for (typename std::list< CGAL::Bbox_2 >::const_iterator bit1 =
                     boxes1.begin(); bit1 != boxes1.end(); bit1++) {
                for (typename std::list< CGAL::Bbox_2 >::const_iterator bit2 =
                         boxes2.begin(); bit2 != boxes2.end(); bit2++) {
                    if (CGAL::do_overlap(*bit1, *bit2)) {
                        return true;
                    }
                }
            }
        }

        return false;
    }

public:

    /*!\brief
     * Constructs for a given \c arc its covering approximation.
     */
    template < class OutputIterator >
    OutputIterator construct_covering_approximation(
            const Arc_2& arc, OutputIterator oi) const {

        CKvA_CERR("\nconstruct_covering_approximation; arc: " << arc
             << ";\n cv:" << arc << "\n");

        // TODO compute more than a single bbox (EB)

        CGAL::Bbox_2 bbox = arc.bbox();

        CKvA_CERR("\nres: " << bbox << "\n");

        *oi++ = bbox;
        return oi;
    }

};


//! checks whether and how two arcs are intersection - with first filtering
template < class CurvedKernelViaAnalysis_2, class FunctorBase >
class Intersect_2 :
        public FunctorBase::Intersect_2 {

public:
    //! this instance template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    CGAL_FILTERED_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;

#if DOXYGEN_RUNNING
    //! type of curve
    typedef typename Curved_kernel_via_analysis_2::Curve_2 Curve_2;

    //! type of point
    typedef typename Curved_kernel_via_analysis_2::Point_2 Point_2;

    //! type of arc
    typedef typename Curved_kernel_via_analysis_2::Arc_2 Arc_2;
#endif

    //! the bae type
    typedef typename FunctorBase::Intersect_2 Base;

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

        CKvA_CERR("\nfiltered_intersect; cv1: " << cv1
             << ";\n cv2:" << cv2 << "");

        if (!Base::_ckva()->may_have_intersection_2_object()(cv1, cv2)) {
            // return no one
            CKvA_CERR("\nfilter: sucessfull\n");

            CGAL_assertion_code(
            {
                std::vector<CGAL::Object> tmp;
                Base::operator()(cv1, cv2, std::back_inserter(tmp));
                CGAL_assertion(tmp.empty());
            });
            return oi;
        }

        // else
        CKvA_CERR("\nfilter: failed\n");

        return Base::operator()(cv1, cv2, oi);
    }
};

template < class CurvedKernelViaAnalysis_2, class FunctorBase >
class Is_on_2 :
        public FunctorBase::Is_on_2 {

public:
    //! this instance template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    CGAL_FILTERED_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;

#if DOXYGEN_RUNNING
    //! type of curve
    typedef typename Curved_kernel_via_analysis_2::Curve_2 Curve_2;

    //! type of point
    typedef typename Curved_kernel_via_analysis_2::Point_2 Point_2;

    //! type of arc
    typedef typename Curved_kernel_via_analysis_2::Arc_2 Arc_2;
#endif

    //! the bae type
    typedef typename FunctorBase::Is_on_2 Base;

    //! the result type
    typedef bool result_type;

    //! standard constructor
    Is_on_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }

    /*!
     * Checks whether \c p lies on \c c
     * \param p The point to test
     * \param c The curve
     * \return (true) if the \c p lies on \c c
     */
    result_type operator()(const Point_2& p, const Curve_2& c) const {

        CKvA_CERR("\nfiltered_is_on; p: " << p << ";\n c:" << c << "");

        return Base::operator()(p, c);
    }

    result_type operator()(const Point_2& p, const Arc_2& arc) const {

        return Base::operator()(p, arc);
    }
};

#undef CGAL_FILTERED_CKvA_2_GRAB_BASE_FUNCTOR_TYPES

} // namespace Filtered_curved_kernel_via_analysis_2_Functors


template <class FCKvA, class BaseCKvA>
struct Filtered_functor_base :
      public BaseCKvA::template rebind< FCKvA >::Functor_base {

    typedef FCKvA Self;

    typedef BaseCKvA Base_ckva;

    typedef typename BaseCKvA::template rebind< Self >::Functor_base
        Functor_base;

// declares curved kernel functors, for each functor defines a member function
// returning an instance of this functor
#define CGAL_FILTERED_CKvA_2_functor_pred(Y, Z) \
    typedef internal::Filtered_curved_kernel_via_analysis_2_Functors:: \
        Y< Self, Functor_base > Y; \
    Y Z() const { return Y(&Self::instance()); }

#define CGAL_FILTERED_CKvA_2_functor_cons(Y, Z) \
    CGAL_FILTERED_CKvA_2_functor_pred(Y, Z)

    CGAL_FILTERED_CKvA_2_functor_pred(Compare_xy_2, compare_xy_2_object);

    CGAL_FILTERED_CKvA_2_functor_pred(Compare_y_near_boundary_2,
                                      compare_y_near_boundary_2_object);

    CGAL_FILTERED_CKvA_2_functor_pred(Compare_y_at_x_2,
                                      compare_y_at_x_2_object);
    CGAL_FILTERED_CKvA_2_functor_pred(Compare_y_at_x_left_2,
                                      compare_y_at_x_left_2_object);
    CGAL_FILTERED_CKvA_2_functor_pred(Compare_y_at_x_right_2,
                                      compare_y_at_x_right_2_object);

    CGAL_FILTERED_CKvA_2_functor_pred(Is_on_2, is_on_2_object);

    CGAL_FILTERED_CKvA_2_functor_pred(
            May_have_intersection_2, may_have_intersection_2_object
    );

    CGAL_FILTERED_CKvA_2_functor_cons(Intersect_2, intersect_2_object);

#undef CGAL_FILTERED_CKvA_2_functor_pred
#undef CGAL_FILTERED_CKvA_2_functor_cons

};

} // namespace internal

/*!\brief
 * Filtered curved kernel, i.e., intersection predicate is filted by first
 * computing a covering approximation. Only if these overlap for two arcs
 * the exact intersection predicate is called.
 */
template < class BaseCKvA_2 >
class Filtered_curved_kernel_via_analysis_2 :
    public internal::Curved_kernel_via_analysis_2_base<
        Filtered_curved_kernel_via_analysis_2< BaseCKvA_2 >,
        BaseCKvA_2, typename BaseCKvA_2::Curve_kernel_2,
        internal::Filtered_functor_base >
{
public:
    //! \name public typedefs
    //!@{

    //! this instance's first template argument
    typedef BaseCKvA_2 Curved_kernel_via_analysis_2;

    //! myself
    typedef Filtered_curved_kernel_via_analysis_2<
            Curved_kernel_via_analysis_2 > Self;

    //! type of curve kernel
    typedef typename
    Curved_kernel_via_analysis_2::Curve_kernel_2 Curve_kernel_2;

    //! type of curve analysis
    typedef typename Curve_kernel_2::Curve_analysis_2 Curve_analysis_2;

    //!@}

public:
    //!\name embedded types  for \c Arrangement_2 package
    //!@{

    //! type of curve_2
    typedef Curve_analysis_2 Curve_2;

    //! type of a point on generic curve
    typedef internal::Point_2< Self > Point_2;

    //! type of an arc on generic curve
    typedef internal::Arc_2< Self > Arc_2;

    //! type of weakly x-monotone arc for \c ArrangementTraits_2
    typedef Arc_2 X_monotone_curve_2;

    //!@}

protected:

    //! base kernel type
    typedef internal::Curved_kernel_via_analysis_2_base<
        Self, Curved_kernel_via_analysis_2, Curve_kernel_2,
        internal::Filtered_functor_base > Base_kernel;

public:
    //! \name Constructors
    //!@{

    //! default constructor
    Filtered_curved_kernel_via_analysis_2() :
        Base_kernel() {
    }

    //! construct using specific \c Curve_kernel_2 instance (for controlling)
    Filtered_curved_kernel_via_analysis_2(const Curve_kernel_2& kernel) :
        Base_kernel(kernel) {
    }

    //!@}

}; // class Filtered_curved_kernel_via_analysis_2

} // namespace CGAL

#endif // CGAL_FILTERED_CURVED_KERNEL_VIA_ANALYSIS_2_IMPL_H
// EOF
