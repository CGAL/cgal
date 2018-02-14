// Copyright (c) 2007,2008,2009,2010,2011 Max-Planck-Institute Saarbruecken (Germany), 
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>

#ifndef CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_SWEEP_CURVES_ADAPTER_2_H
#define CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_SWEEP_CURVES_ADAPTER_2_H 1

/*!\file include/CGAL/Curved_kernel_via_analysis_2/Sweep_curves_adapter_2.h
 * \brief defines class \c Sweep_curves_adapter_2
 *
 * provides a valid model of \c SoX::CurveSweepTraits_2 to use
 * \c Curved_kernel_via_analysis_2 with \c SoX::sweep_curves
 */

#include <CGAL/config.h>

#include <boost/optional.hpp>
#include <boost/none.hpp>
#include <CGAL/iterator.h>
#include <CGAL/Handle_with_policy.h>

#include <CGAL/Arr_enums.h>

#include <CGAL/Curved_kernel_via_analysis_2/Generic_point_2.h>
#include <CGAL/Curved_kernel_via_analysis_2/Generic_arc_2.h>

#ifndef SCA_CERR
//#define SCA_DEBUG_PRINT_CERR
#ifdef SCA_DEBUG_PRINT_CERR
#define SCA_CERR(x) std::cerr << x
#else
#define SCA_CERR(x) static_cast<void>(0)
#endif
#endif

namespace CGAL {

// defines a set of functors required by CurveSweepTraits_2 concept
namespace Sweep_curves_functors {

template < class SweepCurvesAdapter_2 >
class Compare_xy_2
{
    typedef typename SweepCurvesAdapter_2::Generic_point_2 Point_2;

public:
    typedef CGAL::Comparison_result result_type;
    
    //! standard constructor
    Compare_xy_2(SweepCurvesAdapter_2 *adapter) :
        _m_adapter(adapter) {
        CGAL_assertion(adapter != NULL);
    }

    result_type operator()(const Point_2& p1, const Point_2& p2) const {
        result_type res = (*this)(p1, p2, true);
        SCA_CERR("Result: " << res << "\n");
        return res;
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
    result_type operator()(const Point_2& p1, const Point_2& p2, bool) const
    {
        if(p1.is_identical(p2))
            return CGAL::EQUAL;
    
        typedef typename SweepCurvesAdapter_2::Native_arc_2 Native_arc_2;
        typename SweepCurvesAdapter_2::Native_point_2 pt;
        Native_arc_2 arc;
        CGAL::Arr_curve_end end;
        CGAL::Arr_parameter_space loc1, loc2;
        bool inverse = false;

        SCA_CERR("Compare_xy_2: p1: " << p1 << "\n p2: " << p2 << std::endl);

        if(p1.is_finite()) {
            if(p2.is_finite())
                return (_m_adapter->kernel().compare_xy_2_object()(p1.point(),
                    p2.point()));
            pt = p1.point();
            arc = p2.arc();
            end = p2.curve_end();
            loc1 = arc.location(end);
        } else {
            arc = p1.arc();
            end = p1.curve_end();
            loc1 = arc.location(end);
 
            if(!p2.is_finite()) { // both points lie at infinity
                loc2 = p2.arc().location(p2.curve_end());
                if(Native_arc_2::is_on_left_right(loc1)) {
                    if(loc1 != loc2) // cmp + and -oo in x
                        return (loc1 == CGAL::ARR_LEFT_BOUNDARY ?
                                 CGAL::SMALLER : CGAL::LARGER);
                    return (_m_adapter->kernel().
                       compare_y_curve_ends_2_object()(arc, p2.arc(), end));
                }
                // compare curve ends at +/-oo in y
                if(Native_arc_2::is_on_bottom_top(loc2)) {
                    CGAL::Comparison_result res = (_m_adapter->kernel().
                        compare_x_curve_ends_2_object()
                            (arc, end, p2.arc(), p2.curve_end()));
                    if(res == CGAL::EQUAL && end == p2.curve_end() &&
                            loc1 != loc2) {
                        return (loc1 == CGAL::ARR_BOTTOM_BOUNDARY ?
                            CGAL::SMALLER : CGAL::LARGER);        
                    }
                    return res;
                }
                return (loc2 == CGAL::ARR_LEFT_BOUNDARY ? CGAL::LARGER :
                    CGAL::SMALLER);
            }
            pt = p2.point();
            inverse = true; // inverse result since we compare p2 against p1
        }
        CGAL::Comparison_result res;
        if(Native_arc_2::is_on_left_right(loc1)) // p1 (point) against p2 (arc)
            res = (loc1 == CGAL::ARR_LEFT_BOUNDARY ? CGAL::LARGER :
                 CGAL::SMALLER);
        else {
            // point is p1, arc is p2
            // compares a finite point with a curve end at y=+/-oo:
            res = _m_adapter->kernel().kernel().compare_1_object()
                (pt.x(), arc.curve_end_x(end));
                
            if(res == CGAL::EQUAL) // in case of equality use boundary types:
                res = (loc1 == CGAL::ARR_BOTTOM_BOUNDARY ? CGAL::LARGER :
                    CGAL::SMALLER);
        }
        return (inverse ? -res : res);
    }
    
private:
    SweepCurvesAdapter_2 *_m_adapter;
};

template < class SweepCurvesAdapter_2 >
class Less_xy_2
{
    typedef typename SweepCurvesAdapter_2::Generic_point_2 Point_2;

public:
    typedef bool result_type;
    
    //! standard constructor
    Less_xy_2(SweepCurvesAdapter_2 *adapter) :
        _m_adapter(adapter) {
        CGAL_assertion(adapter != NULL);
    }

    /*!
     * returns \c true if p1 \< p2 lexicographical
     */
    result_type operator()(const Point_2& p1, const Point_2& p2) const {
        return (_m_adapter->compare_xy_2_object()(p1, p2) == CGAL::SMALLER);
    }
    
private:
    SweepCurvesAdapter_2 *_m_adapter;
};

template < class SweepCurvesAdapter_2 >
class Compare_y_at_x_2
{
    typedef typename SweepCurvesAdapter_2::Generic_point_2 Point_2;
    typedef typename SweepCurvesAdapter_2::Generic_arc_2 Arc_2;
   
public:
    typedef CGAL::Comparison_result result_type;
    
    //! standard constructor
    Compare_y_at_x_2(SweepCurvesAdapter_2 *adapter) :
        _m_adapter(adapter) {
        CGAL_assertion(adapter != NULL);
    }

    result_type operator()(const Arc_2& cv, const Point_2& p) const {
        result_type res = (*this)(cv, p, true);
        SCA_CERR("Result: " << res << "\n");
        return res;
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
    result_type operator()(const Arc_2& cv, const Point_2& p, bool) const {

        SCA_CERR("Compare_y_at_x_2: cv: " << cv << "\n point: " <<
            p << std::endl);

        typedef typename SweepCurvesAdapter_2::Native_arc_2 Native_arc_2;
        typename SweepCurvesAdapter_2::Native_point_2 pt;
        typename SweepCurvesAdapter_2::Native_point_2::Coordinate_1 x;
        if(cv.is_degenerate()) {

            if(!cv.source().is_finite()) { // degenerate arc at inf
                CGAL_precondition(!p.is_finite()); // p must also lie at inf
                return (_m_adapter->compare_xy_2_object()(p, cv.source()));
            }
            pt = cv.source().point();
            CGAL_precondition_code(
              x = (p.is_finite() ? p.point().x() :
                p.arc().curve_end_x(p.curve_end()));
            );
            // cv.source().x() must be accessible here
            CGAL_precondition(x == pt.x());
            if(p.is_finite())
                return (_m_adapter->kernel().compare_xy_2_object()
                    (p.point(), pt));
            // for infinite curve end: return inversed result
            return -(_m_adapter->kernel().compare_y_at_x_2_object()(pt,
                p.arc()));
        }
        if(p.is_finite())
            return _m_adapter->kernel().compare_y_at_x_2_object()(p.point(),
                cv.arc());

        CGAL::Arr_curve_end end = p.curve_end(), end2;
       // CGAL::Boundary_type bnd_x = p.arc().boundary_in_x(end);
        CGAL::Arr_parameter_space locp = p.arc().location(end);
        
        if(Native_arc_2::is_on_left_right(locp)) {
            CGAL_precondition(locp == cv.arc().location(end));
            // compare two curve ends at +/-oo in x
            return _m_adapter->kernel().compare_y_curve_ends_2_object()
                (p.arc(), cv.arc(), end);
        }
        // p.arc() has vertical asymptote; cases:
        // 1. cv.arc() is vertical => cv.arc().x == p.curve_end_x()
        // 2. cv.arc() has no vertical asymptote at p.curve_end_x()
        // 3. cv.arc() has vertical asymptote at p.curve_end_x()
        // cases 1. and 3. relate to comparison of two inf curve ends
        x = p.arc().curve_end_x(end);
        bool eq_min, eq_max, in_x_range = cv.arc().is_in_x_range(x, &eq_min,
            &eq_max);
        end2 = CGAL::ARR_MIN_END; // relevant cv.arc()'s end for comparison
        (void)in_x_range;
        CGAL_precondition(in_x_range);

        if(!cv.arc().is_vertical()) {
            if(eq_max && Native_arc_2::is_on_bottom_top(
                    cv.arc().location(CGAL::ARR_MAX_END))) {
                end2 = CGAL::ARR_MAX_END;
            } else if(!eq_min || !Native_arc_2::is_on_bottom_top(
                cv.arc().location(CGAL::ARR_MIN_END))) {
              // compare finite point against asymptotic or vertical curve end
               return (locp == CGAL::ARR_BOTTOM_BOUNDARY ? CGAL::SMALLER :
                    CGAL::LARGER);
            }
        } else {
            CGAL_precondition_msg(p.arc().is_vertical(),
                "p.arc is not within vertical arc x-range!!");
            // arc is vertical => infinite end coincides
            return CGAL::EQUAL; // two vertical arcs => coincide
        }
        
        CGAL_precondition_msg(end == end2, "Point is not within the arc's "
            "x-range");

        if(locp != cv.arc().location(end2))
            return (locp == CGAL::ARR_BOTTOM_BOUNDARY ?
                     CGAL::SMALLER : CGAL::LARGER);
        // compare either two asymptotic ends or one vertical arc + asymptote
        // check whether result need to be reversed
        CGAL::Comparison_result res =
             (_m_adapter->kernel().compare_x_curve_ends_2_object()
                  (p.arc(), end, cv.arc(), end2));

        if(locp == cv.arc().location(end2) &&
                !p.arc().is_vertical() && !cv.arc().is_vertical())
        if((end == CGAL::ARR_MAX_END && locp == CGAL::ARR_TOP_BOUNDARY) ||
           (end == CGAL::ARR_MIN_END && locp == CGAL::ARR_BOTTOM_BOUNDARY)) 
            res = -res;
        return res;
    }
    
private:
    SweepCurvesAdapter_2 *_m_adapter;

};

template < class SweepCurvesAdapter_2 >
class Equal_y_at_x_2
{
    typedef typename SweepCurvesAdapter_2::Generic_point_2 Point_2;
    typedef typename SweepCurvesAdapter_2::Generic_arc_2 Arc_2;
   
public:
    typedef bool result_type;
    
    //! standard constructor
    Equal_y_at_x_2(SweepCurvesAdapter_2 *adapter) :
        _m_adapter(adapter) {
        CGAL_assertion(adapter != NULL);
    }

    /*!
     * returns true if \c p lies on curve \c cv
     */
    result_type operator()(const Arc_2& cv, const Point_2& p) const
    {
        return (_m_adapter->compare_y_at_x_2_object()(cv, p) == CGAL::EQUAL);
    }
    
private:
    SweepCurvesAdapter_2 *_m_adapter;

};

template < class SweepCurvesAdapter_2 >
class Multiplicity_of_intersection_2 {

    typedef typename SweepCurvesAdapter_2::Generic_point_2 Point_2;
    typedef typename SweepCurvesAdapter_2::Generic_arc_2 Arc_2;
   
public:
    typedef int result_type;
    
    //! standard constructor
    Multiplicity_of_intersection_2(SweepCurvesAdapter_2 *) {
    }

    /*!\brief 
     * multiplicity of intersection
     * 
     * The intersection multiplicity of \c *this and \c cv2 at point \c p is
     * returned.
     *
     * \pre \c p must be an intersection point.
     * \pre both arcs are not degenerate
     */
    result_type operator()(const Arc_2& cv1, const Arc_2& cv2,
        const Point_2& p) const {

        SCA_CERR("Multiplicity_of_intersection_2: cv1: " << cv1 << "\n cv2: "
            <<  cv2 << std::endl);

        CGAL_precondition(!cv1.is_degenerate());
        CGAL_precondition(!cv2.is_degenerate());
        CGAL_precondition(p.is_finite());
        return cv1.arc().multiplicity_of_intersection(cv2.arc(), p.point());
    }
};

template < class SweepCurvesAdapter_2 >
class Compare_y_right_of_point_2
{
    typedef typename SweepCurvesAdapter_2::Generic_point_2 Point_2;
    typedef typename SweepCurvesAdapter_2::Generic_arc_2 Arc_2;
   
public:
    typedef CGAL::Comparison_result result_type;
    
    //! standard constructor
    Compare_y_right_of_point_2(SweepCurvesAdapter_2 *) {
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

        SCA_CERR("Compare_y_right_of_point_2: cv1: " << cv1 << "\n cv2: " <<
            cv2 << "\n pt: " << p << std::endl);
            
        CGAL_precondition(!cv1.is_degenerate());
        CGAL_precondition(!cv2.is_degenerate());
        CGAL_precondition(p.is_finite());
        return (cv1.arc().compare_y_at_x_right(cv2.arc(), p.point()));
    }
};


template < class SweepCurvesAdapter_2 >
class Source_2
{
    typedef typename SweepCurvesAdapter_2::Generic_point_2 Point_2;
    typedef typename SweepCurvesAdapter_2::Generic_arc_2 Arc_2;
   
public:
    typedef Point_2 result_type;
    
    //! standard constructor
    Source_2(SweepCurvesAdapter_2 *adapter) :
        _m_adapter(adapter) {
        CGAL_assertion(adapter != NULL);
    }

    /*!
     * returns a minimal end of a curve arc
     */
    result_type operator()(const Arc_2& cv) const {
        return cv.source();
    }

private:
    SweepCurvesAdapter_2 *_m_adapter;
    
};

template < class SweepCurvesAdapter_2 >
class Target_2
{
    typedef typename SweepCurvesAdapter_2::Generic_point_2 Point_2;
    typedef typename SweepCurvesAdapter_2::Generic_arc_2 Arc_2;
   
public:
    typedef Point_2 result_type;
    
    //! standard constructor
    Target_2(SweepCurvesAdapter_2 *adapter) :
        _m_adapter(adapter) {
        CGAL_assertion(adapter != NULL);
    }

    /*!
     * returns a maximal end of a curve arc
     */
    result_type operator()(const Arc_2& cv) const {
        return cv.target();
    }
    
private:
    SweepCurvesAdapter_2 *_m_adapter;
    
};

template < class SweepCurvesAdapter_2 >
class Construct_segment_2
{
    typedef typename SweepCurvesAdapter_2::Generic_point_2 Point_2;
    typedef typename SweepCurvesAdapter_2::Generic_arc_2 Arc_2;
   
public:
    typedef Arc_2 result_type;
    
    //! standard constructor
    Construct_segment_2(SweepCurvesAdapter_2 *) {
    }

    /*!
     * constructs a degenerate segment from a given point
     */
    result_type operator()(const Point_2& p) const {

//         SCA_CERR("Construct_segment_2: arc@" << aa.id() <<
//             "; point@" << p.id() << std::endl);
         return Arc_2(p);
    }
};

template < class SweepCurvesAdapter_2 >
class Is_degenerate_2
{
    typedef typename SweepCurvesAdapter_2::Generic_arc_2 Arc_2;
   
public:
    typedef bool result_type;
    
    //! standard constructor
    Is_degenerate_2(SweepCurvesAdapter_2 *) {
    }

    /*!
     * checks whether this arc represents an isolated point (i.e., degenerate)
     */
    result_type operator()(const Arc_2& cv) const {
        
        return cv.is_degenerate();
    }
};

template < class SweepCurvesAdapter_2 >
class Do_overlap_2
{
    typedef typename SweepCurvesAdapter_2::Generic_arc_2 Arc_2;
   
public:
    typedef bool result_type;
    
    //! standard constructor
    Do_overlap_2(SweepCurvesAdapter_2 *) {
    }

    /*!\brief
     * checks whether two curve arcs have infinitely many intersection points,
     * i.e., they overlap
     */
    result_type operator()(const Arc_2& cv1, const Arc_2& cv2) const {
        
        if(cv1.is_degenerate() || cv2.is_degenerate())
            return false;
        return cv1.arc().do_overlap(cv2.arc());
    }
};

template < class SweepCurvesAdapter_2 >
class New_endpoints_2 {

    typedef typename SweepCurvesAdapter_2::Generic_point_2 Point_2;
    typedef typename SweepCurvesAdapter_2::Generic_arc_2 Arc_2;
   
public:
    typedef Arc_2 result_type;
    
    //! standard constructor
    New_endpoints_2(SweepCurvesAdapter_2 *adapter) :
        _m_adapter(adapter) {
        CGAL_assertion(adapter != NULL);
    }

    /*!\brief
     * returns the input segment with new - but equal - endpoints
     */
    result_type operator()(const Arc_2& cv, const Point_2& p,
                                          const Point_2& q) const {

        SCA_CERR("New_endpoints_2: cv: " << cv << "\n p: " <<
            p << "\n q: " << q << std::endl);

        CGAL_precondition(
          _m_adapter->compare_xy_2_object()(cv.source(), p) == CGAL::EQUAL
        );
        CGAL_precondition(
            _m_adapter->compare_xy_2_object()(cv.target(), q) == CGAL::EQUAL
        );
            
        cv.new_endpoints(p, q);
        return cv;
    }
    
private:
    SweepCurvesAdapter_2 *_m_adapter;
  
};

template < class SweepCurvesAdapter_2 >
class New_endpoints_opposite_2 {

    typedef typename SweepCurvesAdapter_2::Generic_point_2 Point_2;
    typedef typename SweepCurvesAdapter_2::Generic_arc_2 Arc_2;
   
public:
    typedef Arc_2 result_type;
    
    //! standard constructor
    New_endpoints_opposite_2(SweepCurvesAdapter_2 *) {
    }

    /*!\brief
     * returns the input segment with new - but equal - endpoints
     * lexicographic order of endpoints is ensured automatically, hence no
     * special handling is required
     */
    result_type operator()(const Arc_2& cv, 
                           const Point_2& /* p */,
                           const Point_2& /* q */) const {
        SCA_CERR("\n\nWARNING!! New_endpoints_opposite_2: cv: " << cv <<
            "\n p: " << p << "\n q: " << q << std::endl);
        CGAL_error_msg("New_endpoints_opposite_2 deprecated and must not be \
             called\n");
        //CGAL_precondition(p.is_finite() && q.is_finite());
        //return cv.new_endpoints(p.point(), q.point());
        return cv;
    }
};

template < class SweepCurvesAdapter_2 >
class Intersect_2 {

    typedef typename SweepCurvesAdapter_2::Generic_arc_2 Arc_2;
   
public:
    //typedef bool result_type;
    
    //! standard constructor
    Intersect_2(SweepCurvesAdapter_2 *) {
    }

    /*!\brief
     * computes intersection points of \c *this and \c cv2, writes the result
     * to the output iterator \c oi
     */
    template <class OutputIterator>
    OutputIterator operator()(const Arc_2& cv1,  const Arc_2& cv2,
            OutputIterator oi) {

        SCA_CERR("Intersect_2: cv1: " << cv1 << "\n cv2: " <<
            cv2 << std::endl);
            
        return cv1.intersect(cv2, oi);
    }
};

template < class SweepCurvesAdapter_2 >
class Intersect_right_of_point_2 {

    typedef typename SweepCurvesAdapter_2::Generic_point_2 Point_2;
    typedef typename SweepCurvesAdapter_2::Generic_arc_2 Arc_2;
   
public:
    typedef bool result_type;
    
    //! standard constructor
    Intersect_right_of_point_2(SweepCurvesAdapter_2 *adapter) :
        _m_adapter(adapter) {
        CGAL_assertion(adapter != NULL);
    }

    /*!\brief
     * computes the next intersection of \c cv1 and \c cv2 right of \c ref
     * in lexicographical order and returns it through \c res argument
     *
     * intersect_right_of_point is not called when using sweep_curves() with 
     * intersection dictionary and without validation of internal structures 
     * (as is standard). Hence we can be lazy here for the moment
     * without losing performance.
     */
    result_type operator()(const Arc_2& cv1, const Arc_2& cv2,
        const Point_2& ref, Point_2& res) const {

        SCA_CERR("Intersect_right_of_point_2: cv1: " << cv1 << "\n cv2: " <<
            cv2 << "\n ref: " << ref << std::endl);
            
        typedef std::vector<Point_2> Point_container;
        Point_container tmp;
        cv1.intersect(cv2, back_inserter(tmp));

        for(typename Point_container::const_iterator it =  tmp.begin();
                it != tmp.end(); it++) {
            // assume points are sorted lexicographical    
            if(_m_adapter->compare_xy_2_object()(*it, ref) == CGAL::LARGER) {
                res = *it;
                SCA_CERR("intersection found: " << res << "\n");
                return true;
            }
        }
        return false;
    }
    
private:
    SweepCurvesAdapter_2 *_m_adapter;

};

template < class SweepCurvesAdapter_2 >
class Make_x_monotone_2 
{
    typedef typename SweepCurvesAdapter_2::Curve_2 Curve_2;
    typedef typename SweepCurvesAdapter_2::Generic_point_2 Point_2;
    typedef typename SweepCurvesAdapter_2::Generic_arc_2 Generic_arc_2;
   
public:
    typedef CGAL::iterator<std::output_iterator_tag, Generic_arc_2> 
      result_type;
    
    //! standard constructor
    Make_x_monotone_2(SweepCurvesAdapter_2 *adapter) :
        _m_adapter(adapter) {
        CGAL_assertion(adapter != NULL);
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
    OutputIterator operator()(const Generic_arc_2& cv,
            OutputIterator oi) const {
        *oi++ = cv;
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
    template<class OutputIterator>
    OutputIterator operator()(const Curve_2& cv, OutputIterator oi) const {

        typedef typename SweepCurvesAdapter_2::Native_arc_2 Native_arc_2;
        typedef typename SweepCurvesAdapter_2::Native_point_2 Native_point_2;
        typedef typename SweepCurvesAdapter_2::Generic_point_2 Generic_point_2;
        
        typedef std::vector<CGAL::Object> Objects;
        Objects objs;
        _m_adapter->kernel().make_x_monotone_2_object()(cv,
            std::back_inserter(objs));
        // sort out normal and degenerate arcs
        for(typename Objects::const_iterator it = objs.begin();
                it != objs.end(); it++) {
            Native_arc_2 arc;
            Native_point_2 pt;
            if(CGAL::assign(arc, *it))
                *oi++ = Generic_arc_2(arc);
            else if(CGAL::assign(pt, *it))
                *oi++ = Generic_arc_2(Generic_point_2(pt));
            else
                CGAL_error_msg("Bogus object..\n");
        }
        return oi;
    }
    
private:
    SweepCurvesAdapter_2 *_m_adapter;
};

} // Sweep_curves_functors

//! \brief a wrapper class for \c Curved_kernel_via_analysis_2
template < class CurvedKernelViaAnalysis_2 >
class Sweep_curves_adapter_2 {
      //: public Curved_kernel_via_analysis_2<CurveKernel_2> {

public:
    //! \name public typedefs
    //!@{

    //! this instance's template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! type of curve kernel
    typedef typename Curved_kernel_via_analysis_2::Curve_kernel_2
    Curve_kernel_2;
    
    //! myself
    typedef Sweep_curves_adapter_2< Curved_kernel_via_analysis_2 > Self;
    
    //!@}

public:
    //! \name Constructors
    //!@{

    //! default constructor
    Sweep_curves_adapter_2() {
    }
    
    //! construct using specific \c CKvA_2 instance (for controlling)
    Sweep_curves_adapter_2(const Curved_kernel_via_analysis_2& kernel) :
        _m_kernel(kernel) {
    }
    
    //!@}
    
    //!\name embedded types and predicates for \c CurveSweepTraits
    //!@{

    //! native Curved_kernel_via_analysis_2 objects
    typedef typename Curved_kernel_via_analysis_2::Curve_2 Curve_2;
    typedef typename Curved_kernel_via_analysis_2::Point_2 Native_point_2;
    typedef typename Curved_kernel_via_analysis_2::Arc_2 Native_arc_2;

    //! generic point (supports infinity)
    typedef internal::Generic_point_2<Self> Generic_point_2;

    //! generic arc (supports isolated points)
    typedef internal::Generic_arc_2<Self> Generic_arc_2;

    //! type of point in model of  \c CurveSweepTraits_2 
    typedef Generic_point_2 Point_2;

    //! type of segment in model of  \c CurveSweepTraits_2 
    typedef Generic_arc_2 Segment_2;

// declares functors, for each functor defines a member function
// returning an instance of this functor
#define CGAL_Sweep_curves_pred(Y, Z) \
    typedef Sweep_curves_functors::Y<Self> Y; \
    Y Z() const { return Y((Sweep_curves_adapter_2 *)this); }
#define CGAL_Sweep_curves_cons(Y, Z) CGAL_Sweep_curves_pred(Y, Z)
    
    CGAL_Sweep_curves_pred(Compare_xy_2, compare_xy_2_object)
    CGAL_Sweep_curves_pred(Less_xy_2, less_xy_2_object)
    CGAL_Sweep_curves_pred(Is_degenerate_2, is_degenerate_2_object)
    
    CGAL_Sweep_curves_pred(Do_overlap_2, do_overlap_2_object)
    CGAL_Sweep_curves_pred(Compare_y_at_x_2, compare_y_at_x_2_object)
    CGAL_Sweep_curves_pred(Equal_y_at_x_2, equal_y_at_x_2_object)
    
    CGAL_Sweep_curves_pred(Multiplicity_of_intersection_2,
            multiplicity_of_intersection_2_object)
    CGAL_Sweep_curves_pred(Compare_y_right_of_point_2,
            compare_y_right_of_point_2_object)
            
    CGAL_Sweep_curves_cons(Source_2, source_2_object)
    CGAL_Sweep_curves_cons(Target_2, target_2_object)
    CGAL_Sweep_curves_cons(Construct_segment_2, construct_segment_2_object)
            
    CGAL_Sweep_curves_cons(New_endpoints_2, new_endpoints_2_object)
    CGAL_Sweep_curves_cons(New_endpoints_opposite_2,
            new_endpoints_opposite_2_object)
            
    CGAL_Sweep_curves_cons(Intersect_2, intersect_2_object)
    CGAL_Sweep_curves_cons(Intersect_right_of_point_2,
        intersect_right_of_point_2_object)
        
    CGAL_Sweep_curves_cons(Make_x_monotone_2, make_x_monotone_2_object)

#undef CGAL_Sweep_curves_pred
#undef CGAL_Sweep_curves_cons

    const Curved_kernel_via_analysis_2& kernel() const {
        return _m_kernel;
    }

    //!@}

protected:
    //!\name private members
    //!@{

    //! reference to Curved_kernel_via_analysis_2 object
    Curved_kernel_via_analysis_2 _m_kernel;
     
    //!@}
}; // class Sweep_curves_adapter_2

} //namespace CGAL

#endif // CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_SWEEP_CURVES_ADAPTER_2
// EOF
