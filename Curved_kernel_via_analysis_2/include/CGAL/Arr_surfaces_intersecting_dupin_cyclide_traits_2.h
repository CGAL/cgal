// Copyright (c) 2007-2008 Max-Planck-Institute Saarbruecken (Germany), 
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://eric@scm.gforge.inria.fr/svn/cgal/trunk/Kernel_23/include/CGAL/enum.h $
// $Id: enum.h 44129 2008-07-12 21:09:38Z spion $
//
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//                 Michael Kerber <mkerber@mpi-inf.mpg.de>

#ifndef CGAL_ARR_SURFACES_INTERSECTING_DUPIN_CYCLIDE_TRAITS_2
#define CGAL_ARR_SURFACES_INTERSECTING_DUPIN_CYCLIDE_TRAITS_2 1

/*!\file include/CGAL/Arr_surfaces_intersecting_dupin_cyclide_traits_2.h
 *\brief Provides traits class to compute arrangement on a Dupin cyclide
 * induced by other intersection curves with other surfaces
 */

#ifndef CGAL_ARRANGEMENT_ON_DUPIN_CYCLIDE
#define CGAL_ARRANGEMENT_ON_DUPIN_CYCLIDE 1
#endif 

#include <CGAL/config.h>

#include <boost/array.hpp>

#include <CGAL/Arr_tags.h>

#include <CGAL/Curved_kernel_via_analysis_2.h>
#include <CGAL/Curved_kernel_via_analysis_2/Curved_kernel_via_analysis_2_functors.h>

#ifdef CGAL_CKVA_COMPILE_RENDERER
#include <CGAL/Curved_kernel_via_analysis_2/Curve_renderer_facade.h>
#endif

#include <CGAL/Algebraic_kernel_d/Dupin_cyclide_3.h>

CGAL_BEGIN_NAMESPACE 

namespace CGALi {

namespace Dupin_cyclide_via_analysis_2_Functors {

#define CGAL_DC_CKvA_2_GRAB_BASE_FUNCTOR_TYPES \
    typedef typename Curved_kernel_via_analysis_2::Point_2 Point_2; \
    typedef typename Curved_kernel_via_analysis_2::Arc_2 Arc_2; \

// end define

/*! A functor that compares the x-coordinates on the 
 * vertical identification
 */
template < class CurvedKernelViaAnalysis_2 >
class Compare_x_on_identification_2 : public 
// not in CKvA_2
Curved_kernel_via_analysis_2_Functors::
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef Curved_kernel_via_analysis_2_Functors::
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;

    CGAL_DC_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef CGAL::Comparison_result result_type;
    
    //! standard constructor
    Compare_x_on_identification_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }
    
public:
    /*! Compare the x-coordinates of 2 given points that lie 
     * on the vertical identification
     * \param point1 the first point
     * \param point2 the second point
     * \return SMALLER - p1 is lexicographically smaller than p2;
     *         EQUAL   - p1 and p2 coincides;
     *         LARGER  - p1 is lexicographically larger than p2;
     * \pre point1 and point2 lie on the horizontal identification arc
     */
    result_type operator()(const Point_2& point1,
                           const Point_2& point2) const {
        CGAL_precondition(
                (point1.location() == CGAL::ARR_BOTTOM_BOUNDARY ||
                 point1.location() == CGAL::ARR_TOP_BOUNDARY)
                && 
                (point2.location() == CGAL::ARR_BOTTOM_BOUNDARY ||
                 point2.location() == CGAL::ARR_TOP_BOUNDARY)
        );

        // Remark: Is filtered
        // both ends lie on the ns-identification
        CGAL::Comparison_result res =
            point1.x().compare(point2.x());
        
        //std::cout << "compare_x_on_ident: " << res << std::endl;
        return res;
    }
};

/*! A functor that compares the y-coordinates on the 
 * horizontal identification
 */
template < class CurvedKernelViaAnalysis_2 >
class Compare_y_on_identification_2 : public 
// not in CKvA_2
Curved_kernel_via_analysis_2_Functors::
Curved_kernel_via_analysis_2_functor_base< CurvedKernelViaAnalysis_2 > {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef Curved_kernel_via_analysis_2_Functors::
    Curved_kernel_via_analysis_2_functor_base< Curved_kernel_via_analysis_2 >
    Base;

    CGAL_DC_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef CGAL::Comparison_result result_type;
    
    //! standard constructor
    Compare_y_on_identification_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }
    
public:
    /*! Compare the y-coordinates of two given points that lie 
     * on the horizontal identification
     * \param point1 the first point
     * \param point2 the second point
     * \return SMALLER - p1 is lexicographically smaller than p2;
     *         EQUAL   - p1 and p2 coincides;
     *         LARGER  - p1 is lexicographically larger than p2;
     * \pre point1 and point2 lie on the horizontal identification arc
     */
    CGAL::Comparison_result operator()(const Point_2& point1,
                                       const Point_2& point2) const {
        CGAL_precondition(
                (point1.location() == CGAL::ARR_LEFT_BOUNDARY ||
                 point1.location() == CGAL::ARR_RIGHT_BOUNDARY)
                && 
                (point2.location() == CGAL::ARR_LEFT_BOUNDARY ||
                 point2.location() == CGAL::ARR_RIGHT_BOUNDARY)
        );
        // used for comparison in top traits!!
        
        // both points lie on the we-identification, i.e., equalx
        CGAL::Object obj1 = 
            point1.curve().asymptotic_value_of_arc(point1.location(),
                                                   point1.arcno());
        CGAL::Object obj2 = 
            point2.curve().asymptotic_value_of_arc(point2.location(),
                                                   point2.arcno());
        
        typename Point_2::Curved_kernel_via_analysis_2::Curve_kernel_2::
            Algebraic_real_1 y1, y2;
        CGAL::Arr_parameter_space ps1, ps2;
        
        if (CGAL::assign(ps1, obj1)) {
            if (CGAL::assign(ps2, obj2)) {
                return CGAL::EQUAL;
            } else {
                CGAL_assertion(CGAL::assign(y2, obj2));
                return (ps1 == CGAL::ARR_BOTTOM_BOUNDARY ?
                        CGAL::SMALLER : CGAL::LARGER);
            }
        } else {
            CGAL_assertion_code(bool check = )
                CGAL::assign(y1, obj1);
            CGAL_assertion(check);
            if (CGAL::assign(ps2, obj2)) {
                return (ps2 == CGAL::ARR_TOP_BOUNDARY ?
                        CGAL::SMALLER : CGAL::LARGER);
            } else {
                CGAL_assertion_code(bool check = )
                    CGAL::assign(y2, obj2);
                CGAL_assertion(check);
                
                // Remark: Is filtered
                // TODO use AK_1
                return Point_2::Curved_kernel_via_analysis_2::instance().
                    kernel().compare_x_2_object()(y1, y2);
            }
        }
    }
};

/*! A functor that compares the equality of points and curves
 */
template < class CurvedKernelViaAnalysis_2 >
class Is_bounded_2 : 
        public CurvedKernelViaAnalysis_2::Functor_base::Is_bounded_2 {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the bae type
    typedef typename Curved_kernel_via_analysis_2::Functor_base::Is_bounded_2
         Base;
    
    CGAL_DC_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef bool result_type;
    
    //! standard constructor
    Is_bounded_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }
    
    /*! Is the end of an x-monotone curve bounded?
     * \param cv The x-monotone curve.
     * \param ce The end of xcv identifier.
     * \return true is the curve end is bounded, and false otherwise
     */
    result_type operator()(const Arc_2& cv, Arr_curve_end ce) const {
        // all points on a dupin cyclide are bounded
        return true;
    }
};



/*! A functor that compares the equality of points and curves
 */
template < class CurvedKernelViaAnalysis_2 >
class Equal_2 : 
        public CurvedKernelViaAnalysis_2::Functor_base::Equal_2 {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the bae type
    typedef typename Curved_kernel_via_analysis_2::Functor_base::Equal_2 Base;
    
    CGAL_DC_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef bool result_type;
    
public:
    //! standard constructor
    Equal_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel),
        _m_compare_x_on_identification(kernel),
        _m_compare_y_on_identification(kernel) {
    }
    
    /*!
     * Check if the two points are the same.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return (true) if the two point are the same; (false) otherwise.
     */
    result_type operator()(const Point_2& p1, const Point_2& p2) const {
        
        CGAL::Arr_parameter_space ps1 = p1.location();
        CGAL::Arr_parameter_space ps2 = p2.location();
        
        if (ps1 != ps2) {
            return false;
        }
        // else
        
        result_type res = false;
        
        switch (ps1) {
        case CGAL::ARR_LEFT_BOUNDARY:
        case CGAL::ARR_RIGHT_BOUNDARY: {
            res = (_m_compare_y_on_identification(p1, p2) == CGAL::EQUAL);
            break;
        }
        case CGAL::ARR_BOTTOM_BOUNDARY:
        case CGAL::ARR_TOP_BOUNDARY: {
            res = (_m_compare_x_on_identification(p1, p2) == CGAL::EQUAL);
            break;
        }
        case CGAL::ARR_INTERIOR:
        default: {
            // Remark: Is filtered
            Base base_equal(this->_ckva());
            
            res = base_equal(p1, p2);
            break;
        }
        }
        
        return res;
    }
    
    /*!
     * Check if the two x-monotone curves are the same 
     * (have the same graph).
     * \param cv1 The first curve
     * \param cv2 The second curve.
     * \return (true) if the two curves are the same; (false) otherwise.
     */
    result_type operator()(const Arc_2& cv1, const Arc_2& cv2) const {
        // Remark: Is filtered
        Base base_equal(this->_ckva());

        result_type res = base_equal(cv1, cv2);
        
        return res;
    }
protected:
    //! comparison instance on ns-identification
    typename Curved_kernel_via_analysis_2::Compare_x_on_identification_2
    _m_compare_x_on_identification;

    //! comparison instance on we-identification
    typename Curved_kernel_via_analysis_2::Compare_y_on_identification_2
    _m_compare_y_on_identification;
    
};

/*!\brief functor that compares vertical alignment of two curves
 * slightly to the right of an intersection
 */
template < class CurvedKernelViaAnalysis_2 >
class Compare_y_at_x_right_2 : 
        public CurvedKernelViaAnalysis_2::Functor_base::Compare_y_at_x_right_2
{

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef typename 
    Curved_kernel_via_analysis_2::Functor_base::Compare_y_at_x_right_2 Base;

    CGAL_DC_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef CGAL::Comparison_result result_type;
    
    //! standard constructor
    Compare_y_at_x_right_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }

public:    
    /*!
     * Compares the y value of two x-monotone curves immediately 
     * to the right of their intersection point. If one of the curves is
     * vertical (emanating upward from p), 
     * it's always considered to be above
     * the other curve.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param p The intersection point.
     * \pre The point p lies on both curves, and both of them must be 
     * also be defined (lexicographically) to its right.
     * \return The relative position of cv1 with respect to 
     * cv2 immdiately to the right of p: SMALLER, LARGER or EQUAL.
     */
    result_type operator()(
            const Arc_2& cv1, const Arc_2& cv2,
            const Point_2& p) const {
        
        // TODO check filtering
        if (p.location() != CGAL::ARR_INTERIOR) {
            
            if (p.location() == CGAL::ARR_LEFT_BOUNDARY ||
                p.location() == CGAL::ARR_RIGHT_BOUNDARY) {
                
                CGAL::Object obj1 = 
                    cv1.curve().asymptotic_value_of_arc(
                            CGAL::ARR_LEFT_BOUNDARY,
                            cv1.arcno()
                    );
                CGAL::Object obj2 = 
                    cv2.curve().asymptotic_value_of_arc(
                            CGAL::ARR_LEFT_BOUNDARY,
                            cv2.arcno()
                    );
                
                CGAL::Arr_parameter_space ps1, ps2;
                
                if (CGAL::assign(ps1, obj1) && CGAL::assign(ps2, obj2)) {
                    if (ps1 != ps2) {
                        CGAL::Comparison_result ret = 
                            (ps1 == CGAL::ARR_TOP_BOUNDARY ? 
                             CGAL::LARGER : CGAL::SMALLER);
                        //std::cout << "pole compare right " 
                        //          << ret << std::endl;
                        return ret;
                    }
                }
                // else
                return cv1.compare_y_near_boundary(cv2, CGAL::ARR_MIN_END);
            }
            
            if (p.location() == CGAL::ARR_BOTTOM_BOUNDARY ||
                p.location() == CGAL::ARR_TOP_BOUNDARY) {
                
                CGAL::Arr_parameter_space ps_y1 = 
                    cv1.location(CGAL::ARR_MIN_END);
                CGAL::Arr_parameter_space ps_y2 = 
                    cv2.location(CGAL::ARR_MIN_END);
                CGAL_assertion(ps_y1 != CGAL::ARR_INTERIOR);
                CGAL_assertion(ps_y2 != CGAL::ARR_INTERIOR);
                if (ps_y1 != ps_y2) {
                    if (ps_y1 == CGAL::ARR_TOP_BOUNDARY) {
                        return CGAL::SMALLER;
                    } else {
                        return CGAL::LARGER;
                    }
                }
                
                CGAL::Comparison_result res = 
                    cv1.compare_x_near_boundary(CGAL::ARR_MIN_END,
                                                cv2, CGAL::ARR_MIN_END);
                if (ps_y1 == CGAL::ARR_BOTTOM_BOUNDARY) {
                    res = -res;
                }
                return res;
            }
        }
        
        CGAL_assertion(p.location() == CGAL::ARR_INTERIOR);
        
        Base base_compare_y_at_x_right(this->_ckva());
        
        return base_compare_y_at_x_right(cv1, cv2, p);
    }
};

/*!\brief functor that compares vertical alignment of two curves
 * slightly to the left of an intersection
 */
template < class CurvedKernelViaAnalysis_2 >
class Compare_y_at_x_left_2 : 
        public CurvedKernelViaAnalysis_2::Functor_base::Compare_y_at_x_left_2 {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef typename 
    Curved_kernel_via_analysis_2::Functor_base::Compare_y_at_x_left_2 Base;

    CGAL_DC_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef CGAL::Comparison_result result_type;
    
    //! standard constructor
    Compare_y_at_x_left_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    }

public:
    /*!
     * Compares the y value of two x-monotone curves immediately 
     * to the left of their intersection point. If one of the curves is
     * vertical (emanating upward from p), 
     * it's always considered to be above
     * the other curve.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param p The intersection point.
     * \pre The point p lies on both curves, and both of them must be 
     * also be defined (lexicographically) to its left.
     * \return The relative position of cv1 with respect to 
     * cv2 immdiately to the left of p: SMALLER, LARGER or EQUAL.
     */
    result_type operator()(
            const Arc_2& cv1, const Arc_2& cv2,
            const Point_2& p) const {
    
        // TODO check filtering
        if (p.location() != CGAL::ARR_INTERIOR) {
            
            if (p.location() == CGAL::ARR_RIGHT_BOUNDARY ||
                p.location() == CGAL::ARR_LEFT_BOUNDARY) {
                
                CGAL::Object obj1 = 
                    cv1.curve().asymptotic_value_of_arc(
                            CGAL::ARR_RIGHT_BOUNDARY,
                            cv1.arcno()
                    );
                CGAL::Object obj2 = 
                    cv2.curve().asymptotic_value_of_arc(
                            CGAL::ARR_RIGHT_BOUNDARY,
                            cv2.arcno()
                    );
                
                CGAL::Arr_parameter_space ps1, ps2;
                
                if (CGAL::assign(ps1, obj1) && CGAL::assign(ps2, obj2)) {
                    if (ps1 != ps2) {
                        CGAL::Comparison_result ret = 
                            (ps1 == CGAL::ARR_TOP_BOUNDARY ? 
                             CGAL::LARGER : CGAL::SMALLER);
                        //std::cout << "pole compare right " 
                        //          << ret << std::endl;
                        return ret;
                    }
                }
                // else
                return cv1.compare_y_near_boundary(cv2, CGAL::ARR_MAX_END);
            }
            
            
            if (p.location() == CGAL::ARR_BOTTOM_BOUNDARY ||
                p.location() == CGAL::ARR_TOP_BOUNDARY) {
                
                CGAL::Arr_parameter_space ps_y1 = 
                    cv1.location(CGAL::ARR_MAX_END);
                CGAL::Arr_parameter_space ps_y2 = 
                    cv2.location(CGAL::ARR_MAX_END);
                CGAL_assertion(ps_y1 != CGAL::ARR_INTERIOR);
                CGAL_assertion(ps_y2 != CGAL::ARR_INTERIOR);
                if (ps_y1 != ps_y2) {
                    if (ps_y1 == CGAL::ARR_TOP_BOUNDARY) {
                        return CGAL::SMALLER;
                    } else {
                        return CGAL::LARGER;
                    }
                }
                
                CGAL::Comparison_result res = 
                    cv1.compare_x_near_boundary(CGAL::ARR_MAX_END,
                                                cv2, CGAL::ARR_MAX_END);
                if (ps_y1 == CGAL::ARR_TOP_BOUNDARY) {
                    res = -res;
                }
                return res;
            }
        }
        
        CGAL_assertion(p.location() == CGAL::ARR_INTERIOR);
        
        Base base_compare_y_at_x_left(this->_ckva());

        return base_compare_y_at_x_left(cv1, cv2, p);
    }
};

/*! A functor that computes the intersections of two arcs
 */
template < class CurvedKernelViaAnalysis_2 >
class Intersect_2 : 
        public CurvedKernelViaAnalysis_2::Functor_base::Intersect_2 {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the bae type
    typedef typename Curved_kernel_via_analysis_2::Functor_base::Intersect_2
         Base;
    
    CGAL_DC_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef std::iterator< output_iterator_tag, CGAL::Object > result_type;
    
public:
    //! standard constructor
    Intersect_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel),
        _m_equal(kernel) {
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
        
        bool do_overlap = cv1.do_overlap(cv2);
        
        std::list< CGAL::Object > intersections;

        if (!do_overlap) {
            CGAL::Arr_parameter_space ps1_min = 
                cv1.location(CGAL::ARR_MIN_END);
            CGAL::Arr_parameter_space ps2_min = 
                cv2.location(CGAL::ARR_MIN_END);
            
            if (ps1_min != CGAL::ARR_INTERIOR && 
                ps1_min == ps2_min) {
                if (_m_equal(cv1.curve_end(CGAL::ARR_MIN_END),
                             cv2.curve_end(CGAL::ARR_MIN_END))) {
                    std::pair< Point_2, unsigned int > 
                        pair(cv1.curve_end(CGAL::ARR_MIN_END),0);
                    *oi++ = CGAL::make_object(pair);
                }
            }
        }
        
        // Remark: Is filtered
        Base base_intersect(this->_ckva());
        base_intersect(cv1, cv2, std::back_inserter(intersections));
        
        for (std::list< CGAL::Object >::const_iterator it = 
                 intersections.begin(); it != intersections.end(); it++) {
            *oi++ = *it;
        }

        if (!do_overlap) {
            CGAL::Arr_parameter_space ps1_max = 
                cv2.location(CGAL::ARR_MAX_END);
            CGAL::Arr_parameter_space ps2_max = 
                cv2.location(CGAL::ARR_MAX_END);
            
            if (ps1_max != CGAL::ARR_INTERIOR && 
                ps1_max == ps2_max) {
                if (_m_equal(cv1.curve_end(CGAL::ARR_MAX_END),
                             cv2.curve_end(CGAL::ARR_MAX_END))) {
                    std::pair< Point_2, unsigned int > 
                        pair(cv1.curve_end(CGAL::ARR_MAX_END),0);
                    *oi++ = CGAL::make_object(pair);
                }
            }
        }

        // TODO Modify Intersect_2 to return points on identifications
        // (can be computed in terms of asymptotic values)
        
        return oi;
    }
protected:
    //! comparison instance on ns-identification
    typename Curved_kernel_via_analysis_2::Equal_2
    _m_equal;
};

template < class CurvedKernelViaAnalysis_2 >
class Make_x_monotone_2 : 
        public CurvedKernelViaAnalysis_2::Functor_base::Make_x_monotone_2 {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef typename 
    Curved_kernel_via_analysis_2::Functor_base::Make_x_monotone_2 Base;
    
    CGAL_DC_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;

    //! the result type 
    typedef std::iterator< output_iterator_tag, CGAL::Object > result_type;
    
    typedef typename Curved_kernel_via_analysis_2::Curve_2 Curve_2;

protected:
    typedef typename Curve_2::Arithmetic_kernel AK;
    typedef typename AK::Integer Integer;
    typedef CGAL::Polynomial< Integer > Polynomial_1;
    typedef CGAL::Polynomial<Polynomial_1> Polynomial_2;
    typedef CGAL::Polynomial<Polynomial_2> Polynomial_3;
    typedef CGAL::Polynomial<Polynomial_3> Polynomial_4;
    typedef typename Curved_kernel_via_analysis_2::Dupin_cyclide_3 
    Dupin_cyclide_3;
    
public:
    //! standard constructor
    Make_x_monotone_2(Curved_kernel_via_analysis_2 *kernel) :
        Base(kernel) {
    } 
    
private:
    Polynomial_3 construct_coefficient(Integer c, int i, int j, int k) {
        std::vector<Integer> coeffs(k+1);
        for (int a = 0; a < k; a++) {
            coeffs[a]=Integer(0);
        }
        coeffs[k] = c;
        Polynomial_1 p_1(coeffs.begin(),coeffs.end());
        
        std::vector< Polynomial_1 > p_1_coeffs(j+1);
        for (int a = 0; a < j; a++) {
            p_1_coeffs[a] = Polynomial_1(0);
        }
        p_1_coeffs[j] = p_1;
        Polynomial_2 p_2(p_1_coeffs.begin(),p_1_coeffs.end());
        
        std::vector<Polynomial_2> p_2_coeffs(i+1);
        for (int a = 0; a < i; a++) {
            p_2_coeffs[a] = Polynomial_2(0);
        }
        p_2_coeffs[i] = p_2;
        Polynomial_3 p_3(p_2_coeffs.begin(),p_2_coeffs.end());
        return p_3;
    }
    
    Polynomial_4 homogenize_trivariate_polynomial(Polynomial_3 p) {
        int n = CGAL::total_degree(p);
        std::vector<Polynomial_3> hat_p_coeffs(n+1);
        for (int i = 0; i <= n; i++) {
            hat_p_coeffs[i] = Polynomial_3(Polynomial_2(Polynomial_1(0)));
        }
        for(int i = 0; i <= p.degree(); i++) {
            if (!p[i].is_zero()) {
                for (int j = 0; j <= p[i].degree(); j++) {
                    if (!p[i][j].is_zero() ) {
                        for (int k = 0; k <= p[i][j].degree(); k++) {
                            if (p[i][j][k] != 0 ) {
                                CGAL_assertion(i + j + k <= n);
                                Polynomial_3 tri_coeff = 
                                    construct_coefficient(p[i][j][k],i,j,k);
                                hat_p_coeffs[n-i-j-k] += tri_coeff;
                            }
                        }
                    }
                }
            }
        }
        Polynomial_4 hat_p(hat_p_coeffs.begin(),hat_p_coeffs.end());
        //CGAL_assertion(hat_p.total_degree() == n);
        return hat_p;
        
    }
    
    Polynomial_2 substitute_xyzw(
            Polynomial_4 p, 
            Polynomial_2 x,
            Polynomial_2 y,
            Polynomial_2 z,
            Polynomial_2 w) {
        
        int i = p.degree();
        typename CGAL::Polynomial_traits_d< Polynomial_3 >::Substitute subst;
        Polynomial_2 pp[3] = {x,y,z};
        Polynomial_2 r = subst(p[i--], pp, pp+3);
        while (i >= 0) { 
            r *= w; 
            r += subst(p[i--], pp, pp+3); 
        }
        return r;
    }
    
public:
    /*!\brief
     * Decomposes the intersection curve defined by \c cv and \c base() 
     * into x-monotone subcurves 
     * and insert them to the given output iterator.
     * \param cv The surface
     * \param oi The output iterator, whose value-type is Object. 
     * The returned objects are all wrappers X_monotone_curve_2 objects.
     * \return The past-the-end iterator.
     */
    template < class OutputIterator >
    OutputIterator operator() (
            const Curve_2& cv, 
            OutputIterator oi) {
        
        Polynomial_4 hom_f = homogenize_trivariate_polynomial(cv.f());
        
#if !NDEBUG            
        std::cout << "cv.f=" << cv.f() << std::endl;
        
        std::cout << "cv.f_hom=" << hom_f  << std::endl;
#endif           
        
        Dupin_cyclide_3 base = this->_ckva()->base();
        
        Polynomial_2 x = base.x_param(),
            y = base.y_param(),
            z = base.z_param(),
            w = base.w_param();
        
        // create intersection curve in parameterization
        Polynomial_2 p = 
            substitute_xyzw(hom_f, x, y, z, w);
        
        // Make it square free
        p = CGAL::CGALi::make_square_free(p);
        // Also, make content square free:
        Polynomial_1 content = p.content();
        Polynomial_1 sf_content = CGAL::CGALi::make_square_free(content);
        if (sf_content != content) {
            p = p/content;
            p = p*sf_content;
        }          
        
        p = CGAL::canonicalize(p);
        
#if QdX_PRINTOUT_INTERSECTION_CURVES
        CGAL::set_ascii_mode(std::cout);
        std::cout << p << std::endl;
        return oi;
#endif
        
#if !NDEBUG
        CGAL::set_ascii_mode(std::cout);
        std::cout << p << std::endl;
        CGAL::set_pretty_mode(std::cout);
        std::cout  << p << std::endl;
#endif
        
        typedef typename 
            Curved_kernel_via_analysis_2::Curve_kernel_2::Curve_analysis_2 
            Curve_analysis_2;
        
        // curve cache
        Curve_analysis_2 curr_curve = 
            Curved_kernel_via_analysis_2::instance().kernel().
            construct_curve_2_object()(p);
        
        // create segments // may contain degenerate segments
        CGAL::CGALi::Make_x_monotone_2< Curved_kernel_via_analysis_2 >
            make_x_monotone(this->_ckva());
        make_x_monotone(curr_curve, oi);
        
        // TODO create objects on identification
        
        return oi;
    }
};

#undef CGAL_CKvA_2_GRAB_BASE_FUNCTOR_TYPES

} // namespace  Dupin_cyclide_via_analysis_2_Functors

template < class ASiDCTraits_2 >
class Dupin_cyclide_point_2 :
        public CGALi::Point_2< ASiDCTraits_2 > {

public:
    
    //! first template argument
    typedef ASiDCTraits_2 ASiDC_traits_2;

    //! base class
    typedef CGALi::Point_2< ASiDC_traits_2 > Base;

    //! itself
    typedef Dupin_cyclide_point_2 Self;

    //! our lovely cyclide
    typedef typename ASiDC_traits_2::Dupin_cyclide_3 Dupin_cyclide_3;

     //! types needed to replicate constructors
    typedef typename Base::Curve_analysis_2 Curve_analysis_2;
    typedef typename Base::X_coordinate_1 X_coordinate_1;
    typedef typename Base::Rep Rep;

       //! for parameterization
    typedef CGAL::Polynomial< double > Poly_double_1;
    typedef CGAL::Polynomial< Poly_double_1 > Poly_double_2;
    
    typedef CGAL::array< double, 4 > Point_double_4;
    

    //!\name replicates all constructors of the base (stupid solution)
    //! see base type constructors for detailed description 
    //!@{

    Dupin_cyclide_point_2() : 
        Base(Rep()) {   
    }

    Dupin_cyclide_point_2(const Self& p) : 
            Base(static_cast<const Base&>(p)) {  
    }

    Dupin_cyclide_point_2(const X_coordinate_1& x, const Curve_analysis_2& c,
         int arcno) :
        Base(x, c, arcno) {
    }

    //!@}
protected:
    //!\name Special constructors for points on the boundary
    //!@{
    
    Dupin_cyclide_point_2(CGAL::Arr_curve_end inf_end, 
            const Curve_analysis_2& c, int arcno) :
        Base(inf_end, c, arcno) {
    }
    
    Dupin_cyclide_point_2(const X_coordinate_1& x, const Curve_analysis_2& c, 
            CGAL::Arr_curve_end inf_end) :
        Base(Rep(x, c, inf_end)) {
    }
    
    Dupin_cyclide_point_2(Rep rep) : 
        Base(rep) {  
    }

    //!@}
public:
    //!\name visualization & approximation
    //!@{

#ifdef CGAL_CKVA_COMPILE_RENDERER
    //! sets up rendering window \c bbox and resolution 
    static void setup_renderer(CGAL::Bbox_2 bbox, int res_w, int res_h) {
        Curve_renderer_facade< ASiDC_traits_2 >::setup(bbox, res_w, res_h);
    }
    
    //! sets up cyclide parameterization equations 
    static void setup_parameterization(const Dupin_cyclide_3& base_surf) {
         
        param_surface(&base_surf);
        param_tube_circle(&base_surf);
        param_outer_circle(&base_surf);
        param_pole(&base_surf);
    }

    //! get / set parameterization polynomial with index \c idx (0..3)
    static const Poly_double_2 *param_surface(
        const Dupin_cyclide_3* base_surf = NULL) {

        static Poly_double_2 _param[4]; // x, y, z, w respectively
        if(base_surf != NULL) 
    // cannot use std::transform because of stupid layout
            _param[0] = CGAL::to_double(base_surf->x_param()),
            _param[1] = CGAL::to_double(base_surf->y_param()),
            _param[2] = CGAL::to_double(base_surf->z_param()),   
            _param[3] = CGAL::to_double(base_surf->w_param());
               
        return _param;
    }

    static const Poly_double_1 *param_tube_circle(
        const Dupin_cyclide_3* base_surf = NULL) {
        
        static Poly_double_1 _param[4]; // see above 
        if(base_surf != NULL) {
            typename Real_embeddable_traits<
                typename Dupin_cyclide_3::Polynomial_1 >::To_double
                    to_double;
            std::transform(base_surf->tube_circle(),
                 base_surf->tube_circle()+4, _param, to_double);
        }
        return _param;
    }    

    static const Poly_double_1 *param_outer_circle(
        const Dupin_cyclide_3* base_surf = NULL) {
        
        static Poly_double_1 _param[4]; // see above 
        if(base_surf != NULL) {
            typename Real_embeddable_traits<
                typename Dupin_cyclide_3::Polynomial_1 >::To_double
                    to_double;
            std::transform(base_surf->outer_circle(),
                 base_surf->outer_circle()+4, _param, to_double);
        }
        
        return _param;
    }    

    static const Point_double_4& param_pole(
        const Dupin_cyclide_3* base_surf = NULL) {

        static Point_double_4 _pole;
        if(base_surf != NULL) {
            typename Real_embeddable_traits<
                typename Dupin_cyclide_3::Point_int_4::value_type >::To_double
                    to_double;
            std::transform(base_surf->pole().begin(),
                base_surf->pole().end(), _pole.begin(), 
                to_double);
        }
        return _pole;
    } 

    /*!\brief
     * computes approximation of a point 
     *
     * returns \c false if the point does not fall within the drawing window
     */
    template < class OutputIterator >
    OutputIterator compute_approximation(OutputIterator oi) const {
    
        double c, invw;
        switch(this->location()) {
        case CGAL::ARR_BOTTOM_BOUNDARY:
        case CGAL::ARR_TOP_BOUNDARY: {

            const Poly_double_1 *outer = param_outer_circle();
            c = CGAL::to_double(this->x());
            invw = 1.0 / outer[3].evaluate(c);
            *oi++ = outer[0].evaluate(c) * invw,
            *oi++ = outer[1].evaluate(c) * invw,
            *oi++ = outer[2].evaluate(c) * invw;
            return oi;
        }
        
        case CGAL::ARR_LEFT_BOUNDARY:
        case CGAL::ARR_RIGHT_BOUNDARY: {
            CGAL::Object o = this->curve().asymptotic_value_of_arc(
                this->location(), this->arcno());
        
            if(const X_coordinate_1 *x =
                 CGAL::object_cast<X_coordinate_1>(&o)) {

                const Poly_double_1 *tube = param_tube_circle();
                c = CGAL::to_double(*x);
                //std::cerr << "horizontal asymptote: \n";// incorrect
                invw = 1.0 / tube[3].evaluate(c);
                *oi++ = tube[0].evaluate(c) * invw,
                *oi++ = tube[1].evaluate(c) * invw,
                *oi++ = tube[2].evaluate(c) * invw;
            } else {
                const Point_double_4& pole = param_pole();
                invw = 1.0 / pole[3];
                *oi++ = pole[0] * invw, *oi++ = pole[1] * invw, 
                *oi++ = pole[2] * invw;
            }
            return oi;
        }
        default:; // make compiler happy
        }
    
        typedef Curve_renderer_facade< ASiDC_traits_2 > Facade;
        
        std::pair< double, double > cc;
        if(!Facade::instance().draw(*this, cc))
            return oi; // bad luck

        const Poly_double_2* params = param_surface();
        typename CGAL::Polynomial_traits_d< Poly_double_2 >::Substitute subst;
        
        double *ptr = (double *)&cc;
        invw = 1.0 / subst(params[3], ptr, ptr+2);
        *oi++ = subst(params[0], ptr, ptr+2) * invw; 
        *oi++ = subst(params[1], ptr, ptr+2) * invw;
        *oi++ = subst(params[2], ptr, ptr+2) * invw;
        return oi;
    }
#endif // CGAL_CKVA_COMPILE_RENDERER
    //!@}

    friend class CGALi::Arc_2<ASiDC_traits_2>;
};

template < class ASiDCTraits_2 >
class Dupin_cyclide_arc_2 :
        public CGALi::Arc_2< ASiDCTraits_2 > {

public:
    
    //! first template argument
    typedef ASiDCTraits_2 ASiDC_traits_2;

    //! base class
    typedef CGALi::Arc_2< ASiDC_traits_2 > Base;

    //! itself
    typedef Dupin_cyclide_arc_2 Self;

    //! point in parametric space
    typedef typename ASiDC_traits_2::Point_2 Point_2;
    
    //! point approximation
    //typedef typename CGAL::Cartesian< double >::Point_3 Approximation_3;

    //! types needed to replicate constructors
    typedef typename Base::Curve_analysis_2 Curve_analysis_2;
    typedef typename Base::X_coordinate_1 X_coordinate_1;
    typedef typename Base::Rep Rep;

    //! for visualization
    typedef CGAL::Polynomial< double > Poly_double_1;
    typedef CGAL::Polynomial< Poly_double_1 > Poly_double_2;

    //!\name replicates all constructors of the base (stupid solution)
    //! see base type constructors for detailed description 
    //!@{

    Dupin_cyclide_arc_2() : 
        Base() {   
    }

    Dupin_cyclide_arc_2(const Self& a) : 
        Base(static_cast<const Base&>(a)) {  
    }

    //!@} 
    //!\name Constructors for non-vertical arcs
    //!@{
       
    Dupin_cyclide_arc_2(const Point_2& p, const Point_2& q, 
        const Curve_analysis_2& c, int arcno, int arcno_p, int arcno_q) : 
        Base(p, q, c, arcno, arcno_p, arcno_q) { 
    }
   
    Dupin_cyclide_arc_2(const Point_2& origin, CGAL::Arr_curve_end inf_end, 
          const Curve_analysis_2& c, int arcno, int arcno_o) :
        Base(origin, inf_end, c, arcno, arcno_o) {
    }
    
    Dupin_cyclide_arc_2(const Point_2& origin, const X_coordinate_1& asympt_x, 
          CGAL::Arr_curve_end inf_end, const Curve_analysis_2& c, int arcno,
          int arcno_o) : Base(origin, asympt_x, inf_end, c, arcno, arcno_o) {
    }

    Dupin_cyclide_arc_2(const Curve_analysis_2& c, int arcno) :
        Base(c, arcno) {
    }
    
    Dupin_cyclide_arc_2(const X_coordinate_1& asympt_x1, 
        CGAL::Arr_curve_end inf_end1, const X_coordinate_1& asympt_x2,
         CGAL::Arr_curve_end inf_end2, const Curve_analysis_2& c, int arcno) :
        Base(asympt_x1, inf_end1, asympt_x2, inf_end2, c, arcno) {
    }

    Dupin_cyclide_arc_2(CGAL::Arr_curve_end inf_endx, 
            const X_coordinate_1& asympt_x, CGAL::Arr_curve_end inf_endy, 
            const Curve_analysis_2& c, int arcno) :
        Base(inf_endx, asympt_x, inf_endy, c, arcno) {
    }
    
    //!@}
    //!\name Constructors for vertical arcs
    //!@{
    
    Dupin_cyclide_arc_2(const Point_2& p, const Point_2& q, 
            const Curve_analysis_2& c) : 
        Base(p, q, c) {          
    }
    
    Dupin_cyclide_arc_2(const Point_2& origin, CGAL::Arr_curve_end inf_end,
          const Curve_analysis_2& c) : Base(origin, inf_end, c) {
    }
    
    Dupin_cyclide_arc_2(const X_coordinate_1& x, const Curve_analysis_2& c) :
        Base(x, c) {
    }
    //!@}
protected:
    //!\name Constructor for replace endpoints + rebind
    //!@{
    
    Dupin_cyclide_arc_2(Rep rep) : 
        Base(rep) { 
    }
    
    //!@}
public:
    //!\name visualization & approximation
    //!@{

#ifdef CGAL_CKVA_COMPILE_RENDERER

#endif    
    //!@}


    // spreading out friends ;)
    friend class CGALi::Arc_2< ASiDC_traits_2 >;
};
        

template <class ASiDC_traits_2, class BaseCKvA>
struct Dupin_cyclide_functor_base :
        public BaseCKvA::template rebind< ASiDC_traits_2 >::Functor_base {

    typedef ASiDC_traits_2 Self;

    typedef typename BaseCKvA::template rebind< Self >::Functor_base
         Functor_base;
        
#if 0 
    // TODO add special construct_point_2 and construct_arc_2 functors
    // -> points/curves on identifications
    CGAL_CKvA_2_functor_cons(Construct_point_2, 
                             construct_point_2_object);

    CGAL_CKvA_2_functor_cons(Construct_arc_2, 
                             construct_arc_2_object);
#endif    

// declares curved kernel functors, for each functor defines a member function
// returning an instance of this functor
#define CGAL_DC_CKvA_2_functor_pred(Y, Z) \
    typedef Dupin_cyclide_via_analysis_2_Functors::Y< Self > Y; \
    Y Z() const { \
        return Y(&Self::instance()); \
    } 
    
#define CGAL_DC_CKvA_2_functor_cons(Y, Z) \
    CGAL_DC_CKvA_2_functor_pred(Y, Z)
    
    CGAL_DC_CKvA_2_functor_pred(Compare_x_on_identification_2, 
                                compare_x_on_identification_2_object);

    CGAL_DC_CKvA_2_functor_pred(Compare_y_on_identification_2, 
                                compare_y_on_identification_2_object);

    CGAL_DC_CKvA_2_functor_pred(Equal_2, equal_2_object);

    CGAL_DC_CKvA_2_functor_pred(Intersect_2, intersect_2_object);

    CGAL_DC_CKvA_2_functor_pred(Compare_y_at_x_right_2, 
                                compare_y_at_x_right_2_object);

    CGAL_DC_CKvA_2_functor_pred(Compare_y_at_x_left_2, 
                                compare_y_at_x_left_2_object);

    CGAL_DC_CKvA_2_functor_cons(Make_x_monotone_2, 
                                make_x_monotone_2_object);

#undef CGAL_DC_CKvA_2_functor_pred
#undef CGAL_DC_CKvA_2_functor_cons

};

} // namespace CGALi

/*!\brief
 * Derived ArrangmentTraits class for points and segments embedded on 
 * a given Dupin cyclide.
 */
template < class BaseCKvA_2 >
class Arr_surfaces_intersecting_dupin_cyclide_traits_2 :
    public CGALi::Curved_kernel_via_analysis_2_base<
        Arr_surfaces_intersecting_dupin_cyclide_traits_2< BaseCKvA_2 >,
        BaseCKvA_2, typename BaseCKvA_2::Curve_kernel_2,
        CGALi::Dupin_cyclide_functor_base >
{
public:
    
    //! this instance's template parameter
    typedef BaseCKvA_2 Curved_kernel_via_analysis_2;
    
    //! the class itself
    typedef Arr_surfaces_intersecting_dupin_cyclide_traits_2 
      < Curved_kernel_via_analysis_2 > 
    Self;
    
    //! type of curve kernel
    typedef typename 
    Curved_kernel_via_analysis_2::Curve_kernel_2 Curve_kernel_2;
    
    // TODO remove depency to Exacus
    //! type of arithmetic traits
    typedef typename CGAL::Get_arithmetic_kernel< 
      typename Curve_kernel_2::Coefficient 
    >::Arithmetic_kernel
    AK;
    
    //! type of Surface
    typedef CGAL::Dupin_cyclide_3< AK > Dupin_cyclide_3;

    typedef CGAL::Algebraic_surface_3< AK > Surface_3;

    //! typedef of Curve_2
    typedef Surface_3 Curve_2;

    //! type of point
    typedef CGALi::Dupin_cyclide_point_2< Self > Point_2;

    //! type of arc
    typedef CGALi::Dupin_cyclide_arc_2< Self > Arc_2;
    
    //! type of x-monotone curve
    typedef Arc_2 X_monotone_curve_2;

        //! Tag to tell that the boundary functors are implemented
    typedef CGAL::Arr_bounded_boundary_tag Boundary_category;


protected:

    //! base kernel type
    typedef CGALi::Curved_kernel_via_analysis_2_base<
        Self, Curved_kernel_via_analysis_2, Curve_kernel_2,
         CGALi::Dupin_cyclide_functor_base >
    Base_kernel;
    
public:
    //!\name Constructors
    //!@{

    // default constructor should be not supported!
    /*\brief 
     * Default constructor 
     */
    Arr_surfaces_intersecting_dupin_cyclide_traits_2() {
    };
    
    /*\brief 
     * Standard constructor that stores \c base as base Dupin cyclide
     */
    Arr_surfaces_intersecting_dupin_cyclide_traits_2(
            const Dupin_cyclide_3& base
    ) :
        _m_base(base) {
        
    };
    
    //@}

    //!\name Accessors
    //!@{
    
    /*!\brief
     * returns the stored base Dupin cyclide
     */
    Dupin_cyclide_3 base() const {
        return this->_m_base;
    }
    
    //!@}
    
protected:
    //! the stored base Dupin cyclide
    mutable Dupin_cyclide_3 _m_base;
    
};

CGAL_END_NAMESPACE

#endif // CGAL_ARR_SURFACES_INTERSECTING_DUPIN_CYCLIDE_TRAITS_2
// EOF
