// TODO add licence
//
// ----------------------------------------------------------------------------
//
// Library       : QdX
//
// File          : CGAL/Arr_surfaces_intersecting_dupin_cyclide_traits_2.h
// QdX_release   : $Name:  $
// Revision      : $Revision: 1.17 $
// Revision_date : $Date: 2008-02-09 16:31:13 $
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//                 Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================

/*! \file Arr_surfaces_intersecting_dupin_cyclide_traits_2.h
 * \brief Provides traits class to compute arrangement on a Dupin cyclide
 * induced by other intersection curves with other surfaces
 */

#ifndef CGAL_ARR_SURFACES_INTERSECTING_DUPIN_CYCLIDE_TRAITS_2
#define CGAL_ARR_SURFACES_INTERSECTING_DUPIN_CYCLIDE_TRAITS_2 1

#include <CGAL/basic.h>
#include <CGAL/Arr_tags.h>

#include <CGAL/Curved_kernel_via_analysis_2.h>
#include <CGAL/Curved_kernel_via_analysis_2/Curved_kernel_via_analysis_2_functors.h>

// TODO remove dependency to Exacus!
#include <QdX/basic.h>
#include <NiX/Get_arithmetic_traits.h>
#include <QdX/SfX/Dupin_cyclide_3.h>

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
    typedef Arity_tag<2>            Arity;
    
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
    typedef Arity_tag<2>            Arity;
    
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
        public CurvedKernelViaAnalysis_2::Base::Is_bounded_2 {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the bae type
    typedef typename Curved_kernel_via_analysis_2::Base::Is_bounded_2 Base;
    
    CGAL_DC_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef bool result_type;
    typedef Arity_tag<2> Arity;
    
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
        public CurvedKernelViaAnalysis_2::Base::Equal_2 {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the bae type
    typedef typename Curved_kernel_via_analysis_2::Base::Equal_2 Base;
    
    CGAL_DC_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef bool result_type;
    typedef Arity_tag<2> Arity;
    
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
        public CurvedKernelViaAnalysis_2::Base::Compare_y_at_x_right_2 {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef typename 
    Curved_kernel_via_analysis_2::Base::Compare_y_at_x_right_2 Base;

    CGAL_DC_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<3>            Arity;
    
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
        public CurvedKernelViaAnalysis_2::Base::Compare_y_at_x_left_2 {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef typename 
    Curved_kernel_via_analysis_2::Base::Compare_y_at_x_left_2 Base;

    CGAL_DC_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<3>            Arity;
    
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
        public CurvedKernelViaAnalysis_2::Base::Intersect_2 {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the bae type
    typedef typename Curved_kernel_via_analysis_2::Base::Intersect_2 Base;
    
    CGAL_DC_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;
    
    //! the result type
    typedef std::iterator< output_iterator_tag, CGAL::Object > result_type;
    typedef Arity_tag<2> Arity;
    
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
        public CurvedKernelViaAnalysis_2::Base::Make_x_monotone_2 {

public:
    //! this instance' first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! the base type
    typedef typename 
    Curved_kernel_via_analysis_2::Base::Make_x_monotone_2 Base;
    
    CGAL_DC_CKvA_2_GRAB_BASE_FUNCTOR_TYPES;

    //! the result type 
    typedef std::iterator< output_iterator_tag, CGAL::Object > result_type;
    typedef Arity_tag<2> Arity;   
    
    typedef typename Curved_kernel_via_analysis_2::Curve_2 Curve_2;

protected:
    typedef typename Curve_2::Arithmetic_traits AT;
    typedef typename AT::Integer Integer;
    typedef typename AT::Poly_int1 Polynomial_1;
    typedef typename AT::Poly_int2 Polynomial_2;
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
        Polynomial_2 r = NiX::substitute_xyz(p[i--],x,y,z);
        while (i >= 0) { r *= w; r += NiX::substitute_xyz(p[i--],x,y,z); }
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
        
        p = CGAL::CGALi::canonicalize_polynomial(p);
        
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

} // namespace CGALi

/*!\brief
 * Derived ArrangmentTraits class for points and segments embedded on 
 * a given Dupin cyclide.
 */
template < class CurvedKernelViaAnalysis_2 >
class Arr_surfaces_intersecting_dupin_cyclide_traits_2 :
    public CurvedKernelViaAnalysis_2::template 
rebind< 
Arr_surfaces_intersecting_dupin_cyclide_traits_2 <
  CurvedKernelViaAnalysis_2 > >::Other
{
public:
    
    //! this instance's template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;
    
    //! the class itself
    typedef Arr_surfaces_intersecting_dupin_cyclide_traits_2 
      < Curved_kernel_via_analysis_2 > 
    Self;
    
    //! type of curve kernel
    typedef typename 
    Curved_kernel_via_analysis_2::Curve_kernel_2 Curve_kernel_2;
    
    // TODO remove depency to Exacus
    //! type of arithmetic traits
    typedef typename NiX::Get_arithmetic_traits< 
      typename Curve_kernel_2::Coefficient 
    >::Arithmetic_traits
    AT;
    
    //! type of Surface
    typedef QdX::Dupin_cyclide_3< AT > Dupin_cyclide_3;

    typedef QdX::Algebraic_surface_3< AT > Surface_3;

    //! typedef of Curve_2
    typedef Surface_3 Curve_2;

    //! type of point
    typedef CGALi::Point_2< Self > Point_2;

    //! type of arc
    typedef CGALi::Arc_2< Self > Arc_2;
    
    //! type of x-monotone curve
    typedef Arc_2 X_monotone_curve_2;

    //! type of Base
    typedef typename Curved_kernel_via_analysis_2::template
    rebind< Self >::Other Base;
    
    //! Tag to tell that the boundary functors are implemented
    typedef CGAL::Arr_bounded_boundary_tag Boundary_category;
    
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
    
    //!\name Predicates and Constructions
    //!@{
    
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
    typedef CGAL::CGALi::Dupin_cyclide_via_analysis_2_Functors::Y< Self > Y; \
    Y Z() const { \
        return Y(&Self::instance()); \
    } \
    
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

    //!@}
    
protected:
    //! the stored base Dupin cyclide
    mutable Dupin_cyclide_3 _m_base;
    
};

CGAL_END_NAMESPACE

#endif // CGAL_ARR_SURFACES_INTERSECTING_DUPIN_CYCLIDE_TRAITS_2
// EOF
