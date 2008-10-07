// Copyright (c) 2006-2008 Max-Planck-Institute Saarbruecken (Germany), 
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
// $URL$
// $Id$
//
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//                 Michael Kerber <mkerber@mpi-inf.mpg.de>

#ifndef CGAL_ALGEBRAIC_KERNEL_D_BITSTREAM_IN_Z_FOR_XY_TRAITS_H
#define CGAL_ALGEBRAIC_KERNEL_D_BITSTREAM_IN_Z_FOR_XY_TRAITS_H 1

/*!\file include/CGAL/Algebraic_kernel_d/Bitstream_in_z_for_xy_traits.h
 * \brief Traits class for CGAL::CGALi::Bitstream_rndl_tree
 */

#warning "Replace Bitstream_in_z_for_xy_traits with Bitstream_coeffient_kernel"

#include <CGAL/config.h>

#include <boost/optional.hpp>

#include <CGAL/Algebraic_kernel_d/Real_embeddable_extension.h>

#include <CGAL/Arrangement_2l/macros.h>

CGAL_BEGIN_NAMESPACE

// predeclaration
template < class SurfaceZAtXyIsolatorTraits >
class Bitstream_in_z_for_xy_traits;

namespace CGALi { 

/*!\brief
 * Representation class for traits of bitstream tree
 */
template < class SurfaceZAtXyIsolatorTraits >
class Bitstream_in_z_for_xy_traits_rep {
public:

    // types
    //! this instance's first template parameter
    typedef SurfaceZAtXyIsolatorTraits Surface_z_at_xy_isolator_traits;

    CGAL_SURFACE_Z_AT_XY_ISOLATOR_TRAITS_SNAP_TYPEDEFS(
            Surface_z_at_xy_isolator_traits
    );
#if DOXYGEN_RUNNING
    //! type of point
    typedef typename Surface_z_at_xy_isolator_traits::Point_2 Point_2;

    //! type of bivariate polynomial
    typedef typename P_curve_2::Polynomial_2 Polynomial_2;
#endif

    //! type of x-coordinate
    typedef typename Curve_kernel_2::X_coordinate_1 X_coordinate_1;
    
    //! type of rational numbers
    typedef typename X_coordinate_1::Rational Rational;

    //! type of curve analysis
    typedef typename Curve_analysis_2::Status_line_1 Status_line_1;

    //! type for intervals
    typedef boost::numeric::interval< Rational > Interval;

    //! type of coefficient for Bitstream_descartes_rndl_tree
    typedef Polynomial_2 Coefficient;
    
    //! type of integer for Bitstream_descartes_rndl_tree
    typedef typename CGAL::Fraction_traits< Rational >::Numerator_type Integer;
    
    //!\name Constructors
    //!@{
    
    /*!\brief
     * default constructor
     */
    Bitstream_in_z_for_xy_traits_rep() {
    };

    
    /*!\brief
     * standard constructor from a given refineable point
     */
    Bitstream_in_z_for_xy_traits_rep(const Point_2& pt,
                                     bool use_artificial_x_interval) :
        _m_point(pt),
        _m_use_artificial_x_interval(use_artificial_x_interval) {
    };
    
    //!@}

private:
    // members
    //! the stored projected point
    mutable Point_2 _m_point;
    
    //! stores whether artificial x intervals is used
    mutable bool _m_use_artificial_x_interval;
    
    //! if x is rational this is its value
    mutable Rational _m_x_rational;
    
    //! current interval radius
    mutable Rational _m_eps;

    //! artificial x-interval
    mutable boost::optional< Interval > _m_x_interval;

    // friends
    friend 
    class 
    CGAL::Bitstream_in_z_for_xy_traits< Surface_z_at_xy_isolator_traits >;
    
}; // Bitstream_in_z_for_xy_traits_rep

} // namespace CGALi


/*!\brief
 * Model of BitstreamDescartesRndlTreeTraits concept to isolate the
 * real roots of a square-free polynomial defined by an algebraic
 * surface in 3d and a line parallel to the z-axis running through a
 * given point, whose coordinates are usally algebraic.
 */
template < class SurfaceZAtXyIsolatorTraits >
class Bitstream_in_z_for_xy_traits : public 
::CGAL::Handle_with_policy< CGAL::CGALi::Bitstream_in_z_for_xy_traits_rep< 
SurfaceZAtXyIsolatorTraits > > {
    
public:
    
    // types
    //! this instance's first template parameter
    typedef SurfaceZAtXyIsolatorTraits Surface_z_at_xy_isolator_traits;

    CGAL_SURFACE_Z_AT_XY_ISOLATOR_TRAITS_SNAP_TYPEDEFS(
            Surface_z_at_xy_isolator_traits
    );
#if DOXYGEN_RUNNING
    //! type of point
    typedef typename Surface_z_at_xy_isolator_traits::Point_2 Point_2;
#endif
    
    //! type of instantiated class
    typedef Bitstream_in_z_for_xy_traits< Surface_z_at_xy_isolator_traits > 
    Self;

    //! type of representation
    typedef 
    CGALi::Bitstream_in_z_for_xy_traits_rep< Surface_z_at_xy_isolator_traits >
    Rep;

    //! type of Base
    typedef CGAL::Handle_with_policy< Rep > Base;
    
private:

    //! type of x-coordinate
    typedef typename Rep::X_coordinate_1 X_coordinate_1;
    
    //! type of status line
    typedef typename Rep::Status_line_1 Status_line_1;

public:
    
    //! type for rational numbers
    typedef typename Rep::Rational Rational;

    //! type for intervals
    typedef typename Rep::Interval Interval;

    //! type for Bitstream_descartes.h
    typedef Polynomial_3 POLY;
    
public:
    // BS-types
    //!\name Types for Bitstream Traits
    //!@{
    
    //! type of coefficient for Bitstream_descartes_rndl_tree
    typedef Polynomial_2 Coefficient;
    
    //! type of integer for Bitstream_descartes_rndl_tree
    typedef typename CGAL::Fraction_traits< Rational >::Numerator_type Integer;
    
    //! type of boundary for Bitstream_descartes_rndl_tree
    typedef Rational Boundary;

    //! type for the box
    typedef std::pair<Interval, Interval> Box;
    
    //!@}

    //!\name Constructors
    //!@{
    
    /*!\ brief
     * Default constructor
     */
    Bitstream_in_z_for_xy_traits() :
        Base(Rep()) {
    }

    /*!\brief
     * Standard constructor to define xy with the help of a point
     */
    Bitstream_in_z_for_xy_traits(
            const Point_2& pt,
            bool use_artificial_x_interval = false) : 
        Base(Rep(pt, use_artificial_x_interval)) {
        CGAL_precondition(pt.arcno() >= 0);
        // checks if rational and if so, sets artifical interval
        _x_iv();
    }

    //!@}
    
    //!\name Access members
    //!@{
    // functions to collaborate with functors

    //! returns stored point
    inline 
    Point_2 point() const {
        return this->ptr()->_m_point;
    }
    
    /*!\brief
     * returns current approximation of p as interval
     *
     * Note that is_zero(p) is not public as it uses internally point_on_curve
     * that uses this class itself -> cyclic dependency
     */
    inline
    Interval approximation(const Coefficient& p) const {
        
        // initialize x_iv
        Interval x_iv(_x_iv());
        
        // initialize y_iv
        Interval y_iv(_y_iv(point()));
        
        return _evaluate_iv_2(p, x_iv, y_iv);
    }
    
private:
    
    //! returns finite x-coordinate of stored point
    inline
    static
    X_coordinate_1 _x(const Point_2& point) {
        return point.x();
    }

    //! returns arcno of stored point
    inline 
    static
    int _arcno(const Point_2& point) {
        return point.arcno();
    }

    /*!\brief
     * returns Event1_info of supporting curve of stored 
     * point at its x-coordinate
     */
    inline
    static 
    Status_line_1 _sl(const Point_2& point) {
        return point.curve().status_line_at_exact_x(_x(point));
    }
    
    //!@}
    
    //!\name Intervals and their evaluations
    //!@{

public:

    //! returns current interval approximation for x
    inline
    Interval x_interval() const {
        return _x_iv();
    }

    //! returns current interval approximation for y
    inline
    Interval y_interval() const {
        return _y_iv(point());
    }
    
private:

    //! returns current interval for x
    inline
    Interval _x_iv(bool try_to_use_exact_rational = false) const {
        bool use_axiv = this->ptr()->_m_use_artificial_x_interval;
        if (use_axiv && this->ptr()->_m_x_interval) {
            if(try_to_use_exact_rational) {
                return Interval(this->ptr()->_m_x_rational,
                                this->ptr()->_m_x_rational);
            }
            return *this->ptr()->_m_x_interval;
        } else {
            if (use_axiv && _x(point()).is_rational()) {
                Interval y_iv = _y_iv(point());
                Rational mid = _x(point()).low();
                this->ptr()->_m_x_rational = mid;
                Rational eps = 
                    (y_iv.upper() - y_iv.lower()) / Rational(2);
                this->ptr()->_m_eps = eps;
                this->ptr()->_m_x_interval  =
                    Interval(mid - eps, mid + eps);
                if(try_to_use_exact_rational) {
                    return Interval(this->ptr()->_m_x_rational,
                                    this->ptr()->_m_x_rational);
                }
                return *this->ptr()->_m_x_interval;
            }
        }
        // else 
        return Interval(_x(point()).low(), _x(point()).high());
    }
    
    //! return current interval for y
    inline
    static
    Interval _y_iv(const Point_2& point) {
        
        typename Arrangement_traits_2::Curve_kernel_2::Lower_boundary_y_2
            lower_boundary_y = 
            Arrangement_traits_2::instance().kernel().
            lower_boundary_y_2_object();
        typename Arrangement_traits_2::Curve_kernel_2::Upper_boundary_y_2
            upper_boundary_y = 
            Arrangement_traits_2::instance().kernel().
            upper_boundary_y_2_object();
        
        return Interval(lower_boundary_y(
                                _sl(point).algebraic_real_2(_arcno(point))
                        ),
                        upper_boundary_y(
                                _sl(point).algebraic_real_2(_arcno(point))
                        )
        );
    }
    
    //! interval evaluation of univariate polynomial 
    inline
    static
    Interval _evaluate_iv_1(const typename Coefficient::NT& p, 
                            const Interval& x_iv) {
        // TASK cache for iv(p)
        int n = p.degree();
        Interval ret(p[n],p[n]);
        for(int i = n - 1; i >= 0; i--) {
            ret *= x_iv;;
            ret += Interval(p[i],p[i]);
        }
        return ret;
    }
    
    //! interval evaluation of bivariate polynomial 
    inline
    static 
    Interval _evaluate_iv_2(const Coefficient& p, 
                            const Interval& x_iv, 
                            const Interval& y_iv) {
        // TASK cache for iv(p)
        int i = p.degree();
        Interval ret = _evaluate_iv_1(p[i--],x_iv);
        while (i >= 0) { 
            ret *= y_iv; 
            ret += _evaluate_iv_1(p[i--],x_iv); 
        }
        return ret;
    }
    
    //!@}
    
    //!\name Representation
    //!@{
public:    
    //! refines stored representation
    inline
    void refine() const {
        
        // Goal: keep xy-approx as quadratic as possible
        
        const Point_2& point = this->ptr()->_m_point;

        bool do_refine_y = true;

        if (this->ptr()->_m_use_artificial_x_interval || !is_x_rational()) {
            
            Interval x_iv = _x_iv();
            Interval y_iv = _y_iv(point);
            
            Rational x_len = x_iv.upper() - x_iv.lower();
            Rational y_len = y_iv.upper() - y_iv.lower();
            
            CGAL_assertion(x_len > 0);
            CGAL_assertion(y_len > 0);
            
            if (x_len > y_len) {
                do_refine_y = false;
            }
        }
        
        // Remark: it is not good to make the x- and y- length similar
        //         long, as this requires an initial lengthy refinement of
        //         y since x is usally very small.
        //         The point is different: One refinement step of x
        //         implies a large number of refinement steps for y
        //         until |x| ~ |y|.
        
        if (do_refine_y) {
            refine_y();
        } else {
            // refine x
            refine_x();
        }
    }

    //! refines x of stored representation
    inline
    void refine_x() const {
        if (this->ptr()->_m_use_artificial_x_interval && 
            this->ptr()->_m_x_interval) {
            this->ptr()->_m_eps /= Rational(2);
            const Rational& mid = this->ptr()->_m_x_rational;
            const Rational& eps = this->ptr()->_m_eps;
            this->ptr()->_m_x_interval = Interval(mid - eps, mid + eps);
        } else {
            _x(point()).refine();
        }
    }

    //! refines y of stored representation
    inline
    void refine_y() const {
        int arc = _arcno(point());
        
        typename Arrangement_traits_2::Curve_kernel_2::Refine_y_2
            refine_y = 
            Arrangement_traits_2::instance().kernel().
            refine_y_2_object();
        
        refine_y(_sl(point()).algebraic_real_2(arc));
    }
    
    //! refines approximation until \c pt does not lie in approximation
    //! \pre pt != stored point
    inline
    void refine_against(const Self& traits) {
        CGAL_precondition(traits.point() != this->point());
        
        Interval this_x_iv = x_interval();
        Interval this_y_iv = y_interval();
        
        Interval that_x_iv = traits.x_interval();
        Interval that_y_iv = traits.y_interval();
        
        while (boost::numeric::overlap(this_x_iv, that_x_iv) && 
               boost::numeric::overlap(this_y_iv, that_y_iv)) {
            this->refine();
            traits.refine();
        }
    }

    Box approximation_square(long prec) {
        
        typename CGAL::Fraction_traits<Rational>::Compose compose;

        Rational bound = (prec < 0) 
            ? compose(CGAL::ipower(Integer(2),-prec),1)
            : compose(1,CGAL::ipower(Integer(2),prec));

        Rational bound_times_2 = (prec < 0) 
            ? compose(CGAL::ipower(Integer(2),-(prec+1)),1)
            : compose(1,CGAL::ipower(Integer(2),prec+1));

        

        // Refine until the approximation is inside the box
        while(y_interval().upper()-y_interval().lower() > bound_times_2 ||
              x_interval().upper()-x_interval().lower() > bound_times_2) {
            refine();
        }
        
        // now, scale to the correct size
        Interval x_iv = _rescale_to(x_interval(),prec);
        CGAL_assertion(x_iv.upper()-x_iv.lower() == bound);
        Interval y_iv = _rescale_to(y_interval(),prec);
        CGAL_assertion(y_iv.upper()-y_iv.lower() == bound);

        return std::make_pair(x_iv,y_iv);

    }

private:

    Interval _rescale_to(Interval iv, long prec) {
        
        // TODO Implement for non-integers

        typedef CGAL::Fraction_traits<Rational> Fraction_traits;
        typename Fraction_traits::Decompose decompose;
        Integer num, denom;
        Integer pow = CGAL::ipower(Integer(2),CGAL::abs(prec)+1);
        Integer low_result, high_result;
        Integer remainder;

        decompose(iv.lower(), num, denom);

        //std::cout << "Lower: " << num << ", " << denom << std::endl;

        if(prec<0) {
            denom = denom * pow;
        } else {
            num = num * pow;
        }
        CGAL::div_mod(num, denom, low_result,remainder);

        if(remainder<0) {
            low_result = low_result - 1;
        }

        decompose(iv.upper(), num, denom);

        //std::cout << "Upper: " << num << ", " << denom << std::endl;

        if(prec<0) {
            denom = denom * pow;
        } else {
            num = num * pow;
        }
        
        CGAL::div_mod(num, denom,high_result,remainder);
        if(remainder>0) {
            high_result = high_result + 1;
        }

        CGAL_assertion(high_result - low_result <=2);
        if(high_result==low_result) {
            
            low_result = low_result - 1;
            high_result = high_result + 1;
            

        } else if(high_result - low_result == 1) {
            if(CGAL::mod(low_result,Integer(2)) == 0) {
                low_result = CGAL::div(low_result,Integer(2));
                high_result = CGAL::div(high_result+1,Integer(2));
            } else {
                low_result = CGAL::div(low_result-1,Integer(2));
                high_result = CGAL::div(high_result,Integer(2));
            }
            pow = pow / 2;
        }

        typename Fraction_traits::Compose compose;

        if(prec < 0) {
            return Interval(compose(pow*low_result,Integer(1)),
                            compose(pow*high_result,1));
        }
        
        return Interval(compose(low_result,pow),
                        compose(high_result,pow));
        
    }
    

    // TODO add refine_y
private:
    /*!\brief
     * tests whether coefficient is zero
     *
     * I.e., whether stored point lies on the given curve
     */
    inline 
    bool is_zero(const Coefficient& f) {
        typename Surface_z_at_xy_isolator_traits::Point_on_curve_2 
            point_on_curve
                    // TODO (point_on_curve_object())
            ;
        bool on_curve = point_on_curve(this->ptr()->_m_point, f);

        return on_curve;
    }

    //!@}

public:
    //!\name Rational neighbors
    //!@{
    
    //! returns \c true iff x is rational
    bool is_x_rational() const {
        return _x(point()).is_rational();
    }
    
    //! return the rational x, if x is rational
    Rational rational_x() const {
        CGAL_precondition(is_x_rational());
        return _x(point()).low();
    }

    // TODO make use of Curve_kernel_2::Boundary_between

    //! returns a rational to the left of x within the current x-range approx
    Rational left_x() const {
        CGAL_precondition(!is_x_rational() || 
                          this->ptr()->_m_use_artificial_x_interval);
        if (is_x_rational()) {
            CGAL_assertion(this->ptr()->_m_use_artificial_x_interval);
            this->_x_iv();
            return 
                (this->ptr()->_m_x_interval->lower() + 
                 this->ptr()->_m_x_rational) / Rational(2);
        }
        return _x(point()).rational_between(X_coordinate_1(_x(point()).low()));
    }

    //! returns a rational to the right of x within the current x-range approx
    Rational right_x() const {
        CGAL_precondition(!is_x_rational() || 
                          this->ptr()->_m_use_artificial_x_interval);
        if (is_x_rational()) {
            CGAL_assertion(this->ptr()->_m_use_artificial_x_interval);
            this->_x_iv();
            return 
                (this->ptr()->_m_x_rational + 
                 this->ptr()->_m_x_interval->upper()) / Rational(2);
        }
        return _x(point()).rational_between(
                X_coordinate_1(_x(point()).high())
        );
    }

    //! returns whether a given \c y is in current y-approximation
    bool is_in_x_interval_interior(Rational x) const {
        Interval x_iv = _x_iv();
        return (x_iv.lower() < x) && (x < x_iv.upper());
    }
    
    
    //! returns mean value of y-approximation
    Rational mean_y() const {
        Interval y_iv = _y_iv(point());
        return (y_iv.upper() - y_iv.lower()) / Rational(2);
    }
    
    //! returns whether a given \c y is in current y-approximation
    bool is_in_y_interval_interior(Rational y) const {
        Interval y_iv = _y_iv(this->ptr()->_m_point);
        return (y_iv.lower() < y) && (y < y_iv.upper());
    }
    
    //!@}


public:

    // TASK document functors
    
    class Approximator {
    public:
        Approximator() {
        }
        
        Approximator(const Self& traits) :
            _m_traits(traits) {
        }
        
    public:
        Integer operator() (Coefficient f, long p) const {
            
            // compute interval approx of f(pt)
            Interval f_eval_iv;
            
            // compute bound
            Rational bound;
            
            if (p < 0) {
                bound = CGAL::ipower((Integer)2,-p);
            } else {
                typename CGAL::Fraction_traits<Rational>::Compose compose;
                bound = compose(1,CGAL::ipower(Integer(2),p));
            }
            
            const Point_2& pt = this->_m_traits.point();

            // initialize x_iv
            Interval x_iv(this->_m_traits._x_iv());
            
            // initialize y_iv
            Interval y_iv(this->_m_traits._y_iv(pt));
            
            // check whether bound is met
            while (true) {

                f_eval_iv = this->_m_traits._evaluate_iv_2(f, x_iv, y_iv);
                
                // TASK precision might be too high if already rational!!!
                if (CGAL::compare(f_eval_iv.upper() - f_eval_iv.lower(), bound)
                    == CGAL::SMALLER) {
                    break;
                }
                
                // if bound not met: refine
                this->_m_traits.refine();
                x_iv = this->_m_traits._x_iv();
                y_iv = this->_m_traits._y_iv(pt);
                
            }

            // until is is refined enough
            
            // now return the correct integer value

            Integer return_value;

            typename CGAL::Fraction_traits<Rational>::Decompose decompose;

            if (CGAL::sign(f_eval_iv.lower()) == CGAL::ZERO) {
                return_value = Integer(0);
            } else if (CGAL::sign(f_eval_iv.lower()) != 
                       CGAL::sign(f_eval_iv.upper())) {
                return_value = Integer(0);
            } else {
                Integer num,denom;
                decompose(f_eval_iv.upper(),num,denom);
                
                if (p < 0) {
                    denom *= CGAL::ipower(Integer(2),-p);
                } else {
                    num *= CGAL::ipower(Integer(2),p);
                }
                return_value = CGAL::div(num,denom);
            }
            
            return return_value;
        }

    private:
        // members
        mutable Self _m_traits;
    };
    
    /*!\brief
     * returns instance of Approximator 
     */
    Approximator approximator_object() const { 
        return Approximator(*this); 
    }
    
    //! befriending Approximator
    friend class Approximator;
    
    
    class Lower_bound_log2_abs {
        
    public:
        Lower_bound_log2_abs() {
        }
        
	Lower_bound_log2_abs(const Self& traits) : 
            _m_traits(traits) {
        }
        
	long operator() (Coefficient f) {
            // called for leading coefficient
            typename CGAL::Fraction_traits<Rational>::Decompose decompose;

            typename CGAL::CGALi::Real_embeddable_extension<Integer>
                ::Floor_log2_abs floor_log2_abs;
            typename CGAL::CGALi::Real_embeddable_extension<Integer>
                ::Ceil_log2_abs ceil_log2_abs;

            Interval f_eval_iv;
            Rational abs_lower,abs_upper;
            Integer lower_num,lower_denom,upper_num,upper_denom;
            
            const Point_2& pt = _m_traits.point();

            // initialize x_iv
            Interval x_iv(this->_m_traits._x_iv(true));
            
            // initialize y_iv
            Interval y_iv(this->_m_traits._y_iv(pt));
            while (true) {
                
                f_eval_iv = this->_m_traits._evaluate_iv_2(f,x_iv, y_iv);
                
                CGAL::Sign lower_sign = CGAL::sign(f_eval_iv.lower());
                
                if (CGAL::sign(f_eval_iv.upper()) == lower_sign) {
                    if(lower_sign == CGAL::POSITIVE) {
                        abs_lower = f_eval_iv.lower();
                        abs_upper = f_eval_iv.upper();
                    } else {
                        abs_lower = CGAL::abs(f_eval_iv.upper());
                        abs_upper = CGAL::abs(f_eval_iv.lower());
                    }
                    decompose(abs_lower,lower_num,lower_denom);
                    decompose(abs_upper,upper_num,upper_denom);

                    long lower_bound = 
                        floor_log2_abs(lower_num) - ceil_log2_abs(lower_denom);
                    long upper_bound = 
                        ceil_log2_abs(upper_num) - floor_log2_abs(upper_denom);

                    CGAL_assertion(upper_bound >= lower_bound);
                    if (upper_bound - lower_bound <= 3) {
                        return lower_bound;
                    }
                }

                this->_m_traits.refine();
                x_iv = this->_m_traits._x_iv(true);
                y_iv = this->_m_traits._y_iv(pt);
                
            }
	}
        
    private:
        // members
        mutable Self _m_traits;
        
    };

    Lower_bound_log2_abs lower_bound_log2_abs_object() const { 
        return Lower_bound_log2_abs(*this); 
    }

    friend class Lower_bound_log2_abs;
    
    
    class Upper_bound_log2_abs_approximator {
        
    private:
        typename CGAL::Fraction_traits<Rational>::Decompose decompose;
        
        typename CGAL::CGALi::Real_embeddable_extension<Integer>
        ::Floor_log2_abs floor_log2_abs;
        typename CGAL::CGALi::Real_embeddable_extension<Integer>
        ::Ceil_log2_abs ceil_log2_abs;
        
        std::vector<long> coefficients_for_point;
        
        int refinements_of_point;
        
        bool zero_test_enabled;
        
        int refinement_limit;
        
        // Stores id of polynomials which are known to vanish (or not to 
        // vanish) at alpha
        std::vector<long> zeroes,non_zeroes;

    public:
        Upper_bound_log2_abs_approximator() {
        }
        
        Upper_bound_log2_abs_approximator(const Self& traits) :
            refinements_of_point(0),
            zero_test_enabled(false),
            refinement_limit(0),
            _m_traits(traits) {
        }
    
        bool initial_upper_bound(
                Coefficient f, long& ub_log2_abs,bool& is_certainly_zero
        ) {
            return improve_upper_bound(f,ub_log2_abs,is_certainly_zero);
        }
        
        bool improve_upper_bound(
                const Coefficient f, long& ub_log2_abs,bool& is_certainly_zero
        ) {

            if (std::find(coefficients_for_point.begin(),
                          coefficients_for_point.end(),
                          f.id()) != coefficients_for_point.end()) {
                coefficients_for_point.clear();
                this->_m_traits.refine();
                ++refinements_of_point;
                if (refinements_of_point >= refinement_limit) {
                    zero_test_enabled = true;
                }
            }
            coefficients_for_point.push_back(f.id());
            if (zero_test_enabled) {
                if (std::find(zeroes.begin(),
                              zeroes.end(),
                              f.id()) != zeroes.end()) {
                    is_certainly_zero = true;
                    return true;
                } else if(std::find(non_zeroes.begin(),
                                    non_zeroes.end(),
                                    f.id()) != non_zeroes.end()) {
                    is_certainly_zero = false;
                } else {
                    is_certainly_zero = this->_m_traits.is_zero(f);
                    if (is_certainly_zero) {
                        zeroes.push_back(f.id());
                        return true;
                    } else {
                        non_zeroes.push_back(f.id());
                    }
                }
            } else {
                is_certainly_zero = false;
            }
            
            Interval f_eval_iv;
            
            const Point_2& pt = this->_m_traits.point();

            // initialize x_iv
            Interval x_iv(this->_m_traits._x_iv(true));
            
            // initialize y_iv
            Interval y_iv(this->_m_traits._y_iv(pt));
            
            f_eval_iv = this->_m_traits._evaluate_iv_2(f,x_iv, y_iv);
            
            CGAL::Sign lower_sign = CGAL::sign(f_eval_iv.lower());
            bool iv_contains_zero = lower_sign != 
                                    CGAL::sign(f_eval_iv.upper());
            
            Rational abs_upper = 
                (CGAL::abs(f_eval_iv.lower()) < CGAL::abs(f_eval_iv.upper())) 
                ?  CGAL::abs(f_eval_iv.upper()) : CGAL::abs(f_eval_iv.lower());
            
            Rational abs_lower(0);
            if (!iv_contains_zero) {
                abs_lower = 
                    (CGAL::abs(f_eval_iv.lower()) < 
                     CGAL::abs(f_eval_iv.upper())) 
                    ?  CGAL::abs(f_eval_iv.lower()) : 
                    CGAL::abs(f_eval_iv.upper());
            }
            
            Integer upper_num, upper_denom;
            decompose(abs_upper,upper_num,upper_denom);
            if(upper_num == Integer(0)) {
                is_certainly_zero = true;
                return true;
            }
            ub_log2_abs = 
                ceil_log2_abs(upper_num) - floor_log2_abs(upper_denom);
            

            if (!iv_contains_zero) {
                Integer lower_num,lower_denom;
                decompose(abs_lower,lower_num,lower_denom);
                
                long lb_log2_abs = 
                    floor_log2_abs(lower_num)  - ceil_log2_abs(lower_denom);
                CGAL_assertion(ub_log2_abs >= lb_log2_abs);

                return ((ub_log2_abs - lb_log2_abs) <= 3);
            } else {
                return false;
            }
        }
        
    private:
        // members
        mutable Self _m_traits;
        
    };
    
    Upper_bound_log2_abs_approximator 
    upper_bound_log2_abs_approximator_object() const {
        return Upper_bound_log2_abs_approximator(*this);
    }        
    
    class Boundary_creator {
    public:
        
        Boundary_creator() {
        }
        
        Boundary operator()(Integer x, long p) {
            Integer num = x, denom;
            if (p < 0) {
                denom = CGAL::ipower(Integer(2), -p);
            } else {
                num *= CGAL::ipower(Integer(2), p);
                denom = 1;
            }
            Rational r = typename CGAL::Fraction_traits< Rational >::Compose()
                (num, denom);
            // TODO check normalization
            CGAL::simplify(r);
            return r;
        } 
    };
    
    // trivial implementations
    //! type of Sign
    typedef typename CGAL::Real_embeddable_traits< Integer >::Sgn Sign;

    //! type of Ceil_log2_abs_Integer
    typedef typename 
        CGAL::CGALi::Real_embeddable_extension<Integer>::Ceil_log2_abs 
        Ceil_log2_abs_Integer;

    //! type of Ceil_log2_abs_long
    typedef typename 
        CGAL::CGALi::Real_embeddable_extension< long >::Ceil_log2_abs 
        Ceil_log2_abs_long;
};

CGAL_END_NAMESPACE

#endif // CGAL_ALGEBRAIC_KERNEL_D_BITSTREAM_IN_Z_FOR_XY_TRAITS_H
// EOF
