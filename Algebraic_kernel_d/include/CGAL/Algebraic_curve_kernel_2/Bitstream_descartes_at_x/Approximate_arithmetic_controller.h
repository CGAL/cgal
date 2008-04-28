// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>
//                 Ralf Schindlbeck
//
// ============================================================================


#ifndef CGAL_APPROXIMATE_ARITHMETIC_CONTROLLER
#define CGAL_APPROXIMATE_ARITHMETIC_CONTROLLER 1

#include <boost/numeric/interval.hpp>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

template<typename Poly_1_,typename AlgebraicReal>
class Approximate_arithmetic_controller {

public:

    typedef Poly_1_ Poly_1;
    typedef AlgebraicReal Algebraic_real;
 
    typedef typename Algebraic_real::Rational Boundary;

    typedef typename Poly_1::NT Coefficient;

    typedef typename boost::numeric::interval<Boundary> Boundary_interval;

    typedef typename CGAL::Coercion_traits<Coefficient, Boundary>::Type 
    Eval_result;

    typedef typename boost::numeric::interval<Eval_result> Eval_interval;

    //Default constructor
    Approximate_arithmetic_controller() {}

    // Constrcutor with Polynomial and Algebraic real
    Approximate_arithmetic_controller(Poly_1 f, Algebraic_real alpha)
        : _m_f(f), _m_alpha(alpha)
    {
#if AcX_USE_DERIVATIVE_OPTION || AcX_USE_BISECTION_OPTION
        f_x = CGAL::diff(_m_f);
        flag = false;
#endif

#if AcX_USE_BISECTION_OPTION
        flag_r = false;
#endif
		iter_refine_steps = 0;
	}

    //! Obains an interval approximation of f(alpha)
    Eval_interval interval_approximation() const {
        Boundary_interval interval(_m_alpha.low(),_m_alpha.high());
	#if AcX_USE_DERIVATIVE_OPTION || AcX_USE_BISECTION_OPTION
		//! AcX_USE_DERIVATIVE_OPTION On or AcX_USE_BISECTION_OPTION On;
		//! flag request, because to save runtime
		if(!flag)
		{
			interval_f_x = evaluate_iv(f_x, interval);
			CGAL::Sign sign_f_x = CGAL::sign(interval_f_x.lower()) * CGAL::sign(interval_f_x.upper());
			//! find the sign of the derived function _m_f on their bounds
			sign_low = f_x.sign_at(_m_alpha.low());

			if(sign_f_x > 0)
			{
				//! f_x have no horizontal tangent => interval good ;-)
				flag = true;
				//std::cout << "function f is strictly monotone" << std::endl;

			#if AcX_USE_BISECTION_OPTION
				// help-variables, needed to save the old/new state of the intervals
				low = _m_alpha.low();
				high = _m_alpha.high();
				low_f = _m_f.evaluate(_m_alpha.low());
				high_f = _m_f.evaluate(_m_alpha.high());
			#endif
			}
			else
			{
				//! f_x have a horizontal tangent => interval was overestimated
				//std::cout << "function f is not monotone" << std::endl;
				interval_f = evaluate_iv( _m_f, interval );
				
			}
		}
		if(flag)
		{
			//! computing direct, without IA => because the interval isn't
			//! overestimated
		#if !AcX_USE_BISECTION_OPTION
			if(sign_low == CGAL::POSITIVE)
			{
				//! slope of _m_f is positive => interval OK
				interval_f = Eval_interval( _m_f.evaluate(_m_alpha.low()), _m_f.evaluate(_m_alpha.high()) );
			}
			else
			{
				//! slope of _m_f is negative => switch interval-boundaries
				interval_f = Eval_interval( _m_f.evaluate(_m_alpha.high()), _m_f.evaluate(_m_alpha.low()) );
			}
		#else
			if(!flag_r)
			{
				if(sign_low == CGAL::POSITIVE)
				{
					//! slope of _m_f is positive => interval OK
					interval_f = Eval_interval( low_f, high_f );
				}
				else
				{
					//! slope of _m_f is negative => switch interval-boundaries
					interval_f = Eval_interval( high_f, low_f );
				}
			}
			else
			{
				//! request if the low boundary hasn't change after the refinement
				if(low == low_r)
				{
					high = high_r;
					high_f = _m_f.evaluate(_m_alpha.high());
					if(sign_low == CGAL::POSITIVE)
					{
						interval_f = Eval_interval( low_f , high_f );
					}
					else
					{
						interval_f = Eval_interval( high_f, low_f  );
					}
				}
				//! request if the high boundary hasn't change after the refinement
				else if(high == high_r)
				{
					low = low_r;
					low_f = _m_f.evaluate(_m_alpha.low());
					if(sign_low == CGAL::POSITIVE)
					{
						interval_f = Eval_interval( low_f , high_f );
					}
					else
					{
						interval_f = Eval_interval( high_f, low_f  );
					}
				}
				else // for quadratic refinement method
				{
					//std::cout << "Refinement method without modified option, step: " << i << std::endl;
					low = low_r;
					high = high_r;
					low_f = _m_f.evaluate(_m_alpha.low());
					high_f = _m_f.evaluate(_m_alpha.high());
					if(sign_low == CGAL::POSITIVE)
					{
						interval_f = Eval_interval( low_f , high_f );
					}
					else
					{
						interval_f = Eval_interval( high_f, low_f  );
					}
				}
			}
		#endif
		}
		return interval_f;
	#else
		//! AcX_USE_DERIVATIVE_OPTION Off or AcX_USE_BISECTION_OPTION Off;
		return evaluate_iv( _m_f, interval );
	#endif
    }

    //! Refines the algebraic real
    int refine_value() const {
        _m_alpha.refine();
	#if AcX_USE_BISECTION_OPTION
		// help-variables, needed to save the new intervals after refinement
		flag_r = true;
		low_r = _m_alpha.low();
		high_r = _m_alpha.high();
	#endif
		iter_refine_steps++;
		return iter_refine_steps;
    }


private:

    Poly_1 _m_f;
    Algebraic_real _m_alpha;
	//! helpvariables for the "derivate"  and "bisection" option
#if AcX_USE_DERIVATIVE_OPTION || AcX_USE_BISECTION_OPTION
	// variable: indicate the first dervative of polynome f
    Poly_1 f_x;
	// variable: denotes the sign of the low boundary of the Algebraic_real alpha 
    mutable CGAL::Sign sign_low;
	// variable: terms the intervals for the poylnome f and its first derivative
    mutable Eval_interval interval_f, interval_f_x;
	// helpvariable: flag for indication, if f is strictly monotone in the boundary
    mutable bool flag;
#endif
	//! helpvariables only for the "bisection" option
#if AcX_USE_BISECTION_OPTION
	// variable: denotes the old state of the maximum borders of alpha
    mutable Boundary low, high;
	// helpvariable: flag_r for indication, if one refinement step is made
    mutable bool flag_r;
	// variable: terms the new state of the maximum borders of alpha
    mutable Boundary low_r, high_r;
	// variable: denotes the state of the maximum borders of the polynome f
    mutable Boundary low_f, high_f;
#endif
	// variable, which counts the refinement steps
    mutable int iter_refine_steps;

}; //class Approximate_arithmetic_controller

} // namespace CGALi

CGAL_END_NAMESPACE

#endif
