// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : 
//
// ============================================================================

#ifndef CGAL_ALGEBRAIC_CURVE_KERNEL_CURVE_ANALYSIS_2_H
#define CGAL_ALGEBRAIC_CURVE_KERNEL_CURVE_ANALYSIS_2_H

#include <CGAL/basic.h>
#include <CGAL/Handle_with_policy.h>

#include <SoX/basic.h>
//#include <CGAL/algorithm.h>
#include <SoX/GAPS/types.h>

CGAL_BEGIN_NAMESPACE

template < class AlgebraicCurveKernel_2, class Rep_ > 
class Curve_analysis_2;

namespace CGALi {

template < class AlgebraicCurveKernel_2 >
class Curve_analysis_2_rep {

public:
	// this first template argument
	typedef AlgebraicCurveKernel_2 Algebraic_curve_kernel_2;

	// myself
    typedef Curve_analysis_2_rep<Algebraic_curve_kernel_2> Self;

	typedef typename Algebraic_curve_kernel_2::Curve_pair_2 Curve_pair_2;

	typedef typename Curve_pair_2::Algebraic_curve_2 Curve_2; 
	
	// constructors
public:
    // default constructor ()
    Curve_analysis_2_rep()  
	{  }
    
    // standard constructor
    Curve_analysis_2_rep(const Curve_2& curve) : curve_(curve)
	{   }

    // data
	// temporarily this implementation uses underlying Curve_2 from SweepX
    mutable Curve_2 curve_;
    
    // befriending the handle
    friend class Curve_analysis_2<Algebraic_curve_kernel_2, Self>;
};
    
//! \brief The class is meant to provide tools to analyse a single curve. 
//! 
//! Analysis describes the curve’s interesting points and how they are 
//! connected. The analysis searches for events. Events only occur at a finite 
//! number of x-coordinates. Each such coordinate defines a 
//! \c CurveVerticalLine_1 of an event. These coordinates also define open
//! intervals on the x-axis. Different \c CurveVerticalLine_1 at values within 
//! one such interval only differ in the values of the \c Algebraic_real_2 
//! entries. Topological information are equal for all x-coordinate inside such
//! an open interval.
template <class AlgebraicCurveKernel_2, 
		  class Rep_ = CGALi::Curve_analysis_2_rep<AlgebraicCurveKernel_2> >
class Curve_analysis_2 : public ::CGAL::Handle_with_policy< Rep_ > 
{
public:
	//!@{
	//! \name typedefs

    //! this instance's first template parameter
	typedef AlgebraicCurveKernel_2 Algebraic_curve_kernel_2;

    //! this instance's second template parameter
	typedef Rep_ Rep;
	
	//! x-coordinate type
	typedef typename Algebraic_curve_kernel_2::X_coordinate_1 X_coordinate_1;

	//! type of a curve point
	typedef typename Algebraic_curve_kernel_2::Xy_coordinate_2 Xy_coordinate_2;

	//! type of a curve
	typedef typename Algebraic_curve_kernel_2::Curve_2 Curve_2;

	//! myself
	typedef Curve_analysis_2<Algebraic_curve_kernel_2, Rep> Self;

	//! type of a vertical line
	typedef CGALi::Curve_vertical_line_1<Self>
		Curve_vertical_line_1;

	//! the handle superclass
    typedef ::CGAL::Handle_with_policy< Rep > Base;

	//!@}
	//!\name Constructors
	//!@{

    //! \brief default constructor
    Curve_analysis_2() : Base(Rep()) 
    {  }

    /*!\brief
     * copy constructor
     */
    Curve_analysis_2(const Self& p) : Base(static_cast<const Base&>(p)) 
	{  }

	//! \brief constructs a curve analysis from a given \c Curve_2 object
	Curve_analysis_2(const Curve_2& c) : Base(Rep(c))
	{  }
           
    /*!\brief
     * constructsa curve analysis from a given represenation
     */
    Curve_analysis_2(Rep rep) : Base(rep) 
    {  }

	//!@}
	//!\name Access functions
	//!@{

	//! \brief returns the defining polynomial of the analysis
	Curve_2 get_polynomial_2() const
	{ 
		return this->ptr()->curve_;
	}

	//! \brief alias for \c get_polynomial_2()
	Curve_2 get_curve_2() const
	{ 
		return get_polynomial_2();
	}

	//! \brief returns number of vertical lines that encode an event
	int number_of_vertical_lines_with_event()
	{
		return this->ptr()->curve_.num_events();
	}

	//! \brief returns an instance of \c CurveVerticalLine_1 at the i-th
	//! event
	//!
	//! \pre 0 <= i < number_of_vertical_lines_with_event()
	Curve_vertical_line_1 vertical_line_at_event(int i)
	{
		CGAL_precondition(i >= 0&&i < number_of_vertical_lines_with_event());
// 		::CGAL::set_pretty_mode(std::cout);
// 		std::cout << "evt info, i = " << i << ": " <<
// 				this->ptr()->curve_.event_info(i) << std::endl;
		return Curve_vertical_line_1(this->ptr()->curve_.event_info(i));
	}

	//! \brief returns an instance of CurveVerticalLine_1 of the i-th 
	//! interval
	//!
 	//! \pre 0 <= i < number_of_vertical_lines_with_event()
	Curve_vertical_line_1 vertical_line_of_interval(int i)
	{
		CGAL_precondition(i >= 0&&i <= number_of_vertical_lines_with_event());
		typedef typename Curve_vertical_line_1::Event1_info Event1_info;
		int n_arcs = this->ptr()->curve_.arcs_over_interval(i);
		// # of arcs to the left and to the right is the same over an interval
		return Curve_vertical_line_1(
			Event1_info(X_coordinate_1(
			this->ptr()->curve_.boundary_value_in_interval(i)),
				n_arcs, n_arcs));
	}

	//! \brief returns vertical_line_at_event(i), if x hits i-th event, 
	//! otherwise vertical_line_of_interval(i), where i is the id of the 
	//! interval \c x lies in
	//!
	//! If \c pertub is \c CGAL::NEGATIVE (CGAL::POSITIVE) and x states an 
	//! event, then \c vertical_line_of_interval(i)
	//! (\c vertical_line_of_interval(i+1)) is returned.
	//! 
	//! \pre \c x is finite
	Curve_vertical_line_1 vertical_line_for_x(X_coordinate_1 x, 
		CGAL::Sign perturb = CGAL::ZERO)
	{
		// CGAL_precondition(x is finite ??);
		int i;
		bool is_evt;
		this->ptr()->curve_.x_to_index(x, i, is_evt);
		if(is_evt) {
			if(perturb == CGAL::ZERO)
				return vertical_line_at_event(i);
			if(perturb == CGAL::POSITIVE)
				i++;
		}
		return vertical_line_of_interval(i);
	}

	//! \brief returns an instance of CurveVerticalLine_1 at a given \c x
	//!
	//! \pre \c x is finite
	Curve_vertical_line_1 vertical_line_at_exact_x(X_coordinate_1 x)
	{
		// CGAL_precondition(x is finite ??);
		int i;
		bool is_evt;
		this->ptr()->curve_.x_to_index(x, i, is_evt);
		if(is_evt) 
			return vertical_line_at_event(i);
		// else construct Event1_info at certain x over the ith interval
		int n_arcs = this->ptr()->curve_.arcs_over_interval(i);
		return Curve_vertical_line_1(Event1_info(x, n_arcs, n_arcs));
	}
	
}; // class Curve_analysis_2

} // namespace CGALi

CGAL_END_NAMESPACE

#endif // CGAL_ALGEBRAIC_CURVE_KERNEL_CURVE_ANALYSIS_1_H
