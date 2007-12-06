// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#ifndef CGAL_ALGEBRAIC_CURVE_KERNEL_CURVE_ANALYSIS_2_H
#define CGAL_ALGEBRAIC_CURVE_KERNEL_CURVE_ANALYSIS_2_H

#include <CGAL/basic.h>
#include <CGAL/Handle_with_policy.h>
#include <CGAL/Algebraic_curve_kernel_2/Xy_coordinate_2.h>
#include <CGAL/Algebraic_curve_kernel_2/Status_line_CA_1.h>

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
    Curve_analysis_2_rep(const Curve_2& curve) : _m_curve(curve)
    {   }

    // data
    // temporarily this implementation uses underlying Curve_2 from SweepX
    mutable Curve_2 _m_curve;
    
    // befriending the handle
    friend class Curve_analysis_2<Algebraic_curve_kernel_2, Self>;
};
    
//! \brief The class is meant to provide tools to analyse a single curve. 
//! 
//! Analysis describes the curve’s interesting points and how they are 
//! connected. The analysis searches for events. Events only occur at a finite 
//! number of x-coordinates. Each such coordinate defines a 
//! \c StatusLine_1 of an event. These coordinates also define open
//! intervals on the x-axis. Different \c StatusLine_1 at values within
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

    //! an instance of a size type
    typedef int size_type;

    //! type of a vertical line
    typedef CGALi::Status_line_CA_1<Self> Status_line_1;
        
    //! type of underlying Event1_info structure
    typedef typename Curve_2::Event1_info Event1_info;

    //! the handle superclass
    typedef ::CGAL::Handle_with_policy<Rep> Base;
    
    //!@}
public:
    //!\name Constructors
    //!@{

    //! \brief default constructor
    Curve_analysis_2() : 
        Base(Rep()) {  
    }

    /*!\brief
     * copy constructor
     */
    Curve_analysis_2(const Self& p) : 
        Base(static_cast<const Base&>(p)) {  
    }

    //! \brief constructs a curve analysis from a given \c Curve_2 object
    //!
    //! for safety purposes implicit conversion from \c Curve_2 is disabled
    explicit Curve_analysis_2(const Curve_2& c) : 
        Base(Rep(c)) {  
    }
           
    /*!\brief
     * constructsa curve analysis from a given represenation
     */
    Curve_analysis_2(Rep rep) : 
        Base(rep) {  
    }

    //!@}
public:
    //!\name Access functions
    //!@{

    //! \brief returns the defining polynomial of the analysis
    Curve_2 polynomial_2() const
    { 
        return this->ptr()->_m_curve;
    }

    //! \brief alias for \c polynomial_2()
    Curve_2 curve_2() const
    { 
        return polynomial_2();
    }

    //! \brief returns number of vertical lines that encode an event
    size_type number_of_status_lines_with_event() const
    {
        return this->ptr()->_m_curve.num_events();
    }

    //! \brief returns an instance of StatusLine_1 at the i-th
    //! event
    //!
    //! \pre 0 <= i < number_of_status_lines_with_event()
    Status_line_1 status_line_at_event(size_type i) const
    {
        CGAL_precondition(i >= 0&&i < number_of_status_lines_with_event());

#ifdef CGAL_ACK_2_USE_STATUS_LINES
        return this->ptr()->_m_curve.status_line_at_event(*this, i);

#else
        typedef typename Curve_2::Event1_info Event1_info;
        Event1_info info = this->ptr()->_m_curve.event_info(i);

        typedef typename Event1_info::Arc_container EArc_container;
        const EArc_container& src = info.get_arcs();

        typename Status_line_1::Arc_container dst(info.num_arcs());
        int k = 0;
        
        for(typename EArc_container::const_iterator eit = src.begin();
                eit != src.end(); eit++, k++) 
            dst[k] = std::make_pair(eit->num_arcs_left(),
                eit->num_arcs_right());
                
        Status_line_1 sline(info.x(), i, *this, info.num_arcs_left(),
            info.num_arcs_right(), dst, info.has_vertical_line());

        typename Status_line_1::Arc_pair minus_inf, plus_inf;
        info.num_arcs_approaching_vertical_asymptote(minus_inf.first,
            minus_inf.second, plus_inf.first, plus_inf.second);

        sline._set_number_of_branches_approaching_infinity(minus_inf,
                plus_inf);
        return sline;
#endif // CGAL_ACK_2_USE_STATUS_LINES
    }

    //! \brief returns an instance of StatusLine_1 of the i-th 
    //! interval
    //!
     //! \pre 0 <= i < number_of_status_lines_with_event()
    Status_line_1 status_line_of_interval(size_type i) const
    {
        CGAL_precondition(i >= 0&&i <= number_of_status_lines_with_event());
        
        size_type n_arcs = this->ptr()->_m_curve.arcs_over_interval(i);
        return Status_line_1(X_coordinate_1(this->ptr()->_m_curve.
              boundary_value_in_interval(i)), i, *this, n_arcs);
    }

    //! \brief returns status_line_at_event(i), if x hits i-th event,
    //! otherwise status_line_of_interval(i), where i is the id of the
    //! interval \c x lies in
    //!
    //! If \c pertub is \c CGAL::NEGATIVE (CGAL::POSITIVE) and x states an 
    //! event, then \c status_line_of_interval(i)
    //! (\c status_line_of_interval(i+1)) is returned.
    //! 
    //! \pre \c x is finite
    Status_line_1 status_line_for_x(X_coordinate_1 x,
        CGAL::Sign perturb = CGAL::ZERO) const
    {
        // CGAL_precondition(x is finite ??);
        size_type i;
        bool is_evt;
        this->ptr()->_m_curve.x_to_index(x, i, is_evt);
        if(is_evt) {
            if(perturb == CGAL::ZERO)
                return status_line_at_event(i);
            if(perturb == CGAL::POSITIVE)
                i++;
        } 
        return status_line_of_interval(i);
    }

    //! \brief returns an instance of StatusLine_1 at a given \c x
    //!
    //! \pre \c x is finite
    Status_line_1 status_line_at_exact_x(X_coordinate_1 x) const
    {
        // CGAL_precondition(x is finite ??);
        size_type i;
        bool is_evt;
        this->ptr()->_m_curve.x_to_index(x, i, is_evt);
        if(is_evt) 
            return status_line_at_event(i);
        
        size_type n_arcs = this->ptr()->_m_curve.arcs_over_interval(i);
        return Status_line_1(x, i, *this, n_arcs);
    }
    
    //!@}
}; // class Curve_analysis_2

} // namespace CGALi

CGAL_END_NAMESPACE

#endif // CGAL_ALGEBRAIC_CURVE_KERNEL_CURVE_ANALYSIS_1_H
