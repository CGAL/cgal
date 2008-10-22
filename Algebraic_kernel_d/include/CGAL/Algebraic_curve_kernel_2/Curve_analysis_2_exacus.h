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

#ifndef CGAL_ACK_CURVE_ANALYSIS_2_EXACUS_H
#define CGAL_ACK_CURVE_ANALYSIS_2_EXACUS_H

#include <CGAL/basic.h>
#include <CGAL/Handle_with_policy.h>
#include <CGAL/Arr_enums.h>

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

    // supporting polynomial type
    typedef typename Algebraic_curve_kernel_2::Polynomial_2
        Polynomial_2;

    // internal curve type (temporary)
    typedef typename Algebraic_curve_kernel_2::Internal_curve_2
        Internal_curve_2;

    // constructors
public:
    // default constructor ()
    Curve_analysis_2_rep()  
    {  }
    
    // standard constructor
    Curve_analysis_2_rep(const Polynomial_2& f)
    {
//////////////////////////////////////////////////////////////////////////////
//// this should be eliminated in the final version once Intenal_curve_2 is
//// replaced by Curve_analysis_2
#ifdef CGAL_USE_CnX_KERNEL
		_m_curve = Internal_curve_2(f);
#else
        _m_curve = Internal_curve_2::get_curve_cache()(f);
#endif
//////////////////////////////////////////////////////////////////////////////
    }

    // data
    // temporarily this implementation uses underlying Curve_2 from SweepX
    mutable Internal_curve_2 _m_curve;
    
    // befriending the handle
    friend class Curve_analysis_2<Algebraic_curve_kernel_2, Self>;
};
    
//! \brief The class is meant to provide tools to analyse a single curve. 
//! 
//! Analysis describes the curves interesting points and how they are 
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
    // type of intenal curve (temporary)
    typedef typename Rep_::Internal_curve_2 Internal_curve_2;

public:
    //!@{
    //! \name public typedefs

    //! this instance's first template parameter
    typedef AlgebraicCurveKernel_2 Algebraic_curve_kernel_2;

    //! this instance's second template parameter
    typedef Rep_ Rep;
    
    //! x-coordinate type
    typedef typename Algebraic_curve_kernel_2::X_coordinate_1 X_coordinate_1;

    //! y-coordinate type
    typedef X_coordinate_1 Y_coordinate_1;

    //! type of a curve point
    typedef typename Algebraic_curve_kernel_2::Xy_coordinate_2 Xy_coordinate_2;

    //! required by Status_line_CA_1
    typedef X_coordinate_1 Algebraic_real_1;

    //! required by Status_line_CA_1
    typedef Xy_coordinate_2 Algebraic_real_2;

    //! supporting polynomial type
    typedef typename Algebraic_curve_kernel_2::Polynomial_2
        Polynomial_2;

    //! myself
    typedef Curve_analysis_2<Algebraic_curve_kernel_2, Rep> Self;

    //! an instance of a size type
    typedef int size_type;

    //! type of a vertical line
    typedef CGALi::Status_line_CA_1<Self> Status_line_1;
        
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

    /*!\brief
     * constructs a curve analysis from a given polynomial
     *
     * for safety purposes implicit conversion from \c Polynomial_2 is disabled
     */
    explicit Curve_analysis_2(const Polynomial_2& f) :
        Base(Rep(f)) {
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

    /*!\brief
     * returns the defining polynomial of the curve analysis
     */
    Polynomial_2 polynomial_2() const
    { 
        return this->ptr()->_m_curve.f();
    }

    // temporary method to access internal curve representation
    Internal_curve_2 _internal_curve() const {
        return this->ptr()->_m_curve;
    }

    //! \brief returns number of vertical lines that encode an event
    size_type number_of_status_lines_with_event() const
    {
        return _internal_curve().num_events();
    }

    //! \brief returns an instance of StatusLine_1 at the i-th
    //! event
    //!
    //! \pre 0 <= i < number_of_status_lines_with_event()
    Status_line_1 status_line_at_event(size_type i) const
    {
        CGAL_precondition(i >= 0&&i < number_of_status_lines_with_event());

#ifdef CGAL_ACK_2_USE_STATUS_LINES
        return _internal_curve().status_line_at_event(*this, i);

#else
        typedef typename Internal_curve_2::Event1_info Event1_info;
        Event1_info info = _internal_curve().event_info(i);

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
        sline.set_isolator(info.refinement());
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
        typedef typename Internal_curve_2::Event1_info Event1_info;
        size_type n_arcs = _internal_curve().arcs_over_interval(i);
        Status_line_1 sline
            (X_coordinate_1(_internal_curve().boundary_value_in_interval(i)),
             i,
             *this,
             n_arcs);
#ifndef CGAL_ACK_2_USE_STATUS_LINES
        Event1_info info 
            = _internal_curve().event_info_at_x
            (_internal_curve().boundary_value_in_interval(i));
        sline.set_isolator(info.refinement());
#endif
        return sline;
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
        
        typedef typename Internal_curve_2::Event1_info Event1_info;
        Event1_info info 
            = _internal_curve().event_info_at_x(x);
        size_type n_arcs = _internal_curve().arcs_over_interval(i);
        Status_line_1 sline(x,i,*this,n_arcs);
        sline.set_isolator(info.refinement());
        return sline;
    }

    /*!\brief
     * returns the index of the event at the status line defined by
     * \c s x-coordinate, or -1 if \c s does not lie on the curve.
     */
    size_type find(const Xy_coordinate_2& s) const {
        return 0;
    }

    /*!\brief
     * returns a \c CGAL::Object that encodes the asymptotic value of a
     * curve-arc approaching the left or the right boundary \c loc of the
     * underlying parameter space.
     *
     * Allowed instantiations of the \c CGAL::Object are \c Algebraic_real_1 ,
     * in case the x-asympote of the arc is finite, or
     * \c CGAL::ARR_BOTTOM_BOUNDARY and \c CGAL::ARR_TOP_BOUNDARY in case
     * the defined arc approaches the respective corners of the parameter
     * space.
     *
     * \pre \c loc is either \c CGAL::ARR_LEFT_BOUNDARY or
     *  \c CGAL::ARR_RIGHT_BOUNDARY
     */  
     CGAL::Object asymptotic_value_of_arc(CGAL::Arr_parameter_space loc,
             size_type arcno) const {

         CGAL_precondition(loc == CGAL::ARR_LEFT_BOUNDARY ||
            loc == CGAL::ARR_RIGHT_BOUNDARY);
         
         return (loc == CGAL::ARR_LEFT_BOUNDARY ?
                 _internal_curve().
                     horizontal_asymptote_for_arc_to_minus_infinity(arcno) :
                 _internal_curve().
                     horizontal_asymptote_for_arc_to_plus_infinity(arcno));
            
     }

    //!@}

    // friend function for id-based hashing
    friend std::size_t hash_value(const Self& x) {
        return static_cast<std::size_t>(x.id());
    }

    friend std::ostream& operator <<(std::ostream& os, const Self&) {
        os << "Unfortunately output operator is not implemented for "
        " Curve_analysis_2. If you feel like doing this - go ahead! \n";
        return os;
    }
  
}; // class Curve_analysis_2

} // namespace CGALi

CGAL_END_NAMESPACE

#endif // CGAL_ALGEBRAIC_CURVE_KERNEL_CURVE_ANALYSIS_1_H
