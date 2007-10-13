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

#ifndef CGAL_ALGEBRAIC_CURVE_KERNEL_CURVE_PAIR_ANALYSIS_2_H
#define CGAL_ALGEBRAIC_CURVE_KERNEL_CURVE_PAIR_ANALYSIS_2_H

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

template < class AlgebraicCurveKernel_2, class Rep_ > 
class Curve_pair_analysis_2;

namespace CGALi {

template < class AlgebraicCurveKernel_2 >
class Curve_pair_analysis_2_rep {

public:
    // this first template argument
    typedef AlgebraicCurveKernel_2 Algebraic_curve_kernel_2;

    // myself
    typedef Curve_pair_analysis_2_rep<Algebraic_curve_kernel_2> Self;

    // type of curve pair
    typedef typename Algebraic_curve_kernel_2::Curve_pair_2 Curve_pair_2;
    
    // type of 1-curve analysis
    typedef typename Algebraic_curve_kernel_2::Curve_analysis_2
        Curve_analysis_2;

    // constructors
public:
    // default constructor ()
    Curve_pair_analysis_2_rep() 
    {   }

     // temporary constructor: directly from curve pair
    Curve_pair_analysis_2_rep(const Curve_pair_2& curve_pair) : 
        _m_ca1(curve_pair.curve1()), _m_ca2(curve_pair.curve2()), 
            _m_curve_pair(curve_pair)
    {   }

    // constructs from two curve analysis instances
    Curve_pair_analysis_2_rep(const Curve_analysis_2& ca1,
        const Curve_analysis_2& ca2) : _m_ca1(ca1), _m_ca2(ca2)
    {  
        _m_curve_pair = Algebraic_curve_kernel_2::get_curve_pair_cache()
            (std::make_pair(ca1.get_polynomial_2(), ca2.get_polynomial_2()));
        //! attention: is it enough to compare ids only ? or for safety
        //! its better to compare polynomials ?
        _m_is_swapped = (_m_curve_pair.curve1().id() !=
             ca1.get_polynomial_2().id());
      //  if(_m_is_swapped)
        //    std::cout << "the content was swapped\n";
    }

    // data
    mutable Curve_analysis_2 _m_ca1, _m_ca2;
    // temporarily this implementation is based on Curve_pair_2 from GAPS
    mutable Curve_pair_2 _m_curve_pair;
    // indicates that the curves in a curve pair were swapped after precaching
    // (this happens when the first curve is defined by a polynomial of higher
    // degree than the second one)
    bool _m_is_swapped;
    
    // befriending the handle
    friend class Curve_pair_analysis_2<Algebraic_curve_kernel_2, Self>;
};
    
//! \brief The class is meant to provide tools to analyse a pair of curves. 
//!
//! Analysis describes the curve pair's interesting points and how they are 
//! connected. The analysis searches for events. Events only occur at a finite 
//! number of x-coordinate. Each such coordinate is covered by a 
//! \c CurvePairVerticalLine_1, originated by the events of a single curve and 
//! also the intersections of two curves. These coordinates also define
//! open intervals on the x-axis. \c CurvePairVerticalLine 1 at values in 
//! between one such interval differ only in the values of the 
//! \c Algebraic_real_2 entries. Topological information are equal.
template <class AlgebraicCurveKernel_2, 
      class Rep_ = CGALi::Curve_pair_analysis_2_rep<AlgebraicCurveKernel_2> >
class Curve_pair_analysis_2 : public ::CGAL::Handle_with_policy< Rep_ > 
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

    //! type of a curve pair
    typedef typename Algebraic_curve_kernel_2::Curve_pair_2 Curve_pair_2;

    //! type of a curve
    typedef typename Algebraic_curve_kernel_2::Curve_2 Curve_2;

    //! type of 1-curve analysis
    typedef typename Algebraic_curve_kernel_2::Curve_analysis_2 
            Curve_analysis_2;

    //! myself
    typedef Curve_pair_analysis_2<Algebraic_curve_kernel_2, Rep> Self;

    //! type of a vertical line
    typedef CGALi::Curve_pair_vertical_line_1<Self>
        Curve_pair_vertical_line_1;
        
    //! the handle superclass
    typedef ::CGAL::Handle_with_policy<Rep> Base;

    //!@}
public:
    //!\name Constructors
    //!@{

    //! \brief default constructor
    Curve_pair_analysis_2() : 
        Base(Rep()) {  
    }

    /*!\brief
     * copy constructor
     */
    Curve_pair_analysis_2(const Self& p) : 
        Base(static_cast<const Base&>(p)) {  
    }

    //! constructs a curve pair analysis defined by analysis
    //! given by \c _m_ca1 and \c _m_ca2. The polynomials defining the analysis
    //! must be squarefree and coprime.
    Curve_pair_analysis_2(const Curve_analysis_2& ca1,
        const Curve_analysis_2& ca2) : 
            Base(Rep(ca1, ca2)) {  
    }
 
    /*!\brief
     * constructs a curve pair analysis from a given represenation
     */
    Curve_pair_analysis_2(Rep rep) : 
        Base(rep) {  
    }

    //!@}
public:
    //!\name Access functions
    //!@{

    //! \brief returns curve analysis for c-"th" curve (0 or 1)
    Curve_analysis_2 get_curve_analysis(bool c) const
    { 
        if(this->ptr()->_m_is_swapped)
            c ^= 1;
        if(c == 0)
            return this->ptr()->_m_ca1;
        return this->ptr()->_m_ca2;
    }

    //! \brief returns number of vertical lines that encode an event
    int number_of_vertical_lines_with_event() const
    {
        return this->ptr()->_m_curve_pair.num_events();
    }

    //! \brief given the i-th event of the curve pair this method returns the
    //! id of the event of the corresponding curve analysis \c c (0 or 1),
    //! or -1, if the curve has no event at this coordinate.
    //!
    //! \pre 0 <= i < number_of_vertical_lines_with_event()
    int event_of_curve_analysis(int i, bool c) const
    {
        CGAL_precondition(i >= 0&&i < number_of_vertical_lines_with_event());
        SoX::Index_triple triple = 
            this->ptr()->_m_curve_pair.event_indices(i);
        if(this->ptr()->_m_is_swapped)
            c ^= 1;
        return (c == 0 ? triple.ffy : triple.ggy);
    }

    //! \brief returns an instance of \c CurvePairVerticalLine_1 at the i-th
    //! event
    //!
    //! \pre 0 <= i < number_of_vertical_lines_with_event()
    Curve_pair_vertical_line_1 vertical_line_at_event(int i) const
    {
        CGAL_precondition(i >= 0&&i < number_of_vertical_lines_with_event());
        return Curve_pair_vertical_line_1(
            this->ptr()->_m_curve_pair.slice_at_event(i), i,
                this->ptr()->_m_is_swapped);
    }

    //! \brief returns an instance of CurvePairVerticalLine_1 of the i-th 
    //! interval between x-events
    //!
    //! \pre 0 <= i < number_of_vertical_lines_with_event() 
    Curve_pair_vertical_line_1 vertical_line_of_interval(int i) const
    {
        CGAL_precondition(i >= 0&&i <= number_of_vertical_lines_with_event());
        
//         std::cout << "interval id: " << 
//             this->ptr()->_m_curve_pair.slice_at_interval(i).id() << "\n";
//         
        return Curve_pair_vertical_line_1(
            this->ptr()->_m_curve_pair.slice_at_interval(i), i,
                this->ptr()->_m_is_swapped); 
    }

    //! \brief returns vertical_line_at_event(i), if x hits i-th event, 
    //! otherwise vertical_line_of_interval(i), where i is the id of the 
    //! interval \c x lies in.
    //!
    //! If \c pertub is \c CGAL::NEGATIVE (\c CGAL::POSITIVE) and x states an 
    //! event, then \c vertical_line_of_interval(i)
    //! (\c vertical_line_of_interval(i+1)) is returned.
    //! 
    //! \pre \c x is finite
    Curve_pair_vertical_line_1 vertical_line_for_x(X_coordinate_1 x, 
        CGAL::Sign perturb = CGAL::ZERO) const
    {
        // CGAL_precondition(x is finite ??);
        int i;
        bool is_evt;
        this->ptr()->_m_curve_pair.x_to_index(x, i, is_evt);
        if(is_evt) {
            if(perturb == CGAL::ZERO)
                return vertical_line_at_event(i);
            if(perturb == CGAL::POSITIVE)
                i++;
        }
        return vertical_line_of_interval(i);
    }

    //! \brief returns an instance of CurvePairVerticalLine_1 at a given \c x
    //!
    //! \pre \c x is finite
    Curve_pair_vertical_line_1 vertical_line_at_exact_x(X_coordinate_1 x) const
    {
        // CGAL_precondition(x is finite ??);
        return vertical_line_for_x(x);
    }

    //!@}
}; // class Curve_pair_analysis_2

} // namespace CGALi

CGAL_END_NAMESPACE

#endif // CGAL_ALGEBRAIC_CURVE_KERNEL_CURVE_PAIR_ANALYSIS_2_H
