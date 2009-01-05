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

#ifndef CGAL_ACK_CURVE_PAIR_ANALYSIS_2_EXACUS_H
#define CGAL_ACK_CURVE_PAIR_ANALYSIS_2_EXACUS_H

#include <CGAL/basic.h>
#include <CGAL/Handle_with_policy.h>
#include <CGAL/Algebraic_curve_kernel_2/Curve_analysis_2_exacus.h>
#include <CGAL/Algebraic_curve_kernel_2/Status_line_CPA_1.h>

#include <SoX/GAPS/types.h>

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

    // internal curve pair type (temporary)
    typedef typename Algebraic_curve_kernel_2::Internal_curve_pair_2
        Internal_curve_pair_2;

    // type of 1-curve analysis
    typedef typename Algebraic_curve_kernel_2::Curve_analysis_2
        Curve_analysis_2;

    // constructors
public:
    // default constructor ()
    Curve_pair_analysis_2_rep() 
    {   }

    // constructs from two curve analysis instances
    Curve_pair_analysis_2_rep(const Curve_analysis_2& ca1,
                              const Curve_analysis_2& ca2) 
        : _m_ca1(ca1), _m_ca2(ca2) {
        
        _m_curve_pair = Internal_curve_pair_2(
            ca1._internal_curve(), ca2._internal_curve());
    }

    // data
    mutable Curve_analysis_2 _m_ca1, _m_ca2;

    // temporarily this implementation is based on Curve_pair_2 from GAPS
    mutable Internal_curve_pair_2 _m_curve_pair;
    
    // befriending the handle
    friend class Curve_pair_analysis_2<Algebraic_curve_kernel_2, Self>;
};
    
//! \brief The class is meant to provide tools to analyse a pair of curves. 
//!
//! Analysis describes the curve pair's interesting points and how they are 
//! connected. The analysis searches for events. Events only occur at a finite 
//! number of x-coordinate. Each such coordinate is covered by a 
//! \c StatusLine_1, originated by the events of a single curve and
//! also the intersections of two curves. These coordinates also define
//! open intervals on the x-axis. \c StatusLine_1 at values in
//! between one such interval differ only in the values of the 
//! \c Algebraic_real_2 entries. Topological information are equal.
template <class AlgebraicCurveKernel_2, 
      class Rep_ = CGALi::Curve_pair_analysis_2_rep<AlgebraicCurveKernel_2> >
class Curve_pair_analysis_2 : public ::CGAL::Handle_with_policy< Rep_ > 
{
    // internal curve pair type (temporary)
    typedef typename Rep_::Internal_curve_pair_2 Internal_curve_pair_2;

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

    //! required by Status_line_CPA_1
    typedef X_coordinate_1 Algebraic_real_1;

    //! required by Status_line_CPA_1
    typedef Xy_coordinate_2 Algebraic_real_2;

    //! type of 1-curve analysis
    typedef typename Algebraic_curve_kernel_2::Curve_analysis_2 
            Curve_analysis_2;

    //! an instance of a size type
    typedef typename Curve_analysis_2::size_type size_type;

    //! myself
    typedef Curve_pair_analysis_2<Algebraic_curve_kernel_2, Rep> Self;

    //! type of a vertical line
    typedef CGALi::Status_line_CPA_1<Self> Status_line_1;
        
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

    /*!\brief
     * constructs a curve pair analysis defined by analyses given by \c _m_ca1
     * and \c _m_ca2.
     *
     * polynomials defining the analysis must be squarefree and coprime.
     */
    explicit Curve_pair_analysis_2(const Curve_analysis_2& ca1,
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

    /*! \brief
     * returns curve analysis for c-"th" curve (0 or 1)
     */
    Curve_analysis_2 curve_analysis(bool c) const
    { 
        if(c == 0)
            return this->ptr()->_m_ca1;
        return this->ptr()->_m_ca2;
    }

    /*! \brief
     * returns number of status lines that encode an event
     */
    size_type number_of_status_lines_with_event() const
    {
        return this->ptr()->_m_curve_pair.num_events();
    }


    /*! \brief
     *  given the i-th event of the curve pair this method returns the
     * id of the event of the corresponding curve analysis \c c (0 or 1),
     * or -1, if the curve has no event at this coordinate.
     *
     * \pre 0 <= i \< number_of_status_lines_with_event()
     */
    size_type event_of_curve_analysis(size_type i, bool c) const
    {
        CGAL_precondition(i >= 0 && i < number_of_status_lines_with_event());
        SoX::Index_triple triple = 
            this->ptr()->_m_curve_pair.event_indices(i);
        return (c == 0 ? triple.ffy : triple.ggy);
    }
    
    size_type event_of_curve_analysis(size_type i, 
                                      const Curve_analysis_2& c) const {
        CGAL_assertion(c.id()==curve_analysis(false).id() ||
                       c.id()==curve_analysis(true).id());
        SoX::Index_triple triple = 
            this->ptr()->_m_curve_pair.event_indices(i);
        return (c.id()==curve_analysis(false).id()) ? triple.ffy : triple.ggy;
    }


    //! Shortcut for coherent tests
    SoX::Index_triple event_indices(size_type i) const {
        return this->ptr()->_m_curve_pair.event_indices(i);
    }

    /*! \brief
     * returns an instance of \c StatusLine_1 at the i-th event
     *
     * \pre 0 <= i \< number_of_status_lines_with_event()
     */
    Status_line_1 status_line_at_event(size_type i) const {
    
        CGAL_precondition(i >= 0 && i < number_of_status_lines_with_event());

#ifdef CGAL_ACK_2_USE_STATUS_LINES
        return this->ptr()->_m_curve_pair.status_line_at_event(*this, i);
          
#else
        typename Internal_curve_pair_2::Event2_slice slice =
            this->ptr()->_m_curve_pair.slice_at_event(i);

        typename Status_line_1::Arc_container arcs(slice.num_arcs());
        std::copy(slice.arcs_at_event().begin(),
            slice.arcs_at_event().end(), arcs.begin());

        Status_line_1 line(i, arcs, *this);

        CGAL_precondition_code(
                int typo;
                typename Status_line_1::Arc_pair pair;
                for (int j = 0; j < slice.num_arcs(); j++) {
                    pair = line.curves_at_event(j);
                    if (pair.first != -1 && pair.second != -1) {
                        typo = 2;
                    } else {
                        typo = (pair.first != -1 ? 0 : 1);
                    }
               /*std::cout << "[" << pair.first << "; " << pair.second << "] and "  <<
                      (int)(slice.arc_at_event(j).first) << "\n";*/
                    CGAL_precondition(typo == slice.arc_at_event(j).first);
                }    
        );
        return line;
#endif // CGAL_ACK_2_USE_STATUS_LINES
    }

    /*! \brief
     * returns an instance of StatusLine_1 of the i-th
     * interval between x-events
     *
     * \pre 0 <= i < number_of_status_lines_with_event()
     */
    Status_line_1 status_line_of_interval(size_type i) const {
    
        CGAL_precondition(i >= 0 && i <= number_of_status_lines_with_event());

#ifdef CGAL_ACK_2_USE_STATUS_LINES

        return this->ptr()->_m_curve_pair.status_line_of_interval(*this, i);
#else

        typename Internal_curve_pair_2::Event2_slice slice =
            this->ptr()->_m_curve_pair.slice_at_interval(i);

        typename Status_line_1::Int_container arcs(slice.num_arcs());
        std::copy(slice.arcs_at_interval().begin(),
            slice.arcs_at_interval().end(), arcs.begin());

        Status_line_1 line(i, arcs, *this);

        CGAL_precondition_code(
                int typo;
                typename Status_line_1::Arc_pair pair;
                for (int j = 0; j < slice.num_arcs(); j++) {
                    pair = line.curves_at_event(j);
                    typo = (pair.first != -1 ? 0 : 1);
                    CGAL_precondition(typo == slice.arc_at_interval(j));
                }
        );
        
        return line;
#endif // CGAL_ACK_2_USE_STATUS_LINES
    }

    /*! \brief
     * returns status_line_at_event(i), if x hits i-th event,
     * otherwise status_line_of_interval(i), where i is the id of the
     * interval \c x lies in.
     *
     * If \c pertub is \c CGAL::NEGATIVE (\c CGAL::POSITIVE) and x states an
     * event, then \c status_line_of_interval(i)
     * (\c status_line_of_interval(i+1)) is returned.
     * 
     * \pre \c x is finite
     */
    Status_line_1 status_line_for_x(X_coordinate_1 x,
        CGAL::Sign perturb = CGAL::ZERO) const {
        
        // CGAL_precondition(x is finite ??);
        size_type i;
        bool is_evt;
        this->ptr()->_m_curve_pair.x_to_index(x, i, is_evt);
        if(is_evt) {
            if(perturb == CGAL::ZERO) {
                Status_line_1 sline = status_line_at_event(i);
                sline._set_x(x);
                return sline;
            }    
            if(perturb == CGAL::POSITIVE)
                i++;
        }
        return status_line_of_interval(i);
    }

    /*! \brief
     * returns an instance of StatusLine_1 at a given \c x
     *
     * \pre \c x is finite
     */
    Status_line_1 status_line_at_exact_x(X_coordinate_1 x) const {
        // CGAL_precondition(x is finite ??);
        size_type i;
        bool is_evt;
        this->ptr()->_m_curve_pair.x_to_index(x, i, is_evt);

        Status_line_1 sline;
        if(is_evt)
            sline = status_line_at_event(i);
        else
            sline = status_line_of_interval(i);
        sline._set_x(x);
        return sline;
    }

    //!@}
protected:

    // temporary method returns curve pair representaion
    Internal_curve_pair_2 _internal_curve_pair() const {
        return this->ptr()->_m_curve_pair;
    }

    // befriending status line class
    friend class CGALi::Status_line_CPA_1<Self>;
}; // class Curve_pair_analysis_2

} // namespace CGALi

CGAL_END_NAMESPACE

#endif // CGAL_ALGEBRAIC_CURVE_KERNEL_CURVE_PAIR_ANALYSIS_2_H
