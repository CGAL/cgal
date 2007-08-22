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

#ifndef CGAL_ALGEBRAIC_CURVE_KERNEL_CURVE_VERTICAL_LINE_1_H
#define CGAL_ALGEBRAIC_CURVE_KERNEL_CURVE_VERTICAL_LINE_1_H

#include <CGAL/basic.h>
#include <CGAL/Handle_with_policy.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

template < class CurveAnalysis_2, class Rep_ > 
class Curve_vertical_line_1;

template <class CurveAnalysis_2, class Rep>
std::ostream& operator<< (std::ostream&, 
    const Curve_vertical_line_1<CurveAnalysis_2, Rep>&);

template < class CurveAnalysis_2 >
class Curve_vertical_line_1_rep {

    // this template argument
    typedef CurveAnalysis_2 Curve_analysis_2;

    // myself
    typedef Curve_vertical_line_1_rep<Curve_analysis_2> Self;

    // type of x-coordinate
    typedef typename Curve_analysis_2::X_coordinate_1
                X_coordinate_1; 

    // type of a curve point
    typedef typename Curve_analysis_2::Xy_coordinate_2
                Xy_coordinate_2; 

    //! type of underlying \c Event1_info
    typedef typename Curve_analysis_2::Curve_2::Event1_info Event1_info;

    // constructors
public:
    // default constructor ()
    Curve_vertical_line_1_rep() 
    {   }

    // standard constructor
    Curve_vertical_line_1_rep(const Event1_info& info,
            const Curve_analysis_2& ca_2, int index) : 
        _m_event_info(info), _m_ca_2(ca_2), _m_index(index)
    {   }
    
    // data
    // temporary added underlying Event1_info object
    mutable Event1_info _m_event_info;
    
    mutable Curve_analysis_2 _m_ca_2; // supporting curve analysis
    
    mutable int _m_index; // this vertical line id (# of event or # of
            // intervaldepending on whether or not this vertical line encodes
            // an event
        
    // befriending the handle
    friend class Curve_vertical_line_1<Curve_analysis_2, Self>;
};

//! \brief The class provides information about the intersections of a curve 
//! with a vertical line at a given finite x-coordinate. 
//!
//! Note that a curve can have a vertical line component at this coordinate
//! and non-vertical components may intersect the vertical line respectively. 
//! With the help of this class' methods one is able to compute the local 
//! topology of the curve at the given vertical line. Note that vertical lines 
//! at x = +/-oo are not allowed, since different events (curve ends going to 
//! infinity with different non-horizontal asymptotes) would have equal 
//! y-coordinate (+/-oo), which confuses more than it helps. Note in addition 
//! that curve ends approaching the vertical asymptote introduce an event 
//! (depending on whether approaching +oo or -oo - but the event with 
//! coordinates (x,?oo), resp. (x,+oo), occur only once, if they occur, and 
//! they imply not to be associated with a instance of \c Algebraic_real_2.
template <class CurveAnalysis_2, 
          class Rep_ = CGALi::Curve_vertical_line_1_rep<CurveAnalysis_2> >
class Curve_vertical_line_1
      : public ::CGAL::Handle_with_policy< Rep_ > {
public:
    //!@{
    //!\name typedefs

    //! this instance's first template parameter
    //! model of AlgebraicKernel_d_2::CurveAnalysis_2
    typedef CurveAnalysis_2 Curve_analysis_2;
    
    //! this instance's second template parameter
    typedef Rep_ Rep;

    //! this instance itself
    typedef Curve_vertical_line_1<Curve_analysis_2, Rep> Self;

    //! type of x-coordinate
    typedef typename Curve_analysis_2::X_coordinate_1 X_coordinate_1; 

    //! type of a curve point
    typedef typename Curve_analysis_2::Xy_coordinate_2 Xy_coordinate_2;

    //! type of underlying \c Event1_info
    typedef typename Curve_analysis_2::Curve_2::Event1_info Event1_info;
    
     //! the handle superclass
    typedef ::CGAL::Handle_with_policy< Rep > Base;
    
    //!@}
public:
    //!\name constructors
    //!@{

    /*!\brief
     * Default constructor
     */
    Curve_vertical_line_1() : 
        Base(Rep()) {   
    }

    /*!\brief
     * copy constructor
     */
    Curve_vertical_line_1(const Self& p) : 
            Base(static_cast<const Base&>(p)) {  
    }

    /*!\brief
     * Constructs a new instance from \c Event1_info object, supporting
     * \c Curve_analysis_2 and vertical line \c index
     *
     * \c index encodes the index of event or interval depending on whether
     * this vertical line is an event-line or lies in the interval
     */
    Curve_vertical_line_1(const Event1_info& info,
            const Curve_analysis_2& ca_2, int index) : 
        Base(Rep(info, ca_2, index)) {   
        
        CGAL_precondition(id >= 0 && index <
             this->ptr()->_m_ca_2.
             number_of_vertical_lines_with_event() + (is_event()? 0: 1));
    }
        
    /*!\brief
     * constructs from a given represenation
     */
    Curve_vertical_line_1(Rep rep) : 
        Base(rep) {  
    }
    
    //!@}
public:
    //!\name access functions
    //!@{
    
    //! \brief returns the x-coordinate of the vertical line (always a finite
    //! value).
    X_coordinate_1 x() const
    {
        return this->ptr()->_m_event_info.x();
    }
    
    //! \brief returns this vertical line supporting curve analysis
    Curve_analysis_2 get_curve_analysis_2() const
    {
        return this->ptr()->_m_ca_2;
    }
    
    //! \brief returns this vertical line index (event or interval index)
    int get_index() const
    {
        return this->ptr()->_m_index;
    }
        
    //! \brief returns true in case the given curve contains the vertical line
    //! as a component
    bool covers_line() const
    {    
        return this->ptr()->_m_event_info.has_vertical_line();
    }

    //! \brief returns number of distinct and finite intersections of a curve
    //! with a (intended) vertical line ignoring a real vertical line
    //! component of the curve at the given x-coordinate.
    int number_of_events() const
    {
        return (this->ptr()->_m_event_info.num_arcs());
    }

    //!\brief  returns an object of type \c Xy_coordinate_2 for the j-th event
    //!
    //! \pre 0 <= j < num_of_events()
    // TODO add get_xy_coordinate() with the same functionality
    Xy_coordinate_2 get_algebraic_real_2(int j) const
    {
        CGAL_precondition(0 <= j&&j < number_of_events());
        // how to get the pointer to the curve ?
        // we have to store such a pointer for vertical line
        // TODO isn't it a good idea to cache the construction?
        // this way we get filter failures, e.g., the same point 
        // is represented twice when accessing it twice
        return Xy_coordinate_2(x(), 
            this->ptr()->_m_ca_2.get_polynomial_2(), j);
    }

    //!\brief returns the number of branches of the curve connected to j-th
    //! event immediately to the left, to the right, respectively, as a pair of
    //! unsigned int ignoring vertical curve components at the given 
    //! x-coordinate.
    //!
    //! \pre 0 <= j < num_of_events()
    std::pair<int, int> get_number_of_incident_branches(int j) const
    {
        CGAL_precondition(0 <= j&&j < number_of_events());
        return std::make_pair(
            this->ptr()->_m_event_info.num_arcs_left(j),
            this->ptr()->_m_event_info.num_arcs_right(j));
    }

    //!\brief returns the number of vertical asymptotes at x of the curve
    //! approaching y=-oo from left and right. A vertical line being component 
    //! of the curve is ignored.
    std::pair<int, int> 
        get_number_of_branches_approaching_minus_infinity() const
    {
        int left, right, dummy;
        this->ptr()->_m_event_info.num_arcs_approaching_vertical_asymptote(
            left, right, dummy, dummy);
        return std::make_pair(left, right);
    }

    //!\brief returns the number of vertical asymptotes at x of the curve
    //! approaching y=+oo from left and right. A vertical line being component 
    //! of the curve is ignored.
    std::pair<int, int> 
        get_number_of_branches_approaching_plus_infinity() const
    {
        int left, right, dummy;
        this->ptr()->_m_event_info.num_arcs_approaching_vertical_asymptote(
            dummy, dummy, left, right);
        return std::make_pair(left, right);
    }

    //! \brief returns true if curve has vertical line component or curve \c f
    //! has intersection with f_y at \c x
    bool is_event() const
    {
        return this->ptr()->_m_event_info.is_event();
    }
    
    // temporary access function (for testing)
    Event1_info get_info() const
    {
        return this->ptr()->_m_event_info;
    }

    //!@}
}; // class Curve_vertical_line_1

template <class CurveAnalysis_2, class Rep>
std::ostream& operator<< (
        std::ostream& os, 
        const CGALi::Curve_vertical_line_1<CurveAnalysis_2, Rep>& cp_line) {
    os << (cp_line.get_info());
    return os;
}

} // namespace CGALi

CGAL_END_NAMESPACE

#endif // CGAL_ALGEBRAIC_CURVE_KERNEL_CURVE_VERTICAL_LINE_1_H

