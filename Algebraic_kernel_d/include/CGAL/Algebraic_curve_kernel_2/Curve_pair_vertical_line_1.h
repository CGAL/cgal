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

#ifndef CGAL_ALGEBRAIC_CURVE_KERNEL_CURVE_PAIR_VERTICAL_LINE_1_H
#define CGAL_ALGEBRAIC_CURVE_KERNEL_CURVE_PAIR_VERTICAL_LINE_1_H

#include <CGAL/basic.h>
#include <CGAL/Handle_with_policy.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

template < class CurvePairAnalysis_2, class Rep_ > 
class Curve_pair_vertical_line_1;

template <class CurvePairAnalysis_2, class Rep>
std::ostream& operator<< (std::ostream&, 
    const Curve_pair_vertical_line_1<CurvePairAnalysis_2, Rep>&);

template < class CurvePairAnalysis_2 >
class Curve_pair_vertical_line_1_rep {

    // this template argument
    typedef CurvePairAnalysis_2 Curve_pair_analysis_2;

    // myself
    typedef Curve_pair_vertical_line_1_rep<Curve_pair_analysis_2> Self;

    // type of x-coordinate
    typedef typename Curve_pair_analysis_2::X_coordinate_1
                X_coordinate_1; 

    // type of a curve point
    typedef typename Curve_pair_analysis_2::Xy_coordinate_2
                Xy_coordinate_2; 

    // type of internal algebraic curve pair
    typedef typename Curve_pair_analysis_2::Curve_pair_2 Curve_pair_2;

    // type of underlying \c Event2_slice
    typedef typename Curve_pair_2::Event2_slice Event2_slice;
    
    // constructors
public:
    // default constructor ()
    Curve_pair_vertical_line_1_rep()  
    {   }
    
    // standard constructor
    Curve_pair_vertical_line_1_rep(const Event2_slice& event_slice) : 
        _m_event_slice(event_slice)
    {   }
    
    // data
    // underlying Event2_slice object (temporary)
    mutable Event2_slice _m_event_slice;

    // befriending the handle
    friend class Curve_pair_vertical_line_1<Curve_pair_analysis_2, Self>;
};

//! \brief The class provides information about the intersections of a pair of 
//! curves with a (intended) vertical line (ignoring vertical lines of the 
//! curves themselves). 
//! 
//! Each intersection of a curve with the vertical line defined by some given x
//! induces an event. An event can be asked for its coordinates 
//! (\c Algebraic_real_2) and the involved curve(s). Note that the involvement 
//! also holds for curve ends approaching the vertical asymptote. 
//! Curve_pair_vertical_line_1 at x = +/-oo are not allowed.
template <class CurvePairAnalysis_2, 
      class Rep_ = CGALi::Curve_pair_vertical_line_1_rep<CurvePairAnalysis_2> >
class Curve_pair_vertical_line_1 : 
    public ::CGAL::Handle_with_policy< Rep_ > 
{
public:
    //!@{
    //!\name typedefs

    //! this instance's first template parameter
    typedef CurvePairAnalysis_2 Curve_pair_analysis_2;
    
    //! this instance's second template parameter
    typedef Rep_ Rep;

    //! this instance itself
    typedef Curve_pair_vertical_line_1<Curve_pair_analysis_2, Rep> Self;

    //! type of internal algebraic curve pair
    typedef typename Curve_pair_analysis_2::Curve_pair_2 Curve_pair_2;

    //! type of x-coordinate
    typedef typename Curve_pair_analysis_2::X_coordinate_1 X_coordinate_1; 

    //! type of a curve point
    typedef typename Curve_pair_analysis_2::Xy_coordinate_2 Xy_coordinate_2;

    //! type of underlying \c Event2_slice
    typedef typename Curve_pair_2::Event2_slice Event2_slice;

     //! the handle superclass
    typedef ::CGAL::Handle_with_policy< Rep > Base;
    
    //!@}
public:
    //!\name constructors
    //!@{

    /*!\brief
     * Default constructor
     */
    Curve_pair_vertical_line_1() : 
        Base(Rep()) {   
    }

    /*!\brief
     * copy constructor
     */
    Curve_pair_vertical_line_1(const Self& p) : 
            Base(static_cast<const Base&>(p)) {  
    }

    /*!\brief
     * Constructs a new instance from \c Event2_slice object
     */
    Curve_pair_vertical_line_1(const Event2_slice& slice) : 
        Base(Rep(slice)) {   
    }
        
    /*!\brief
     * constructs from a given represenation
     */
    Curve_pair_vertical_line_1(Rep rep) : 
        Base(rep) {  
    }
    
    //!@}
public:
    //!@}
    //!\name access functions
    //!@{
    
    // overriden member function while the class doesn't have it's own
    // representation
    bool is_identical(const Self& line) const
    {
        return this->ptr()->_m_event_slice.is_identical(line.get_slice());
    }
    
    //! \brief returns the x-coordinate of the vertical line (always a finite
    //! value).
    X_coordinate_1 x() const
    {
        int id = this->ptr()->_m_event_slice.id();
        if(is_event())
            return this->ptr()->_m_event_slice.curve_pair().event_x(id);
        return X_coordinate_1(this->ptr()->_m_event_slice.curve_pair().
                boundary_value_in_interval(id));
    }
        
    //! \brief returns number of distinct and finite intersections of a pair 
    //! of curves  with a (intended) vertical line ignoring a real vertical 
    //! line component of the curve at the given x-coordinate.
    int number_of_events() const
    {
        return (this->ptr()->_m_event_slice.num_arcs());
    }

    //! \brief returns the y-position of the k-th event of the c-th (0 or 1)
    //! curve in the sequence of events. 
    //!
    //! Note that each event is formed by the first, the second, or both curves
    //!
    //! \pre 0 <= k < "number of arcs defined for curve[c] at x()"
    int get_event_of_curve(int k, bool c) const
    {
        typename Event2_slice::Int_container ic = 
            this->ptr()->_m_event_slice.arcno_to_pos
                    ((c == 0 ? SoX::CURVE1 : SoX::CURVE2));
        CGAL_precondition(0 <= k && k < static_cast<int>(ic.size()));
        return (ic[k]);
    }

    //! \brief returns the multiplicity of intersection defined at event with
    //! position \c j. May return 0 in case multiplicity is unknown.
    //!
    //! \pre There is an intersection of both curves at j-th event
    //! \pre 0 <= j < number_of_events()
    int get_multiplicity_of_intersection(int j) const
    {
        CGAL_precondition(0 <= j && j < number_of_events());
        CGAL_precondition(is_event());
        const typename Event2_slice::Arc_at_event_2& arc = 
            this->ptr()->_m_event_slice.arc_at_event(j);
        CGAL_precondition(arc.first == SoX::CURVE_BOTH);
        return arc.second;
    }

    //! \brief returns a pair of \c int indicating whether event \c j is formed
    //! by which arc numbers of the first and the second curve, or -1, if the 
    //! corresponding curve is not involved.
    //!
    //! \pre 0 <= j < number_of_events()
    std::pair<int, int> get_curves_at_event(int j) const
    {
        CGAL_precondition(0 <= j && j < number_of_events());
        typedef std::vector<typename Event2_slice::Arc_at_event_2> 
            Arc_vector;
        Arc_vector res(2); // maximum two elements
        typename Arc_vector::iterator end = 
            this->ptr()->_m_event_slice.pos_to_arc(j, res.begin()), it;
        int idx[2] = {-1, -1};
        it = res.begin();
        idx[((*it).first == SoX::CURVE1 ? 0 : 1)] = (*it).second;
        if(++it != end) 
            idx[((*it).first == SoX::CURVE1 ? 0 : 1)] = (*it).second;
        return std::make_pair(idx[0], idx[1]);
    }

    //! \brief returns true if a curve has an event or in case there is an
    //! intersection of both curves.
    bool is_event() const
    {    
        return this->ptr()->_m_event_slice.is_event();
    }

    //! \brief returns true if there is an intersection of both curves.
    bool is_intersection() const
    {
        CGAL_precondition(is_event());
        return this->ptr()->_m_event_slice.is_intersection();
    }
    
    // temporary access function (for testing)
    Event2_slice get_slice() const
    {
        return this->ptr()->_m_event_slice;
    }
    
    //!@}
}; // class Curve_pair_vertical_line_1

template <class CurvePairAnalysis_2, class Rep>
std::ostream& operator<< (
        std::ostream& os, 
        const CGALi::Curve_pair_vertical_line_1<CurvePairAnalysis_2, Rep>&
             cpv_line) {
    os << (cpv_line.get_slice());
    return os;
}

} // namespace CGALi

CGAL_END_NAMESPACE

#endif // CGAL_ALGEBRAIC_CURVE_KERNEL_CURVE_PAIR_VERTICAL_LINE_1_H

