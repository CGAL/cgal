// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#ifndef CGAL_ALGEBRAIC_CURVE_KERNEL_STATUS_LINE_CA_1_H
#define CGAL_ALGEBRAIC_CURVE_KERNEL_STATUS_LINE_CA_1_H

#include <CGAL/basic.h>
#include <CGAL/Handle_with_policy.h>

#include <CGAL/Algebraic_kernel_d/Bitstream_descartes.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes_rndl_tree_traits.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_coefficient_kernel_at_alpha.h>

namespace CGAL {

namespace internal {

template < class CurveAnalysis_2, class Rep_ >
class Status_line_CA_1;

template <class CurveAnalysis_2, class Rep>
std::ostream& operator<< (std::ostream&,
    const Status_line_CA_1<CurveAnalysis_2, Rep>&);

#if !CGAL_ACK_USE_EXACUS
template < typename AlgebraicCurveKernel_2 >
class Event_line_builder;

template < typename AlgebraicCurveKernel_2 >
class Shear_transformation;
#endif


template < class AlgebraicCurveKernel_2 >
class Status_line_CA_1_rep {

public:

    // this template argument
    typedef AlgebraicCurveKernel_2 Algebraic_curve_kernel_2;

    typedef typename Algebraic_curve_kernel_2::Curve_analysis_2
      Curve_analysis_2;

    // myself
    typedef Status_line_CA_1_rep<Algebraic_curve_kernel_2> Self;

    // type of x-coordinate
    typedef typename Curve_analysis_2::Algebraic_real_1
                Algebraic_real_1;

    // type of a curve point
    typedef typename Curve_analysis_2::Algebraic_real_2
                Algebraic_real_2;

    // type of bivariate Polynomial
    typedef typename Curve_analysis_2::Polynomial_2
                Polynomial_2;

    // an instance of a size type
    typedef typename Curve_analysis_2::size_type size_type;

    // encodes number of arcs to the left and to the right
    typedef std::pair<size_type, size_type> Arc_pair;

    // container of arcs
    typedef std::vector<Arc_pair> Arc_container;

    // Isolator type
    typedef typename Curve_analysis_2::Bitstream_descartes Bitstream_descartes;

    // constructors

    // default constructor ()
    Status_line_CA_1_rep()
    {   }

    // constructs status line over interval
    Status_line_CA_1_rep(
            Algebraic_real_1 x, size_type i,
            const Curve_analysis_2& ca, size_type n_arcs) :
            _m_kernel(ca.kernel()),
            _m_x(x), _m_index(i), _m_ca(ca),/*_m_num_arcs(n_arcs, n_arcs),*/
            _m_total_arcs(n_arcs), _m_vertical_line(false), _m_event(false),
            _m_num_arcs_minus_inf(0, 0), _m_num_arcs_plus_inf(0, 0),
                _m_xy_coords(n_arcs)  {
    }

    // constructs status line at events
    Status_line_CA_1_rep(
        Algebraic_real_1 x, size_type i,
        const Curve_analysis_2& ca,
        size_type , size_type ) :
      _m_kernel(ca.kernel()),
             _m_x(x), _m_index(i), _m_ca(ca),
             /*_m_num_arcs(n_arcs_left, n_arcs_right),*/ _m_total_arcs(0),
             _m_vertical_line(false), _m_event(true),
             _m_num_arcs_minus_inf(0, 0), _m_num_arcs_plus_inf(0, 0) {
    };

    //! kernel instance
    // TODO remove kernel?
    const Algebraic_curve_kernel_2 *_m_kernel;

    //! x-coordinate of event info
    mutable Algebraic_real_1 _m_x;

    //! this status line id (# of event or # of interval depending on whether
    //! or not this status line encodes an event)
    size_type _m_index;

    //! underlying curve analysis
    Curve_analysis_2 _m_ca;

    //! number of incident arcs to the left and to the right
    //Arc_pair _m_num_arcs;

    //! sequence of arcs crossing this status line (valid only event lines)
    mutable boost::optional<Arc_container> _m_arcs;

    //! number of arcs intersecting this status line
    mutable int _m_total_arcs;

    //! curve has vertical line at this x-coordinate
    mutable bool _m_vertical_line;

    //! decsribes an event
    mutable bool _m_event;

    //! number of arcs running down the pole
    Arc_pair _m_num_arcs_minus_inf;

    //! number of arcs running up the pole
    Arc_pair _m_num_arcs_plus_inf;

    /*// matchings valid?
    mutable bool matching_valid_;

    // match side arcs to event arcs
    mutable std::vector< int > matching_[2];

    // total number of arcs
    int numarcs_at_;

    // arc number of lowest event
    int lowest_event_;

    // stores multiplicities
    std::vector< int > multiplicities_;*/

    // stores algebraic real over the vertical line
    mutable std::vector<boost::optional< Algebraic_real_2 > >_m_xy_coords;

    // stores the isolator instance
    mutable boost::optional<Bitstream_descartes> isolator;

     // befriending the handle
    friend class Status_line_CA_1<Curve_analysis_2, Self>;
    //friend class Curve_analysis_2;
    //friend class Event_line_builder<Curve_analysis_2>;
    //friend class Shear_transformation<Curve_analysis_2>;

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
template <class AlgebraicCurveKernel_2,
          class Rep_ = internal::Status_line_CA_1_rep<AlgebraicCurveKernel_2> >
class Status_line_CA_1
      : public ::CGAL::Handle_with_policy< Rep_ > {
public:
    //!@{
    //!\name typedefs

    //! this instance's first template parameter
    //! model of AlgebraicKernel_d_2
    typedef AlgebraicCurveKernel_2 Algebraic_curve_kernel_2;

    typedef typename Algebraic_curve_kernel_2::Curve_analysis_2 Curve_analysis_2;

    //! this instance's second template parameter
    typedef Rep_ Rep;

    //! this instance itself
    typedef Status_line_CA_1<Algebraic_curve_kernel_2, Rep> Self;

    //! type of x-coordinate
    typedef typename Curve_analysis_2::Algebraic_real_1 Algebraic_real_1;

    //! type of a curve point
    typedef typename Curve_analysis_2::Algebraic_real_2 Algebraic_real_2;

    //! an instance of a size type
    typedef typename Curve_analysis_2::size_type size_type;

    //! encodes number of arcs to the left and to the right
    typedef std::pair<size_type, size_type> Arc_pair;

    //! container of arcs
    typedef std::vector<Arc_pair> Arc_container;

    //! Local isolator type
    typedef typename Rep::Bitstream_descartes Bitstream_descartes;

     //! the handle superclass
    typedef ::CGAL::Handle_with_policy< Rep > Base;

    //!@}
public:
    //!\name constructors
    //!@{

    /*!\brief
     * Default constructor
     */
    Status_line_CA_1() :
        Base(Rep()) {
    }

    /*!\brief
     * copy constructor
     */
#ifdef DOXYGEN_RUNNING
    Status_line_CA_1(const Self& p) :
            Base(static_cast<const Base&>(p)) {
    }
#endif

    /*!\brief
     * constructs a status line over the \c i-th interval with x-coordinate
     * \c x
     *
     * \c n_arcs defines # of curve arcs over the interval
     *
     * \pre specified x-coordinate belongs to \c i-th interval
     */
    Status_line_CA_1(
        Algebraic_real_1 x, size_type i, const Curve_analysis_2& ca,
        size_type n_arcs) :
        Base(Rep(x, i, ca, n_arcs)) {

        CGAL_precondition(n_arcs >= 0);
        CGAL_precondition_code(
            bool is_event;
            size_type idx;
            ca.x_to_index(x, idx, is_event);
            CGAL_precondition(!is_event && idx == i);
        );
    }


    /*!\brief
     * constructs a status line at the event with x-coorinate \c x
     *
     * \c arcs container defines # of incident arcs to the left and to the
     * right of each intersection point of the curve \c ca with this status
     * line, sorted by y-s in ascending order.
     * \c n_arcs_left and \c n_arcs_right specify total number of incident
     * arcs to the left and to the right respectively.
     * \c has_v_line specifies whether the curve has a vertical line as a
     * component at this x-coordinate
     *
     * \pre there is a curve event at specified x-coordinate
     */
    Status_line_CA_1(Algebraic_real_1 x,
        size_type i, const Curve_analysis_2& ca,
        size_type n_arcs_left, size_type n_arcs_right,
        Arc_container arcs, bool has_v_line = false) :
        Base(Rep(x, i, ca, n_arcs_left, n_arcs_right)) {

        CGAL_precondition(n_arcs_left >= 0 && n_arcs_right >= 0);
        CGAL_precondition_code(
            bool is_event;
            size_type idx;
            ca.x_to_index(x, idx, is_event);
            CGAL_precondition(idx == i);
        );
        _set_arcs(arcs);
        if(has_v_line)
            _set_v_line();
    }

    /*!\brief
     * constructs a status line at the event with x-coorinate \c x
     *
     * arcs and vertical line flag can be set later
     */
    Status_line_CA_1(Algebraic_real_1 x,
        size_type i, const Curve_analysis_2& ca,
        size_type n_arcs_left, size_type n_arcs_right) :
        Base(Rep(x, i, ca, n_arcs_left, n_arcs_right)) {

        CGAL_precondition(n_arcs_left >= 0 && n_arcs_right >= 0);
        CGAL_precondition_code(
            bool is_event;
            size_type idx;
            ca.x_to_index(x, idx, is_event);
            CGAL_precondition(idx == i);
        );
    }

    /*!\brief
     * constructs from a given represenation
     */
    Status_line_CA_1(Rep rep) :
        Base(rep) {
    }

    //!@}
public:
    //!\name access functions
    //!@{

    /*! \brief
     *  returns the x-coordinate of the status line (always a finite value)
     */
    Algebraic_real_1 x() const {
        return this->ptr()->_m_x;
    }

    /*! \brief
     *  returns this status line CurveAnalysis_2 object
     */
    Curve_analysis_2 curve_analysis_2() const {
        return this->ptr()->_m_ca;
    }

    /*! \brief
     *  returns this status line index (event or interval index)
     */
    size_type index() const {
        return this->ptr()->_m_index;
    }

    /*! \brief
     *  returns \c true in case the given curve contains the status line
     *  as a component
     */
    bool covers_line() const {
        return this->ptr()->_m_vertical_line;
    }

    /*!\brief
     *  returns \c true if the curve f has intersection with f_y at x
     */
    bool has_f_fy_intersection() const {
        return (is_event() && !covers_line());
    }

    /*! \brief
     *  returns \c true if one of \c covers_line of \c has_f_fy_intersection
     *  evaluates to \c true
     */
    bool is_event() const {
        return this->ptr()->_m_event;
    }

    /*! \brief
     * returns \c true if the ith point is an event point
     */
    bool is_event(size_type i) const {
        // TODO: Make it possible to detect singularities as well
        Arc_pair branches = number_of_incident_branches(i);
        return branches.first!=1 || branches.second!=1;
    }


    /*! \brief
     * returns number of distinct and finite intersections of a curve with a
     * (intended) vertical line ignoring a real vertical line component of the
     * curve at the given x-coordinate.
     */
    size_type number_of_events() const {
        return this->ptr()->_m_total_arcs;
    }

    /*!\brief
     *  returns \c Algebraic_real_2 for j-th event over this vert line
     *
     * \pre 0 <= j < num_of_events()
     */
    Algebraic_real_2 algebraic_real_2(size_type j) const
    {
        CGAL_precondition(0 <= j&&j < number_of_events());
        if(!this->ptr()->_m_xy_coords[j])
          this->ptr()->_m_xy_coords[j] = Algebraic_real_2(x(),
                this->ptr()->_m_ca, j);
        return *(this->ptr()->_m_xy_coords[j]);
    }

    /*!\brief
     * alias for \c get_algebraic_real_2()
     */
    Algebraic_real_2 xy_coordinate_2(size_type j) const {
        return algebraic_real_2(j);
    }

    /*!\brief
     *  returns the number of branches of the curve connected to j-th
     *  event immediately to the left, to the right, respectively, as a pair of
     *  unsigned int ignoring vertical curve components at the given
     *  x-coordinate.
     *
     * \pre 0 <= j < num_of_events()
     */
    Arc_pair number_of_incident_branches(int j) const {

        CGAL_precondition(0 <= j&&j < number_of_events());
        if(!is_event())
            return Arc_pair(1, 1);

        return (*(this->ptr()->_m_arcs))[j];
    }

    /*! \brief
     * returns the number of vertical asymptotes at x of the curve
     * approaching y=-oo from left and right. A vertical line being component
     * of the curve is ignored.
     */
    const Arc_pair& number_of_branches_approaching_minus_infinity() const {
        return this->ptr()->_m_num_arcs_minus_inf;
    }

    /*! \brief
     *  returns the number of vertical asymptotes at x of the curve
     *  approaching y=+oo from left and right. A vertical line being component
     *  of the curve is ignored.
     */
    const Arc_pair& number_of_branches_approaching_plus_infinity() const {
        return this->ptr()->_m_num_arcs_plus_inf;
    }
protected:
    Algebraic_curve_kernel_2* kernel() const {
        return this->ptr()->_m_kernel;
    }


    //!@}
public:
    //!@{

    /*!\brief
     * sets # of arcs running from left and right to -inf and +inf at vertical
     * asymptote
     */
    void _set_number_of_branches_approaching_infinity(
            const Arc_pair& minus_inf, const Arc_pair& plus_inf) {

        CGAL_precondition(minus_inf.first >= 0 && minus_inf.second >= 0);
        CGAL_precondition(plus_inf.first >= 0 && plus_inf.second >= 0);

        this->ptr()->_m_num_arcs_minus_inf = minus_inf;
        this->ptr()->_m_num_arcs_plus_inf = plus_inf;

        if(!this->ptr()->_m_event)
            this->ptr()->_m_event = (minus_inf.first + minus_inf.second +
                plus_inf.first + plus_inf.second > 0);
    }

    void _set_arcs(const Arc_container& arcs) const {
        CGAL_precondition(is_event());
        this->ptr()->_m_arcs = arcs;
        this->ptr()->_m_total_arcs = static_cast<int>(arcs.size());
        this->ptr()->_m_xy_coords.resize(arcs.size());
    }

    void _set_v_line() const {
        CGAL_precondition(is_event());
        this->ptr()->_m_vertical_line = true;
    }

    //!@}
    //!\name IO
    //!@{

    void write(std::ostream& os) const {

        os << "status_line [CA@" << this->ptr()->_m_ca.id() << std::flush;
#if CGAL_ACK_USE_EXACUS
        os << "; x = " << x() << "; #events: " << number_of_events() << "; "
           << std::flush;
#else
        os << "; x = " << CGAL::to_double(x()) << "; #events: "
           << number_of_events() << "; " << std::flush;
#endif


        if(is_event()) {
            os << "incident branches: {" << std::flush;
//            typename Arc_container::const_iterator ait =
//                (*this->ptr()->_m_arcs).begin();
            for(int i = 0; i < number_of_events(); i++) {

                Arc_pair arc_pair = number_of_incident_branches(i);

                if(i!=0) {
                    os << ", " << std::flush;
                }
                Algebraic_real_2 xy = algebraic_real_2(i);
                typedef typename Bitstream_descartes::Bound Bound;
                Bound th = CGAL::ipower(Bound(1,2),53);
                std::pair<double,double> d_pair
                   = xy.to_double();
                os << "y=" << d_pair.second << ", " << std::flush;
                os << "(" << arc_pair.first << ", " << arc_pair.second << ")"
                   << std::flush;
            }
            os << "}";
            Arc_pair p = number_of_branches_approaching_minus_infinity();
            if(p.first + p.second > 0)
                os << "; approaching -oo: (" << p.first << "; " <<
                    p.second << ")" << std::flush;
            p = number_of_branches_approaching_plus_infinity();
            if(p.first + p.second > 0)
                os << "; approaching +oo: (" << p.first << "; " <<
                    p.second << ")" << std::flush;
            if(covers_line())
                os << "; covers line" << std::flush;
        } else
            os << "interval line" << std::flush;

        os << "]" << std::flush;
    }

    //!@}

    //! Sets the isolator instance
    void set_isolator (const Bitstream_descartes& isolator) const {
        this->ptr()->isolator = isolator;
    }

    //! Returns the isolator instance
    Bitstream_descartes& isolator() const {
        CGAL_assertion(bool(this->ptr()->isolator));
        return this->ptr()->isolator.get();
    }

    //! Returns whether an isolator has been given for that status line
    bool has_isolator() const {
        return this->ptr()->isolator;
    }

    typename Bitstream_descartes::Bound lower_bound(int index) const {
        return isolator().left_bound(index);
    }

    typename Bitstream_descartes::Bound upper_bound(int index) const {
        return isolator().right_bound(index);
    }

    typename Bitstream_descartes::Bound interval_length(int index) const {
        return isolator().right_bound(index)-
               isolator().left_bound(index);
    }

    int get_upper_bound_for_multiplicity(int index) const {
        return isolator().get_upper_bound_for_multiplicity(index);
    }

    void refine(int index) const {
        return isolator().refine_interval(index);
    }

    void refine_to(int index, typename Bitstream_descartes::Bound b) {
            while(upper_bound(index) - lower_bound(index) > b) {
                refine(index);
            }
    }


    //! these are our friends
    //friend class Curve_analysis_2;

}; // class Status_line_CA_1

template <class CurveAnalysis_2, class Rep>
std::ostream& operator<< (
        std::ostream& os,
        const internal::Status_line_CA_1<CurveAnalysis_2, Rep>& line) {

    line.write(os);
    return os;
}

} // namespace internal

} //namespace CGAL

#endif // CGAL_ALGEBRAIC_CURVE_KERNEL_STATUS_LINE_CA_1_H

