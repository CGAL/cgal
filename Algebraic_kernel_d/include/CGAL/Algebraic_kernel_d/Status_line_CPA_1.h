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

#ifndef CGAL_ALGEBRAIC_CURVE_KERNEL_STATUS_LINE_CPA_1_H
#define CGAL_ALGEBRAIC_CURVE_KERNEL_STATUS_LINE_CPA_1_H

#include <CGAL/basic.h>
#include <CGAL/Handle_with_policy.h>
#include <CGAL/algorithm.h>

namespace CGAL {

namespace internal {

template < class CurvePairAnalysis_2, class Rep_ >
class Status_line_CPA_1;

template <class CurvePairAnalysis_2, class Rep>
std::ostream& operator<< (std::ostream&,
    const Status_line_CPA_1<CurvePairAnalysis_2, Rep>&);

template < class CurvePairAnalysis_2 >
class Status_line_CPA_1_rep {

    // this template argument
    typedef CurvePairAnalysis_2 Curve_pair_analysis_2;

    // myself
    typedef Status_line_CPA_1_rep<Curve_pair_analysis_2> Self;

    // type of x-coordinate
    typedef typename Curve_pair_analysis_2::Algebraic_real_1
                Algebraic_real_1;

    // type of a curve point
    typedef typename Curve_pair_analysis_2::Algebraic_real_2
                Algebraic_real_2;

    // an instance of a size type
    typedef typename Curve_pair_analysis_2::size_type size_type;

     // encodes number of arcs to the left and to the right
    typedef std::pair<size_type, size_type> Arc_pair;

    // container of arcs
    typedef std::vector<Arc_pair> Arc_container;

    // container of integers ?
    typedef std::vector<size_type> Int_container;

    // constructors
public:
    // default constructor ()
    Status_line_CPA_1_rep()
    {   }

    // constructs an empty status line object
    Status_line_CPA_1_rep(size_type i, Curve_pair_analysis_2 cpa) :
        _m_index(i), _m_cpa(cpa), _m_event(false), _m_intersection(false) {
    }

    // stores this status line interval or event index of a curve pair
    size_type _m_index;

    // represents x-coordinate of event of rational value over interval
    // computed only by demand
    mutable boost::optional<Algebraic_real_1> _m_x;

    // for each event point stores a pair of arcnos of the 1st and 2nd curve
    // or -1 if respective curve is not involved
    mutable Arc_container _m_arcs;

    // inverse mapping from arcnos of the 1st and 2nd curve to respective
    // y-position
    mutable Int_container _m_arcno_to_pos[2];

    // stores multiplicities of intersection points (-1 if there is no 2-curve
    // intersection)
    mutable Int_container _m_mults;

    // underlying curve pair analysis
    Curve_pair_analysis_2 _m_cpa;

    // is there an event
    mutable bool _m_event;

    // is there is an intersection of both curves
    mutable bool _m_intersection;

    // befriending the handle
    friend class Status_line_CPA_1<Curve_pair_analysis_2, Self>;
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
      class Rep_ = internal::Status_line_CPA_1_rep<CurvePairAnalysis_2> >
class Status_line_CPA_1 :
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
    typedef Status_line_CPA_1<Curve_pair_analysis_2, Rep> Self;

    //! type of x-coordinate
    typedef typename Curve_pair_analysis_2::Algebraic_real_1 Algebraic_real_1;

    //! type of a curve point
    typedef typename Curve_pair_analysis_2::Algebraic_real_2 Algebraic_real_2;

    //! an instance of a size type
    typedef typename Curve_pair_analysis_2::size_type size_type;

    //! encodes number of arcs to the left and to the right
    typedef std::pair<size_type, size_type> Arc_pair;

    //! container of arcs
    typedef std::vector<Arc_pair> Arc_container;

    //! container of integers ?
    typedef std::vector<size_type> Int_container;

     //! the handle superclass
    typedef ::CGAL::Handle_with_policy< Rep > Base;

    //!@}
public:
    //!\name constructors
    //!@{

    /*!\brief
     * Default constructor
     */
    Status_line_CPA_1() :
        Base(Rep()) {
    }

    /*!\brief
     * copy constructor
     */
#ifdef DOXYGEN_RUNNING
    Status_line_CPA_1(const Self& p) :
            Base(static_cast<const Base&>(p)) {
    }
#endif

    /*!\brief
     * constructs undefined status line
     */
    Status_line_CPA_1(size_type i, Curve_pair_analysis_2 cpa) :
        Base(Rep(i, cpa)) {
    }

    /*!\brief
     * constructs a status line at the \c i -th event of a curve pair
     *
     * each element of \c arcs is a pair with the first item specifying the
     * type of event (0 - event of the 1st curve, 1 - of the second curve,
     * 2 - of both curves), and the second item - multiplicity of intersection
     * or -1 if not available
     */
    Status_line_CPA_1(size_type i, const Arc_container& arcs,
            Curve_pair_analysis_2 cpa) :
        Base(Rep(i, cpa)) {
        _set_event_arcs(arcs);
    }

     /*!\brief
     * constructs a status line over the \c i -th interval of a curve pair
     *
     * each element of \c arcs specifies to which curve a respective arc
     * belongs to (0 - arc of the 1st curve, 1 - arc of the 2nd curve)
     * \c is_swapped defines that the curves in targeting curve pair analysis
     * were swapped during precaching
     */
    Status_line_CPA_1(size_type i, const Int_container& arcs,
            Curve_pair_analysis_2 cpa) :
        Base(Rep(i, cpa)) {
        _set_interval_arcs(arcs);
    }

protected:
    /*!\brief
     * constructs from a given represenation
     */
    Status_line_CPA_1(Rep rep) :
        Base(rep) {
    }

    //!@}
public:
    //!\name access functions
    //!@{

    /*! \brief
     * returns the x-coordinate of the vertical line (always a finite value).
     */
    Algebraic_real_1 x() const {
        // unless x-coordiate was explicitly set with _set_x: compute its value
        if(!this->ptr()->_m_x) {
            this->ptr()->_m_x = (is_event() ?
#if CGAL_ACK_USE_EXACUS
                this->ptr()->_m_cpa._internal_curve_pair().event_x(index()) :
                Algebraic_real_1(this->ptr()->_m_cpa._internal_curve_pair().
                               bound_value_in_interval(index())));
#else
                this->ptr()->_m_cpa.event_x(index()) :
                Algebraic_real_1(this->ptr()->_m_cpa.
                    bound_value_in_interval(index())));
#endif
        }
        return *(this->ptr()->_m_x);
    }

    //! returns this vertical line's index (event or interval index)
    size_type index() const {
        CGAL_precondition(this->ptr()->_m_index>=0);
        return this->ptr()->_m_index;
    }

    /*! \brief
     *  returns number of distinct and finite intersections of a pair
     *  of curves  with a (intended) vertical line ignoring a real vertical
     *  line component of the curve at the given x-coordinate.
     */
    size_type number_of_events() const {
        return static_cast<size_type>(this->ptr()->_m_arcs.size());
    }


    /*! \brief
     *  returns the y-position of the k-th event of
     *  the curve in the sequence of events.
     *
     * Note that each event is formed by the 1st, 2nd, or both curves
     *
     * \pre 0 <= k < "number of arcs defined for curve c at x()"
     */
    size_type event_of_curve(size_type k,
                             const typename Curve_pair_analysis_2
                                 ::Curve_analysis_2& c) const {
        CGAL_assertion(c.id()==this->ptr()->_m_cpa.curve_analysis(false).id()||
                       c.id()==this->ptr()->_m_cpa.curve_analysis(true).id());
        bool b = (c.id()==this->ptr()->_m_cpa.curve_analysis(true).id());
        return event_of_curve(k,b);
    }




    /*! \brief
     *  returns the y-position of the k-th event of the c-th (0 or 1)
     * curve in the sequence of events.
     *
     * Note that each event is formed by the 1st, 2nd, or both curves
     *
     * \pre 0 <= k < "number of arcs defined for curve[c] at x()"
     */
    size_type event_of_curve(size_type k, bool c) const {

        CGAL_precondition_msg(0 <= k &&
            k < static_cast<size_type>(this->ptr()->_m_arcno_to_pos[c].size()),
                "Invalid arc number of the c-th curve specified");
        return this->ptr()->_m_arcno_to_pos[c][k];
    }

    /*! \brief
     *  returns the multiplicity of intersection defined at event with
     * position \c j. May return -1 in case multiplicity is unknown.
     *
     * \pre There is an intersection of both curves at j-th event
     * \pre 0 <= j < number_of_events()
     */
    size_type multiplicity_of_intersection(size_type j) const
    {
        CGAL_precondition(0 <= j && j < number_of_events());
        CGAL_precondition(is_intersection());
        CGAL_precondition(this->ptr()->_m_arcs[j].first != -1 &&
                          this->ptr()->_m_arcs[j].second != -1);

        return this->ptr()->_m_mults[j];
    }

    /*! \brief
     * returns a pair of \c int indicating whether event \c j is formed
     * by which arc numbers of the first and the second curve, or -1, if the
     * corresponding curve is not involved.
     *
     * \pre 0 <= j < number_of_events()
     */
    Arc_pair curves_at_event(size_type j) const
    {
        CGAL_precondition(0 <= j && j < number_of_events());
        const Arc_pair& arc = this->ptr()->_m_arcs[j];
        return arc;
    }

    /*!
     * returns an index indicating whether event \c j is formed
     * by which arc numbers of the curve \c ca, or -1, if the
     * corresponding curve is not involved.
     */
    Arc_pair curves_at_event(size_type j,
                             const typename Curve_pair_analysis_2
                                 ::Curve_analysis_2& c1,
                             const typename Curve_pair_analysis_2
                                 ::Curve_analysis_2& CGAL_precondition_code(c2)) const
    {

        CGAL_precondition(0 <= j && j < number_of_events());
        CGAL_assertion(c1.id()!=c2.id());
        CGAL_assertion
            (c1.id()==this->ptr()->_m_cpa.curve_analysis(false).id()||
             c1.id()==this->ptr()->_m_cpa.curve_analysis(true).id());
        CGAL_assertion
            (c2.id()==this->ptr()->_m_cpa.curve_analysis(false).id()||
             c2.id()==this->ptr()->_m_cpa.curve_analysis(true).id());
        bool b = (c1.id()==this->ptr()->_m_cpa.curve_analysis(false).id());
        const Arc_pair& arc_pair = curves_at_event(j);
        return b ? arc_pair : std::make_pair(arc_pair.second, arc_pair.first);
    }

    /*! \brief
     *  returns true if a curve has an event or in case there is an
     *  intersection of both curves.
     */
    bool is_event() const {
        return this->ptr()->_m_event;
    }

    /*! \brief
     * returns true if there is an intersection of both curves.
     */
    bool is_intersection() const {
        return this->ptr()->_m_intersection;
    }

    //!@}
public:
    //!@{

    /*!\brief
     * sets x-coordinate of a status line
     */
    void _set_x(Algebraic_real_1 x) const {
        this->ptr()->_m_x = x;
    }

    /*!\brief
     * sets arcs at event (use at your own risk!)
     */
    void _set_event_arcs(const Arc_container& arcs) const {

        size_type k = 0, arcf = 0, arcg = 0;
        this->ptr()->_m_arcs.resize(arcs.size());
        this->ptr()->_m_mults.resize(arcs.size());
        this->ptr()->_m_event = true;

        for(typename Arc_container::const_iterator ait = arcs.begin();
                ait != arcs.end(); ait++, k++) {

            if(ait->first == 0) { // 1st curve
                this->ptr()->_m_arcs[k].first = arcf++;
                this->ptr()->_m_arcs[k].second = -1;
                this->ptr()->_m_arcno_to_pos[0].push_back(k);

            } else if(ait->first == 1) { // 2nd curve
                this->ptr()->_m_arcs[k].first = -1;
                this->ptr()->_m_arcs[k].second = arcg++;
                this->ptr()->_m_arcno_to_pos[1].push_back(k);

            } else if(ait->first == 2) { // intersection
                this->ptr()->_m_arcs[k].first = arcf++;
                this->ptr()->_m_arcs[k].second = arcg++;
                this->ptr()->_m_arcno_to_pos[0].push_back(k);
                this->ptr()->_m_arcno_to_pos[1].push_back(k);
                this->ptr()->_m_intersection = true;

            } else
                CGAL_error_msg("Bogus curve index..");
            this->ptr()->_m_mults[k] = ait->second;
        }
    }

    /*!\brief
     * sets arcs over interval (use at your own risk!)
     */
    void _set_interval_arcs(const Int_container& arcs) const {

        this->ptr()->_m_arcs.resize(arcs.size());
        this->ptr()->_m_event = false;
        this->ptr()->_m_intersection = false;

        size_type k = 0, arcf = 0, arcg = 0;
        for(typename Int_container::const_iterator ait = arcs.begin();
                ait != arcs.end(); ait++, k++) {
            if(*ait == 0) { // 1st curve
                this->ptr()->_m_arcs[k].first = arcf++;
                this->ptr()->_m_arcs[k].second = -1;
                this->ptr()->_m_arcno_to_pos[0].push_back(k);

            } else if(*ait == 1) { // 2nd curve
                this->ptr()->_m_arcs[k].first = -1;
                this->ptr()->_m_arcs[k].second = arcg++;
                this->ptr()->_m_arcno_to_pos[1].push_back(k);

            } else
                CGAL_error_msg("Bogus curve index..");
        }
    }

    //!@}
public:
    //!\name IO
    //!@{

    void write(std::ostream& os) const {

        os << "status_line [CPA@" << this->ptr()->_m_cpa.id();

        os << "; x = " << (index()==-1 ? 999.999 : CGAL::to_double(x())) << "; #events: " << number_of_events() << "; ";

        typename Arc_container::const_iterator ait =
                this->ptr()->_m_arcs.begin();
        if(is_event())
            os << "arcs at event: {";
        else
            os << "arcs of interval: {";

        for(; ait != this->ptr()->_m_arcs.end(); ait++) {
            if(ait != this->ptr()->_m_arcs.begin())
                os << ", ";
            os << "(" << ait->first << "; " << ait->second << ")";
        }
        os << "}, arcno2pos: (";

        CGAL::output_range(os, this->ptr()->_m_arcno_to_pos[0].begin(),
                this->ptr()->_m_arcno_to_pos[0].end(), ",");
        os << "), (";
        CGAL::output_range(os, this->ptr()->_m_arcno_to_pos[1].begin(),
                this->ptr()->_m_arcno_to_pos[1].end(), ",");
        os << ")]";
    }

    //!@}
}; // class Status_line_CPA_1

template <class CurvePairAnalysis_2, class Rep>
std::ostream& operator<< (std::ostream& os,
        const internal::Status_line_CPA_1<CurvePairAnalysis_2, Rep>& sline) {

    sline.write(os);
    return os;
}

} // namespace internal

} //namespace CGAL

#endif // CGAL_ALGEBRAIC_CURVE_KERNEL_STATUS_LINE_CPA_1_H
