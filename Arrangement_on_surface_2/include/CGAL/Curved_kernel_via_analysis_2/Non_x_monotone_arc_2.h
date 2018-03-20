// Copyright (c) 2008,2009,2010,2011 Max-Planck-Institute Saarbruecken (Germany), 
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>

#ifndef CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_NON_X_MONOTONE_ARC_2_H
#define CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_NON_X_MONOTONE_ARC_2_H

/*!\file include/CGAL/Curved_kernel_via_analysis_2/Non_x_monotone_arc_2.h
 * \brief defines class \c Non_x_monotone_arc_2
 *  
 * non x-monotone arc of a generic curve
 */

#include <CGAL/config.h>
#include <CGAL/Handle_with_policy.h>

#include <CGAL/Curved_kernel_via_analysis_2/Arc_2.h>


namespace CGAL {

namespace internal {

template < class CurvedKernelViaAnalysis_2 >
class Non_x_monotone_arc_2_rep { 

public:

    // this instance's first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;
    
    // myself
    typedef Non_x_monotone_arc_2_rep< Curved_kernel_via_analysis_2 > Self;

    // type of an arc on generic curve
    typedef typename Curved_kernel_via_analysis_2::Arc_2 Arc_2;

    // a list of x-monotone arcs
    typedef std::vector< Arc_2 > Arc_vector;

public:    

    /*!\brief
     * Default constructor
     */
    Non_x_monotone_arc_2_rep() {
    }

    /*!\brief 
     * constructs from one x-monotone arc
     */
    Non_x_monotone_arc_2_rep(const Arc_2& arc) {
        _m_x_monotone_arcs.push_back(arc);
    }

    /*!\brief 
     * constructs a non x-monotone arc from the list of x-monotone pieces
     */
    template <class InputIterator>
    Non_x_monotone_arc_2_rep(InputIterator start, InputIterator end) :
        _m_x_monotone_arcs(std::distance(start, end)) {

        std::copy(start, end, _m_x_monotone_arcs.begin());
    }

    //! a list of x-monotone pieces this arc consists of
    mutable Arc_vector _m_x_monotone_arcs;

}; // Non_x_monotone_arc_2_rep

/*! \brief
 * Class representing a (not necessary x-monotone) curve arc
 *
 * This class represents a not necessarily x-monotone curve arc. The arc is
 * given as a list of connected x-monotone pieces.
 *
 * By constructing a new arc, its validity is checked to ensure that its
 * x-monotone pieces form a single chain
 */
template < class CurvedKernelViaAnalysis_2,
        class Rep_  = Non_x_monotone_arc_2_rep< CurvedKernelViaAnalysis_2 > >
class Non_x_monotone_arc_2 :
    public CGAL::Handle_with_policy< Rep_ > {
    
public:
    //!\name public typedefs
    //!@{
    
    //! this instance's first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! this instance's second template parameter
    typedef Rep_ Rep;

    //! this instance itself
    typedef Non_x_monotone_arc_2< Curved_kernel_via_analysis_2, Rep > Self;
    
    //! type of curve kernel
    typedef typename Curved_kernel_via_analysis_2::Curve_kernel_2
        Curve_kernel_2;
       
    //! type of analysis of a pair of curves
    typedef typename Curve_kernel_2::Curve_analysis_2 Curve_analysis_2;
    
    //! type of a point on generic curve
    typedef typename Curved_kernel_via_analysis_2::Point_2 Point_2;

    //! type of an x-monotone arc on generic curve
    typedef typename Curved_kernel_via_analysis_2::Arc_2 Arc_2;

    //! iterator type to range through the list of x-monotone arcs
    typedef typename Rep::Arc_vector::const_iterator Arc_const_iterator;
    
    //! the handle superclass
    typedef ::CGAL::Handle_with_policy< Rep > Base;

    //!@}
public:
    //!\name basic constructors
    //!@{

    /*!\brief
     * Default constructor
     */
    Non_x_monotone_arc_2() : 
        Base(Rep()) {   
    }

    /*!\brief
     * copy constructor
     */
    Non_x_monotone_arc_2(const Self& a) :
        Base(static_cast<const Base&>(a)) {  
    }

    /*! \brief
     * constructs an arc from one x-monotone piece
     */
    Non_x_monotone_arc_2(const Arc_2& arc) :
        Base(Rep(arc)) {
    }

    /*! \brief
     * constructs an arc from a list x-monotone pieces specified by
     * iterator range <tt>[start; end)</tt>
     *
     * template argument type of \c InputIterator is \c Arc_2
     * 
     * \pre the x-monotone arcs must be connected into a single chain
     * \pre either all x-monotone arcs must be vertical or non-vertical
     */
    template <class InputIterator>
    Non_x_monotone_arc_2(InputIterator start, InputIterator end) :
        Base(Rep(start, end)) {

        CGAL_precondition(start != end);
        CGAL_precondition_code(
            InputIterator oi = start;
            InputIterator next = ++oi;
            bool vertical = oi->is_vertical();
            Curve_analysis_2 curve = oi->curve();
            for(; next != end; next++) {
                // ensure that supporting curves are identical for all arcs
                CGAL_precondition(next->curve().is_identical(curve));
                // either all arcs are vertical or not
                CGAL_precondition(vertical == next->is_vertical());
                /* this condition is hard to check because end-points can be
                 reordered in Arc_2 constructor
                CGAL_precondition(oi->curve_end(CGAL::ARR_MAX_END).
                    compare_xy(next->curve_end(CGAL::ARR_MIN_END)) ==
                            CGAL::EQUAL);*/
                oi = next;
            }
        );
    }

    //!@}
public:
    //! \name Access functions
    //!@{
    
    /*!\brief
     * returns the number of x-monotone arcs this object consists of
     */
    int number_of_x_monotone_arcs() const {
        return static_cast<int>(this->ptr()->_m_x_monotone_arcs.size());
    }
    
    /*!\brief
     * returns iterator pointing to the first x-monotone arc in the list
     */
    Arc_const_iterator begin() const {
        return this->ptr()->_m_x_monotone_arcs.begin(); 
    }

    /*!\brief
     * returns iterator pointing to beyond the last x-monotone arc in the list
     */
    Arc_const_iterator end() const {
        return this->ptr()->_m_x_monotone_arcs.end();
    }

#if 1 // not needed
    /*!\brief
     * returns a distinct \c ith x-monotone piece of the arc
     */
    const Arc_2& x_monotone_arc(int i) const {
        CGAL_precondition(i >= 0);
        CGAL_precondition(i < number_of_x_monotone_arcs());
        return this->ptr()->_m_x_monotone_arcs[i];
    }
#endif	

    /*!\brief
     * returns the supporting curve
     */
    Curve_analysis_2 curve() const {
        CGAL_precondition(number_of_x_monotone_arcs() > 0);
        return this->ptr()->_m_x_monotone_arcs[0].curve();
    }

    /*!\brief
     * returns \c true if this arc consists of vertical segments
     */
    bool is_vertical() const {
        CGAL_precondition(number_of_x_monotone_arcs() > 0);
        return this->ptr()->_m_x_monotone_arcs[0].is_vertical();
    } 
    
    //!@}
}; // Non_x_monotone_arc_2

/*!\relates Non_x_monotone_arc_2
 * \brief 
 * output operator
 */
template < class CurvedKernelViaAnalysis_2, class Rep_>
std::ostream& operator <<(std::ostream& os,
    const Non_x_monotone_arc_2<CurvedKernelViaAnalysis_2, Rep_>& arc) {

    os << "List of x-monotone arcs: [\n";
    int i = 0;

    typename Non_x_monotone_arc_2<CurvedKernelViaAnalysis_2, Rep_>::
        Arc_const_iterator ait;
    for(ait = arc.begin(); ait != arc.end(); ait++, i++)
        os << i << ": " << *ait << std::endl;
    os << "]\n\n";
    return os;
}

#if DOXYGEN_RUNNING

/*!\relates Non_x_monotone_arc_2
 * \brief
 * draws an arc to Qt window \c qt
 */
template < class CurvedKernelViaAnalysis_2, class Rep_>
CGAL::Qt_widget& operator <<(CGAL::Qt_widget& ws,
        const Non_x_monotone_arc_2<CurvedKernelViaAnalysis_2, Rep_>& arc) {

    QPainter *ppnt = &ws.get_painter();

    int n_arcs = arc.number_of_x_monotone_arcs();

    QBrush defbrush = ppnt->brush();
    // forbid drawing connection points between x-monotone pieces
    if(n_arcs > 1)
        ppnt->setBrush(Qt::NoBrush);

    typename Non_x_monotone_arc_2<CurvedKernelViaAnalysis_2, Rep_>::
        Arc_const_iterator ait;
    for(ait = arc.begin(); ait != arc.end(); ait++)
        ws << *ait;

    if(n_arcs > 1) {
        ppnt->setBrush(defbrush);

        typeof(arc.x_monotone_arc(0)) const *xarc = &arc.x_monotone_arc(0);
        if(xarc->location(CGAL::ARR_MIN_END) == CGAL::ARR_INTERIOR)
            ws << xarc->curve_end(CGAL::ARR_MIN_END);

        xarc = &arc.x_monotone_arc(n_arcs-1);
        if(xarc->location(CGAL::ARR_MAX_END) == CGAL::ARR_INTERIOR)
            ws << xarc->curve_end(CGAL::ARR_MAX_END);
    }
    return ws;
};

#endif // DOXYGEN_RUNNING

} // namespace internal

} //namespace CGAL

#endif // CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_NON_X_MONOTONE_ARC_2_H
// EOF
