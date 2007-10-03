// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//                 Arno Eigenwillig <arno@mpi-inf.mpg.de>
//                 Pavel Emeliyanenko <asm@mpi-sb.mpg.de> 
//
// ============================================================================

#ifndef CGAL_GPA_MAKE_X_MONOTONE_H
#define CGAL_GPA_MAKE_X_MONOTONE_H

/*! \file GPA_2/Make_x_monotone.h
 *  defines \c Make_x_monotone functor
 */

#include <CGAL/basic.h>
#include <CGAL/Handle_with_policy.h>

//#include <CGAL/GPA_2/Curve_interval_arcno_cache.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

/*!\brief 
 * splits a complete curve into x-monotone sweepable arcs
 *
 * The given curve is split into sweepable arcs by cutting it into
 * connected, either x-monotone pieces of constant interior arc number 
 * at every event x-coordinate or vertical (\c GPA_2::Arc_2 objects). Isolated
 * points are stored as \c GPA_2::Point_2 objects.
 * 
 * The resulting arcs and points are written to the output iterator as
 * polymorph \c CGAL::Object. Past-the-end value of the iterator is returned.
 *
 * Arcs extending to infinity get an endpoint with x-coordinate
 * plus/minus infinity, as appropriate. Note that infinity here is
 * understood in an "infimaximal" sense: smaller/larger than any
 * event, modelling the behaviour for all sufficiently small/large
 * x-coordinates. The same for arcs running up/down a vertical asymtote
 * of the curve.
 */
template <class GPA_2>
struct Make_x_monotone_2 :
    public Binary_function< typename GPA_2::Curve_2, 
            std::iterator<output_iterator_tag, CGAL::Object>,
            std::iterator<output_iterator_tag, CGAL::Object> > {
            
    //!\name public typedefs            
    //!@{
    
    //! type of x-coordinate
    typedef typename GPA_2::X_coordinate_1 X_coordinate_1;
    
    //! type of a finite point on curve
    typedef typename GPA_2::Xy_coordinate_2 Xy_coordinate_2;
    
    //! type of generic curve
    typedef typename GPA_2::Curve_2 Curve_2;
    
    //! type of 1-curve analysis
    typedef typename GPA_2::Curve_analysis_2 Curve_analysis_2;
    
    //! type of vertical line
    typedef typename Curve_analysis_2::Curve_vertical_line_1
        Curve_vertical_line_1;
    
    //! type of point on curve
    typedef typename GPA_2::Point_2 Point_2;
    
    //! type of curve arc
    typedef typename GPA_2::Arc_2 Arc_2;
    
    //!@}
    //!\name standard constuctor functor invokation
    //!@{

    //! standard constructor
    Make_x_monotone_2(GPA_2 *gpa) :
        _m_gpa(gpa) {
        CGAL_assertion(gpa != NULL);
    }
    
    //! function-call operator
    template <class OutputIterator>
    OutputIterator operator()(Curve_2 curve, OutputIterator oi) {

        if(NiX::total_degree(curve.f()) < 1) // use CGAL::Total_degree ? 
            return oi;
        
        std::cout << "processing curve: " << curve.f() << "\n";    
            
        Curve_analysis_2 ca_2(curve);
        Curve_vertical_line_1 evt_line1, evt_line2, 
            int_line = ca_2.vertical_line_of_interval(0);
        int total_events = ca_2.number_of_vertical_lines_with_event();
        // handle special case of a curve without any events
        if(total_events == 0) {
            for(int k = 0; k < int_line.number_of_events(); k++) 
                *oi++ = CGAL::make_object(Arc_2(curve, k));
            return oi;
        }
        _m_curve = curve;
        const typename GPA_2::Curve_interval_arcno_cache& map_interval_arcno = 
            _m_gpa->get_interval_arcno_cache();
        
        std::pair<int, CGAL::Boundary_type> info1, info2;
        std::vector<Point_2> min_pts, max_pts;
        X_coordinate_1 min_x, max_x;
        int i, j, k, n;
        Arc_2 arc;
        // first handle segments before first event
        evt_line1 = ca_2.vertical_line_at_event(0);
        max_x = evt_line1.x();
        
        for(k = 0; k < evt_line1.number_of_events(); k++) 
            max_pts.push_back(Point_2(max_x, curve, k));
                        
        std::cout << "handling events over the 1st interval\n";    
        for(k = 0; k < int_line.number_of_events(); k++) {
            
            info1 = map_interval_arcno(evt_line1, 1, k); 
            std::cout << "k = " << k << "; evt arcno: " <<
                info1.first << "; tendency: " << info1.second << "\n";
            
            if(info1.second != CGAL::NO_BOUNDARY) 
                arc = Arc_2(CGAL::MIN_END, max_x, (info1.second < 0 ?
                    CGAL::MIN_END : CGAL::MAX_END), curve, k);
            else
                arc = Arc_2(max_pts[info1.first], CGAL::MIN_END, curve, k,
                    info1.first);
            *oi++ = CGAL::make_object(arc);
        }
        min_pts = max_pts;
        max_pts.clear();
        min_x = max_x;
               
        // next handle arcs between events, including isolated points
        for(i = 0; i < total_events-1; i++) {
            evt_line1 = ca_2.vertical_line_at_event(i);
            evt_line2 = ca_2.vertical_line_at_event(i+1);
            max_x = evt_line2.x();
            
            std::cout << "is_fixed flag: " << evt_line1.get_info().is_fixed()
                 << "\n";
            oi = _handle_vertical_and_isolated(evt_line1, min_x, min_pts, oi);
                                
            n = evt_line2.number_of_events();
            for(k = 0; k < n; k++) 
                max_pts.push_back(Point_2(max_x, curve, k));
            
            n = ca_2.vertical_line_of_interval(i+1).number_of_events();
            CGAL::Curve_end inf1_end, inf2_end;
            for(k = 0; k < n; k++) {
                std::cout << "getting first event arcno: k = " << k << "\n";
                info1 = map_interval_arcno(evt_line1, 0, k); 
                std::cout << "getting second event arcno: k = " << k << "\n";
                info2 = map_interval_arcno(evt_line2, 1, k); 
                inf2_end = (info2.second < 0 ? CGAL::MIN_END : 
                    CGAL::MAX_END);
                if(info1.second != CGAL::NO_BOUNDARY) {
                    inf1_end = (info1.second < 0 ? CGAL::MIN_END : 
                        CGAL::MAX_END);
                    if(info2.second != CGAL::NO_BOUNDARY) 
                        arc = Arc_2(min_x, max_x, inf1_end, inf2_end,curve,k);
                    else 
                        arc = Arc_2(max_pts[info2.first], min_x, inf1_end,
                            curve, k, info2.first);
                } else if(info2.second != CGAL::NO_BOUNDARY) {
                    arc = Arc_2(min_pts[info1.first],  max_x, inf2_end, curve,
                        k, info1.first);
                } else 
                    arc = Arc_2(min_pts[info1.first], max_pts[info2.first],
                        curve, k, info1.first, info2.first);
                *oi++ = CGAL::make_object(arc);
            }
            min_pts = max_pts;
            max_pts.clear();
            min_x = max_x;
        }
        
        // here: min_x/min_pts hold information about the last event line
        // event_line2 - points to the last event line
        // vertical line or isolated points at last event?
        evt_line2 = ca_2.vertical_line_at_event(total_events-1);
        min_x = evt_line2.x();
                
        oi = _handle_vertical_and_isolated(evt_line2, min_x, min_pts, oi);
        
        std::cout << "handling events over the last interval\n";
        n = ca_2.vertical_line_of_interval(total_events).number_of_events();
        for(k = 0; k < n; k++) {
        
            info1 = map_interval_arcno(evt_line2, 0, k); 
            std::cout << "k = " << k << "; evt arcno: " <<
                info1.first << "; tendency: " << info1.second << "\n";
            if(info1.second != CGAL::NO_BOUNDARY) 
                arc = Arc_2(CGAL::MAX_END, min_x, (info1.second < 0 ?
                    CGAL::MIN_END : CGAL::MAX_END), curve, k);
            else
                arc = Arc_2(min_pts[info1.first], CGAL::MAX_END, curve, k,
                    info1.first);
            *oi++ = CGAL::make_object(arc);
        }
        return oi;
    }
    
    //!@}
private:        
    //!\name private members
    //!@{
    
    template <class OutputIterator>
    OutputIterator _handle_vertical_and_isolated(Curve_vertical_line_1 cv_line,
        X_coordinate_1 x, std::vector<Point_2> pts, OutputIterator oi) const {
        
               
        int n = cv_line.number_of_events(), j;
        if(cv_line.covers_line()) { // look for vertical arcs
            if(n > 0) {
                // the first vertical ray
                *oi++ = CGAL::make_object(Arc_2(pts[0], CGAL::MIN_END,
                    _m_curve));
                for(j = 1; j < n-1; j++)  // interior bounded arcs
                    *oi++ = CGAL::make_object(Arc_2(pts[j], pts[j+1],
                        _m_curve));
                    // the last vertical ray
                    *oi++ = CGAL::make_object(Arc_2(pts[n-1], CGAL::MAX_END,
                        _m_curve));
            } else // unbounded vertical line
                *oi++ = CGAL::make_object(Arc_2(x, _m_curve));
            return oi;
        } 
        // look for isolated points
        std::pair<int, int> ipair;
        for(j = 0; j < n; j++) {
            ipair = cv_line.get_number_of_incident_branches(j);
        ///////////////////////////////////////////////////////////
        /// ATTENTION: isolated points lying on regular arcs are not handled!!
        ///////////////////////////////////////////////////////////
            if(ipair.first == 0&&ipair.second == 0) 
                *oi++ = CGAL::make_object(Point_2(x, _m_curve, j));
        }
        return oi;
    }
            
    //! pointer to \c GPA_2 
    GPA_2 *_m_gpa;
    //! to avoid passing curve as a parameter
    Curve_2 _m_curve;
    
    //!@}
};  // struct Make_x_monotone

/*!
 * \brief 
 * short-hand for \c Curve_to_segments<AlgebraicSegment2,...>(c,oi)
 */
/*template <class AlgebraicSegment2, class AlgebraicCurve2, class OutputIterator>
OutputIterator curve_to_segments(AlgebraicCurve2 c, OutputIterator oi) {
    Curve_to_segments<AlgebraicSegment2, AlgebraicCurve2, OutputIterator> 
        c2s;
    return c2s(c, oi);
}*/

} // namespace CGALi

CGAL_END_NAMESPACE

#endif // CGAL_GPA_MAKE_X_MONOTONE_H
