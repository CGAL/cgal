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

#ifndef CGAL_CURVED_KERNEL_CURVE_INTERVAL_ARCNO_CACHE_H
#define CGAL_CURVED_KERNEL_CURVE_INTERVAL_ARCNO_CACHE_H

/*! \file Curved_kernel_via_analysis_2/Curve_interval_arcno_cache.h
 *  defines \c Curve_interval_arcno_cache functor
 */

#include <CGAL/basic.h>
#include <CGAL/Handle_with_policy.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

/*!\brief 
 * given an arc number over an interval \c i this functor computes a
 * corresponding event arcno (and asymptotic tendency) over certain vertical
 * line which lies on the given side w.r.t. the interval \c i. 
 *
 * For caching issues for each accessed \c Curve_2 object we store
 * \c Interval_arcno_map structure. In its turn, \c Interval_arcno_map
 * for each \c Status_line_1 object stores the precomputed mapping
 * from interval arcnos (on left and right sides) to event arcnos 
 */
template <class CurvedKernelViaAnalysis_2>
struct Curve_interval_arcno_cache {

    //!\name public typedefs            
    //!@{

    //! this instance's first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;
    
    //! type of generic curve
    typedef typename Curved_kernel_via_analysis_2::Curve_2 Curve_2;
    
    //! type of 1-curve analysis
    typedef typename Curved_kernel_via_analysis_2::Curve_analysis_2
        Curve_analysis_2;
    
    //! type of x-coordinate
    typedef typename Curved_kernel_via_analysis_2::X_coordinate_1
        X_coordinate_1;
    
    //! type of status line
    typedef typename Curve_analysis_2::Status_line_1
        Status_line_1;
    
    //! type of point on curve
    typedef typename Curved_kernel_via_analysis_2::Point_2 Point_2;
    
    //! type of curve arc
    typedef typename Curved_kernel_via_analysis_2::Arc_2 Arc_2;
    
    //! event arc number descriptor: stores an arc number along with curve's
    //! end type (+/-oo or \c NO_BOUNDARY )
    typedef std::pair<int, CGAL::Boundary_type> Arcno_desc;
    
    typedef Status_line_1 first_argument_type;
    typedef bool                  second_argument_type;
    typedef int                   thirg_argument_type;
    
    typedef Arcno_desc            result_type;
    typedef Arity_tag<3> Arity;   
    
    //!@}
    //!\name standard constuctor functor invokation
    //!@{

    //! standard constructor
    Curve_interval_arcno_cache(Curved_kernel_via_analysis_2 *kernel) :
        _m_curved_kernel_2(kernel), _m_last_curve_id(-1) {
        
        CGAL_assertion(kernel != NULL);
    }
    
    //!\brief given arcno over an interval, this computes a corresponding event
    //! arcno on status line \c cv_line lying on the given \c side w.r.t. the
    //! interval
    //!
    //! \c side = 0: left side; \c side = 1: right side
    result_type operator()(const Status_line_1& cv_line, bool side, 
        int interval_arcno) const {
        
        CGAL_precondition(interval_arcno >= 0);
        if(!cv_line.is_event()) // # of arcs over interval is constant
            return std::make_pair(interval_arcno, CGAL::NO_BOUNDARY);
        
        Curve_analysis_2 ca_2 = cv_line.curve_analysis_2();
        int curve_id = ca_2.curve_2().id();
        if(_m_last_curve_id != curve_id) {
            typename Curve_to_interval_arcno_map::iterator it;
            if(_m_last_curve_id != -1) { 
            // commit the last changes to the global curve arcno map
                it = _m_curve_arcno_map.find(_m_last_curve_id);
                CGAL_precondition(it != _m_curve_arcno_map.end());
                it->second = _m_last_interval_map; 
            }        
            // modify _m_last_interval_map to be pointing to the
            // Interval_arcno_map  corresponding to the given curve_id
            it = _m_curve_arcno_map.insert(std::make_pair(curve_id,
                    Interval_arcno_map())).first;
            _m_last_interval_map = it->second;
            _m_last_curve_id = curve_id;
        }
        // cache status lines by id(); shall we use x-coordinate instead of
        // id to ensure uniqueness ?
        typename Interval_arcno_map::const_iterator it = 
            _m_last_interval_map.find(cv_line.x().id());
        if(it != _m_last_interval_map.end()) {
            const Arcno_vector_pair& tmp = it->second;
            // note: side index must be reversed 
            CGAL_precondition(interval_arcno < static_cast<int>(side == 1 ?
                tmp.first.size() : tmp.second.size()));
            return (side == 1 ? tmp.first[interval_arcno] :
                tmp.second[interval_arcno]);
        }
        
        int n_left = ca_2.status_line_of_interval(
                cv_line.index()).number_of_events(), // # arcs to the left
            n_right = ca_2.status_line_of_interval( // # arcs to the right
                cv_line.index()+1).number_of_events();
        std::pair<int, int> n_mininf, n_maxinf, ipair;
        n_mininf = cv_line.number_of_branches_approaching_minus_infinity();
        n_maxinf = cv_line.number_of_branches_approaching_plus_infinity();
        // we must also account for the number of asymptotic arcs approaching
        // this status line
        Arcno_vector_pair vpair;
        vpair.first.resize(n_left); //+ n_mininf.first + n_maxinf.first);
        vpair.second.resize(n_right); //+ n_mininf.second + n_maxinf.second);
        
        int i, _arcno, left_i = 0, right_i = 0;
        // process arcs approaching -oo from left and right
        for(; left_i < n_mininf.first; left_i++)
            vpair.first[left_i] = std::make_pair(left_i, CGAL::MINUS_INFINITY);
            
        for(; right_i < n_mininf.second; right_i++)
            vpair.second[right_i] = std::make_pair(right_i,
                CGAL::MINUS_INFINITY);
                
        // process all finite events over this vertical line
        for(_arcno = 0; _arcno < cv_line.number_of_events(); _arcno++) {
            ipair = cv_line.number_of_incident_branches(_arcno);
            for(i = 0; i < ipair.first; i++)
                vpair.first[left_i++] = std::make_pair(_arcno,
                    CGAL::NO_BOUNDARY);
                    
            for(i = 0; i < ipair.second; i++)
            vpair.second[right_i++] = std::make_pair(_arcno,
                CGAL::NO_BOUNDARY);
        }
        // process arcs approaching +oo from left and right
        // this is not clear.. why arcnos at +oo are mapped to the same
        // interval arcnos ?
        for(i = 0; i < n_maxinf.first; i++, left_i++)
            vpair.first[left_i] = std::make_pair(left_i, CGAL::PLUS_INFINITY); 
                
        for(i = 0; i < n_maxinf.second; i++, right_i++)
            vpair.second[right_i] = std::make_pair(right_i,
                CGAL::PLUS_INFINITY);
        
        CGAL_precondition(left_i == static_cast<int>(vpair.first.size()) && 
            right_i == static_cast<int>(vpair.second.size()));
        _m_last_interval_map.insert(std::make_pair(cv_line.x().id(), vpair));
        CGAL_precondition(interval_arcno < static_cast<int>(side == 1 ?
                vpair.first.size() : vpair.second.size()));
        return (side == 1 ? vpair.first[interval_arcno] :
               vpair.second[interval_arcno]);
    }

    //!@}
private: 
    //!\name private members
    //!@{

    //! pointer to \c Curved_kernel_via_analysis_2 
    Curved_kernel_via_analysis_2 *_m_curved_kernel_2;
    
    //! a pair of vectors (right and left side of event-line respectively)
    typedef std::pair<std::vector<Arcno_desc>, std::vector<Arcno_desc> >
        Arcno_vector_pair;
        
    //! maps from \c Status_line_1 id to a container of interval ->
    //! event arcnos 
    // TODO: use hash instead ?? cache status lines by x-coordinate ?
    typedef std::map<int, Arcno_vector_pair> Interval_arcno_map;
    
    //! maps from \c Curve_2 id to \c Interval_arcno_map
    typedef std::map<int, Interval_arcno_map> Curve_to_interval_arcno_map;
    
    //! arcno mapping container instance
    mutable Curve_to_interval_arcno_map _m_curve_arcno_map;
    
    //! an id of the last queried \c Curve_2 object
    mutable int _m_last_curve_id;
    
    //! stores an address of the last accessed map 
    mutable Interval_arcno_map _m_last_interval_map;
    
    //!@}
};  // struct Curve_interval_arcno_cache

} // namespace CGALi

CGAL_END_NAMESPACE

#endif // CGAL_CURVED_KERNEL_CURVE_INTERVAL_ARCNO_CACHE_H
