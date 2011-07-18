// Copyright (c) 2008,2009,2010,2011 Max-Planck-Institute Saarbruecken (Germany), 
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s): Eric Berberich       <eric@mpi-inf.mpg.de>

#ifndef CGAL_ARR_TOPOLOGY_TRAITS_SIGN_OF_PATH_H
#define CGAL_ARR_TOPOLOGY_TRAITS_SIGN_OF_PATH_H

/*!\file include/CGAL/Arr_topology_traits/Sign_of_path.h
 * \brief Definition of the Sign_of_path class template
 */

#include <CGAL/config.h>

namespace CGAL {

namespace internal {

template < class GeometryTraits_2, class TopologyTraits_2 >
class Sign_of_path {
    
public:
    //!\name Public types
    //!@{
    
    //! this instance's first template parameter
    typedef GeometryTraits_2 Geometry_traits_2;

    //! this instance's second template parameter
    typedef TopologyTraits_2 Topology_traits_2;
    
    //! type of point
    typedef typename Geometry_traits_2::Point_2 Point_2;
    
    //! type of x-monotone curve 
    typedef typename Geometry_traits_2::X_monotone_curve_2  
    X_monotone_curve_2;
    
    //! type of dcel
    typedef typename Topology_traits_2::Dcel Dcel;
    
    //! type of vertex
    typedef typename Dcel::Vertex Vertex;
    
    //! type of halfedge
    typedef typename Dcel::Halfedge Halfedge;

    //!@}

    //!\name Constructors
    //!@{

    //! Constructs instance from \c traits
    Sign_of_path(const Topology_traits_2 *traits) :
        _m_topology_traits(traits) {
        CGAL_precondition(traits != NULL);
    };

    //!@}

    //!\name function-call operators
    //!@{

    //! computes sign of closed path 
    CGAL::Sign operator()(const Halfedge* he1, const Halfedge* he2) const {
    
        // status: move to arr
        
#if CGAL_ARR_TOPOLOGY_TRAITS_VERBOSE
        std::cout << "Sign_of_path(he1->he2)?" << std::endl;
#endif

        CGAL::Sign result = CGAL::ZERO;
        
        if (he1->next() == he2 && he2->next() == he1) {
            return result;
        }
        
        // Start with the next of prev1:
        const Halfedge* curr = he1;
        const Halfedge* next = curr->next();
        
        CGAL_assertion(next != he1);
        
        do {
            
            if (!curr->has_null_curve() && !next->has_null_curve()) {
                
#if 0
      //           std::cout << "curr: " << curr->curve() << std::endl;
//                 std::cout << "currd: " << curr->direction() << std::endl;
//                 std::cout << "currf: " << (curr->is_on_inner_ccb() ?
//                                            curr->inner_ccb()->face() :
//                                            curr->outer_ccb()->face()) 
//                           << std::endl;
                
//                 int k = 0;
//                 const Halfedge* circ = curr;
//                 const Halfedge* first = circ;
//                 do {
//                     k++;
//                     circ = circ->next()->opposite();
//                 } while (circ != first);
                
//                 if (k > 2) {
//                     std::cout << k << " = degree > 2" << std::endl;
//                     std::cout << "PT: " 
//                               << curr->vertex()->point() << std::endl;
//                     circ = curr;
//                     first = circ;
//                     std::cout << "******************************************" 
//                               << std::endl;    
//                     do {
//                         if (!circ->has_null_curve()) {
//                             std::cout << "ccurve: " << circ->curve() 
//                                       << std::endl;
//                         }
//                         std::cout << "cdir: " << circ->direction() 
//                                   << std::endl;
// #if 0
//                         if (!curr->vertex()->has_null_point()) {
//                             std::cout << "cface: " << (circ->is_on_inner_ccb() ? 
//                                                        circ->inner_ccb()->face() : 
//                                                        circ->outer_ccb()->face())
//                                       << std::endl;
//                         }
// #endif
//                         std::cout << "******************************************" << std::endl;
//                         circ = circ->next()->opposite();
//                     } while (circ != first);
//                 }
#endif                


                CGAL::Sign tmp = 
                    _m_topology_traits->_sign_of_subpath(curr, next);
                
                if (tmp != CGAL::ZERO) {
                    switch (result) {
                    case ZERO:
                        result = tmp;
                        break;
                    default:
                        CGAL_assertion(result == -tmp || result == tmp);
                        result = CGAL::ZERO;
                    }
                }
            }
            
            curr = next;
            next = curr->next();
        } while (curr != he2);

#if CGAL_ARR_TOPOLOGY_TRAITS_VERBOSE  
        std::cout << "Sign_of_path(he1->he2) = " << result << std::endl;
#endif        
        return result;
    }
    
    //! computes sign of open path 
    // TODO: move to arr
    CGAL::Sign operator()(const Halfedge* he1, const Halfedge* he2, 
                          const X_monotone_curve_2& cv) const {
        
#if CGAL_ARR_TOPOLOGY_TRAITS_VERBOSE
        std::cout << "Sign_of_path(he1->he2->cv[->he1])?" << std::endl;
#endif
        
        //   source(he1)->target(he1)->....
        // ->source(he2)->target(he2)->cv[->source(he1)]

        CGAL::Sign result = this->operator()(he1, he2);

        typename Geometry_traits_2::Parameter_space_in_x_2 
            parameter_space_in_x =
            _m_topology_traits->geometry_traits()->
            parameter_space_in_x_2_object();
        typename Geometry_traits_2::Parameter_space_in_y_2 
            parameter_space_in_y =
            _m_topology_traits->geometry_traits()->
            parameter_space_in_y_2_object();
        
        // check whether cv can influence the counters
        
        CGAL::Arr_parameter_space ps_min_x = 
            parameter_space_in_x(cv, CGAL::ARR_MIN_END);
        CGAL::Arr_parameter_space ps_min_y = 
            parameter_space_in_y(cv, CGAL::ARR_MIN_END);
        
        CGAL::Arr_parameter_space ps_max_x = 
            parameter_space_in_x(cv, CGAL::ARR_MAX_END);  
        CGAL::Arr_parameter_space ps_max_y = 
            parameter_space_in_y(cv, CGAL::ARR_MAX_END);  
        
        if (ps_min_x != CGAL::ARR_INTERIOR || 
            ps_min_y != CGAL::ARR_INTERIOR || 
            ps_max_x != CGAL::ARR_INTERIOR || 
            ps_max_y != CGAL::ARR_INTERIOR) {
            
            // sign can change!
            
            bool equalmin = false;
            
            Point_2 minp = 
                _m_topology_traits->geometry_traits()->
                construct_min_vertex_2_object()(cv);
            
            CGAL::Arr_curve_end he2_tgt_ind =
                (he2->direction() == CGAL::ARR_LEFT_TO_RIGHT ? 
                 CGAL::ARR_MAX_END : CGAL::ARR_MIN_END
                );
            
            CGAL::Arr_parameter_space he2_tgt_psx = 
                parameter_space_in_x(he2->curve(), he2_tgt_ind);
            CGAL::Arr_parameter_space he2_tgt_psy = 
                parameter_space_in_y(he2->curve(), he2_tgt_ind);
            
            bool he2_tgt_on_boundary = 
                (he2_tgt_psx != CGAL::ARR_INTERIOR ||
                 he2_tgt_psy != CGAL::ARR_INTERIOR);

            bool min_on_boundary = 
                (ps_min_x != CGAL::ARR_INTERIOR || 
                 ps_min_y != CGAL::ARR_INTERIOR);
            
            if (he2_tgt_on_boundary == min_on_boundary) {
                if (he2_tgt_on_boundary) {
                    // compare at boundary
                    equalmin = 
                        _m_topology_traits->are_equal(he2->vertex(), 
                                                      cv, CGAL::ARR_MIN_END, 
                                                      ps_min_x, ps_min_y);
                } else {
                    equalmin = (_m_topology_traits->geometry_traits()->
                                compare_xy_2_object()(
                                        he2->vertex()->point(), minp
                                ) == CGAL::EQUAL);
                }
            }
            
#if CGAL_ARR_TOPOLOGY_TRAITS_VERBOSE
            std::cout << "equalmin: " << equalmin << std::endl;
#endif
            if (he2_tgt_on_boundary) {
                
                CGAL::Sign tmp1 = 
                    _m_topology_traits->_sign_of_subpath(
                            he2, true, // target!
                            cv, 
                            equalmin ? CGAL::ARR_MIN_END : CGAL::ARR_MAX_END
                    );
                
#if CGAL_ARR_TOPOLOGY_TRAITS_VERBOSE
                std::cout << "target(he2)->cv: " << tmp1 << std::endl;
#endif
                if (tmp1 != CGAL::ZERO) {
                    switch (result) {
                    case ZERO:
                        result = tmp1;
                        break;
                    default:
                        CGAL_assertion(result == -tmp1 || result == tmp1);
                        result = CGAL::ZERO;
                    }
                }
            }
            
            CGAL::Arr_curve_end he1_src_ind =
                (he1->direction() == CGAL::ARR_LEFT_TO_RIGHT ? 
                 CGAL::ARR_MIN_END : CGAL::ARR_MAX_END
                );
            
            CGAL::Arr_parameter_space he1_src_psx = 
                parameter_space_in_x(he1->curve(), he1_src_ind);
            CGAL::Arr_parameter_space he1_src_psy = 
                parameter_space_in_y(he1->curve(), he1_src_ind);
            
            bool he1_src_on_boundary = 
                (he1_src_psx != CGAL::ARR_INTERIOR ||
                 he1_src_psy != CGAL::ARR_INTERIOR);

            if (he1_src_on_boundary) {
                
                CGAL::Sign tmp2 = 
                    _m_topology_traits->_sign_of_subpath(
                            he1, false, // source
                            cv,
                            equalmin ? CGAL::ARR_MAX_END : CGAL::ARR_MIN_END
                    );
                
#if CGAL_ARR_TOPOLOGY_TRAITS_VERBOSE          
                std::cout << "cv->source(he1): " << tmp2 << std::endl;
#endif
                
                if (tmp2 != CGAL::ZERO) {
                    switch (result) {
                    case ZERO:
                        result = tmp2;
                        break;
                    default:
                        CGAL_assertion(result == -tmp2 || result == tmp2);
                        result = CGAL::ZERO;
                    }
                }
            }
        }
        
#if CGAL_ARR_TOPOLOGY_TRAITS_VERBOSE
        std::cout << "Sign_of_path(he1->he2->cv[->he1]) = " 
                  << result << std::endl;
#endif
        return result;
    }

    //!@}

protected:
    //! stores pointer to geometry traits
    const Topology_traits_2 *_m_topology_traits;
};

} // namespace internal

} //namespace CGAL

#endif // CGAL_ARR_TOPOLOGY_TRAITS_SIGN_OF_PATH_H
// EOF
