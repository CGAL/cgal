// Copyright (c) 1997  Tel-Aviv University (Israel).
// All rights reserved.
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
// Author(s)     : Eric Berberich <erric@mpi-inf.mpg.de>
//                 Ron Wein <wein@post.tau.ac.il>

#ifndef CGAL_ARR_QdX_OVERLAY_HELPER_H
#define CGAL_ARR_QdX_OVERLAY_HELPER_H

/*!
 * Definition of the Arr_qdx_overlay_helper class-template.
 */

#include <CGAL/Arr_topology_traits/Arr_qdx_construction_helper.h>

CGAL_BEGIN_NAMESPACE

/*! \class Arr_qdx_overlay_helper
 * A helper class for the overlay sweep-line visitor, suitable for the overlay
 * of Arrangement_on_surface_2 objects instantiated with a topology-traits
 * class for unbounded curves in the plane.
 */
template <class Traits_,
          class ArrangementRed_,
          class ArrangementBlue_,
          class Arrangement_,
          class Event_,
          class Subcurve_>
class Arr_qdx_overlay_helper
{
public:
    
    typedef Traits_                                        Traits_2; 
    typedef Arrangement_                                   Arrangement_2;
    typedef Event_                                         Event;
    typedef Subcurve_                                      Subcurve;
    
    typedef typename Traits_2::X_monotone_curve_2          X_monotone_curve_2;
    typedef typename Traits_2::Point_2                     Point_2;
    
    // The input arrangements (the "red" and the "blue" one):
    typedef ArrangementRed_                                 Arrangement_red_2;
    typedef typename Arrangement_red_2::Halfedge_const_handle
    Halfedge_handle_red;
    typedef typename Arrangement_red_2::Face_const_handle   Face_handle_red;
    typedef typename Arrangement_red_2::Vertex_const_handle Vertex_handle_red;
    
    typedef ArrangementBlue_                                Arrangement_blue_2;
    typedef typename Arrangement_blue_2::Halfedge_const_handle
    Halfedge_handle_blue;
    typedef typename Arrangement_blue_2::Face_const_handle  Face_handle_blue;
    typedef typename Arrangement_blue_2::Vertex_const_handle
    Vertex_handle_blue;
    
    // Define the helper class for the construction visitor.
    typedef Arr_qdx_construction_helper<Traits_2,
                                             Arrangement_2,
                                             Event,
                                             Subcurve>    Construction_helper;
    
protected:
    // Data members:
    const typename Arrangement_red_2::Topology_traits   *m_red_top_traits;
    const typename Arrangement_blue_2::Topology_traits  *m_blue_top_traits;

    //! Red quadric face
    Face_handle_red m_red_tf;

    //! Blue quadric face
    Face_handle_blue m_blue_tf;
    
public:
    
    /*! Constructor, given the input red and blue arrangements. */
    Arr_qdx_overlay_helper (const Arrangement_red_2 *red_arr,
                            const Arrangement_blue_2 *blue_arr) :
        m_red_top_traits (red_arr->topology_traits()),
        m_blue_top_traits (blue_arr->topology_traits())
    {}
    
    /// \name Notification functions.
    //@{
    
    /* A notification issued before the sweep process starts. */
    void before_sweep() {
        // Get the faces in both arrangements.
        // TODO set top_face in QdX_TOP_TRAITS
        m_red_tf = Face_handle_red(m_red_top_traits->top_face());
        m_blue_tf = Face_handle_blue(m_blue_top_traits->top_face());
    }
    
    /*!
     * A notification invoked before the sweep-line starts handling the given
     * event.
     */  
    void before_handle_event (Event* event) {
        // similar to torus
        Arr_curve_end ind = (event->number_of_left_curves() == 0 &&
                         event->number_of_right_curves() != 0) ?
            CGAL::ARR_MIN_END : CGAL::ARR_MAX_END;
        
        const Subcurve  *xc = (ind == CGAL::ARR_MIN_END) ?
            (*(event->right_curves_begin())) :
            (*(event->left_curves_begin()));
        
        if (event->parameter_space_in_x() == CGAL::ARR_LEFT_BOUNDARY) {
            // left side
            switch (xc->color()) {
            case Traits_2::RED :
                m_red_tf = xc->red_halfedge_handle()->twin()->face();
                break;
                
            case Traits_2::BLUE :
                m_blue_tf = xc->blue_halfedge_handle()->twin()->face();
                break;
                
            case Traits_2::RB_OVERLAP :
                m_red_tf = xc->red_halfedge_handle()->twin()->face();
                m_blue_tf = xc->blue_halfedge_handle()->twin()->face();
                break;
            }
        } else { 
            
            if (event->parameter_space_in_y() == CGAL::ARR_TOP_BOUNDARY) {
                
                Halfedge_handle_red red_he = xc->red_halfedge_handle();
                Halfedge_handle_blue blue_he = xc->blue_halfedge_handle();
                
                // top side
                switch (xc->color()) {
                case Traits_2::RED :
                    
                    CGAL_assertion(red_he != Halfedge_handle_red());
                    
                    if (ind == CGAL::ARR_MIN_END) {
                        this->m_red_tf = 
                            (red_he->direction() == CGAL::ARR_RIGHT_TO_LEFT ? 
                             red_he->twin() : red_he)->face();
                    } else {
                        this->m_red_tf = 
                            (red_he->direction() == CGAL::ARR_LEFT_TO_RIGHT ? 
                             red_he->twin() : red_he)->face();
                    }
                    
                    break;
                case Traits_2::BLUE :
                    
                    CGAL_assertion(blue_he != Halfedge_handle_blue());
                    
                    if (ind == CGAL::ARR_MIN_END) {
                        this->m_blue_tf = 
                            (blue_he->direction() == CGAL::ARR_RIGHT_TO_LEFT ?
                             blue_he->twin() : blue_he)->face();
                    } else {
                        this->m_blue_tf = 
                            (blue_he->direction() == CGAL::ARR_LEFT_TO_RIGHT ?
                             blue_he->twin() : blue_he)->face();
                    }
                    
                    break;
                case Traits_2::RB_OVERLAP :
                    
                    CGAL_assertion(red_he != Halfedge_handle_red());
                    CGAL_assertion(blue_he != Halfedge_handle_blue());

                    if (ind == CGAL::ARR_MIN_END) {
                        this->m_red_tf = 
                            (red_he->direction() == CGAL::ARR_RIGHT_TO_LEFT ? 
                             red_he->twin() : blue_he)->face();
                        this->m_blue_tf = 
                            (blue_he->direction() == CGAL::ARR_RIGHT_TO_LEFT ?
                             blue_he->twin() : blue_he)->face();
                    } else {
                        this->m_red_tf = 
                            (red_he->direction() == CGAL::ARR_LEFT_TO_RIGHT ? 
                             red_he->twin() : red_he)->face();
                        this->m_blue_tf = 
                            (blue_he->direction() == CGAL::ARR_LEFT_TO_RIGHT ?
                             blue_he->twin() : blue_he)->face();
                    }
                    break;
                }
            } 
        }
    }
    
    //@}
    
    /*! Get the current red top face. */
    Face_handle_red red_top_face () const
    {
        return (m_red_tf);
    }
    
    /*! Get the current blue top face. */
    Face_handle_blue blue_top_face () const
    {
        return (m_blue_tf);
    }
};

CGAL_END_NAMESPACE

#endif // ARR_QdX_OVERLAY_HELPER
