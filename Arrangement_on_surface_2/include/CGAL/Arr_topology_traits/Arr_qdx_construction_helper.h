// Copyright (c) 2006, 2007  Tel-Aviv University (Israel), Max-Planck-Institut
// fuer Informatik (Germany).
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
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//                 Ron Wein <wein@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>
//                 

#ifndef CGAL_ARR_QdX_CONSTRUCTION_HELPER_H
#define CGAL_ARR_QdX_CONSTRUCTION_HELPER_H

/*! \file
 * Definition of the Arr_qdx_construction_helper class-template.
 */

#include <CGAL/Sweep_line_empty_visitor.h>
#include <CGAL/Unique_hash_map.h>

CGAL_BEGIN_NAMESPACE

/*! \class Arr_qdx_construction_helper
 * A helper class for the construction sweep-line visitor, suitable
 * for an Arrangement_on_surface_2 instantiated with a topology-traits class
 * for curves on an ellipsoid
 */
template <class Traits_, class Arrangement_, class Event_, class Subcurve_> 
class Arr_qdx_construction_helper
{
public:
    typedef Traits_                                         Traits_2;
    typedef Arrangement_                                    Arrangement_2;
    typedef Event_                                          Event;
    typedef Subcurve_                                       Subcurve;
    
    typedef typename Traits_2::X_monotone_curve_2           X_monotone_curve_2;
    typedef typename Traits_2::Point_2                      Point_2;
    
    typedef Sweep_line_empty_visitor<Traits_2, Subcurve, Event>
    Base_visitor;
    
    typedef typename Arrangement_2::Vertex_handle           Vertex_handle;
    typedef typename Arrangement_2::Halfedge_handle         Halfedge_handle;
    typedef typename Arrangement_2::Face_handle             Face_handle;
    
    typedef typename Subcurve::Halfedge_indices_list        Indices_list;
    typedef CGAL::Unique_hash_map<Halfedge_handle, Indices_list> 
                                                          Halfedge_indices_map;
    
protected:
    typedef typename Arrangement_2::Topology_traits         Topology_traits;
    
    typedef typename Topology_traits::Vertex                DVertex;
    typedef typename Topology_traits::Halfedge              DHalfedge;
    
    // Data members:
    
    //! The topology-traits class
    Topology_traits * m_top_traits;
    Arr_accessor<Arrangement_2>  m_arr_access;  // An arrangement accessor.
    
    //! The unbounded arrangement face
    Face_handle m_top_face;
    
    //! Indices of the curves that "see" the rightmost face
    Indices_list m_subcurves_at_tf;
    
    //! A pointer to a map of halfedges to indices lists
    // (stored in the visitor class)
    Halfedge_indices_map * m_he_ind_map_p;
    
public:
    /*! Constructor. */
    Arr_qdx_construction_helper(Arrangement_2 * arr) :
        m_top_traits(arr->topology_traits()),
        m_arr_access(*arr),
        m_he_ind_map_p(NULL)
    {}
    
    /*! Destructor. */
    virtual ~Arr_qdx_construction_helper() {}
    
    /// \name Notification functions.
    //@{
    
    /* A notification issued before the sweep process starts. */
    virtual void before_sweep()
    {
        // Get the unbounded face.
        m_top_face = Face_handle(m_top_traits->top_face());
    }
    
    /*! A notification invoked before the sweep-line starts handling the given
     * event.
     */
    virtual void before_handle_event(Event * event)
    {
        const CGAL::Arr_parameter_space ps_x = event->parameter_space_in_x();
        const CGAL::Arr_parameter_space ps_y = event->parameter_space_in_y();
        
        if (ps_x == CGAL::ARR_LEFT_BOUNDARY) {
            CGAL_assertion(ps_y == CGAL::ARR_INTERIOR);
            // The event has only one right curve.
            CGAL_assertion((event->number_of_left_curves() == 0) &&
                           (event->number_of_right_curves() == 1));
            
            const X_monotone_curve_2 & xc =
                (*(event->right_curves_begin()))->last_curve();
            
            if (m_top_traits->leftmost_vertex() == NULL)
            {
                m_arr_access.create_boundary_vertex (xc, CGAL::ARR_MIN_END,
                                                     ps_x, ps_y);
            }

            event->set_vertex_handle(
                    Vertex_handle(m_top_traits->leftmost_vertex())
            );
            return;
        }
        
        if (ps_x == CGAL::ARR_RIGHT_BOUNDARY) {
            CGAL_assertion(ps_y == CGAL::ARR_INTERIOR);
            // The event has only one left curve.
            CGAL_assertion((event->number_of_left_curves() == 1) &&
                           (event->number_of_right_curves() == 0));
            
            const X_monotone_curve_2 & xc =
                (*(event->left_curves_begin()))->last_curve();
            
            if (m_top_traits->rightmost_vertex() == NULL)
            {
                m_arr_access.create_boundary_vertex (xc, CGAL::ARR_MAX_END,
                                                     ps_x, ps_y);
            } 

            event->set_vertex_handle(
                    Vertex_handle(m_top_traits->rightmost_vertex())
            );  
            return;
        }
        
        if (ps_y != CGAL::ARR_INTERIOR) {
            CGAL_assertion(ps_y == CGAL::ARR_TOP_BOUNDARY ||
                           ps_y == CGAL::ARR_BOTTOM_BOUNDARY);
            CGAL_assertion(ps_x == CGAL::ARR_INTERIOR);
            // there is exactly one event
            CGAL_assertion(event->number_of_left_curves() +
                           event->number_of_right_curves() >= 1);
            const X_monotone_curve_2 & xc =
                (ps_y == CGAL::ARR_BOTTOM_BOUNDARY ? 
                 // bottom
                 (event->number_of_right_curves() > 0 ? 
                  (*(event->right_curves_begin()))->last_curve() :
                  (*(event->left_curves_rbegin()))->last_curve()) :
                 // top
                 (event->number_of_left_curves() > 0 ? 
                  (*(event->left_curves_rbegin()))->last_curve() :
                  (*(event->right_curves_begin()))->last_curve()) 
                );
            
            CGAL::Arr_curve_end ind =  
                (ps_y == CGAL::ARR_BOTTOM_BOUNDARY ? 
                 // bottom
                 (event->number_of_right_curves() > 0 ? 
                  CGAL::ARR_MIN_END : CGAL::ARR_MAX_END) : 
                 // top
                 (event->number_of_left_curves() > 0 ? 
                  CGAL::ARR_MAX_END : CGAL::ARR_MIN_END)
                );
            
            DVertex * v = 
                m_top_traits->vertex_on_identification(event->point());
            
            Vertex_handle vh(v);
            
            if (v == NULL) {
                vh = m_arr_access.create_boundary_vertex (xc, ind,
                                                     ps_x, ps_y);
            } 
            
            event->set_vertex_handle(vh);
            return;
        }
    }
    
    /*! A notification invoked when a new subcurve is created. */
    virtual void add_subcurve(Halfedge_handle he, Subcurve * sc) { 
        
        // Check whether the halfedge (or its twin) lie on the top face.
        Halfedge_handle     he_on_top_face;
        
        CGAL::Arr_parameter_space ps_y_min = 
            this->m_top_traits->geometry_traits()->
            parameter_space_in_y_2_object()(he->curve(), CGAL::ARR_MIN_END);
        CGAL::Arr_parameter_space ps_y_max = 
            this->m_top_traits->geometry_traits()->
            parameter_space_in_y_2_object()(he->curve(), CGAL::ARR_MAX_END);
        
        if (ps_y_min == CGAL::ARR_BOTTOM_BOUNDARY) {
            he_on_top_face = 
                (he->direction() == CGAL::ARR_RIGHT_TO_LEFT ? he : he->twin());
        } else if (ps_y_max == CGAL::ARR_TOP_BOUNDARY) {
            he_on_top_face = 
                (he->direction() == CGAL::ARR_LEFT_TO_RIGHT ? he : he->twin());
        } else {
            return;
        }
        
        // Associate all curve indices of subcurves that lie in the current
        // top face with a halfedge that see this halfedge.
        if (m_he_ind_map_p != NULL) {
            Indices_list& list_ref = (*m_he_ind_map_p)[he_on_top_face];
            list_ref.splice (list_ref.end(), m_subcurves_at_tf);
        } else {
            m_subcurves_at_tf.clear();
        }
        CGAL_assertion (m_subcurves_at_tf.empty());
        
        return;
    }
    
    /*! Collect a subcurve index that does not see any status-line from below.
     */
    void add_subcurve_in_top_face(unsigned int index)
    {
        m_subcurves_at_tf.push_back(index);
        return;
    }
    
    /*! A notification invoked before the given event it deallocated. */
    void before_deallocate_event(Event * event) { 
        return; 
    }
    //@} 
    
    /*! Set the map that maps each halfedge to the list of subcurve indices
     * that "see" the halfedge from below.
     */
    void set_halfedge_indices_map(Halfedge_indices_map & table)
    {
        m_he_ind_map_p = &table;
        return;
    }
    
    /*! Determine if we should swap the order of predecessor halfedges when
     * calling insert_at_vertices_ex() .
     */
    bool swap_predecessors(Event * event) const
    {
        // If we insert an edge whose right end has boundary condition
        // before the curve of identification
        // have to flip the order of predecessor halfegdes.
        return (event->parameter_space_in_y() == CGAL::ARR_TOP_BOUNDARY &&
                event->parameter_space_in_x() == CGAL::ARR_INTERIOR);
    }
    
    /*! Get the current top face. */
    Face_handle top_face() const { 
        return m_top_face; 
    }
};

CGAL_END_NAMESPACE

#endif // CGAL_ARR_QdX_CONSTRUCTION_HELPER
