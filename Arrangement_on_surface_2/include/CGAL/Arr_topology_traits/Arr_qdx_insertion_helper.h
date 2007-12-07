// Copyright (c) 2007  Tel-Aviv University (Israel), Max-Planck-Institute for
// Computer Science (Germany)
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
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ron Wein <wein@post.tau.ac.il>
//                 Eric Berberich <eric@mpi-inf.mpg.de>

#ifndef CGAL_ARR_QDX_INSERTION_HELPER_H
#define CGAL_ARR_QDX_INSERTION_HELPER_H

/*!
 * Definition of the Arr_qdx_insertion_helper class-template.
 */

#include <CGAL/Sweep_line_2/Arr_construction_sl_visitor.h>
#include <CGAL/Arr_topology_traits/Arr_qdx_construction_helper.h>

CGAL_BEGIN_NAMESPACE

/*! \class Arr_qdx_insertion_helper
 * A helper class for the insertion sweep-line visitors, suitable
 * for an Arrangement_on_surface_2 instantiated with a topology-traits class
 * for unbounded curves in the plane.
 */
template <class Traits_, class Arrangement_, class Event_, class Subcurve_> 
class Arr_qdx_insertion_helper :
    public Arr_qdx_construction_helper<Traits_, Arrangement_,
                                            Event_, Subcurve_>
{
public:
    
    typedef Traits_                                      Traits_2;
    typedef Arrangement_                                 Arrangement_2;
    typedef Event_                                       Event;
    typedef Subcurve_                                    Subcurve;

    typedef typename Traits_2::X_monotone_curve_2        X_monotone_curve_2;
    typedef typename Traits_2::Point_2                   Point_2;
    
    
    typedef Arr_qdx_construction_helper<Traits_2,
                                             Arrangement_2,
                                             Event,
                                             Subcurve> Base;
    
    typedef Sweep_line_empty_visitor<Traits_2,
                                   Subcurve,
                                   Event>              Base_visitor;
    
    typedef Arr_qdx_insertion_helper<Traits_2,
                                          Arrangement_2,
                                          Event,
                                          Subcurve>    Self;
    
    typedef Arr_construction_sl_visitor<Self>            Parent_visitor;
    
    typedef typename Arrangement_2::Face_handle          Face_handle;
    
    typedef typename Base::Indices_list                  Indices_list;
    typedef typename Base::Halfedge_indices_map          Halfedge_indices_map;
    
protected:
    
    typedef typename Base::Topology_traits               Topology_traits;
    typedef typename Base::Vertex_handle                 Vertex_handle;
    typedef typename Base::Halfedge_handle               Halfedge_handle;
    
public:
    
    /*! Constructor. */
    Arr_qdx_insertion_helper (Arrangement_2 *arr) :
        Base (arr)
    {}
    
    /*! Destructor. */
    virtual ~Arr_qdx_insertion_helper(){}
    
    /// \name Notification functions.
    //@{
    
    /* A notification issued before the sweep process starts. */
    virtual void before_sweep ();
    
    /*!
     * A notification invoked before the sweep-line starts handling the given
     * event.
     */
    virtual void before_handle_event (Event* event);
    //@}
};

//-----------------------------------------------------------------------------
// Memeber-function definitions:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// A notification issued before the sweep process starts.
//
template <class Tr, class Arr, class Evnt, class Sbcv> 
void Arr_qdx_insertion_helper<Tr,Arr,Evnt,Sbcv>::before_sweep ()
{
    this->m_top_face = Face_handle(this->m_top_traits->bottom_face());
    return;
}

//-----------------------------------------------------------------------------
// A notification invoked before the sweep-line starts handling the given
// event.
//
template <class Tr, class Arr, class Evnt, class Sbcv> 
void Arr_qdx_insertion_helper<Tr,Arr,Evnt,Sbcv>::before_handle_event
    (Event* event)
{
    // Ignore events that do not have boundary conditions.
    const Arr_parameter_space ps_x = event->parameter_space_in_x();
    const Arr_parameter_space ps_y = event->parameter_space_in_y();
    
    if (ps_x == ARR_INTERIOR && ps_y == ARR_INTERIOR) {
        return;
    }
    
    // The is exactly one curve incident to an event with boundary conditions.
    // Obtain this curve and check whether it 
    // already exists in the arrangement.
    CGAL_assertion(((event->number_of_left_curves() == 0) &&
                    (event->number_of_right_curves() == 1)) ||
                   ((event->number_of_left_curves() == 1) &&
                    (event->number_of_right_curves() == 0)));
    
    const Arr_curve_end ind = ((event->number_of_left_curves() == 0 &&
                            event->number_of_right_curves() == 1) ? 
                           ARR_MIN_END :  ARR_MAX_END);
    const X_monotone_curve_2& xc = (ind == ARR_MIN_END) ?
        (*(event->right_curves_begin()))->last_curve() :
        (*(event->left_curves_begin()))->last_curve();
    
    Halfedge_handle he = xc.halfedge_handle();

    if (he == Halfedge_handle())
    {
        // The curve is not in the arrangement, 
        // use the base construction helper
        // to handle the event:
        Base::before_handle_event (event);
        return;
    }
    
    // In case we encounter an existing curve incident to the left or 
    // top face, we have to update the current top face
    if (ps_x < 0) {
        // left side
        CGAL_assertion (ind == ARR_MIN_END);
        this->m_top_face = he->twin()->face();
    } else if (ps_y == BEFORE_DISCONTINUITY) {
        // top side (taken from torus)

        if (ind == CGAL::ARR_MIN_END) {
            this->m_top_face = 
                (he->direction() == CGAL::ARR_RIGHT_TO_LEFT ? 
                 he->twin() : he)->face();
        } else {
            this->m_top_face = 
                (he->direction() == CGAL::ARR_LEFT_TO_RIGHT ? 
                 he->twin() : he)->face();
        }
    }
    return;
}

CGAL_END_NAMESPACE

#endif
