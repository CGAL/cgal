// Copyright (c) 2005,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Ron Wein         <wein@post.tau.ac.il>
//                 Efi Fogel        <efif@post.tau.ac.il>
//                 Baruch Zukerman  <baruchzu@post.tau.ac.il>

#ifndef CGAL_ARR_ON_SURFACE_WITH_HISTORY_2_FUNCTIONS_H
#define CGAL_ARR_ON_SURFACE_WITH_HISTORY_2_FUNCTIONS_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
 * Member-function definitions for the Arrangement_on_surface_with_history_2
 * class.
 */

namespace CGAL {

//-----------------------------------------------------------------------------
// Default constructor.
//
template<class GeomTr, class TopTr>
Arrangement_on_surface_with_history_2<GeomTr, TopTr>::
Arrangement_on_surface_with_history_2 () :
  Base_arr_2 ()
{
  m_observer.attach (*this);
}

//-----------------------------------------------------------------------------
// Copy constructor.
//
template<class GeomTr, class TopTr>
Arrangement_on_surface_with_history_2<GeomTr,TopTr>::
Arrangement_on_surface_with_history_2 (const Self& arr) :
  Base_arr_2 ()
{
  assign (arr);
  m_observer.attach (*this);
}

//-----------------------------------------------------------------------------
// Constructor given a traits object.
//
template<class GeomTr, class TopTr>
Arrangement_on_surface_with_history_2<GeomTr,TopTr>::
Arrangement_on_surface_with_history_2 (const Geometry_traits_2 * tr) :
  Base_arr_2 (static_cast<const Data_traits_2*> (tr))
{
  m_observer.attach (*this);
}

//-----------------------------------------------------------------------------
// Assignment operator.
//
template<class GeomTr, class TopTr>
Arrangement_on_surface_with_history_2<GeomTr,TopTr>&
Arrangement_on_surface_with_history_2<GeomTr,TopTr>::
operator=(const Self& arr)
{
  // Check for self-assignment.
  if (this == &arr)
    return (*this);
  
  assign (arr);
  return (*this);
}

//-----------------------------------------------------------------------------
// Assign an arrangement with history.
//
template<class GeomTr, class TopTr>
void Arrangement_on_surface_with_history_2<GeomTr,TopTr>::
assign(const Self& arr)
{
  // Clear the current contents of the arrangement.
  clear();
  
  // Assign the base arrangement.
  Base_arr_2::assign (arr);
  
  // Create duplicates of the stored curves and map the curves of the
  // original arrangement to their corresponding duplicates.
  typedef std::map<const Curve_halfedges*, Curve_halfedges*>  Curve_map;
  typedef typename Curve_map::value_type                      Curve_map_entry;
 
  Curve_map                 cv_map;
  Curve_const_iterator      ocit;
  const Curve_2            *p_cv;
  Curve_halfedges          *dup_c;

  for (ocit = arr.curves_begin(); ocit != arr.curves_end(); ++ocit)
  {
    // Create a duplicate of the current curve.
    dup_c = m_curves_alloc.allocate (1);
    
    p_cv = &(*ocit);
    m_curves_alloc.construct (dup_c, *p_cv);
    m_curves.push_back (*dup_c);
    
    // Assign a map entry.
    cv_map.insert (Curve_map_entry (&(*ocit), dup_c));
  }

  // Go over the list of halfedges in our arrangement. The curves associated
  // with these edges sotre pointers to the curves in the original
  // arrangement, so we now have to modify these pointers, according to the
  // mapping we have just created. While doing so, we also construct the set
  // of edges associated with each (duplicated) curve in our arrangement.
  Data_iterator                           dit;
  std::list<Curve_2*>                     dup_curves;
  typename std::list<Curve_2*>::iterator  iter;
  Edge_iterator                           eit;
  Halfedge_handle                         e;
  const Curve_halfedges                  *org_c;

  for (eit = this->edges_begin(); eit != this->edges_end(); ++eit)
  {
    e = eit;
    dup_curves.clear();
    for (dit = e->curve().data().begin(); 
         dit != e->curve().data().end(); ++dit)
    {
      org_c = static_cast<Curve_halfedges*>(*dit);
      dup_c = (cv_map.find (org_c))->second;
      
      dup_curves.push_back (dup_c);
      dup_c->_insert (e);
    }

    // Replace the curve pointers associated with the edge.
    e->curve().data().clear();
    for (iter = dup_curves.begin(); iter != dup_curves.end(); ++iter)
      e->curve().data().insert (*iter);
  }

  return;
}

//-----------------------------------------------------------------------------
// Destructor.
//
template<class GeomTr, class TopTr>
Arrangement_on_surface_with_history_2<GeomTr,TopTr>::
~Arrangement_on_surface_with_history_2 ()
{
  clear();
}

//-----------------------------------------------------------------------------
// Clear the arrangement.
//
template<class GeomTr, class TopTr>
void Arrangement_on_surface_with_history_2<GeomTr,TopTr>::clear ()
{
  // Free all stored curves.
  Curve_iterator         cit = m_curves.begin();
  Curve_halfedges       *p_cv;
  
  while (cit != m_curves.end())
  {
    p_cv = &(*cit);
    ++cit;
    
    m_curves.erase (p_cv);
    m_curves_alloc.destroy (p_cv);
    m_curves_alloc.deallocate (p_cv, 1);
  }
  m_curves.destroy();
  
  // Clear the base arrangement.
  Base_arr_2::clear();
  
  return;
}

//-----------------------------------------------------------------------------
// Split a given edge into two at the given split point.
//
template<class GeomTr, class TopTr>
typename Arrangement_on_surface_with_history_2<GeomTr,TopTr>::Halfedge_handle
Arrangement_on_surface_with_history_2<GeomTr,TopTr>::
split_edge(Halfedge_handle e, const Point_2& p)
{
  // Split the curve associated with the halfedge e at the given point p.
  Data_x_curve_2       cv1, cv2;
  
  this->m_geom_traits->split_2_object() (e->curve(), p,
                                       cv1, cv2);

  // cv1 always lies to the left of cv2. If e is directed from left to right,
  // we should split and return the halfedge associated with cv1, and
  // otherwise we should return the halfedge associated with cv2 after the
  // split.
  if (e->direction() == ARR_LEFT_TO_RIGHT)
  {
    return (Base_arr_2::split_edge (e, cv1, cv2));
  }
  else
  {
    return (Base_arr_2::split_edge (e, cv2, cv1));
  }
}

//-----------------------------------------------------------------------------
// Merge two edges to form a single edge.
//
template<class GeomTr, class TopTr>
typename Arrangement_on_surface_with_history_2<GeomTr,TopTr>::Halfedge_handle
Arrangement_on_surface_with_history_2<GeomTr,TopTr>::
merge_edge(Halfedge_handle e1, Halfedge_handle e2)
{
  CGAL_precondition_msg (are_mergeable(e1, e2), 
                         "Edges are not mergeable.");

  // Merge the two curves.
  Data_x_curve_2       cv;
  
  this->m_geom_traits->merge_2_object()(e1->curve(), e2->curve(), cv);
  
  return (Base_arr_2::merge_edge (e1, e2, cv));
}

//-----------------------------------------------------------------------------
// Check if two edges can be merged to a single edge.
//
template<class GeomTr, class TopTr>
bool Arrangement_on_surface_with_history_2<GeomTr,TopTr>::are_mergeable
    (Halfedge_const_handle e1,
     Halfedge_const_handle e2) const
{
  // Both halfedges must be non-fictitious.
  if (e1->is_fictitious() || e2->is_fictitious())
    return (false);

  // In order to be mergeable, the two halfedges must share a common
  // end-vertex. We assign vh to be this vertex.
  Vertex_const_handle      vh;
  
  if (e1->target() == e2->source() || e1->target() == e2->target())
  {
    vh = e1->target();
  }
  else
  {
    if (e1->source() == e2->source() || e1->source() == e2->target())
    {
      vh = e1->source();
    }
    else
    {
      // No common end-vertex: the edges are not mergeable.
      return (false);
    }
  }
  
  // If there are other edges incident to vh, it is impossible to remove it
  // and merge the two edges.
  if (vh->degree() != 2)
    return (false);
  
  // Check whether the curves associated with the two edges are mergeable.
  return (this->m_geom_traits->are_mergeable_2_object()(e1->curve(),
                                                        e2->curve()));
}

//-----------------------------------------------------------------------------
// Register a new observer (so it starts receiving notifications).
//
template<class GeomTr, class TopTr>
void Arrangement_on_surface_with_history_2<GeomTr,TopTr>::
_register_observer(Arr_observer<Self> *p_obs)
{
  Base_arr_2::_register_observer
    (reinterpret_cast<Arr_observer<Base_arr_2>*>(p_obs));
  return;
}

//-----------------------------------------------------------------------------
// Unregister an observer (so it stops receiving notifications).
//
template<class GeomTr, class TopTr>
bool Arrangement_on_surface_with_history_2<GeomTr,TopTr>::
_unregister_observer(Arr_observer<Self> *p_obs)
{
  return (Base_arr_2::_unregister_observer 
          (reinterpret_cast<Arr_observer<Base_arr_2>*>(p_obs)));
}

} //namespace CGAL

#endif
