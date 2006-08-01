// Copyright (c) 2006  Tel-Aviv University (Israel).
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
// $Source: $
// $Revision$ $Date$
// $Name:  $
//
// Author(s)     : Ron Wein   <wein@post.tau.ac.il>

#ifndef CGAL_ENVELOPE_DIVIDE_AND_CONQUER_2_IMPL_H
#define CGAL_ENVELOPE_DIVIDE_AND_CONQUER_2_IMPL_H

/*! \file
 * Definitions of the functions of the Envelope_divide_and_conquer_2 class.
 */

CGAL_BEGIN_NAMESPACE

// ---------------------------------------------------------------------------
// Construct the lower/upper envelope of the given list of non-vertical curves.
//
template <class Traits, class Diagram>
void Envelope_divide_and_conquer_2<Traits,Diagram>::
_construct_envelope_non_vertical (Curve_pointer_iterator begin,
                                  Curve_pointer_iterator end,
                                  Envelope_diagram_1& out_d)
{
  out_d.clear();
  
  if (begin == end)
  {
    return;
  }
  
  // Check if the range contains just a single curve.
  Curve_pointer_iterator    iter = begin;
  ++iter;
  
  if (iter == end)
  {
    // Construct a singleton diagram, which matches a single curve.
    _construct_singleton_diagram (*(*begin), out_d);
  }
  else
  {
    // Divide the given range of curves into two.
    Curve_pointer_iterator  div_it = begin;
    unsigned int            count = 0;
    
    for (iter = begin; iter != end; ++iter)
    {
      if (count % 2 == 0)
        ++div_it;
      
      count++;
    }
    
    // Construct the diagrams (envelopes) for the two sub-ranges recursively 
    // and then merge the two diagrams to obtain the result.
    Envelope_diagram_1   d1;
    Envelope_diagram_1   d2;
    
    _construct_envelope_non_vertical (begin, div_it,
                                      d1);
    
    _construct_envelope_non_vertical (div_it, end,
                                      d2);

    _merge_envelopes (d1, d2, out_d);
  }
  
  return;
}

// ---------------------------------------------------------------------------
// Construct a singleton diagram, which matches a single curve.
//
template <class Traits, class Diagram>
void Envelope_divide_and_conquer_2<Traits,Diagram>::
_construct_singleton_diagram (const X_monotone_curve_2& cv,
                              Envelope_diagram_1& out_d)
{
  CGAL_assertion (out_d.leftmost() == out_d.rightmost());
  CGAL_assertion (out_d.leftmost()->is_empty());
  
  // Check if the given curve is bounded from the left and from the right.
  if (traits->infinite_in_x_2_object() (cv, MIN_END) != FINITE)
  {
    if (traits->infinite_in_x_2_object() (cv, MAX_END) != FINITE)
    {
      // The curve is defined over (-oo, oo), so its diagram contains
      // only a single edge.
      out_d.leftmost()->add_curve (cv);
      
      return;
    }

    // The curve is defined over (-oo, x], where x is finite.
    // Create a vertex and associate it with the right endpoint of cv.
    CGAL_precondition
      (traits->infinite_in_y_2_object() (cv, MAX_END) == FINITE);
    
    Vertex_handle  v = 
      out_d.new_vertex (traits->construct_max_vertex_2_object() (cv));
    Edge_handle    e_right = out_d.new_edge();
    
    v->add_curve (cv);
    v->set_left (out_d.leftmost());
    v->set_right (e_right);
    
    // The leftmost edge is associated with cv, and the rightmost is empty.
    out_d.leftmost()->add_curve (cv);
    out_d.leftmost()->set_right (v);
    
    e_right->set_left (v);
    out_d.set_rightmost (e_right);
    
    return;
  }
  
  if (traits->infinite_in_x_2_object() (cv, MAX_END) != FINITE)
  {
    // The curve is defined over [x, +oo), where x is finite.
    // Create a vertex and associate it with the left endpoint of cv.
    CGAL_precondition
      (traits->infinite_in_y_2_object() (cv, MIN_END) == FINITE);
    
    Vertex_handle  v = 
      out_d.new_vertex (traits->construct_min_vertex_2_object() (cv));
    Edge_handle    e_left = out_d.new_edge();
    
    v->add_curve (cv);
    v->set_left (e_left);
    v->set_right (out_d.rightmost());
    
    // The rightmost edge is associated with cv, and the leftmost is empty.
    out_d.rightmost()->add_curve (cv);
    out_d.rightmost()->set_left (v);
    
    e_left->set_right (v);
    out_d.set_leftmost (e_left);
    
    return;
  }
  
  // If we reached here, the curve is defined over a bounded x-range.
  // We therefore create the following diagram:
  //
  //             (empty)    v1     e       v2   (empty)
  //      -oo -------------(+)============(+)------------ +oo
  //
  CGAL_precondition
    (traits->infinite_in_y_2_object() (cv, MIN_END) == FINITE);
  CGAL_precondition
    (traits->infinite_in_y_2_object() (cv, MAX_END) == FINITE);
  
  Vertex_handle  v1 = 
    out_d.new_vertex (traits->construct_min_vertex_2_object() (cv));
  Vertex_handle  v2 = 
    out_d.new_vertex (traits->construct_max_vertex_2_object() (cv));
  Edge_handle    e_left = out_d.new_edge();
  Edge_handle    e_right = out_d.new_edge();
  Edge_handle    e = out_d.leftmost();
  
  v1->add_curve (cv);
  v1->set_left (e_left);
  v1->set_right (e);
  
  v2->add_curve (cv);
  v2->set_left (e);
  v2->set_right (e_right);
  
  e->add_curve (cv);
  e->set_left (v1);
  e->set_right (v2);
  
  e_left->set_right (v1);
  e_right->set_left (v2);
  
  out_d.set_leftmost (e_left);
  out_d.set_rightmost (e_right);
  
  return;
}

// ---------------------------------------------------------------------------
// Merge two minimization (or maximization) diagrams.
//
template <class Traits, class Diagram>
void Envelope_divide_and_conquer_2<Traits,Diagram>::
_merge_envelopes (const Envelope_diagram_1& d1,
                  const Envelope_diagram_1& d2,
                  Envelope_diagram_1& out_d)
{
  Edge_const_handle    e1 = d1.leftmost();
  bool                 is_leftmost1 = true;
  Vertex_const_handle  v1;
  Edge_const_handle    e2 = d2.leftmost();
  bool                 is_leftmost2 = true;
  Vertex_const_handle  v2;
  Vertex_const_handle  next_v;
  bool                 next_exists = true;
  Comparison_result    res_v;
  bool                 same_x = false;
  
  do
  {
    // Locate the vertex that has smaller x-coordinate between v1 and v2.
    // If both have the same x-ccordinate, find the one that should be in
    // the envelope.
    same_x = false;
    
    if (e1 == d1.rightmost())
    {
      if (e2 == d2.rightmost())
      {
        // Both current edges do not have a vertex to their right.
        next_exists = false;
      }
      else
      {
        // e1 is not bounded from the right while e2 is.
        v2 = e2->right();
        next_v = v2;
        res_v = LARGER;
      }
    }
    else if (e2 == d2.rightmost())
    {
      // e2 is not bounded from the right while e1 is.
      v1 = e1->right();
      next_v = v1;
      res_v = SMALLER;
    }
    else
    {
      v1 = e1->right();
      v2 = e2->right();
      res_v = _compare_vertices (v1, v2, same_x);
      next_v = (res_v == SMALLER) ? v1 : v2;
    }
    
    // Check if the current edges represent empty intervals or not.
    if (! e1->is_empty() && ! e2->is_empty())
    {
      // Both edges are not empty, and there are curves defined on them.
      _merge_two_intervals (e1, is_leftmost1,
                            e2, is_leftmost2,
                            next_v, next_exists,
                            (res_v == SMALLER) ? 1 : 2,
                            out_d);
    }
    else if (! e1->is_empty() && e2->is_empty())
    {
      // e1 is not empty but e2 is empty:
      _merge_single_interval (e1,
                              next_v, next_exists,
                              (res_v == SMALLER),
                              out_d);
    }
    else if (e1->is_empty() && ! e2->is_empty())
    {
      // e1 is empty and e2 is not empty:
      _merge_single_interval (e2,
                              next_v, next_exists,
                              (res_v != SMALLER),
                              out_d);
    }
    else
    {
      // Both edges are empty: append an empty edge to out_d:
      if (next_exists)
      {
        Vertex_handle  new_v = _append_vertex (out_d, next_v->point(), e1);
        new_v->add_curves (next_v->curves_begin(), next_v->curves_end());
      }
    }
    
    // Proceed to the next diagram edge(s), if possible.
    if (next_exists)
    {
      // Check if we should proceed on d1 or on d2.
      if (res_v == SMALLER)
      {
        e1 = v1->right();
        is_leftmost1 = false;

        if (same_x)
        {
          e2 = v2->right();
          is_leftmost2 = false;
        }
      }
      else if (res_v == LARGER)
      {
        e2 = v2->right();
        is_leftmost2 = false;

        if (same_x)
        {
          e1 = v1->right();
          is_leftmost1 = false;
        }
      }
      else
      {
        e1 = v1->right();
        is_leftmost1 = false;
                
        e2 = v2->right();
        is_leftmost2 = false;
      }
    }
    
  } while (next_exists);
  
  return;
}

// ---------------------------------------------------------------------------
// Compare two diagram vertices.
//
template <class Traits, class Diagram>
Comparison_result Envelope_divide_and_conquer_2<Traits,Diagram>::
_compare_vertices (Vertex_const_handle v1,
                   Vertex_const_handle v2,
                   bool& same_x) const
{
  Comparison_result   res = traits->compare_x_2_object() (v1->point(),
                                                          v2->point());
  
  if (res != EQUAL)
  {
    same_x = false;
    return (res);
  }
  else
  {
    same_x = true;
  }

  // In case the x-coordinates of the two vertices are equal:
  res = traits->compare_xy_2_object() (v1->point(),
                                       v2->point());
  
  if ((env_type == LOWER && res == SMALLER) ||
      (env_type == UPPER && res == LARGER))
    return (SMALLER);
  else if ((env_type == LOWER && res == LARGER) ||
           (env_type == UPPER && res == SMALLER))
    return (LARGER);
  
  // The two vertices represent equal points:
  return (EQUAL);
}

// ---------------------------------------------------------------------------
// Deal with an interval which is non-empty in one of the merged diagrams and
// empty in the other.
//
template <class Traits, class Diagram>
void Envelope_divide_and_conquer_2<Traits,Diagram>::
_merge_single_interval (Edge_const_handle e,
                        Vertex_const_handle v, bool v_exists,
                        bool same_org,
                        Envelope_diagram_1& out_d)
{
  if (! v_exists)
  {
    // The non-empty edge e is unbounded from the right, so we simply have
    // to update the rightmost edge in out_d.
    out_d.rightmost()->add_curves (e->curves_begin(), e->curves_end());
    return;
  }
  
  Vertex_handle      new_v;
  
  if (same_org)
  {
    // The non-empty edge ends at v, so we simply insert it to out_d.
    new_v = _append_vertex (out_d, v->point(), e);
    new_v->add_curves (v->curves_begin(), v->curves_end());
    
    return;
  }
  
  // If v is not on e, we should insert it to the merged diagram only if it
  // is below (or above, in case of an upper envelope) the curves of e.
  Comparison_result  res = traits->compare_y_at_x_2_object() (v->point(),
                                                              e->curve());
  
  if ((res == EQUAL) ||
      (env_type == LOWER && res == SMALLER) ||
      (env_type == UPPER && res == LARGER))
  {
    new_v = _append_vertex (out_d, v->point(), e);
    new_v->add_curves (v->curves_begin(), v->curves_end());
  }
  
  if (res == EQUAL)
  {
    // In case of equality, append e's curves to those of the new vertex.
    new_v->add_curves (e->curves_begin(), e->curves_end());
  }
  
  return;
}

// ---------------------------------------------------------------------------
// Merge two non-empty intervals into the merged diagram.
//
template <class Traits, class Diagram>
void Envelope_divide_and_conquer_2<Traits,Diagram>::
_merge_two_intervals (Edge_const_handle e1, bool is_leftmost1,
                      Edge_const_handle e2, bool is_leftmost2,
                      Vertex_const_handle v, bool v_exists,
                      int org_v,
                      Envelope_diagram_1& out_d)
{
  // Get the relative position of two curves associated with e1 and e2
  // at the rightmost of the left endpoints of e1 and e2.
  Comparison_result      u_res;
  bool                   equal_at_v = false;
  Point_2                pu;
  bool                   pu_exists = false;
 
  u_res = traits->compare_y_position_2_object() (e1->curve(),
                                                 e2->curve());
  
  // Flip the result in case of an upper envelope.
  if (env_type == UPPER)
    u_res = CGAL::opposite (u_res);
  
  // Use the current rightmost of the diagram as a reference point.
  bool               v_rm_exists = (out_d.leftmost() != out_d.rightmost());
  Vertex_handle      v_rm;

  if (v_rm_exists)
    v_rm = out_d.rightmost()->left();

  // Find the next intersection of the envelopes to the right of the current
  // rightmost point in the merged diagram.
  std::list<CGAL::Object>           objects;
  CGAL::Object                      obj;
  X_monotone_curve_2                icv;
  std::pair<Point_2, unsigned int>  ipt;
  
  traits->intersect_2_object() (e1->curve(), e2->curve(),
                                std::back_inserter(objects));
  
  while (! objects.empty())
  {
    // Pop the xy-lexicographically smallest intersection object.
    obj = objects.front();
    objects.pop_front();
    
    if (CGAL::assign(ipt, obj))
    {
      // We have a simple intersection point.
      if (v_rm_exists &&
          traits->compare_xy_2_object() (ipt.first, 
                                         v_rm->point()) != LARGER)
      {
        // The point is to the left of the current rightmost vertex in out_d,
        // so we skip it and continue examining the next intersections.
        // However, we update the last intersection point observed.
        pu = ipt.first;
        pu_exists = true;
        continue;
      }
      
      if (v_exists)
      {
        Comparison_result res = traits->compare_xy_2_object() (ipt.first, 
                                                               v->point());
        
        if (res == EQUAL)
          // v is an intersection points, so both curves are equal there:
          equal_at_v = true;
        
        if (res != SMALLER)
        {
          // We passed the next vertex, so we can stop here.
          break;
        }
      }
      
      // Create a new vertex in the output diagram that corrsponds to the
      // current intersection point.
      Vertex_handle        new_v;
      
      if (pu_exists)
      {
        // Update the relative position of the two curves, which is their
        // order immediately to the right of their last observed intersection
        // point pu.
        u_res = traits->compare_y_at_x_right_2_object() (e1->curve(),
                                                         e2->curve(),
                                                         pu);
 
        if (env_type == UPPER)
          u_res = CGAL::opposite (u_res);
      }
      CGAL_assertion (u_res != EQUAL);
      
      if (u_res == SMALLER)
        new_v = _append_vertex (out_d, ipt.first, e1);
      else
        new_v = _append_vertex (out_d, ipt.first, e2);
      
      // Note that the new vertex is incident to all curves in e1 and in e2.
      new_v->add_curves (e1->curves_begin(), e1->curves_end());
      new_v->add_curves (e2->curves_begin(), e2->curves_end());
      
      // Update the handle to the rightmost vertex in the output diagram.
      v_rm = new_v;
      v_rm_exists = true;

      // Get the curve order immediately to the right of the intersection
      // point. Note that in case of even (non-zero) multiplicity the order
      // remains the same.
      if (ipt.second % 2 == 1)
      {
        // Odd multiplicity: flip the current comparison result.
        u_res = CGAL::opposite (u_res);
      }
      else if (ipt.second == 0)
      {
        // The multiplicity is unknown, so we have to compare the curves to
        // the right of their intersection point.
        u_res = traits->compare_y_at_x_right_2_object() (e1->curve(),
                                                         e2->curve(),
                                                         ipt.first);
        CGAL_assertion (u_res != EQUAL);
      }
    }
    else
    {
      // We have an x-monotone curve representing an overlap of the two
      // curves.
      bool     assign_success = CGAL::assign (icv, obj);
      
      CGAL_assertion (assign_success);
      if (! assign_success)
        continue;

      // Get the endpoints of the overlapping curves.
      const bool  has_left = 
        (traits->infinite_in_x_2_object() (icv, MIN_END) == FINITE);
      const bool  has_right = 
        (traits->infinite_in_x_2_object() (icv, MAX_END) == FINITE);
      Point_2     p1, p2;
      
      if (has_left)
        p1 = traits->construct_min_vertex_2_object() (icv);
      
      if (has_right)
        p2 = traits->construct_max_vertex_2_object() (icv);
      
      // Check if the overlapping curve is not relevant to our range.
      if (v_rm_exists && has_right &&
          traits->compare_xy_2_object() (p2, 
                                         v_rm->point()) != LARGER)
      {
        // The right point of the overlappinf curve is to the left of the
        // current rightmost vertex in out_d, so we skip it and continue
        // examining the next intersections.
        // However, we update the last intersection point observed.
        pu = p2;
        pu_exists = true;
        continue;
      }
      
      if (v_exists && has_right)
      {
        Comparison_result res = traits->compare_xy_2_object() (p1, 
                                                               v->point());
        
        if (res == EQUAL)
          // v is an intersection points, so both curves are equal there:
          equal_at_v = true;
        
        if (res != SMALLER)
        {
          // We passed the next vertex, so we can stop here.
          break;
        }
      }
      
      // There is an overlap between the range [u, v] and icv.
      if (has_left && 
          (! v_rm_exists ||
           traits->compare_xy_2_object() (p1, 
                                          v_rm->point()) == LARGER))
      {
        // Create an output edge that represent the portion of [u, v] to the
        // left of the overlapping curve.
        Vertex_handle        new_v;

        if (pu_exists)
        {
          // Update the relative position of the two curves, which is their
          // order immediately to the right of their last observed intersection
          // point pu.
          u_res = traits->compare_y_at_x_right_2_object() (e1->curve(),
                                                           e2->curve(),
                                                           pu);
 
          if (env_type == UPPER)
            u_res = CGAL::opposite (u_res);
        }        
        CGAL_assertion (u_res != EQUAL);

        if (u_res == SMALLER)
          new_v = _append_vertex (out_d, ipt.first, e1);
        else
          new_v = _append_vertex (out_d, ipt.first, e2);
        
        // Note that the new vertex is incident to all curves in e1 and
        // in e2.
        new_v->add_curves (e1->curves_begin(), e1->curves_end());
        new_v->add_curves (e2->curves_begin(), e2->curves_end());

        // Update the handle to the rightmost vertex in the output diagram.
        v_rm = new_v;
        v_rm_exists = true;
      }
      
      if (has_right &&
          (! v_exists ||
           traits->compare_xy_2_object() (p2, 
                                          v->point()) == SMALLER))
      {
        // Create an edge that represents the overlapping curve.
        Vertex_handle        new_v;
        
        new_v = _append_vertex (out_d, ipt.first, e1);
        new_v->left()->add_curves (e2->curves_begin(), e2->curves_end());
        
        new_v->add_curves (e1->curves_begin(), e1->curves_end());
        new_v->add_curves (e2->curves_begin(), e2->curves_end());
        
        // Update the handle to the rightmost vertex in the output diagram.
        v_rm = new_v;
        v_rm_exists = true;

        // Compare the curves to the right of p2.
        u_res = traits->compare_y_at_x_right_2_object() (e1->curve(),
                                                         e2->curve(),
                                                         p2);
        CGAL_assertion (u_res != EQUAL);
      }
      else
      {
        // The overlapping curves reaches v.
        equal_at_v = true;
        u_res = EQUAL;
        break;
      }
      
    }
    
  } // End of the traversal over the intersection objects.
  
  // Handle the portion after the intersection objects.
  if (equal_at_v)
  {
    CGAL_assertion (v_exists);
    
    // In this case the two urves intersect (or overlap) at v.
    Vertex_handle        new_v;
    
    if (u_res == SMALLER)
      new_v = _append_vertex (out_d, v->point(), e1);
    else
      new_v = _append_vertex (out_d, v->point(), e2);
    
    if (u_res == EQUAL)
      new_v->left()->add_curves (e1->curves_begin(), e1->curves_end());
    
    // Note that the new vertex is incident to all curves in e1 and in e2.
    new_v->add_curves (e1->curves_begin(), e1->curves_end());
    new_v->add_curves (e2->curves_begin(), e2->curves_end());
    
    return;
  }
  
  if (! v_exists)
  {
    // Both edges are unbounded from the right, so we simply have
    // to update the rightmost edge in out_d.
    CGAL_assertion (u_res != EQUAL);
    
    if (u_res == SMALLER)
      out_d.rightmost()->add_curves (e1->curves_begin(), e1->curves_end());
    else
      out_d.rightmost()->add_curves (e2->curves_begin(), e2->curves_end());
    
    return;
  }
  
  // Check if we need to insert v into the diagram.
  CGAL_assertion (u_res != EQUAL);
  
  if (u_res == SMALLER)
  {    
    // The final part of the interval is taken from e1.
    Vertex_handle        new_v;

    if (org_v == 1)
    {
      // In case v is also from e1, append it to the merged diagram.
      new_v = _append_vertex (out_d, v->point(), e1);
      new_v->add_curves (v->curves_begin(), v->curves_end());
    }
    else
    {
      // If v is from e2, check if it below (or above, in case of an upper
      // envelope) cv1 to insert it.
      const Comparison_result  res = 
        traits->compare_y_at_x_2_object() (v->point(),
                                           e1->curve());
      
      if (res == EQUAL ||
          (env_type == LOWER && res == SMALLER) ||
          (env_type == UPPER && res == LARGER))
      {
        new_v = _append_vertex (out_d, v->point(), e1);
        new_v->add_curves (v->curves_begin(), v->curves_end());
        
        if (res == EQUAL)
          new_v->add_curves (e1->curves_begin(), e1->curves_end());
      }
    }
  }
  else
  {
    // The final part of the interval is taken from e2.
    Vertex_handle        new_v;

    if (org_v == 2)
    {
      // In case v is also from e2, append it to the merged diagram.
      new_v = _append_vertex (out_d, v->point(), e2);
      new_v->add_curves (v->curves_begin(), v->curves_end());
    }
    else
    {
      // If v is from e1, check if it below (or above, in case of an upper
      // envelope) cv2 to insert it.
      const Comparison_result  res = 
        traits->compare_y_at_x_2_object() (v->point(),
                                           e2->curve());
      
      if (res == EQUAL ||
          (env_type == LOWER && res == SMALLER) ||
          (env_type == UPPER && res == LARGER))
      {
        new_v = _append_vertex (out_d, v->point(), e2);
        new_v->add_curves (v->curves_begin(), v->curves_end());
        
        if (res == EQUAL)
          new_v->add_curves (e2->curves_begin(), e2->curves_end());
      }
    }
  }

  return;
}

// ---------------------------------------------------------------------------
// Append a vertex to the given diagram.
//
template <class Traits, class Diagram>
typename Envelope_divide_and_conquer_2<Traits,Diagram>::Vertex_handle
Envelope_divide_and_conquer_2<Traits,Diagram>::_append_vertex
        (Envelope_diagram_1& diag,
         const Point_2& p, Edge_const_handle e)
{
  // Create the new vertex and the new edge.
  Vertex_handle   new_v = diag.new_vertex (p);
  Edge_handle     new_e = diag.new_edge();
  
  if (! e->is_empty())
    new_e->add_curves (e->curves_begin(), e->curves_end());
  
  // Connect the new vertex.
  new_v->set_left (new_e);
  new_v->set_right (diag.rightmost());
  
  if (diag.leftmost() != diag.rightmost())
  {
    // The diagram is not empty. Connect the new edge to the left of the
    // rightmost edge of the diagram.
    new_e->set_right (new_v);
    new_e->set_left (diag.rightmost()->left());
    diag.rightmost()->left()->set_right (new_e);
    diag.rightmost()->set_left (new_v);
  }
  else
  {
    // The diagram is empty: Make the new edge the leftmost.
    new_e->set_right (new_v);
    diag.set_leftmost (new_e);
    diag.rightmost()->set_left (new_v);      
  }
  
  return (new_v);
}    

// ---------------------------------------------------------------------------
// Merge the vertical segments into the envelope given as a minimization
// (or maximization) diagram.
//
template <class Traits, class Diagram>
void Envelope_divide_and_conquer_2<Traits,Diagram>::
_merge_vertical_segments (Curve_pointer_list& vert_list,
                          Envelope_diagram_1& out_d)
{
    /* RWRW!
    // Sort the vertical segments by their increasing x-coordinate.
    Vertical_strict_weak_ordering vert_order(traits);

    // RWRW: Perform bubble-sort, as the following does not compile under Windows:
    // vert_list.sort (vert_order);
    typename Curve_pointer_list::iterator  i_iter, j_iter;

    for (i_iter = vert_list.begin(); i_iter != vert_list.end(); i_iter++)
    {
      j_iter = i_iter;
      j_iter++;
      for (; j_iter != vert_list.end(); j_iter++)
      {
        if (! vert_order (*i_iter, *j_iter))
        {
          // Swap the pointers.
          const M_curve_2 *temp = *i_iter;
          *i_iter = *j_iter;
          *j_iter = temp;
        }
      }
    }

    // Go over all vertical segments that are to the left of the leftmost
    // vertex of the diagram.
    Curve_pointer_iterator    iter = vert_list.begin();
    M_diagram_vertex_1 *u = NULL;
    M_diagram_vertex_1 *v = out_d.leftmostP;
    Comparison_result  res;
    Point_2            q;

    while (v != NULL)
    {
      while (iter != vert_list.end() &&
             traits->compare_x_2_object() (traits->construct_min_vertex_2_object()((*iter)->xcv),
                                v->p) == SMALLER)
      {
        // Get the lower (or the upper) point of the vertical segment.
        res = traits->compare_xy_2_object() (traits->construct_min_vertex_2_object()((*iter)->xcv),
          traits->construct_max_vertex_2_object()((*iter)->xcv));

        if ((env_type == LOWER && res == SMALLER) ||
          (env_type == UPPER && res == LARGER))
          q = traits->construct_min_vertex_2_object()((*iter)->xcv);
        else
          q = traits->construct_max_vertex_2_object()((*iter)->xcv);

        // Act according to the previous vertex u.
        if (u == NULL)
        {
          // The vertical segment is to the left of the leftmost diagram
          // vertex.
          M_diagram_vertex_1 *new_v1 = vert_alloc.Allocate();
          M_diagram_vertex_1 *new_v2 = vert_alloc.Allocate();
          M_diagram_edge_1   *new_e = edge_alloc.Allocate();
          M_diagram_edge_1   *empty_e = edge_alloc.Allocate();

          new_v1->p = q;
          new_v2->p = q;

          new_v1->leftP = NULL;
          new_v1->rightP = new_e;

          new_e->mcvP = *iter;
          new_e->leftP = new_v1;
          new_e->rightP = new_v2;

          new_v2->leftP = new_e;
          new_v2->rightP = empty_e;

          empty_e->mcvP = NULL;
          empty_e->leftP = new_v2;
          empty_e->rightP = out_d.leftmostP;

          out_d.leftmostP->leftP = empty_e;
          out_d.leftmostP = new_v1;

          // Update the pointer to diagram vertex immediately to the left of v.
          u = new_v2;
        }
        else if (traits->compare_x_2_object() (q, u->p) == EQUAL)
        {
          // The vertical segment has the same x-coordinate as u.
          res = traits->compare_xy_2_object() (q, u->p);

          // Insert a new curve only if it is below (or above, in case of an
          // upper envelope) u->p.
          if ((env_type == LOWER && res == SMALLER) ||
            (env_type == UPPER && res == LARGER))
          {
            M_diagram_edge_1   *new_e = edge_alloc.Allocate();
            M_diagram_vertex_1 *new_v = vert_alloc.Allocate();

            new_v->p = q;

            new_e->mcvP = *iter;
            new_e->leftP = u;
            new_e->rightP = new_v;

            new_v->leftP = new_e;
            new_v->rightP = u->rightP;

            u->rightP->leftP = new_v;
            u->p = q;
            u->rightP = new_e;

            // Update the pointer to diagram vertex immediately to the 
            // left of v.
            u = new_v;
          }
        }
        else
        {
          // The vertical segment is placed in between u and v.
          bool                     add_q = false;
          const X_monotone_curve_2 *cvP = _curve_to_left(v);

          if (cvP == NULL)
          {
            // The edge between u and v is empty:
            add_q = true;
          }
          else
          {
            // Check whether q lies below (or above, in case of an upper
            // envelope) the curves of the edge to the left of u.
            res = traits->compare_y_at_x_2_object() (q, *cvP);

            add_q = (res == EQUAL) ||
                    (env_type == LOWER && res == SMALLER) ||
                    (env_type == UPPER && res == LARGER);
          }

          if (add_q)
          {
            // Cut the edge to the left of v and insert the vertical segment.
            M_diagram_vertex_1 *new_v1 = vert_alloc.Allocate();
            M_diagram_vertex_1 *new_v2 = vert_alloc.Allocate();
            M_diagram_edge_1   *new_e = edge_alloc.Allocate();
            M_diagram_edge_1   *dup_e = edge_alloc.Allocate();

            new_v1->p = q;
            new_v2->p = q;
            *dup_e = *(v->leftP);

            new_v1->leftP = v->leftP;
            new_v1->rightP = new_e;

            new_e->mcvP = *iter;
            new_e->leftP = new_v1;
            new_e->rightP = new_v2;

            new_v2->leftP = new_e;
            new_v2->rightP = dup_e;

            dup_e->leftP = new_v2;
            dup_e->rightP = v;

            v->leftP->rightP = new_v1;
            v->leftP = dup_e;

            // Update the pointer to diagram vertex immediately to the 
            // left of v.
            u = new_v2;
          }
        }

        // Move to the next vertical segment.
        iter++;
      }

      // Move to the next diagram vertex.
      u = v;
      if (v->rightP != NULL)
        v = v->rightP->rightP;
      else
        v = NULL;
    }

    // Deal with all segments located to the right of the diagram.
    while (iter != vert_list.end())
    {
      // Get the lower (or the upper) point of the vertical segment.
      res = traits->compare_xy_2_object() (traits->construct_min_vertex_2_object()((*iter)->xcv),
        traits->construct_max_vertex_2_object()((*iter)->xcv));

      if ((env_type == LOWER && res == SMALLER) ||
        (env_type == UPPER && res == LARGER))
        q = traits->construct_min_vertex_2_object()((*iter)->xcv);
      else
        q = traits->construct_max_vertex_2_object()((*iter)->xcv);

      // The vertical segment is to the right of the rightmost diagram vertex.
      M_diagram_vertex_1 *new_v1 = vert_alloc.Allocate();
      M_diagram_vertex_1 *new_v2 = vert_alloc.Allocate();
      M_diagram_edge_1   *new_e = edge_alloc.Allocate();
      M_diagram_edge_1   *empty_e = edge_alloc.Allocate();

      new_v1->p = q;
      new_v2->p = q;

      empty_e->mcvP = NULL;
      empty_e->leftP = out_d.rightmostP;
      empty_e->rightP = new_v1;

      new_v1->leftP = empty_e;
      new_v1->rightP = new_e;

      new_e->mcvP = *iter;
      new_e->leftP = new_v1;
      new_e->rightP = new_v2;

      new_v2->leftP = new_e;
      new_v2->rightP = NULL;

      out_d.rightmostP->rightP = empty_e;
      out_d.rightmostP = new_v2;

      // Move to the next vertical segment.
      iter++;
    }
    */

  return;
}

CGAL_END_NAMESPACE

#endif

