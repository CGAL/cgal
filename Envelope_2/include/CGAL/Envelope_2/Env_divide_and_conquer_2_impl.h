// Copyright (c) 2006  Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein   <wein@post.tau.ac.il>

#ifndef CGAL_ENVELOPE_DIVIDE_AND_CONQUER_2_IMPL_H
#define CGAL_ENVELOPE_DIVIDE_AND_CONQUER_2_IMPL_H

/*! \file
 * Definitions of the functions of the Envelope_divide_and_conquer_2 class.
 */

#include <boost/optional.hpp>

namespace CGAL {

// ---------------------------------------------------------------------------
// Construct the lower/upper envelope of the given list of non-vertical curves.
//
template <class Traits, class Diagram>
void Envelope_divide_and_conquer_2<Traits,Diagram>::
_construct_envelope_non_vertical(Curve_pointer_iterator begin,
                                 Curve_pointer_iterator end,
                                 Envelope_diagram_1& out_d)
{
  out_d.clear();
  
  if (begin == end)
    return;
  
  // Check if the range contains just a single curve.
  Curve_pointer_iterator    iter = begin;
  ++iter;
  
  if (iter == end)
  {
    // Construct a singleton diagram, which matches a single curve.
    _construct_singleton_diagram(*(*begin), out_d);
  }
  else
  {
    // Divide the given range of curves into two.
    std::size_t size = std::distance(begin, end);
    Curve_pointer_iterator  div_it = begin;
    std::advance(div_it, size / 2);
    
    // Construct the diagrams (envelopes) for the two sub-ranges recursively 
    // and then merge the two diagrams to obtain the result.
    Envelope_diagram_1   d1;
    Envelope_diagram_1   d2;
    
    _construct_envelope_non_vertical(begin, div_it, d1);
    
    _construct_envelope_non_vertical(div_it, end, d2);

    _merge_envelopes(d1, d2, out_d);

    // Print the minimization diagram.
    /* RWRW:
    Edge_const_handle    e = out_d.leftmost();
    Vertex_const_handle  v;

    std::cout << "The diagram: ";
    while (e != out_d.rightmost())
    {
      if (! e->is_empty())
        std::cout << e->curve() << "  ";
      else
        std::cout << "[empty]" << "  ";

      v = e->right();
      std::cout << "(" << v->point() << ")  ";
      
      e = v->right();
    }
    std::cout << "[empty]" << std::endl;
    */
  }
  
  return;
}

// ---------------------------------------------------------------------------
// Construct a singleton diagram, which matches a single curve.
//
template <class Traits, class Diagram>
void Envelope_divide_and_conquer_2<Traits,Diagram>::
_construct_singleton_diagram(const X_monotone_curve_2& cv,
                             Envelope_diagram_1& out_d)
{
  CGAL_assertion(out_d.leftmost() == out_d.rightmost());
  CGAL_assertion(out_d.leftmost()->is_empty());
  
  // Check if the given curve is bounded from the left and from the right.
  if (traits->parameter_space_in_x_2_object()(cv, ARR_MIN_END) != ARR_INTERIOR)
  {
    if (traits->parameter_space_in_x_2_object()(cv, ARR_MAX_END) != ARR_INTERIOR)
    {
      // The curve is defined over (-oo, oo), so its diagram contains
      // only a single edge.
      out_d.leftmost()->add_curve(cv);
      
      return;
    }

    // The curve is defined over (-oo, x], where x is finite.
    // Create a vertex and associate it with the right endpoint of cv.
    CGAL_precondition
      (traits->parameter_space_in_y_2_object()(cv, ARR_MAX_END) == ARR_INTERIOR);
    
    Vertex_handle v =
      out_d.new_vertex(traits->construct_max_vertex_2_object()(cv));
    Edge_handle   e_right = out_d.new_edge();
    
    v->add_curve(cv);
    v->set_left(out_d.leftmost());
    v->set_right(e_right);
    
    // The leftmost edge is associated with cv, and the rightmost is empty.
    out_d.leftmost()->add_curve(cv);
    out_d.leftmost()->set_right(v);
    
    e_right->set_left(v);
    out_d.set_rightmost(e_right);
    
    return;
  }
  
  if (traits->parameter_space_in_x_2_object()(cv, ARR_MAX_END) != ARR_INTERIOR)
  {
    // The curve is defined over [x, +oo), where x is finite.
    // Create a vertex and associate it with the left endpoint of cv.
    CGAL_precondition
      (traits->parameter_space_in_y_2_object()(cv, ARR_MIN_END) == ARR_INTERIOR);
    
    Vertex_handle  v = 
      out_d.new_vertex(traits->construct_min_vertex_2_object()(cv));
    Edge_handle    e_left = out_d.new_edge();
    
    v->add_curve(cv);
    v->set_left(e_left);
    v->set_right(out_d.rightmost());
    
    // The rightmost edge is associated with cv, and the leftmost is empty.
    out_d.rightmost()->add_curve(cv);
    out_d.rightmost()->set_left(v);
    
    e_left->set_right(v);
    out_d.set_leftmost(e_left);
    
    return;
  }
  
  // If we reached here, the curve is defined over a bounded x-range.
  // We therefore create the following diagram:
  //
  //             (empty)    v1     e       v2   (empty)
  //      -oo -------------(+)============(+)------------ +oo
  //
  CGAL_precondition
    (traits->parameter_space_in_y_2_object()(cv, ARR_MIN_END) == ARR_INTERIOR);
  CGAL_precondition
    (traits->parameter_space_in_y_2_object()(cv, ARR_MAX_END) == ARR_INTERIOR);
  
  Vertex_handle  v1 = 
    out_d.new_vertex(traits->construct_min_vertex_2_object()(cv));
  Vertex_handle  v2 = 
    out_d.new_vertex(traits->construct_max_vertex_2_object()(cv));
  Edge_handle    e_left = out_d.new_edge();
  Edge_handle    e_right = out_d.new_edge();
  Edge_handle    e = out_d.leftmost();
  
  v1->add_curve(cv);
  v1->set_left(e_left);
  v1->set_right(e);
  
  v2->add_curve(cv);
  v2->set_left(e);
  v2->set_right(e_right);
  
  e->add_curve(cv);
  e->set_left(v1);
  e->set_right(v2);
  
  e_left->set_right(v1);
  e_right->set_left(v2);
  
  out_d.set_leftmost(e_left);
  out_d.set_rightmost(e_right);
  
  return;
}

// ---------------------------------------------------------------------------
// Merge two minimization (or maximization) diagrams.
//
template <class Traits, class Diagram>
void Envelope_divide_and_conquer_2<Traits,Diagram>::
_merge_envelopes(const Envelope_diagram_1& d1,
                 const Envelope_diagram_1& d2,
                 Envelope_diagram_1& out_d)
{
  Edge_const_handle    e1 = d1.leftmost();
  bool                 is_leftmost1 = true;
  Vertex_const_handle  v1 = Vertex_const_handle();
  Edge_const_handle    e2 = d2.leftmost();
  bool                 is_leftmost2 = true;
  Vertex_const_handle  v2 = Vertex_const_handle();
  Vertex_const_handle  next_v = Vertex_const_handle();
  bool                 next_exists = true;
  Comparison_result    res_v = EQUAL;
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
      res_v = _compare_vertices(v1, v2, same_x);
      next_v = (res_v == SMALLER) ? v1 : v2;
    }
    
    // Check if the current edges represent empty intervals or not.
    if (! e1->is_empty() && ! e2->is_empty())
    {
      // Both edges are not empty, and there are curves defined on them.
      _merge_two_intervals(e1, is_leftmost1, e2, is_leftmost2,
                           next_v, next_exists, res_v, out_d);
      
    }
    else if (! e1->is_empty() && e2->is_empty())
    {
      // e1 is not empty but e2 is empty:
      _merge_single_interval(e1, e2, next_v, next_exists, res_v, out_d);
    }
    else if (e1->is_empty() && ! e2->is_empty())
    {
      // e1 is empty and e2 is not empty:
      _merge_single_interval(e2, e1, next_v, next_exists,
                             CGAL::opposite(res_v), out_d);
    }
    else
    {
      // Both edges are empty: append an empty edge to out_d:
      if (next_exists)
      {
        Vertex_handle  new_v = _append_vertex(out_d, next_v->point(), e1);
        switch(res_v)
        {
        case SMALLER:
          new_v->add_curves(v1->curves_begin(), v1->curves_end()); break;
        case LARGER:
          new_v->add_curves(v2->curves_begin(), v2->curves_end()); break;
        case EQUAL:
          new_v->add_curves(v1->curves_begin(), v1->curves_end());
          new_v->add_curves(v2->curves_begin(), v2->curves_end());
          break;
        }
      }
    }
    
    // Proceed to the next diagram edge(s), if possible.
    if (next_exists)
    {
      // Check if we should proceed on d1 or on d2.
      // \todo: we do not need 3 cases, only two.
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
_compare_vertices(Vertex_const_handle v1,
                  Vertex_const_handle v2,
                  bool& same_x) const
{
  Comparison_result res =
    traits->compare_x_2_object()(v1->point(), v2->point());
  
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
  res = traits->compare_xy_2_object()(v1->point(), v2->point());

  // In case of upper envlope we take the opposite result
  if (env_type == UPPER)
    return CGAL::opposite(res);
  return res;
}

// ---------------------------------------------------------------------------
// Deal with an interval which is non-empty in one of the merged diagrams and
// empty in the other.
//
template <class Traits, class Diagram>
void Envelope_divide_and_conquer_2<Traits,Diagram>::
_merge_single_interval(Edge_const_handle e, Edge_const_handle other_edge,
                       Vertex_const_handle v, bool v_exists,
                       Comparison_result origin_of_v,
                       Envelope_diagram_1& out_d)
{
  if (! v_exists)
  {
    // The non-empty edge e is unbounded from the right, so we simply have
    // to update the rightmost edge in out_d.
    out_d.rightmost()->add_curves(e->curves_begin(), e->curves_end());
    return;
  }
  
  Vertex_handle      new_v;
  
  if (origin_of_v == SMALLER)
  {
    // The non-empty edge ends at v, so we simply insert it to out_d.
    new_v = _append_vertex(out_d, v->point(), e);
    new_v->add_curves(v->curves_begin(), v->curves_end());
    
    return;
  }

  if (origin_of_v == EQUAL) // the edges have vertices at the same place.
  {
    new_v = _append_vertex(out_d, v->point(), e);
    new_v->add_curves(e->right()->curves_begin(), e->right()->curves_end());
    new_v->add_curves(other_edge->right()->curves_begin(), 
                      other_edge->right()->curves_end());
    return;
  }

  // If v is not on e, we should insert it to the merged diagram only if it
  // is below (or above, in case of an upper envelope) the curves of e.
  Comparison_result res =
    traits->compare_y_at_x_2_object()(v->point(), e->curve());
  
  if ((res == EQUAL) ||
      (env_type == LOWER && res == SMALLER) ||
      (env_type == UPPER && res == LARGER))
  {
    new_v = _append_vertex(out_d, v->point(), e);
    new_v->add_curves(v->curves_begin(), v->curves_end());

    if (res == EQUAL)
      {
        // In case of equality, append e's curves to those of the new vertex.
        new_v->add_curves(e->curves_begin(), e->curves_end());
      }
  }
}

//! \brief Functions that should be on Arr_traits_adaptor.
/*@{*/

//! Compare the $y$-coordinates of two curves at their endpoints
/*! The function compares the $y$ values of two curves with a joint 
  range of $x$ values, at the end of the joint range.
  \param xcv1 The first curve
  \param xcv2 The second curve
  \param curve_end ARR_MIN_END - compare the $y$ value of the smaller 
  endpoint, ARR_MAX_END - compare the $y$ value of the larger endpoint.
  \pre The two $x$-monotone curves need to have a partially overlapping 
  $x$-ranges.
  \return 
  \todo Move it to Arr_traits_adaptor ?
*/
template <class Traits, class Diagram>
Comparison_result Envelope_divide_and_conquer_2<Traits,Diagram>::
compare_y_at_end(const X_monotone_curve_2& xcv1,
                 const X_monotone_curve_2& xcv2,
                 Arr_curve_end curve_end) const
{
  CGAL_precondition(traits->is_in_x_range_2_object()(xcv1, xcv2));

  typedef typename Traits::Compare_xy_2               Compare_xy_2;
  typedef typename Traits::Compare_y_at_x_2           Compare_y_at_x_2;
  typedef typename Traits::Construct_min_vertex_2     Construct_min_vertex_2;
  typedef typename Traits::Construct_max_vertex_2     Construct_max_vertex_2;

  Compare_y_at_x_2  compare_y_at_x = traits->compare_y_at_x_2_object();
  Construct_min_vertex_2 min_vertex =
    traits->construct_min_vertex_2_object();
  Construct_max_vertex_2 max_vertex =
    traits->construct_max_vertex_2_object();
  
  // First check whether any of the curves is defined at x boundary.
  const Arr_parameter_space ps_x1 =
    traits->parameter_space_in_x_2_object()(xcv1, curve_end);
  const Arr_parameter_space ps_x2 =
    traits->parameter_space_in_x_2_object()(xcv2, curve_end);
  Comparison_result         res;
  
  if (ps_x1 != ARR_INTERIOR) {
    if (ps_x2 != ARR_INTERIOR) {
      // Compare the relative position of the curves at x boundary.
      return (traits->compare_y_near_boundary_2_object()(xcv1, xcv2,
                                                         curve_end));
    }
    
    // Check if the left end of xcv2 lies at y boundary.
    const Arr_parameter_space ps_y2 =
      traits->parameter_space_in_y_2_object()(xcv2, curve_end);
    
    if (ps_y2 == ARR_BOTTOM_BOUNDARY)
      return (LARGER);          // xcv2 is obviously below xcv1.
    else if (ps_y2 == ARR_TOP_BOUNDARY)
      return (SMALLER);         // xcv2 is obviously above xcv1.
          
    // Compare the position of the left end of xcv2 (which is a normal
    // point) to xcv1.
    res = (curve_end == ARR_MIN_END) ?
      compare_y_at_x(min_vertex(xcv2), xcv1) :
      compare_y_at_x(max_vertex(xcv2), xcv1);

    // Swap the result.
    return CGAL::opposite(res);
  }
  else if (ps_x2 != ARR_INTERIOR) {
    // Check if the left end of xcv1 lies at y boundary.
    const Arr_parameter_space ps_y1 =  traits->parameter_space_in_y_2_object()
      (xcv1, curve_end);
    
    if (ps_y1 == ARR_BOTTOM_BOUNDARY)
      return (SMALLER);         // xcv1 is obviously below xcv2.
    else if (ps_y1 == ARR_TOP_BOUNDARY)
      return (LARGER);          // xcv1 is obviously above xcv2.
    
    // Compare the position of the left end of xcv1 (which is a normal
    // point) to xcv2.
    res = (curve_end == ARR_MIN_END) ?
      compare_y_at_x(min_vertex(xcv1), xcv2) :
      compare_y_at_x(max_vertex(xcv1), xcv2);
    return (res);
  }
  
  // Check if the left curve end lies at y = +/- oo.
  const Arr_parameter_space ps_y1 =
    traits->parameter_space_in_y_2_object()(xcv1, curve_end);
  const Arr_parameter_space ps_y2 =
    traits->parameter_space_in_y_2_object()(xcv2, curve_end);
  Comparison_result         l_res;
  
  if (ps_y1 != ARR_INTERIOR) {
    if (ps_y2 != ARR_INTERIOR) {
      // The curve ends have boundary conditions with oposite signs in y,
      // we readily know their relative position (recall that they do not
      // instersect).
      if ((ps_y1 == ARR_BOTTOM_BOUNDARY) && (ps_y2 == ARR_TOP_BOUNDARY))
        return (SMALLER);
      else if ((ps_y1 == ARR_TOP_BOUNDARY) && (ps_y2 == ARR_BOTTOM_BOUNDARY))
        return (LARGER);

      // Both curves have vertical asymptotes with the same sign in y.
      // Check which asymptote is the rightmost. Note that in this case
      // the vertical asymptotes cannot be equal.
      l_res = traits->compare_x_curve_ends_2_object()(xcv1, curve_end,
                                                      xcv2, curve_end);
      CGAL_assertion(l_res != EQUAL);
      
      if (ps_y1 == ARR_TOP_BOUNDARY)
        return (l_res);
      else
        return CGAL::opposite(l_res);
    }

    // xcv1 has a vertical asymptote and xcv2 has a normal left endpoint.
    // Compare the x-positions of this endpoint and the asymptote.
    const Point_2& left2 =
      (curve_end == ARR_MIN_END) ? min_vertex(xcv2) : max_vertex(xcv2);
        
    l_res =
      traits->compare_x_point_curve_end_2_object()(left2, xcv1, curve_end);
    
    if (l_res == LARGER) {
      // left2 lies in the x-range of xcv1, so it is safe to compare:
      res = compare_y_at_x(left2, xcv1);
      return CGAL::opposite(res);
    }
    else
      // xcv1 is below or above xcv2.
      return ((ps_y1 == ARR_BOTTOM_BOUNDARY) ? SMALLER : LARGER);
  }
  else if (ps_y2 != ARR_INTERIOR) {
    // xcv2 has a vertical asymptote and xcv1 has a normal left endpoint.
    // Compare the x-positions of this endpoint and the asymptote.
    const Point_2&  left1 = 
      (curve_end == ARR_MIN_END) ? min_vertex(xcv1) : max_vertex(xcv1);
        
    l_res =
      traits->compare_x_point_curve_end_2_object()(left1, xcv2, curve_end);
    
    return ((l_res == LARGER) ?
            // left1 lies in the x-range of xcv2, so it is safe to compare:
            (compare_y_at_x(left1, xcv2)) :
            ((ps_y2 == ARR_BOTTOM_BOUNDARY) ? LARGER : SMALLER));
  }

  // In this case we compare two normal points.
  Compare_xy_2            compare_xy = traits->compare_xy_2_object();

  // Get the left endpoints of xcv1 and xcv2.
  const Point_2&  left1 = 
    (curve_end == ARR_MIN_END) ? min_vertex(xcv1) : max_vertex(xcv1);
  const Point_2&  left2 = 
    (curve_end == ARR_MIN_END) ? min_vertex(xcv2) : max_vertex(xcv2);

  // Locate the rightmost point of left1 and left2 and compare its position
  // to the other curve.
  l_res = compare_xy(left1, left2);
  
  return ((l_res != SMALLER) ?
    // left1 is in the x-range of xcv2:
          compare_y_at_x(left1, xcv2) : 
    // left2 is in the x-range of xcv1:
          CGAL::opposite(compare_y_at_x(left2, xcv1)));
}
/*@}*/

// ---------------------------------------------------------------------------
// Merge two non-empty intervals into the merged diagram.
//
template <class Traits, class Diagram>
void Envelope_divide_and_conquer_2<Traits,Diagram>::
_merge_two_intervals(Edge_const_handle e1, bool is_leftmost1, 
                     Edge_const_handle e2, bool is_leftmost2,
                     Vertex_const_handle v, bool v_exists,
                     Comparison_result origin_of_v,
                     Envelope_diagram_1& out_d)
{
  typedef std::pair<Point_2, typename Traits::Multiplicity>  Intersection_point;

  Comparison_result                current_res;
  bool                             equal_at_v = false;
 
  // Get the relative position of two curves associated with e1 and e2
  // at the rightmost of the left endpoints of e1 and e2.
  current_res = compare_y_at_end(e1->curve(), e2->curve(), ARR_MIN_END);
  // Flip the result in case of an upper envelope.
  if (env_type == UPPER)
    current_res = CGAL::opposite(current_res);

  // Use the current rightmost of the two left vertices as a reference point.
  // This is the rightmost vertex in the current minimization diagram (out_d).
  // The intersection points/curves that interest us are the ones in
  // [v_leftmost, v].
  boost::optional<Vertex_const_handle>  v_leftmost =
    boost::optional<Vertex_const_handle>();

  if (is_leftmost1 == true) {
    if (is_leftmost2 == false)
      v_leftmost = e2->left();
  }
  else {
    if (is_leftmost2 == true)
      v_leftmost = e1->left();
    else
    {
      if ((traits->compare_xy_2_object()(e1->left()->point(),
                                         e2->left()->point()) == LARGER))
        v_leftmost = e1->left();
      else
        v_leftmost = e2->left();
    }
  }

  // Find the next intersection of the envelopes to the right of the current
  // rightmost point in the merged diagram.
  // \todo Use the faster object_cast.
  std::list<CGAL::Object>           objects;
  CGAL::Object                      obj;
  const X_monotone_curve_2*         intersection_curve;
  const Intersection_point*         intersection_point;
  
  traits->intersect_2_object()(e1->curve(), e2->curve(),
                               std::back_inserter(objects));
  
  while (! objects.empty()) {
    // Pop the xy-lexicographically smallest intersection object.
    obj = objects.front();
    objects.pop_front();
    
    if ((intersection_point = CGAL::object_cast<Intersection_point>(&obj)) !=
        NULL)
    {
      // We have a simple intersection point.
      bool is_in_x_range = true; // true if the intersection point is to the 
                                 // right of v_leftmost.
      // check if we are before the leftmost point.
      if (v_leftmost &&
          traits->compare_xy_2_object() (intersection_point->first, 
                                         (*v_leftmost)->point()) != LARGER)
      {
        // The point is to the left of the current rightmost vertex in out_d,
        // so we skip it and continue examining the next intersections.
        // However, we update the last intersection point observed.
        is_in_x_range = false;
      }
      
      // check if we arrived at the rightmost point (stop if indeed we are
      // there).
      if (is_in_x_range && v_exists) {
        Comparison_result res = traits->compare_xy_2_object() 
          (intersection_point->first, v->point());
        
        // v is an intersection points, so both curves are equal there:
        if (res == EQUAL)
          equal_at_v = true;
        
        // We passed the next vertex, so we can stop here.
        if (res == LARGER)
          break;
      }
      
      // Create a new vertex in the output diagram that corrsponds to the
      // current intersection point.
      if (is_in_x_range) {
        CGAL_assertion(current_res != EQUAL);

        Vertex_handle new_v = (current_res == SMALLER) ?
          _append_vertex(out_d, intersection_point->first, e1) :
          _append_vertex(out_d, intersection_point->first, e2);
        
        // if we are at v, then this is a special case that is handled after
        // the loop. We need to add the curves from the original vertices.
        if (equal_at_v == false) {
          // Note that the new vertex is incident to all curves in e1 and in e2.
          new_v->add_curves(e1->curves_begin(), e1->curves_end());
          new_v->add_curves(e2->curves_begin(), e2->curves_end());
        }
        
        // Update the handle to the rightmost vertex in the output diagram.
        v_leftmost = new_v;
      }

      // Update the relative position of the two curves, which is their
      // order immediately to the right of their last observed intersection
      // point.
      
      // Get the curve order immediately to the right of the intersection
      // point. Note that in case of even (non-zero) multiplicity the order
      // remains the same.
      if ((current_res != EQUAL) && (intersection_point->second % 2 == 1)) {
        // Odd multiplicity: flip the current comparison result.
        current_res = CGAL::opposite(current_res);
      }
      // if we are at v, then we may not be able to call compare_y_at_x_right_2.
      else if (((intersection_point->second == 0) || (current_res == EQUAL)) &&
               (equal_at_v == false))
      {
        // The multiplicity is unknown, so we have to compare the curves to
        // the right of their intersection point.
        current_res =
          traits->compare_y_at_x_right_2_object()(e1->curve(), e2->curve(),
                                                  intersection_point->first);
        // Flip the result in case of an upper envelope.
        if (env_type == UPPER)
          current_res = CGAL::opposite (current_res);
      }
    }
    else {
      // We have an x-monotone curve representing an overlap of the two
      // curves.
      intersection_curve = CGAL::object_cast<X_monotone_curve_2>(&obj);
      
      if (intersection_curve == NULL)
        CGAL_error_msg("unrecognized intersection object.");

      // Get the endpoints of the overlapping curves.
      const bool  has_left = 
        (traits->parameter_space_in_x_2_object() 
         (*intersection_curve, ARR_MIN_END) == ARR_INTERIOR);
      const bool  has_right = 
        (traits->parameter_space_in_x_2_object() 
         (*intersection_curve, ARR_MAX_END) == ARR_INTERIOR);
      Point_2     p_left, p_right;
      
      if (has_left)
        p_left = traits->construct_min_vertex_2_object()(*intersection_curve);
      
      if (has_right)
        p_right = traits->construct_max_vertex_2_object()(*intersection_curve);

      bool is_in_x_range = true;
      // Check if the overlapping curve is not relevant to our range.
      if (v_leftmost && has_right &&
          (traits->compare_xy_2_object()(p_right, (*v_leftmost)->point()) !=
           LARGER))
      {
        // The right point of the overlappinf curve is to the left of the
        // current rightmost vertex in out_d, so we skip it and continue
        // examining the next intersections.
        // However, we update the last intersection point observed.
        is_in_x_range = false;
      }
      
      if (is_in_x_range && v_exists && has_left) {
        Comparison_result res =
          traits->compare_xy_2_object()(p_left, v->point());
        
        // v is an intersection points, so both curves are equal there:
        if (res == EQUAL)
          equal_at_v = true;
        
        // We passed the next vertex, so we can stop here.
        if (res == LARGER)
          break;
      }
      
      // There is an overlap between the range [u, v] and intersection_curve.
      if (is_in_x_range && has_left && 
          (! v_leftmost ||
           (traits->compare_xy_2_object()(p_left, (*v_leftmost)->point()) ==
            LARGER)))
      {
        // Create an output edge that represent the portion of [u, v] to the
        // left of the overlapping curve.
        CGAL_assertion(current_res != EQUAL);

        Vertex_handle new_v = (current_res == SMALLER) ?
          _append_vertex(out_d, p_left, e1) :
          _append_vertex(out_d, p_left, e2);
        
        // if we are at v, then this is a special case that is handled after
        // the loop. We need to add the curves from the original vertices.
        if (equal_at_v == false) {
          // Note that the new vertex is incident to all curves in e1 and
          // in e2.
          new_v->add_curves(e1->curves_begin(), e1->curves_end());
          new_v->add_curves(e2->curves_begin(), e2->curves_end());
        }
        
        // Update the handle to the rightmost vertex in the output diagram.
        v_leftmost = new_v;
      }
      
      if (is_in_x_range && has_right &&
          (! v_exists ||
           (traits->compare_xy_2_object()(p_right, v->point()) == SMALLER)))
      {
        // Create an edge that represents the overlapping curve.
        Vertex_handle new_v = _append_vertex(out_d, p_right, e1);
        new_v->left()->add_curves(e2->curves_begin(), e2->curves_end());
        
        // We are not at v becuase p_right is smaller than v.
        // The special case that we are at v is handled in the next
        // condition.
        // If we were at v, then this was a special case that is handled
        // later.
        CGAL_assertion(equal_at_v == false);
        new_v->add_curves(e1->curves_begin(), e1->curves_end());
        new_v->add_curves(e2->curves_begin(), e2->curves_end());
        
        // Update the handle to the rightmost vertex in the output diagram.
        v_leftmost = new_v;
      }
      
      if (has_right == false || 
          (v_exists && traits->compare_xy_2_object()(p_right, v->point()) !=
           SMALLER))
      {
        // The overlapping curves reaches v.
        if (v_exists) {
          Vertex_handle new_v = _append_vertex(out_d, v->point(), e1);
          new_v->left()->add_curves(e2->curves_begin(), e2->curves_end());
        }
        
        equal_at_v = true;
        current_res = EQUAL;
        break;
      }

      // We arrive here only if we are not at v and the overlap has a right
      // endpoint.
      // Compare the curves to the right of p_right.
      current_res =
        traits->compare_y_at_x_right_2_object()(e1->curve(), e2->curve(),
                                                p_right);
      // Flip the result in case of an upper envelope.
      if (env_type == UPPER)
        current_res = CGAL::opposite (current_res);
      
    }
  } // End of the traversal over the intersection objects.


  // Handle the portion after the intersection objects.
  if (equal_at_v) {
    CGAL_assertion (v_exists);

    // v_leftmost should be our vertex at v.
    // In this case the two curves intersect (or overlap) at v.
    // We need to add the correct curves to v_leftmost.

    Vertex_handle v_to_be_updated = out_d.rightmost()->left();

    if (origin_of_v == EQUAL) {
      // If the vertices of the edge are the same, we have to get the 
      // curves from there:
      v_to_be_updated->add_curves(e1->right()->curves_begin(), 
                                  e1->right()->curves_end());
      v_to_be_updated->add_curves(e2->right()->curves_begin(), 
                                  e2->right()->curves_end());
    }
    else {
      // We add the curves from the original vertex and from the edge of the
      // other diagram.
      Edge_const_handle e = (origin_of_v == SMALLER) ? e2 : e1;

      v_to_be_updated->add_curves (v->curves_begin(), v->curves_end());
      v_to_be_updated->add_curves (e->curves_begin(), e->curves_end());
    }
    
    return;
  }
  
  if (! v_exists) {
    // Both edges are unbounded from the right, so we simply have
    // to update the rightmost edge in out_d.
    switch (current_res)
    {
    case SMALLER:
      out_d.rightmost()->add_curves(e1->curves_begin(), e1->curves_end());
      return;
    case LARGER:
      out_d.rightmost()->add_curves(e2->curves_begin(), e2->curves_end());
      return;
    case EQUAL:
      out_d.rightmost()->add_curves(e1->curves_begin(), e1->curves_end());
      out_d.rightmost()->add_curves(e2->curves_begin(), e2->curves_end());
      return;
    default:
      CGAL_error_msg("should not reach here.");
      return;
    }
  }

  // origin_of_v could be EQUAL but the curves do not intersect.
  // This is because of the fact that v could be the endpoint of the NEXT
  // curve (which is lower than the currrent curve. The second diagram, however,
  // has a curve that ends at v.
  // For example:
  // First diagram is the segment: [(0, -1), (1, 0)]
  // Second diagram of the two segments: [(0, 0), (1, 1)], [(1, 0), (2, 1)]
  
  // Check if we need to insert v into the diagram.
  if (current_res == SMALLER) {    
    // The final part of the interval is taken from e1.
    Vertex_handle        new_v;

    if (origin_of_v == SMALLER) {
      // In case v is also from e1, append it to the merged diagram.
      new_v = _append_vertex(out_d, v->point(), e1);
      new_v->add_curves(v->curves_begin(), v->curves_end());
    }
    else {
      // if origin_of_v is EQUAL then the two diagram have a vertex at
      // exact same place.
      if (origin_of_v == EQUAL) {
        new_v = _append_vertex(out_d, v->point(), e1);
        new_v->add_curves(v->curves_begin(), v->curves_end());
        
        // adding the curves of the vertex of the first diagram (vertices are
        // equal...)
        new_v->add_curves(e1->right()->curves_begin(), 
                          e1->right()->curves_end());
      }
      else {
        // If v is from e2, check if it below (or above, in case of an upper
        // envelope) cv1 to insert it.
        const Comparison_result  res = 
          traits->compare_y_at_x_2_object()(v->point(), e1->curve());
        
        if (res == EQUAL ||
            ((env_type == LOWER) && (res == SMALLER)) ||
            ((env_type == UPPER) && (res == LARGER)))
          {
            new_v = _append_vertex(out_d, v->point(), e1);
            new_v->add_curves(v->curves_begin(), v->curves_end());
            
            if (res == EQUAL)
              new_v->add_curves(e1->curves_begin(), e1->curves_end());
          }
      }
    }
  }
  else {
    // The final part of the interval is taken from e2.
    Vertex_handle        new_v;

    if (origin_of_v != SMALLER) {
      // In case v is also from e2, append it to the merged diagram.
      new_v = _append_vertex(out_d, v->point(), e2);
      new_v->add_curves(v->curves_begin(), v->curves_end());
      
      // if origin_of_v is EQUAL then the two diagram have a vertex at
      // exact same place.
      if (origin_of_v == EQUAL) {
        // adding the curves of the vertex of the first diagram (vertices are
        // equal...)
        new_v->add_curves(e1->right()->curves_begin(),
                          e1->right()->curves_end());
      }
    }
    else {
      // If v is from e1, check if it below (or above, in case of an upper
      // envelope) cv2 to insert it.
      const Comparison_result  res = 
        traits->compare_y_at_x_2_object()(v->point(), e2->curve());
      
      if (res == EQUAL ||
          ((env_type == LOWER) && (res == SMALLER)) ||
          ((env_type == UPPER) && (res == LARGER)))
      {
        new_v = _append_vertex(out_d, v->point(), e2);
        new_v->add_curves(v->curves_begin(), v->curves_end());
        
        if (res == EQUAL)
          new_v->add_curves(e2->curves_begin(), e2->curves_end());
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
Envelope_divide_and_conquer_2<Traits,Diagram>::
_append_vertex(Envelope_diagram_1& diag,
               const Point_2& p, Edge_const_handle e)
{
  // Create the new vertex and the new edge.
  Vertex_handle   new_v = diag.new_vertex(p);
  Edge_handle     new_e = diag.new_edge();
  
  if (! e->is_empty())
    new_e->add_curves(e->curves_begin(), e->curves_end());
  
  // Connect the new vertex.
  new_v->set_left(new_e);
  new_v->set_right(diag.rightmost());
  
  if (diag.leftmost() != diag.rightmost()) {
    // The diagram is not empty. Connect the new edge to the left of the
    // rightmost edge of the diagram.
    new_e->set_right(new_v);
    new_e->set_left(diag.rightmost()->left());
    diag.rightmost()->left()->set_right(new_e);
    diag.rightmost()->set_left(new_v);
  }
  else {
    // The diagram is empty: Make the new edge the leftmost.
    new_e->set_right(new_v);
    diag.set_leftmost(new_e);
    diag.rightmost()->set_left(new_v);      
  }
  
  return (new_v);
}    

// ---------------------------------------------------------------------------
// Merge the vertical segments into the envelope given as a minimization
// (or maximization) diagram.
//
template <class Traits, class Diagram>
void Envelope_divide_and_conquer_2<Traits,Diagram>::
_merge_vertical_segments(Curve_pointer_vector& vert_vec,
                         Envelope_diagram_1& out_d)
{
  // Sort the vertical segments by their increasing x-coordinate.
  Less_vertical_segment  les_vert(traits);

  std::sort(vert_vec.begin(), vert_vec.end(), les_vert);

  // Proceed on the diagram and on the sorted sequence of vertical segments
  // and merge them into the diagram.
  typename Traits_adaptor_2::Compare_x_2             comp_x =
    traits->compare_x_2_object();
  typename Traits_adaptor_2::Compare_xy_2            comp_xy =
    traits->compare_xy_2_object();
  typename Traits_adaptor_2::Compare_y_at_x_2        comp_y_at_x =
    traits->compare_y_at_x_2_object();
  typename Traits_adaptor_2::Construct_min_vertex_2  min_vertex =
    traits->construct_min_vertex_2_object();
  typename Traits_adaptor_2::Construct_max_vertex_2  max_vertex =
    traits->construct_max_vertex_2_object();

  Edge_handle             e = out_d.leftmost();
  Vertex_handle           v = Vertex_handle();
  Curve_pointer_iterator  iter = vert_vec.begin();
  Curve_pointer_iterator  next;
  Comparison_result       res;
  bool                    in_e_range;
  bool                    on_v;
  Point_2                 p;

  while (iter != vert_vec.end()) {
    // Check if the current vertical segment is on the x-range of the current
    // edge.
    if (e != out_d.rightmost()) {
      // The current edge is not the rightmost one: we compare the x-coordinate
      // of the vertical segment to its right vertex.
      v = e->right();

      res = comp_x(min_vertex(**iter), v->point());
      in_e_range = (res != LARGER);
      on_v = (res == EQUAL);
    }
    else {
      // This is the rightmost edge, so the vertical segment must lie on its
      // x-range.
      in_e_range = true;
      on_v = false;
    }
    
    // If the current vertical segment is not in the x-range of the current
    // edge, we proceed to the next edge.
    if (! in_e_range) {
      e = v->right();
      continue;
    }

    // Go over all vertical segments that share the same x-coordinate and
    // find the one(s) with the smallest endpoint (or largest endpoint, if
    // we construct an upper envelope). 
    std::list<X_monotone_curve_2>    env_cvs;

    env_cvs.push_back(**iter);
    next = iter;
    ++next;
    while ((next != vert_vec.end()) &&
           (comp_x(min_vertex(**iter), min_vertex(**next)) == EQUAL))
    {
      if (env_type == LOWER) {
        // Compare the lower endpoints of both curves.
        res = comp_xy(min_vertex(env_cvs.front()), min_vertex(**next));

        // Update the list of vertical segments with minimal endpoints as
        // necessary.
        if (res == EQUAL) {
          env_cvs.push_back(**next);
        }
        if (res == LARGER) {
          env_cvs.clear();
          env_cvs.push_back(**next);
        }
      }
      else {
        // Compare the upper endpoints of both curves.
        res = comp_xy(max_vertex(env_cvs.front()), max_vertex(**next));

        // Update the list of vertical segments with maximal endpoints as
        // necessary.
        if (res == EQUAL) {
          env_cvs.push_back(**next);
        }
        if (res == SMALLER) {
          env_cvs.clear();
          env_cvs.push_back(**next);
        }
      }

      ++next;
    }

    // Compare the endpoint to the diagram feature.
    p = (env_type == LOWER) ?
      min_vertex(env_cvs.front()) : max_vertex(env_cvs.front());

    if (on_v) {
      // Compare p to the current vertex.
      res = comp_xy(p, v->point());

      if (res == EQUAL) {
        // Add curves to the current vertex.
        v->add_curves(env_cvs.begin(), env_cvs.end());
      }
      else if ((env_type == LOWER && res == SMALLER) ||
               (env_type == UPPER && res == LARGER))
      {
        // Replace the list of curves associated with the vertex.
        v->clear_curves();
        v->add_curves(env_cvs.begin(), env_cvs.end());
      }
    }
    else
    {
      // p lies in the interior of the current edge.
      Vertex_handle   new_v;

      if (e->is_empty()) {
        // Split the empty edge and associate the new vertex with the
        // vertical segments.
        new_v = _split_edge(out_d, p, e);
        new_v->add_curves(env_cvs.begin(), env_cvs.end());
      }
      else {
        // Compare p with the current curve.
        res = comp_y_at_x(p, e->curve());
        
        if (((env_type == LOWER) && (res != LARGER)) ||
            ((env_type == UPPER) && (res != SMALLER)))
        {
          new_v = _split_edge(out_d, p, e);
          new_v->add_curves(env_cvs.begin(), env_cvs.end());

          if (res == EQUAL)
            new_v->add_curve(e->curve());
        }
      }
    }

    // Proceed to the next vertical segment with larger x-coordinate.
    iter = next;
  }

  return;
}

// ---------------------------------------------------------------------------
// Split a given diagram edge by inserting a vertex in its interior.
//
template <class Traits, class Diagram>
typename Envelope_divide_and_conquer_2<Traits,Diagram>::Vertex_handle
Envelope_divide_and_conquer_2<Traits,Diagram>::
_split_edge(Envelope_diagram_1& diag, const Point_2& p, Edge_handle e)
{
  // Create the new vertex and the new edge.
  Vertex_handle   new_v = diag.new_vertex(p);
  Edge_handle     new_e = diag.new_edge();
  
  // Duplicate the curves container associated with e.
  if (! e->is_empty())
    new_e->add_curves(e->curves_begin(), e->curves_end());
  
  // Connect the new vertex between e and new_e.
  new_v->set_left(e);
  new_v->set_right(new_e);
  
  new_e->set_left(new_v);
  if (e != diag.rightmost())
    new_e->set_right(e->right());
  else
    diag.set_rightmost(new_e);

  e->set_right(new_v);

  // Return the new vertex.
  return (new_v);
}

} //namespace CGAL

#endif
