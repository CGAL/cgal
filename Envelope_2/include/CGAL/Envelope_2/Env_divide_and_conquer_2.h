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

#ifndef CGAL_ENVELOPE_DIVIDE_AND_CONQUER_2_H
#define CGAL_ENVELOPE_DIVIDE_AND_CONQUER_2_H

#include <CGAL/Env_default_diagram_1.h>

CGAL_BEGIN_NAMESPACE

/*! \class
 * A class implementing the divide-and-conquer algorithm for computing the
 * lower (or upper) envelope of a set of curves.
 */
template <class Traits_,
          class Diagram_ = Env_default_diagram_1<Traits_> >
class Envelope_divide_and_conquer_2
{
public:

  typedef Traits_                                  Traits_2;
  typedef typename Traits_2::Point_2               Point_2;
  typedef typename Traits_2::X_monotone_curve_2    X_monotone_curve_2;
  typedef typename Traits_2::Curve_2               Curve_2;

  typedef Diagram_                                 Envelope_diagram_1;
  
  typedef Envelope_divide_and_conquer_2<Traits_2, Envelope_diagram_1>  Self;

  enum Envelope_type
  {
    LOWER,
    UPPER
  };

protected:

  typedef typename Envelope_diagram_1::Vertex_const_handle Vertex_const_handle;
  typedef typename Envelope_diagram_1::Vertex_handle       Vertex_handle;
  typedef typename Envelope_diagram_1::Edge_const_handle   Edge_const_handle;
  typedef typename Envelope_diagram_1::Edge_handle         Edge_handle;

  typedef std::list<X_monotone_curve_2 *>          Curve_pointer_list;
  typedef typename Curve_pointer_list::iterator    Curve_pointer_iterator;

  // Data members:
  Traits_2        *traits;        // The traits object.
  bool             own_traits;    // Whether we own the traits object.
  Envelope_type    env_type;      // Either LOWER or UPPER.

  // Copy constructor and assignment operator - not supported.
  Envelope_divide_and_conquer_2 (const Self& );
  Self& operator= (const Self& );

public:

  /*!
   * Constructor.
   */
  Envelope_divide_and_conquer_2 () :
    own_traits(true),
    env_type(LOWER)
  {
    traits = new Traits_2;
  }

  /*!
   * Constructor with a traits object.
   * \param _traits The traits object.
   */
  Envelope_divide_and_conquer_2 (const Traits* _traits) :
    traits (_traits),
    own_traits(false),
    env_type(LOWER)
  {}

  /*!
   * Destructor.
   */
  ~Envelope_divide_and_conquer_2 ()
  {
    if (own_traits)
      delete traits;
  }

  /*!
   * Construct the lower (or upper) envelope to the given range of curves.
   * \param begin An iterator pointing at the beginning of the curves range. 
   * \param end A past-the-end iterator for the curves range.
   * \param type The envelope type (LOWER or UPPER).
   * \param diagram Output: The minimization (or maximization) diagram.
   */
  template <class CurvesIterator>
  void insert_curves (const CurvesIterator& begin,
                      const CurvesIterator& end,
                      const Envelope_type& type,
                      Envelope_diagram_1& diagram)
  {
    // Subdivide the curves into x-monotone subcurves.
    CurvesIterator                     it;
    std::list<CGAL::Object>            objects;
    std::list<CGAL::Object>::iterator  obj_it;
    X_monotone_curve_2                 xcv;
    std::list<X_monotone_curve_2>      x_curves;

    for (it = begin; it != end; it++)
    {
      // Split the current curve to x-monotone subcurves.
      objects.clear();
      traits->make_x_monotone_2_object()(*it, std::back_inserter(objects));

      for (obj_it = objects.begin(); obj_it != objects.end(); ++obj_it)
      {
        if(CGAL::assign (xcv, *obj_itr))
          x_curves.push_back (xcv);
      }
    }

    // Construct the envelope of the x-monotone curves.
    insert_x_monotone_curves (x_curves.begin(), x_curves.end(),
                              type,
                              diagram);
    return;
  }

  /*!
   * Construct the lower (or upper) envelope to the given range of
   * x-monotone curves.
   * \param begin An iterator pointing at the beginning of the curves range. 
   * \param end A past-the-end iterator for the curves range.
   * \param type The envelope type (LOWER or UPPER).
   * \param diagram Output: The minimization (or maximization) diagram.
   */
  template <class XCurvesIterator>
  void insert_x_monotone_curves (const XCurvesIterator& begin,
                                 const XCurvesIterator& end,
                                 const Envelope_type& type,
                                 Envelope_diagram_1& diagram)
  {
    // Set the envelope type.
    env_type = type;

    // Separate the regular curves from the vertical ones.
    typename Traits_2::Is_vertical_2  is_vertical = 
                                              traits->is_vertical_2_object();

    Curve_pointer_list    reg_list;
    Curve_pointer_list    vert_list;
    XCurvesIterator       iter;

    for (iter = begin; iter != end; ++iter)
    {
      if (is_vertical (*iter))
        vert_list.push_back (&(*iter));
      else
        reg_list.push_back (&(*iter));
    }

    // Construct the envelope for the non-vertical curves.
    _construct_envelope_non_vertical (reg_list.begin(), reg_list.end(),
                                      diagram);

    // Merge the vertical segments.
    if (vert_list.size() > 0)
      _merge_vertical_segments (vert_list,
                                diagram);

    return;
  }

  /*!
   * Get the traits object.
   * \return A pointer to the traits object.
   */
  Traits* get_traits () const
  {
    return (traits);
  }

protected:

  /*!
   * Construct the lower/upper envelope of the given list of non-vertical
   * curves.
   * \param begin The first x-monotone curve.
   * \param end A past-the-end iterator for the curves.
   * \param out_d Output: The minimization (or maximization) diagram.
   */
  void _construct_envelope_non_vertical (Curve_pointer_iterator begin,
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
      _construct_singleton_digram (*begin, out_d);
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

  /*!
   * Construct a singleton diagram, which matches a single curve.
   * \param cv The x-monotone curve.
   * \param out_d Output: The minimization (or maximization) diagram.
   */
  void _construct_singleton_diagram (const X_monotone_curve_2& cv,
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
      Edge_hanfle    e_right = out_d.new_edge();

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
      Edge_hanfle    e_left = out_d.new_edge();

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

  /*
   * Merge two minimization (or maximization) diagrams.
   * \param d1 The first diagram, 
   *           representing the envelope of the curve set C1.
   * \param d1 The first diagram,
   *           representing the envelope of the curve set C1.
   * \param out_d Output: The merged diagram, representing the envelope of
   *                      the union of C1 and C2.
   */
  void _merge_envelopes (const Envelope_diagram_1& d1,
                         const Envelope_diagram_1& d2,
                         Envelope_diagram_1& out_d)
  {
    Edge_const_handle    e1 = d1.leftmost();
    Vertex_const_handle  v1;
    Edge_const_handle    e2 = d2.leftmost();
    Vertex_const_handle  v2;
    bool                 prev_exists = false;
    Vertex_const_handle  prev_v;
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
        _merge_two_intervals (prev_v, prev_exists,
                              e1, e2,
                              next_v, next_exists,
                              (res_v == SMALLER) ? 1 : 2,
                              out_d);
      }
      else if (! e1->is_empty() && e2->is_empty())
      {
        // e1 is not empty but e2 is empty:
        _merge_single_interval (e1,
                                next_v, next_exists,
                                (res == SMALLER),
                                out_d);
      }
      else if (e1->is_empty() && ! e2->is_empty())
      {
        // e1 is empty and e2 is not empty:
        _merge_single_interval (e2,
                                next_v, next_exists,
                                (res != SMALLER),
                                out_d);
      }
      else
      {
        // Both edges are empty: append an empty edge to out_d:
        if (next_exists)
          _append_vertex (out_d, next_v, e1);
      }

        // Proceed to the next diagram edge(s), if possible.
      if (next_exists)
      {
        // Check if we should proceed on d1 or on d2.
        if (res_v == SMALLER)
        {
          prev_v = v1;
          prev_exists= true;
          
          e1 = v1->right();
          if (same_x)
            e2 = v2->right();
        }
        else if (res_v == LARGER)
        {
          prev_v = v2;
          prev_exists= true;
          
          e2 = v2->right();
          if (same_x)
            e1 = v1->right();
        }
        else
        {
          prev_v = v1;
          prev_exists= true;

          e1 = v1->right();
          e2 = v2->right();
        }
      }

      /* RWRW - is this needed?
      // Make sure that u is not to the left of the rightmost vertex in the
      // output diagram.
      if (traits->compare_x_2_object() (u->p, out_d.rightmostP->p) == SMALLER)
        u = out_d.rightmostP;
      */
    } while (next_exists);

    return;
  }

  /*!
   * Compare two vertices.
   * \param v1 The first vertex.
   * \param v2 The second vertex.
   * \param same_x Output parameter: TRUE iff x(v1->p) = x(v2->p).
   * \return SMALLER if x(v1->p) < x(v2->p). Or, in case x(v1->p) = x(v2->p):
   *                 if we compute the lower envelope and y(v1->p) < y(v2->p),
   *                 if we compute the upper envelope and y(v1->p) > y(v2->p).
   *         LARGER if x(v1->p) > x(v2->p). Or, in case x(v1->p) = x(v2->p):
   *                if we compute the lower envelope and y(v1->p) > y(v2->p),
   *                if we compute the upper envelope and y(v1->p) < y(v2->p).
   *         EQUAL if v1->p = v2->0.
   */
  Comparison_result _compare_vertices (const M_diagram_vertex_1 *v1,
                                       const M_diagram_vertex_1 *v2,
                                       bool& same_x) const
  {
    Comparison_result   res = traits->compare_x_2_object() (v1->p, v2->p);

    if (res != EQUAL)
    {
      same_x = false;
      return (res);
    }
    else
    {
      same_x = true;
    }

    // In case x(v1->p) = x(v2->p):
    res = traits->compare_xy_2_object() (v1->p, v2->p);

    if ((env_type == LOWER && res == SMALLER) ||
        (env_type == UPPER && res == LARGER))
      return (SMALLER);
    else if ((env_type == LOWER && res == LARGER) ||
             (env_type == UPPER && res == SMALLER))
      return (LARGER);

    return (EQUAL);
  }

  /*!
   * Get a pointer to an x-monotone curve that lies to the left of the given 
   * vertex, or return NULL if there is no such curve.
   * \param v The vertex.
   * \return A pointer to a curve that lies to the left of v.
   */
  const X_monotone_curve_2 *_curve_to_left (const M_diagram_vertex_1 *v) const
  {
    if (v->leftP == NULL)
      // No edge to the left of v:
      return (NULL);
    else if (v->leftP->mcvP == NULL)
      // The edge to the left is empty (contains to curves)
      return (NULL);

    // Return a pointer to the first curve in the list.
    return (&(v->leftP->mcvP->xcv));
  }

  /*!
   * Deal with an interval which is non-empty in one of the merged diagram and
   * empty in the other.
   * \param e The non-empty edge.
   * \param v The next vertex (which is right of out_d.rightmostP).
   * \param same_org Whether e and v originate from the same diagram.
   * \param out_d The merged diagram.
   */
  void _merge_single_interval (const M_diagram_edge_1* e,
                               const M_diagram_vertex_1* v,
                               const bool& same_org,
                               M_diagram_1& out_d) const
  {
    if (same_org)
    {
      // The non-empty edge ends at v, so we insert it to out_d.
      out_d.append_vertex (v->p, e);
      return;
    }

    // If v is not on e, we should insert it to the mered diagram only if it
    // is below (or above, in case of an upper envelope) the curves of e.
    const X_monotone_curve_2& cv = e->mcvP->xcv;
    Comparison_result         res = traits->compare_y_at_x_2_object() (v->p, cv);

    if ((env_type == LOWER && res == SMALLER) ||
        (env_type == UPPER && res == LARGER))
    {
      out_d.append_vertex (v->p, e);
    }

    return;
  }

  /*!
   * Merge two non-empty intervals into the merged diagram.
   * \param u The previous vertex of the merged diagrams we have visited.
   * \param e1 The first non-empty edge.
   * \param cv1 A curve that belongs to the first edge in this interval.
   * \param e2 The second non-empty edge.
   * \param cv2 A curve that belongs to the second edge in this interval.
   * \param v The next vertex (which is right of out_d.rightmostP).
   * \param org_v The origin of v: 1 if it is from e1, 2 if it is from e2.
   * \param out_d The merged diagram.
   */
  void _merge_two_intervals (const M_diagram_vertex_1* u,
                             const M_diagram_edge_1* e1,
                             const X_monotone_curve_2& cv1,
                             const M_diagram_edge_1* e2,
                             const X_monotone_curve_2& cv2,
                             const M_diagram_vertex_1* v,
                             const int& org_v,
                             M_diagram_1& out_d) const
  {
    // Compare the two curves to the right of (out_d.rightmostP).
    Comparison_result res = traits->curves_compare_y_at_x (cv1, cv2,
                                                           u->p);

    
   /* if (org_v == 1)
    {
      res = traits->compare_y_at_x_2_object() (v->p, cv2);
    }
    else
    {
      res = traits->compare_y_at_x_2_object() (v->p, cv1);
      if(res == SMALLER)
        res = LARGER;
      else
        if(res == LARGER)
          res = SMALLER;
    }*/

    
    if (res == EQUAL)
      res = traits->compare_y_at_x_left_2_object() (cv1, cv2,
                                                    v->p);
    if (env_type == UPPER)
    {
      // Flip the result in case of an upper envelope.
      if (res == SMALLER)
        res = LARGER;
      else if (res == LARGER)
        res = SMALLER;
    }

    // Find the next intersection of the envelopes to the right of the current
    // rightmost point in the merged diagram.
    CGAL::Object       obj;
    X_monotone_curve_2 icv;
    Point_2            p1, p2;

    while (true)
    {
      std::list<CGAL::Object>    objects;
      traits->intersect_2_object()(cv1, cv2, std::back_inserter(objects));
      CGAL::Object obj;
      
      // Stop if no intersection has been found.
      if (objects.empty())
        break;

      for(std::list<CGAL::Object>::iterator itr = objects.begin();
          itr != objects.end();
          ++itr)
      {
        std::pair<Point_2, unsigned int>  x_pt; 
        if(CGAL::assign(x_pt, *itr))
        {
          if(traits->compare_xy_2_object()(x_pt.first, out_d.rightmostP->p) == LARGER)
          {
            obj = make_object(x_pt.first);
            break;
          }
        }
        else
        {
          X_monotone_curve_2 cv;
          CGAL::assign(cv, *itr);
          CGAL_assertion(CGAL::assign(cv, *itr));
          Point_2 pt = traits->construct_min_vertex_2_object()(cv);
          if(traits->compare_xy_2_object()(pt, out_d.rightmostP->p) == LARGER)
          {
            obj = *itr;
            break;
          }
        }
      }

      if(obj.is_empty())
        break;
     /* obj = traits->nearest_intersection_to_right (cv1, cv2,
                                                   out_d.rightmostP->p);*/

      // Check for overlaps (when the returned object is a curve).
      if (CGAL::assign (icv, obj))
      {
        // Assign the leftmost endpoint in the overlapping curve to p1,
        // and the rightmost endpoint to p2.
        p1 = traits->construct_min_vertex_2_object() (icv);
        p2 = traits->construct_max_vertex_2_object() (icv);

        if (traits->compare_x_2_object() (p1, p2) != SMALLER)
        {
          Point_2 tmp = p1;
          p1 = p2;
          p2 = tmp;
        }

        // Disregard intersection points that are not to the left of v.
        if (traits->compare_x_2_object() (p1, v->p) != SMALLER)
          break;

        // Deal with overlaps:
        // If the overlap does not start at the rightmost vertex of the merged
        // diagram, insert p1 as a vertex.
        if (traits->compare_x_2_object() (p1, out_d.rightmostP->p) == LARGER)
        {
          // RWRW: bug fix!
          if (res == EQUAL)
            res = SMALLER;

          if (res == SMALLER)
            out_d.append_vertex (p1, e1);
          else if (res == LARGER)
            out_d.append_vertex (p1, e2);
          else
            // This case should never occur:
            CGAL_assertion (res != EQUAL);
        }

        // Make the leftmost point of p2 and v->p a vertex. The edge to its
        // left contains all curves in e1 and in e2.
        bool  reached_v = (traits->compare_x_2_object() (p2, v->p) != SMALLER);

        if (reached_v)
          out_d.append_vertex (v->p, e1);
        else
          out_d.append_vertex (p2, e1);

        M_diagram_edge_1 *left_e = out_d.rightmostP->leftP;
        left_e->mcvP = e2->mcvP;

        if (reached_v)
          return;

        // Compare the two curves to the right of the overlap.
        res = traits->compare_y_at_x_right_2_object() (cv1, cv2,
                                                   p2);

        if (env_type == UPPER)
        {
          // Flip the result in case of an upper envelope.
          if (res == SMALLER)
            res = LARGER;
          else if (res == LARGER)
            res = SMALLER;
        }
      }
      else if (CGAL::assign (p1, obj))
      {
        // Disregard intersection points that are not to the left of v.
        if (traits->compare_x_2_object() (p1, v->p) != SMALLER)
          break;

        // RWRW: bug fix!
        if (res == EQUAL)
            res = SMALLER;

        // Make p1 the new rightmost vertex of the merged diagram.
        if (res == SMALLER)
          out_d.append_vertex (p1, e1);
        else if (res == LARGER)
          out_d.append_vertex (p1, e2);
        else
          // This case should never occur (should be treated as an overlap):
          CGAL_assertion(res != EQUAL);

        // Compare the two curves to the right of the intersection point.
        res = traits->compare_y_at_x_right_2_object() (cv1, cv2,
                                                   p1);

        // This case should never occur (should be treated as an overlap):
        // RWRW: bug fix!
        if (res == EQUAL) res = SMALLER;
        // CGAL_assertion(res != EQUAL);

        if (env_type == UPPER)
        {
          // Flip the result in case of an upper envelope.
          if (res == SMALLER)
            res = LARGER;
          else if (res == LARGER)
            res = SMALLER;
        }
      }
      else
      {
        // This case should never occur:
        CGAL_assertion (false);
      }
    }

    // Check if v should also be inserted to the merged diagram.
    if (res == SMALLER)
    {
      // The final part of the interval is taken from e1.
      if (org_v == 1)
      {
        // In case v is also from e1, append it to the merged diagram.
        out_d.append_vertex (v->p, e1);
      }
      else
      {
        // If v is from e2, check if it below (or above, in case of an upper
        // envelope) cv1 to insert it.
        res = traits->compare_y_at_x_2_object() (v->p, cv1);

        if (res == EQUAL ||
            (env_type == LOWER && res == SMALLER) ||
            (env_type == UPPER && res == LARGER))
        {
          out_d.append_vertex (v->p, e1);
        }
      }
    }
    else if (res == LARGER)
    {
      // The final part of the interval is taken from e2.
      if (org_v == 2)
      {
        // In case v is also from e2, append it to the merged diagram.
        out_d.append_vertex (v->p, e2);
      }
      else
      {
        // If v is from e1, check if it below (or above, in case of an upper
        // envelope) cv2 to insert it.
        res = traits->compare_y_at_x_2_object() (v->p, cv2);

        if (res == EQUAL ||
            (env_type == LOWER && res == SMALLER) ||
            (env_type == UPPER && res == LARGER))
        {
          out_d.append_vertex (v->p, e2);
        }
      }
    }

    return;
  }

  /*!
   * Append a vertex to the given diagram: The new vertex that represents the 
   * given point as the new rightmost vertex of the diagram. The edge 
   * between the current rightmost vertex and the new one contains the same 
   * curves as the input edge.
   * \param diag The diagram.
   * \param p The point that the new vertex is associated with.
   * \param e The input edge.
   * \return A handle for the vertex.
   */
  Vertex_handle _append_vertex (Envelope_diagram_1& diag,
                                const Point_2& p, Edge_handle e)
  {
    // Create the new vertex and the new edge.
    Vertex_handle   new_v = diag.new_vertex (p);
    Edge_handle     new_e = diag.new_edge();
    
    if (! e->is_empty())
      new_e->add_curves (e->curves_begin(), e->curves_end());
    
    // Connect them to the right of the current rightmost vertex of the
    // diagram.
    new_v->set_left (new_e);
    new_v->set_right (diag.rightmost());
    new_e->set_right (new_v);
    new_e->set_left (diag.rightmost()->left());
    diag.rightmost()->left()->set_right (new_e);
    diag.rightmost()->set_left (new_v);
    
    return (new_v);
  }    


  /*! \struct
   * A functor used to sort vertical segments by their x-coordinate.
   */
  class Vertical_strict_weak_ordering
  {
  private:    
    const Traits         *traits;

  public:
    Vertical_strict_weak_ordering (const Traits *_traits) :
      traits(_traits)
    {}

    bool operator() (const M_curve_2 *mcv1, const M_curve_2 *mcv2) const
    {
      return (traits->compare_x_2_object() (traits->construct_min_vertex_2_object()(mcv1->xcv),
              traits->construct_min_vertex_2_object()(mcv1->xcv)) == SMALLER);
    }
  };

  /*!
   * Merge the vertical segments into the lower/upper envelope given as a
   * minimization (or maximization) diagram.
   * \param vert_list The list of vertical segments.
   * \param out_d The input minimization (or maximization) diagram.
   *             The function merges the vertical segments into this diagram.
   */
  void _merge_vertical_segments (Curve_pointer_list& vert_list,
                                 M_diagram_1& out_d)
  {
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

    return;
  }

};

CGAL_END_NAMESPACE

#endif

