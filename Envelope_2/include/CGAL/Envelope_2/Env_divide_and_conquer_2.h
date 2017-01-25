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

#ifndef CGAL_ENVELOPE_DIVIDE_AND_CONQUER_2_H
#define CGAL_ENVELOPE_DIVIDE_AND_CONQUER_2_H

#include <CGAL/license/Envelope_2.h>


#include <CGAL/Arr_enums.h>
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>

#include <boost/optional.hpp>

#include <vector>


namespace CGAL {

/*! \class
 * A class implementing the divide-and-conquer algorithm for computing the
 * lower (or upper) envelope of a set of curves.
 */
template <class Traits_, class Diagram_>
class Envelope_divide_and_conquer_2
{
public:

  typedef Traits_                                  Traits_2;
  typedef typename Traits_2::Point_2               Point_2;
  typedef typename Traits_2::X_monotone_curve_2    X_monotone_curve_2;
  typedef typename Traits_2::Curve_2               Curve_2;

  typedef Diagram_                                 Envelope_diagram_1;

protected:
  
  typedef Envelope_divide_and_conquer_2<Traits_2, Envelope_diagram_1>  Self;

  enum Envelope_type
  {
    LOWER,
    UPPER
  };

  typedef typename Envelope_diagram_1::Vertex_const_handle Vertex_const_handle;
  typedef typename Envelope_diagram_1::Vertex_handle       Vertex_handle;
  typedef typename Envelope_diagram_1::Edge_const_handle   Edge_const_handle;
  typedef typename Envelope_diagram_1::Edge_handle         Edge_handle;

  typedef std::vector<X_monotone_curve_2 *>        Curve_pointer_vector;
  typedef typename Curve_pointer_vector::iterator  Curve_pointer_iterator;

  typedef Arr_traits_adaptor_2<Traits_2>           Traits_adaptor_2;

  // Data members:
  Traits_adaptor_2  *traits;        // The traits object.
  bool               own_traits;    // Whether we own the traits object.
  Envelope_type      env_type;      // Either LOWER or UPPER.

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
    traits = new Traits_adaptor_2;
  }

  /*!
   * Constructor with a traits object.
   * \param _traits The traits object.
   */
  Envelope_divide_and_conquer_2 (const Traits_2* _traits) :
    own_traits(false),
    env_type(LOWER)
  {
    traits = static_cast<const Traits_adaptor_2*> (_traits);
  }

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
   * \param type The envelope type (true for lower, false of upper).
   * \param diagram Output: The minimization (or maximization) diagram.
   */
  template <class CurvesIterator>
  void insert_curves (CurvesIterator begin, CurvesIterator end,
                      bool type,
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
        if(CGAL::assign (xcv, *obj_it))
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
   * \param type The envelope type (true for lower, false for upper).
   * \param diagram Output: The minimization (or maximization) diagram.
   */
  template <class XCurvesIterator>
  void insert_x_monotone_curves (XCurvesIterator begin, XCurvesIterator end,
                                 bool type,
                                 Envelope_diagram_1& diagram)
  {
    // Set the envelope type.
    env_type = (type ? LOWER : UPPER);

    // Separate the regular curves from the vertical ones.
    typename Traits_2::Is_vertical_2  is_vertical = 
                                              traits->is_vertical_2_object();

    Curve_pointer_vector  reg_vec;
    Curve_pointer_vector  vert_vec;
    XCurvesIterator       iter;

    for (iter = begin; iter != end; ++iter)
    {
      if (is_vertical (*iter))
        vert_vec.push_back (&(*iter));
      else
        reg_vec.push_back (&(*iter));
    }

    // Construct the envelope for the non-vertical curves.
    _construct_envelope_non_vertical (reg_vec.begin(), reg_vec.end(),
                                      diagram);

    // Merge the vertical segments.
    if (vert_vec.size() > 0)
      _merge_vertical_segments (vert_vec,
                                diagram);

    return;
  }

  /*!
   * Get the traits object.
   * \return A pointer to the traits object.
   */
  Traits_2* get_traits () const
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
                                         Envelope_diagram_1& out_d);

  /*!
   * Construct a singleton diagram, which matches a single curve.
   * \param cv The x-monotone curve.
   * \param out_d Output: The minimization (or maximization) diagram.
   */
  void _construct_singleton_diagram (const X_monotone_curve_2& cv,
                                     Envelope_diagram_1& out_d);

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
                         Envelope_diagram_1& out_d);

  /*!
   * Compare two vertices.
   * \param v1 The first vertex.
   * \param v2 The second vertex.
   * \param same_x Output parameter: TRUE iff x(v1) = x(v2).
   * \return SMALLER if x(v1) < x(v2). Or, in case x(v1) = x(v2), and
   *                 - we compute the lower envelope, and y(v1) < y(v2),
   *                 - we compute the upper envelope, and y(v1) > y(v2).
   *         LARGER if x(v1) > x(v2). Or, in case x(v1) = x(v2), and
   *                - we compute the lower envelope, and y(v1) > y(v2),
   *                - we compute the upper envelope, and y(v1) < y(v2).
   *         EQUAL if v1 = v2.
   */
  Comparison_result _compare_vertices (Vertex_const_handle v1,
                                       Vertex_const_handle v2,
                                       bool& same_x) const;

  /*!
   * Deal with an interval which is non-empty in one of the merged diagrams and
   * empty in the other.
   * \param e The non-empty edge.
   * \param other_edge The empty edge.
   * \param v The next vertex to the right.
   * \param v_exists Whether the next vertex exists.
   * \param origin_of_v The origin of v: SMALLER if it is from e, 
   *                    LARGER if it is from other_edge. 
   *                    EQUAL result means that both edges have vertex at 
   *                    the same place.
   * \param out_d The merged diagram.
   */
  void _merge_single_interval (Edge_const_handle e, 
                               Edge_const_handle other_edge,
                               Vertex_const_handle v, bool v_exists,
                               Comparison_result origin_of_v,
                               Envelope_diagram_1& out_d);
  
  
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
  Comparison_result compare_y_at_end(const X_monotone_curve_2& xcv1,
                                     const X_monotone_curve_2& xcv2,
                                     Arr_curve_end curve_end) const;



  /*!
   * Merge two non-empty intervals into the merged diagram.
   * \param e1 The first non-empty edge.
   * \param is_leftmost1 Is it the leftmost edge in its diagram.
   * \param e2 The second non-empty edge.
   * \param is_leftmost2 Is it the leftmost edge in its diagram.
   * \param v The next vertex.
   * \param v_exists Whether such a vertex exists.
   * \param origin_of_v The origin of v: SMALLER if it is from e1, 
   *                    otherwise it is from e2. EQUAL result means that
   *                    both diagram have vertex at the same place (but v
   *                    is still taken from e2.
   * \param out_d The merged diagram.
   */
  void _merge_two_intervals (Edge_const_handle e1, bool is_leftmost1,
                             Edge_const_handle e2, bool is_leftmost2,
                             Vertex_const_handle v, bool v_exists,
                             Comparison_result origin_of_v,
                             Envelope_diagram_1& out_d);

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
                                const Point_2& p, Edge_const_handle e);

  /*! \struct
   * A functor used to sort vertical segments by their x-coordinate.
   */
  class Less_vertical_segment
  {
  private:    

    typename Traits_2::Compare_x_2             comp_x;
    typename Traits_2::Construct_min_vertex_2  min_vertex;

  public:

    Less_vertical_segment (const Traits_2 *traits) :
        comp_x(traits->compare_x_2_object()),
        min_vertex(traits->construct_min_vertex_2_object())
    {}

    bool operator() (const X_monotone_curve_2 *cv1,
                     const X_monotone_curve_2 *cv2) const
    {
      return (comp_x (min_vertex (*cv1),
                      min_vertex (*cv2)) == SMALLER);
    }
  };

  /*!
   * Merge the vertical segments into the lower/upper envelope given as a
   * minimization (or maximization) diagram.
   * \param vert_vec The list of vertical segments.
   * \param out_d The input minimization (or maximization) diagram.
   *             The function merges the vertical segments into this diagram.
   */
  void _merge_vertical_segments (Curve_pointer_vector& vert_vec,
                                 Envelope_diagram_1& out_d);

  /*!
   * Split a given diagram edge by inserting a vertex in its interior.
   * \param diag The diagram.
   * \param p The point that the new vertex is associated with.
   * \param e The edge to split.
   * \return A handle for the vertex.
   */
  Vertex_handle _split_edge (Envelope_diagram_1& diag,
                             const Point_2& p, Edge_handle e);

};

} //namespace CGAL

#include <CGAL/Envelope_2/Env_divide_and_conquer_2_impl.h>

#endif
