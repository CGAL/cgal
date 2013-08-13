// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s)     : Tali Zvi <talizvi@post.tau.ac.il>,
//                 Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ron Wein <wein@post.tau.ac.il>

#ifndef CGAL_SWEEP_LINE_SUBCURVE_H
#define CGAL_SWEEP_LINE_SUBCURVE_H

/*! \file
 * Defintion of the Sweep_line_subcurve class.
 */

#include <CGAL/Sweep_line_2/Sweep_line_functors.h>
#include <CGAL/Sweep_line_2/Sweep_line_event.h>
#include <CGAL/Multiset.h>
#include <CGAL/assertions.h>

namespace CGAL {

/*! \class Sweep_line_subcurve
 *
 * This is a wrapper class to X_monotone_curve_2 in the traits class, that
 * contains data that is used when applying the sweep algorithm on a set of
 * x-monotone curves.
 *
 * The information contained in this class is:
 * - the remaining x-monotone curve that is to the right of the current sweep
 *   line.
 * - two event points which are associated with the left and right end of the
 *   curve.
 * - an iterator that points to the location of the subcurve at the status line.
 * - two pointers to subcurves that are the originating subcurves in case of
 *   an overlap, otherwise thay are both NULL.
 *
 */
template <typename Traits_>
class Sweep_line_subcurve {
public:

  typedef Traits_                                    Traits_2;
  typedef typename Traits_2::Point_2                 Point_2;
  typedef typename Traits_2::X_monotone_curve_2      X_monotone_curve_2;

  typedef Sweep_line_subcurve<Traits_2>              Self;
  typedef Curve_comparer<Traits_2, Self>             Compare_curves;
  typedef Multiset<Self*,
                   Compare_curves,
                   CGAL_ALLOCATOR(int)>              Status_line;
  typedef typename Status_line::iterator             Status_line_iterator;

  typedef Sweep_line_event<Traits_2, Self>           Event;

protected:
  // Data members:
  X_monotone_curve_2 m_lastCurve;   // The portion of the curve that lies to
                                    // the right of the last event point 
                                    // that occured on the curve.

  Event* m_left_event;              // The event associated with the left end.
  Event* m_right_event;             // The event associated with the right end
  
  Status_line_iterator m_hint;      // The location of the subcurve in the
                                    // status line (the Y-structure).

  Self* m_orig_subcurve1;           // The overlapping hierarchy
  Self* m_orig_subcurve2;           // (relevant only in case of overlaps).

public:

  /*! Default constructor. */
  Sweep_line_subcurve() :
    m_orig_subcurve1(NULL),
    m_orig_subcurve2(NULL)
  {}

  /*! Constructor given a curve. */
  Sweep_line_subcurve(const X_monotone_curve_2& curve) :
    m_lastCurve(curve),
    m_orig_subcurve1(NULL),
    m_orig_subcurve2(NULL)
  {}

  /*! Initialize the subcurves by setting the curve. */
  void init(const X_monotone_curve_2& curve) { m_lastCurve = curve; }

  /*! Destructor. */
  ~Sweep_line_subcurve() {}

  /*! Get the last intersecing curve so far (const version). */
  const X_monotone_curve_2& last_curve() const { return m_lastCurve; }

  /*! Get the last intersecing curve so far (non-const version). */
  X_monotone_curve_2& last_curve() { return m_lastCurve; }

  /*! Set the last intersecing curve so far. */
  void set_last_curve(const X_monotone_curve_2& cv) { m_lastCurve = cv; }

  /*! Check if the given event is the matches the right-end event. */  
  template <typename SweepEvent>
  bool is_end_point(const SweepEvent* event) const
  { return (m_right_event == (Event*)event); }
 
  /*! Get the event that corresponds to the left end of the subcurve. */
  Event* left_event() const { return m_left_event; }

  /*! Get the event that corresponds to the right end of the subcurve. */
  Event* right_event() const { return m_right_event; }

  /*! Set the event that corresponds to the left end of the subcurve. */
  template<class SweepEvent>
  void set_left_event(SweepEvent* event) { m_left_event =(Event*)event; }

  /*! Set the event that corresponds to the right end of the subcurve. */
  template<class SweepEvent>
  void set_right_event(SweepEvent* event) { m_right_event = (Event*)event; }

  /*! Get the location of the subcurve in the status line .*/
  Status_line_iterator hint() const { return m_hint; }

  /*! Set the location of the subcurve in the status line .*/
  void set_hint(Status_line_iterator hint) { m_hint = hint; }

  /*! Get the subcurves that originate an overlap. */
  Self* originating_subcurve1() { return m_orig_subcurve1; }

  Self* originating_subcurve2() { return m_orig_subcurve2; }

  /*! Set the subcurves that originate an overlap. */
  void set_originating_subcurve1(Self* orig_subcurve1)
  { m_orig_subcurve1 = orig_subcurve1; }

  void set_originating_subcurve2(Self* orig_subcurve2)
  { m_orig_subcurve2 = orig_subcurve2; }

  /*! Get all the leaf-nodes in the hierarchy of overlapping subcurves. */
  template <typename OutputIterator>
  OutputIterator all_leaves(OutputIterator oi)
  {
    if (m_orig_subcurve1 == NULL) {
      *oi = this;
      ++oi;
      return oi;
    }

    oi = m_orig_subcurve1->all_leaves(oi);
    oi = m_orig_subcurve2->all_leaves(oi);
    return oi;
  }

  /*! Check if the given subcurve is a node in the overlapping hierarchy. */
  bool is_inner_node(Self *s)
  {
    if (this == s) return true;
    if (m_orig_subcurve1 == NULL) return false;
    return (m_orig_subcurve1->is_inner_node(s) ||
            m_orig_subcurve2->is_inner_node(s));
  }

  /*! Check if the given subcurve is a leaf in the overlapping hierarchy. */
  bool is_leaf(Self* s)
  {
    if (m_orig_subcurve1 == NULL) return (this == s);
    return (m_orig_subcurve1->is_leaf(s) ||
            m_orig_subcurve2->is_leaf(s));
  }

  /*! Check if the two hierarchies contain the same leaf nodes. */
  bool has_same_leaves(Self *s)
  {
    std::list<Self*> my_leaves;
    std::list<Self*> other_leaves;
    
    this->all_leaves (std::back_inserter(my_leaves));
    s->all_leaves (std::back_inserter(other_leaves));

    typename std::list<Self*>::iterator  iter;
    for (iter = my_leaves.begin(); iter != my_leaves.end(); ++iter) {
      if (std::find(other_leaves.begin(), other_leaves.end(), *iter) ==
          other_leaves.end())
        return false;
    }

    for (iter = other_leaves.begin(); iter != other_leaves.end(); ++iter) {
      if (std::find(my_leaves.begin(), my_leaves.end(), *iter) ==
          my_leaves.end())
        return false;
    }

    return true;
  }

  /*! Check if the two hierarchies contain a common leaf node. */
  bool has_common_leaf(Self *s)
  {
    std::list<Self*> my_leaves;
    std::list<Self*> other_leaves;
    
    this->all_leaves(std::back_inserter(my_leaves));
    s->all_leaves(std::back_inserter(other_leaves));

    typename std::list<Self*>::iterator iter;
    for (iter = my_leaves.begin(); iter != my_leaves.end(); ++iter) {
      if (std::find(other_leaves.begin(), other_leaves.end(), *iter) !=
          other_leaves.end())
        return true;
    }
    return false;
  }

  /*! Get all distinct nodes from the two hierarchies. */
  template <class OutputIterator>
  OutputIterator distinct_nodes(Self* s, OutputIterator oi)
  {
    if (m_orig_subcurve1 == NULL) {
      if (s->is_leaf(this)) {
        *oi = this;
        ++oi;
      }
      return oi;
    }

    if (! s->is_inner_node (m_orig_subcurve1)) {
      *oi = m_orig_subcurve1;
      ++oi;
    }
    else {
      oi = m_orig_subcurve1->distinct_nodes(s, oi);
    }

    if (! s->is_inner_node (m_orig_subcurve2)) {
      *oi = m_orig_subcurve2;
      ++oi;
    }
    else {
      oi = m_orig_subcurve2->distinct_nodes(s, oi);
    }

    return oi;
  }

  /*! Get the depth of the overlap hierarchy. */
  unsigned int overlap_depth()
  {
    if (m_orig_subcurve1 == NULL) return (1);

    unsigned int depth1 = m_orig_subcurve1->overlap_depth();
    unsigned int depth2 = m_orig_subcurve2->overlap_depth();
    if (depth1 > depth2) return (depth1 + 1);
    else return (depth2 + 1);
  }
 
#ifdef CGAL_SL_VERBOSE
  void Print() const;
#endif
};

#ifdef CGAL_SL_VERBOSE
  template<class Traits>
  void Sweep_line_subcurve<Traits>::Print() const
  {
    std::cout << "Curve " << this 
              << "  (" << m_lastCurve << ") " 
              << " [sc1: " << m_orig_subcurve1
              << ", sc2: " << m_orig_subcurve2 << "]"
              << std::endl;
  }
#endif

} //namespace CGAL

#endif
