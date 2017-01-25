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
// Author(s) : Tali Zvi        <talizvi@post.tau.ac.il>,
//             Baruch Zukerman <baruchzu@post.tau.ac.il>
//             Ron Wein        <wein@post.tau.ac.il>
//             Efi Fogel       <efif@gmail.com>

#ifndef CGAL_SWEEP_LINE_EVENT_H
#define CGAL_SWEEP_LINE_EVENT_H

#include <CGAL/license/Sweep_line_2.h>


/*! \file
 * Defintion of the Sweep_line_event class.
 */

#include <list>

#include <CGAL/Sweep_line_2/Sweep_line_subcurve.h>

namespace CGAL {

/*! \class Sweep_line_event
 *
 * A class associated with an event in a sweep line algorithm.
 * An intersection point in the sweep line algorithm is refered to as an event.
 * This class contains the information that is associated with any given
 * event point. This information contains the following:
 * - the actual point
 * - a list of curves that pass through the event point and defined to
 *   the left of the event point.
 * - a list of curves that pass through the event point and defined to
 *   the right of the event point.
 *
 * The class mostly exists to store information and does not have any
 * significant functionality otherwise.
 *
 */
template <typename Traits_, typename Subcurve_>
class Sweep_line_event {
public:
  typedef Traits_                                       Traits_2;
  typedef typename Traits_2::X_monotone_curve_2         X_monotone_curve_2;
  typedef typename Traits_2::Point_2                    Point_2;

  // should be ok, as Traits_ has already extended by Basic_sweep_line
  typedef typename internal::Arr_complete_left_side_category<Traits_2>::Category
                                                        Left_side_category;
  typedef typename internal::Arr_complete_bottom_side_category<Traits_2>::Category
                                                        Bottom_side_category;
  typedef typename internal::Arr_complete_top_side_category<Traits_2>::Category
                                                        Top_side_category;
  typedef typename internal::Arr_complete_right_side_category<Traits_2>::Category
                                                        Right_side_category;

  typedef Subcurve_                                     Subcurve;
  //template<typename SC>
  //struct SC_container { typedef std::list<SC> other; };
  typedef std::list<Subcurve*>                          Subcurve_container;
  typedef typename Subcurve_container::iterator         Subcurve_iterator;
  typedef typename Subcurve_container::const_iterator   Subcurve_const_iterator;
  typedef typename Subcurve_container::reverse_iterator
    Subcurve_reverse_iterator;

  /*! \enum The event type (with other information bits). */
  enum Attribute {
    DEFAULT = 0,
    LEFT_END = 1,            // A curve's left-end is on the event point.
    RIGHT_END = 2,           // A curve's right-end is on the event point.
    ACTION = 4,              // An action point.
    QUERY = 8,               // A query point.
    INTERSECTION = 16,       // Two curves intersects at their interior.
    WEAK_INTERSECTION = 32,  // A curve's end-point is on the interior
                             // of another curve (also may indicate overlap).
    OVERLAP = 64             // Endpoint of an overlap subcurve
  };

protected:
  // Data members:
  Point_2 m_point;                  // The point associated with the event.

  Subcurve_container m_leftCurves;  // The curves lying to the left of the
                                    // event and incident to it.

  Subcurve_container m_rightCurves; // The curves lying to the right of the
                                    // event and incident to it, sorted by
                                    // their y-value to the right of the point.

  char m_type;                      // The event type.
  char m_ps_x;                      // The boundary condition in x.
  char m_ps_y;                      // The boundary condition in y.
  char m_closed;                    // Is the event closed (associated with
                                    // a valid point.

public:
  /*! Default constructor. */
  Sweep_line_event() :
    m_type(0),
    m_ps_x(static_cast<char>(ARR_INTERIOR)),
    m_ps_y(static_cast<char>(ARR_INTERIOR)),
    m_closed(1)
  {}

  /*! Initialize an event that is associated with a valid point. */
  void init(const Point_2& point, Attribute type,
            Arr_parameter_space ps_x, Arr_parameter_space ps_y)
  {
    m_point = point;
    m_type = type;
    m_ps_x = static_cast<char>(ps_x);
    m_ps_y = static_cast<char>(ps_y);
    m_closed = 1;
  }

  /*! Initialize an event associates with an open curve end. */
  void init_at_open_boundary(Attribute type,
                             Arr_parameter_space ps_x,
                             Arr_parameter_space ps_y)
  {
    m_type = type;
    m_ps_x = ps_x;
    m_ps_y = ps_y;
    m_closed = 0;
  }

  /*! Add a subcurve to the container of left curves. */
  void add_curve_to_left(Subcurve* curve)
  {
    // Look for the subcurve.
    Subcurve_iterator iter;

    //std::cout << "add_curve_to_left, curve: ";
    //curve->Print();

    for (iter = m_leftCurves.begin(); iter != m_leftCurves.end(); ++iter) {
      //std::cout << "add_curve_to_left, iter: ";
      //(*iter)->Print();

      // Do nothing if the curve exists.
      if ((curve == *iter) || (*iter)->is_inner_node(curve)) {
        //std::cout << "add_curve_to_left, curve exists" << std::endl;
        return;
      }

      // Replace the existing curve in case of overlap.
      // EBEB 2011-10-27: Fixed to detect overlaps correctly
      if ((curve != *iter) && (curve->has_common_leaf(*iter))) {
        //std::cout << "add_curve_to_left, curve overlaps" << std::endl;
        *iter = curve;
        return;
      }
    }

    // The curve does not exist - insert it to the container.
    m_leftCurves.push_back(curve);
    // std::cout << "add_curve_to_left, pushed back" << std::endl;

    //this->Print();
  }

  /*! Add a subcurve to the container of left curves (without checks). */
  void push_back_curve_to_left(Subcurve* curve)
  { m_leftCurves.push_back(curve); }

  /*! Add a subcurve to the container of right curves. */
  std::pair<bool, Subcurve_iterator>
  add_curve_to_right(Subcurve* curve, const Traits_2* tr)
  {
    if (m_rightCurves.empty()) {
      m_rightCurves.push_back(curve);
      return (std::make_pair(false, m_rightCurves.begin()));
    }

    // Check if its an event at open boundary,
    // and if so then there is no overlap
    //(there cannot be two non-overlap curves at the same event at open
    // boundary).
    if (!this->is_closed())
      return (std::make_pair(true, m_rightCurves.begin()));

    Subcurve_iterator iter = m_rightCurves.begin();
    Comparison_result res;

    while ((res = tr->compare_y_at_x_right_2_object()
            (curve->last_curve(), (*iter)->last_curve(), m_point)) == LARGER)
    {
      ++iter;
      if (iter == m_rightCurves.end()) {
        m_rightCurves.insert(iter, curve);
        return std::make_pair(false, --iter);
      }
    }

    if (res == EQUAL) //overlap !!
    {
      return std::make_pair(true, iter);
    }

    m_rightCurves.insert(iter, curve);
    return std::make_pair(false, --iter);
  }

  /*! Add two Subcurves to the right of the event.
   * \pre The event does not contain any right curves, and the order of sc1
   *      and sc2 is correct.
   */
  std::pair<bool, Subcurve_iterator>
  add_curve_pair_to_right(Subcurve* sc1, Subcurve* sc2)
  {
    m_rightCurves.push_back(sc1);
    m_rightCurves.push_back(sc2);

    Subcurve_iterator iter = m_rightCurves.end();
    --iter;
    return (std::make_pair(false, iter));
  }

  /*! Remove a curve from the set of left curves. */
  void remove_curve_from_left(Subcurve* curve)
  {
    Subcurve_iterator iter;
    for (iter = m_leftCurves.begin(); iter!= m_leftCurves.end(); ++iter) {
      if (curve->has_common_leaf(*iter)) {
        m_leftCurves.erase(iter);
        return;
      }
    }
  }

  /*! Returns an iterator to the first curve to the left of the event. */
  Subcurve_iterator left_curves_begin() { return (m_leftCurves.begin()); }

  /*! Returns an iterator to the one past the last curve to the left
      of the event. */
  Subcurve_iterator left_curves_end() { return (m_leftCurves.end()); }

  /*! Returns an iterator to the first curve to the right of the event. */
  Subcurve_iterator right_curves_begin() { return (m_rightCurves.begin()); }

  /*! Returns an iterator to the one past the last curve to the right
      of the event. */
  Subcurve_iterator right_curves_end() { return (m_rightCurves.end()); }

  /*! Returns a const iterator to the first curve to the left of the event. */
  Subcurve_const_iterator left_curves_begin() const
  { return m_leftCurves.begin(); }

  /*! Returns a const iterator to the past the end curve to the left
    of the event. */
  Subcurve_const_iterator left_curves_end() const
  { return m_leftCurves.end(); }

  /*! Returns a const iterator to the first curve to the right of the event. */
  Subcurve_const_iterator right_curves_begin() const
  { return m_rightCurves.begin(); }

  /*! Returns a const iterator to the past the end curve to the right
    of the event. */
  Subcurve_const_iterator right_curves_end() const
  { return m_rightCurves.end(); }

  /*! Returns a reverse_iterator to the first curve of the reversed list
      of the right curves of the event. */
  Subcurve_reverse_iterator right_curves_rbegin()
  { return m_rightCurves.rbegin(); }

  /*! Returns a reverse_iterator to the past-end curve of the reversed list
      of the right curves of the event. */
  Subcurve_reverse_iterator right_curves_rend()
  { return m_rightCurves.rend(); }

  /*! Returns a reverse_iterator to the first curve of the reversed list
      of the left curves of the event. */
  Subcurve_reverse_iterator left_curves_rbegin()
  { return m_leftCurves.rbegin(); }

  /*! Returns a reverse_iterator to the past-end curve of the reversed list
      of the left curves of the event. */
  Subcurve_reverse_iterator left_curves_rend()
  { return m_leftCurves.rend(); }

  /*! Returns the number of curves defined to the left of the event. */
  size_t number_of_left_curves() { return m_leftCurves.size(); }

  /*! Returns the number of curves defined to the right of the event. */
  size_t number_of_right_curves() { return (m_rightCurves.size()); }

  /*! Checks if at least one curve is defined to the left of the event. */
  bool has_left_curves() const { return (! m_leftCurves.empty()); }

  /*! Checks if at least one curve is defined to the right of the event. */
  bool has_right_curves() const { return (! m_rightCurves.empty()); }

  /*!
   * Get the actual event point (const version).
   * \pre The event is associated with a valid point.
   */
  const Point_2& point() const
  {
    CGAL_precondition(is_closed());
    return m_point;
  }

  /*!
   * Get the actual event point (non-const version).
   * \pre The event is associated with a valid point.
   */
  Point_2& point()
  {
    CGAL_precondition(is_closed());
    return m_point;
  }

  /*!
   * Get a curve associated with the event (const version).
   * \pre The event has incident curves.
   */
  const X_monotone_curve_2& curve() const
  {
    if (has_left_curves())
      return (m_leftCurves.front()->last_curve());

    CGAL_assertion (has_right_curves());
    return (m_rightCurves.front()->last_curve());
  }

  /*! Set the event point. */
  void set_point(const Point_2& pt) { m_point = pt; }

  /// \name Get the event attributes.
  //@{
  bool is_left_end() const { return ((m_type & LEFT_END) != 0); }

  bool is_right_end() const { return ((m_type & RIGHT_END) != 0); }

  bool is_intersection() const { return ((m_type & INTERSECTION ) != 0); }

  bool is_action() const { return ((m_type & ACTION ) != 0); }

  bool is_query() const { return ((m_type & QUERY) != 0); }

  bool is_weak_intersection() const
  { return((m_type & WEAK_INTERSECTION) != 0); }

  bool is_overlap() const { return ((m_type & OVERLAP) != 0); }
  //@}

  /// \name Set the event attributes.
  //@{
  void set_left_end() { m_type |= LEFT_END; }

  void set_right_end() { m_type |= RIGHT_END; }

  void set_intersection() { m_type |= INTERSECTION; }

  void set_action() { m_type |= ACTION; }

  void set_query() { m_type |= QUERY; }

  void set_weak_intersection() { m_type |= WEAK_INTERSECTION; }

  void set_overlap() { m_type |= OVERLAP; }

  void set_attribute(Attribute type) { m_type |= type; }
  //@}

  /// \name Get the boundary conditions of the event.
  //@{
  inline bool is_closed() const { return (m_closed != 0); }

  inline bool is_on_boundary () const
  {
    return ((m_ps_x != static_cast<char>(ARR_INTERIOR)) ||
            (m_ps_y != static_cast<char>(ARR_INTERIOR)));
  }

  inline Arr_parameter_space parameter_space_in_x() const
  { return (Arr_parameter_space(m_ps_x)); }

  inline Arr_parameter_space parameter_space_in_y() const
  { return (Arr_parameter_space(m_ps_y)); }
  //@}

  /*! Replace the set of left subcurves. */
  template <typename InputIterator>
  void replace_left_curves(InputIterator begin, InputIterator end)
  {
    Subcurve_iterator left_iter = m_leftCurves.begin();
    for (InputIterator iter = begin; iter != end; ++iter, ++left_iter)
      *left_iter = static_cast<Subcurve*>(*iter);
    m_leftCurves.erase(left_iter, m_leftCurves.end());
  }

  bool is_right_curve_bigger(Subcurve* c1, Subcurve* c2)
  {
    Subcurve_iterator   iter;
    for (iter = m_rightCurves.begin(); iter != m_rightCurves.end(); ++iter) {
      if ((*iter == c1) ||
          (static_cast<Subcurve*>((*iter)->originating_subcurve1()) == c1) ||
          (static_cast<Subcurve*>((*iter)->originating_subcurve2()) == c1))
        return false;

      if ((*iter == c2) ||
          (static_cast<Subcurve*>((*iter)->originating_subcurve1()) == c2) ||
          (static_cast<Subcurve*>((*iter)->originating_subcurve2()) == c2))
        return true;
    }

    return true;
  }

  /*! Check if the two curves are negihbors to the left of the event. */
  bool are_left_neighbours(Subcurve* c1, Subcurve* c2)
  {
    Subcurve_iterator left_iter = m_leftCurves.begin();
    for (; left_iter != m_leftCurves.end(); ++left_iter) {
      if (*left_iter == c1) {
        Subcurve_iterator temp = left_iter;
        ++temp;
        if (temp != m_leftCurves.end()) return (*temp == c2);
        return false;
      }

      if (*left_iter == c2) {
        Subcurve_iterator temp = left_iter;
        ++temp;
        if (temp!=m_leftCurves.end()) return (*temp == c1);
        return false;
      }
    }

    return false;
  }

#ifdef CGAL_SL_VERBOSE
  void Print();
#endif
};

#ifdef CGAL_SL_VERBOSE
  template<typename Traits, typename Subcurve>
  void Sweep_line_event<Traits, Subcurve>::Print()
  {
    std::cout << "\tEvent info: "  << "\n" ;
    if (this->is_closed())
      std::cout << "\t" << m_point << "\n" ;
    else {
      std::cout << "\t";
      Arr_parameter_space ps_x = this->parameter_space_in_x();
      Arr_parameter_space ps_y = this->parameter_space_in_y();

      switch (ps_x) {
       case ARR_LEFT_BOUNDARY:  std::cout << "left boundary"; break;
       case ARR_RIGHT_BOUNDARY: std::cout << "right boundary"; break;
       case ARR_INTERIOR:
       default:
        switch (ps_y) {
         case ARR_BOTTOM_BOUNDARY: std::cout << "bottom boundary"; break;
         case ARR_TOP_BOUNDARY:    std::cout << "top boundary"; break;
         case ARR_INTERIOR:
         default:
          CGAL_error();
        }
      }
    }
    std::cout << "\n";

    std::cout << "\tLeft curves: \n" ;
    Subcurve_iterator iter;
    for (iter = m_leftCurves.begin(); iter != m_leftCurves.end(); ++iter) {
      std::cout << "\t";
      (*iter)->Print();
      std::cout << "\n";
    }
    std::cout << std::endl;
    std::cout << "\tRight curves: \n" ;
    for (iter = m_rightCurves.begin(); iter != m_rightCurves.end(); ++iter) {
      std::cout << "\t";
      (*iter)->Print();
      std::cout << "\n";
    }

    std::cout << std::endl;
  }
#endif // CGAL_SL_VERBOSE

} //namespace CGAL

#endif
