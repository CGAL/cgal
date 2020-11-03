// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Tali Zvi        <talizvi@post.tau.ac.il>,
//             Baruch Zukerman <baruchzu@post.tau.ac.il>
//             Ron Wein <wein@post.tau.ac.il>
//             Efi Fogel <efif@gmail.com>

#ifndef CGAL_SURFACE_SWEEP_2_NO_OVERLAP_EVENT_BASE_H
#define CGAL_SURFACE_SWEEP_2_NO_OVERLAP_EVENT_BASE_H

#include <CGAL/license/Surface_sweep_2.h>

/*! \file
 *
 * Defintion of the No_overlap_event_base class.
 */

#include <list>

namespace CGAL {
namespace Surface_sweep_2 {

// This class is used to test if `void* T::for_compact_container()`
// exists, to avoid adding a void* pointer to the Event_base structure
// if Point_2 can already be used as a handle for this
template<typename T>
struct has_for_compact_container
{
private:

  template<typename U> static
  auto test(void*) -> decltype(std::declval<U>().for_compact_container() == nullptr, Tag_true());

  template<typename> static Tag_false test(...);

public:

  static constexpr bool value = std::is_same<decltype(test<T>(nullptr)),Tag_true>::value;
};

template <typename Point_2, bool HasFor>
class Event_base_for_compact_container { };

template <typename Point_2>
class Event_base_for_compact_container<Point_2, false>
{
  void* p = nullptr;
public:

  void* operator()(const Point_2&) const
  {
    return p;
  }
  void operator() (Point_2&, void* ptr)
  {
    p = ptr;
  }
};

template <typename Point_2>
class Event_base_for_compact_container<Point_2, true>
{

public:

  void* operator()(const Point_2& p) const
  {
    return p.for_compact_container();
  }
  void operator() (Point_2& p, void* ptr)
  {
    p.for_compact_container(ptr);
  }
};


/*! \class No_overlap_event_base
 *
 * A class associated with an event in a surface-sweep algorithm.
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
template <typename GeometryTraits_2, typename Subcurve_>
class No_overlap_event_base {
public:
  typedef GeometryTraits_2                              Geometry_traits_2;
  typedef Subcurve_                                     Subcurve;

private:
  typedef Geometry_traits_2                             Gt2;

public:
  typedef typename Gt2::X_monotone_curve_2              X_monotone_curve_2;
  typedef typename Gt2::Point_2                         Point_2;

  typedef typename internal::Arr_complete_left_side_category<Gt2>::Category
                                                        Left_side_category;
  typedef typename internal::Arr_complete_bottom_side_category<Gt2>::Category
                                                        Bottom_side_category;
  typedef typename internal::Arr_complete_top_side_category<Gt2>::Category
                                                        Top_side_category;
  typedef typename internal::Arr_complete_right_side_category<Gt2>::Category
                                                        Right_side_category;

  typedef std::list<Subcurve*>                          Subcurve_container;
  typedef typename Subcurve_container::iterator         Subcurve_iterator;
  typedef typename Subcurve_container::const_iterator   Subcurve_const_iterator;
  typedef typename Subcurve_container::reverse_iterator
    Subcurve_reverse_iterator;

  typedef Event_base_for_compact_container
  <Point_2, has_for_compact_container<Point_2>::value>    For_compact_container;

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

  Subcurve_container m_left_curves;  // The curves lying to the left of the
                                    // event and incident to it.

  Subcurve_container m_right_curves; // The curves lying to the right of the
                                    // event and incident to it, sorted by
                                    // their y-value to the right of the point.

  char m_type;                      // The event type.
  char m_ps_x;                      // The boundary condition in x.
  char m_ps_y;                      // The boundary condition in y.
  char m_closed;                    // Is the event closed (associated with
                                    // a valid point.

  // A handle for the compact container (either using the functions of
  // `m_point` if available, or an additional pointer)
  For_compact_container m_for_compact_container;

public:
  /*! Default constructor. */
  No_overlap_event_base() :
    m_type(0),
    m_ps_x(static_cast<char>(ARR_INTERIOR)),
    m_ps_y(static_cast<char>(ARR_INTERIOR)),
    m_closed(1)
  {}

  /*! Squat the content of Point_2 for the pointer of Compact Container */
  void* for_compact_container() const
  {
    return m_for_compact_container(m_point);
  }
  void for_compact_container (void* p)
  {
    m_for_compact_container(m_point, p);
  }

  /*! Initialize an event that is associated with a valid point. */
  void init(const Point_2& point, Attribute type,
            Arr_parameter_space ps_x, Arr_parameter_space ps_y)
  {
    m_point = point;
    m_type = static_cast<char>(type);
    m_ps_x = static_cast<char>(ps_x);
    m_ps_y = static_cast<char>(ps_y);
    m_closed = 1;
  }

  /*! Initialize an event associates with an open curve end. */
  void init_at_open_boundary(Attribute type,
                             Arr_parameter_space ps_x,
                             Arr_parameter_space ps_y)
  {
    m_type = static_cast<char>(type);
    m_ps_x = static_cast<char>(ps_x);
    m_ps_y = static_cast<char>(ps_y);
    m_closed = 0;
  }

  /*! Add a subcurve to the container of left curves. */
  void add_curve_to_left(Subcurve* curve)
  {
    // Look for the subcurve.
    for (Subcurve_iterator iter = m_left_curves.begin();
         iter != m_left_curves.end(); ++iter)
      if (curve == *iter) return;       // do nothing if the curve exists.

    // The curve does not exist - insert it to the container.
    m_left_curves.push_back(curve);
  }

  /*! Add a subcurve to the container of left curves (without checks). */
  void push_back_curve_to_left(Subcurve* curve)
  { m_left_curves.push_back(curve); }

  /*! Add a subcurve to the container of right curves (without checks). */
  void push_back_curve_to_right(Subcurve* curve)
  { m_right_curves.push_back(curve); }

  /*! Add a subcurve to the container of right curves. */
  std::pair<bool, Subcurve_iterator>
  add_curve_to_right(Subcurve* curve, const Gt2* tr)
  {
    if (m_right_curves.empty()) {
      m_right_curves.push_back(curve);
      return (std::make_pair(false, m_right_curves.begin()));
    }

    // Check if its an event at open boundary,
    // and if so then there is no overlap
    //(there cannot be two non-overlap curves at the same event at open
    // boundary).
    if (!this->is_closed())
      return (std::make_pair(true, m_right_curves.begin()));

    Subcurve_iterator iter = m_right_curves.begin();
    Comparison_result res;

    while ((res = tr->compare_y_at_x_right_2_object()
            (curve->last_curve(), (*iter)->last_curve(), m_point)) == LARGER)
    {
      ++iter;
      if (iter == m_right_curves.end()) {
        m_right_curves.insert(iter, curve);
        return std::make_pair(false, --iter);
      }
    }

    // Cannot overlap!
    CGAL_assertion(res != EQUAL);

    m_right_curves.insert(iter, curve);
    return std::make_pair(false, --iter);
  }

  /*! Add two Subcurves to the right of the event.
   * \pre The event does not contain any right curves, and the order of sc1
   *      and sc2 is correct.
   */
  std::pair<bool, Subcurve_iterator>
  add_curve_pair_to_right(Subcurve* sc1, Subcurve* sc2)
  {
    m_right_curves.push_back(sc1);
    m_right_curves.push_back(sc2);

    Subcurve_iterator iter = m_right_curves.end();
    --iter;
    return (std::make_pair(false, iter));
  }

  /*! Returns an iterator to the first curve to the left of the event. */
  Subcurve_iterator left_curves_begin() { return (m_left_curves.begin()); }

  /*! Return an iterator to the one past the last curve to the left
   * of the event.
   */
  Subcurve_iterator left_curves_end() { return (m_left_curves.end()); }

  /*! Remove from the left-curves container a single element.
   */
  Subcurve_iterator left_curves_erase(Subcurve_iterator it)
  { return m_left_curves.erase(it); }

  /*! Return an iterator to the first curve to the right of the event. */
  Subcurve_iterator right_curves_begin() { return (m_right_curves.begin()); }

  /*! Return an iterator to the one past the last curve to the right
   * of the event.
   */
  Subcurve_iterator right_curves_end() { return (m_right_curves.end()); }

  /*! Remove from the right-curves container a single element.
   */
  Subcurve_iterator right_curves_erase(Subcurve_iterator it)
  { return m_right_curves.erase(it); }

  /*! Returns a const iterator to the first curve to the left of the event. */
  Subcurve_const_iterator left_curves_begin() const
  { return m_left_curves.begin(); }

  /*! Return a const iterator to the past the end curve to the left
   * of the event.
   */
  Subcurve_const_iterator left_curves_end() const
  { return m_left_curves.end(); }

  /*! Return a const iterator to the first curve to the right of the event. */
  Subcurve_const_iterator right_curves_begin() const
  { return m_right_curves.begin(); }

  /*! Return a const iterator to the past the end curve to the right
   * of the event.
   */
  Subcurve_const_iterator right_curves_end() const
  { return m_right_curves.end(); }

  /*! Return a reverse_iterator to the first curve of the reversed list
   * of the right curves of the event.
   */
  Subcurve_reverse_iterator right_curves_rbegin()
  { return m_right_curves.rbegin(); }

  /*! Return a reverse_iterator to the past-end curve of the reversed list
   * of the right curves of the event.
   */
  Subcurve_reverse_iterator right_curves_rend()
  { return m_right_curves.rend(); }

  /*! Return a reverse_iterator to the first curve of the reversed list
   * of the left curves of the event.
   */
  Subcurve_reverse_iterator left_curves_rbegin()
  { return m_left_curves.rbegin(); }

  /*! Return a reverse_iterator to the past-end curve of the reversed list
   * of the left curves of the event.
   */
  Subcurve_reverse_iterator left_curves_rend()
  { return m_left_curves.rend(); }

  /*! Return the number of curves defined to the left of the event. */
  size_t number_of_left_curves() { return m_left_curves.size(); }

  /*! Return the number of curves defined to the right of the event. */
  size_t number_of_right_curves() { return (m_right_curves.size()); }

  /*! Check whether at least one curve is defined to the left of the event. */
  bool has_left_curves() const { return (! m_left_curves.empty()); }

  /*! Checks if at least one curve is defined to the right of the event. */
  bool has_right_curves() const { return (! m_right_curves.empty()); }

  /*! Obtain the actual event point (const version).
   * \pre The event is associated with a valid point.
   */
  const Point_2& point() const
  {
    CGAL_precondition(is_closed());
    return m_point;
  }

  /*! Obtain the actual event point (non-const version).
   * \pre The event is associated with a valid point.
   */
  Point_2& point()
  {
    CGAL_precondition(is_closed());
    return m_point;
  }

  /*! Obtain a curve associated with the event (const version).
   * \pre The event has incident curves.
   */
  const X_monotone_curve_2& curve() const
  {
    if (has_left_curves()) return (m_left_curves.front()->last_curve());

    CGAL_assertion(has_right_curves());
    return (m_right_curves.front()->last_curve());
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

  inline bool is_on_boundary() const
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
    Subcurve_iterator left_iter = m_left_curves.begin();
    for (InputIterator iter = begin; iter != end; ++iter, ++left_iter)
      *left_iter = *iter;
    m_left_curves.erase(left_iter, m_left_curves.end());
  }

  /*! Check if the two curves are negihbors to the left of the event. */
  bool are_left_neighbours(Subcurve* c1, Subcurve* c2)
  {
    Subcurve_iterator left_iter = m_left_curves.begin();
    for (; left_iter != m_left_curves.end(); ++left_iter) {
      if (*left_iter == c1) {
        Subcurve_iterator temp = left_iter;
        ++temp;
        if (temp != m_left_curves.end()) return (*temp == c2);
        return false;
      }

      if (*left_iter == c2) {
        Subcurve_iterator temp = left_iter;
        ++temp;
        if (temp != m_left_curves.end()) return (*temp == c1);
        return false;
      }
    }

    return false;
  }

#ifdef CGAL_SS_VERBOSE
  void Print();
#endif
};

#ifdef CGAL_SS_VERBOSE
template <typename Traits, typename Subcurve>
void No_overlap_event_base<Traits, Subcurve>::Print()
{
  std::cout << "\tEvent info: "  << "\n" ;
  if (this->is_closed()) std::cout << "\t" << m_point << "\n" ;
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
       default: CGAL_error();
      }
    }
  }
  std::cout << "\n";

  std::cout << "\tLeft curves: \n" ;
  Subcurve_iterator iter;
  for (iter = m_left_curves.begin(); iter != m_left_curves.end(); ++iter) {
    std::cout << "\t";
    (*iter)->Print();
    std::cout << "\n";
  }
  std::cout << std::endl;
  std::cout << "\tRight curves: \n" ;
  for (iter = m_right_curves.begin(); iter != m_right_curves.end(); ++iter) {
    std::cout << "\t";
    (*iter)->Print();
    std::cout << "\n";
  }

  std::cout << std::endl;
}
#endif // CGAL_SS_VERBOSE

} // namespace Surface_sweep_2
} // namespace CGAL

#endif
