// Copyright (c) 2025 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Efi Fogel <efifogel@gmail.com>

#ifndef CGAL_DO_INTERSECT_SURFACE_SWEEP_2_H
#define CGAL_DO_INTERSECT_SURFACE_SWEEP_2_H

#include <CGAL/license/Surface_sweep_2.h>

/*! \file
 *
 * Definition of the Do_intersect_surface_sweep_2 class.
 */

#include <list>
#include <vector>

#include <CGAL/No_intersection_surface_sweep_2.h>
#include <CGAL/Surface_sweep_2/Random_access_output_iterator.h>
#include <CGAL/algorithm.h>

namespace CGAL {
namespace Surface_sweep_2 {

/*! \class
 *
 * Do_intersect_surface_sweep_2 is a class that implements the Bentley and Ottmann sweep line algorithm.
 * It detects intersections.
 */

template <typename Visitor_>
class Do_intersect_surface_sweep_2 : public No_intersection_surface_sweep_2<Visitor_> {
public:
  using Visitor = Visitor_;

private:
  using Base = No_intersection_surface_sweep_2<Visitor>;

public:
  using Geometry_traits_2 = typename Base::Geometry_traits_2;
  using Event = typename Base::Event;
  using Subcurve = typename Base::Subcurve;
  using Allocator = typename Base::Allocator;

  using Traits_adaptor_2 = typename Base::Traits_adaptor_2;
  using Point_2 = typename Traits_adaptor_2::Point_2;
  using X_monotone_curve_2 = typename Traits_adaptor_2::X_monotone_curve_2;
  using Multiplicity = typename Traits_adaptor_2::Multiplicity;

  using Event_queue_iterator = typename Base::Event_queue_iterator;
  using Event_subcurve_iterator = typename Event::Subcurve_iterator;

  using Attribute = typename Event::Attribute;

  using Subcurve_container = std::list<Subcurve*>;
  using Subcurve_iterator = typename Subcurve_container::iterator;

  using Status_line_iterator = typename Base::Status_line_iterator;

  using Intersection_point = std::pair<Point_2, Multiplicity>;
  using Intersection_result = std::variant<Intersection_point, X_monotone_curve_2>;
  using Intersection_vector = std::vector<Intersection_result>;
  using vector_inserter = Random_access_output_iterator<Intersection_vector>;
  using Subcurve_alloc = typename Base::Subcurve_alloc;

protected:
  using All_sides_oblivious_category = typename Base::All_sides_oblivious_category;

  // Parameter spaces that are either oblivious or open cannot have points
  // on the boundary, and in particular intersection points. In other words,
  // intersection points of an oblivious or open parameter space is interior
  // by definition.
  using Sides_category = typename Base::Sides_category;

public:
  /*! constructs.
   * \param visitor A pointer to a sweep-line visitor object.
   */
  Do_intersect_surface_sweep_2(Visitor* visitor, bool closed = true) : Base(visitor), m_closed(closed) {}

  /*! constructs.
   * \param traits A pointer to a sweep-line traits object.
   * \param visitor A pointer to a sweep-line visitor object.
   */
  Do_intersect_surface_sweep_2(const Geometry_traits_2* traits, Visitor* visitor, bool closed = true) :
    Base(traits, visitor),
    m_closed(closed)
  {}

  /*! destructs. */
  virtual ~Do_intersect_surface_sweep_2() {}

  /*! runs the sweep-line algorithm on a given range of x-monotone curves.
   * \param curves_begin An iterator for the first curve in the range.
   * \param curves_end A past-the-end iterator for the range.
   * \pre The value-type of CurveInputIterator is X_monotone_curve_2.
   */
  template <typename CurveInputIterator>
  void do_intersect_sweep(CurveInputIterator curves_begin, CurveInputIterator curves_end) {
    this->m_visitor->before_sweep();
    bool overlap = _init_do_intersect_sweep(curves_begin, curves_end);
    if (! overlap) this->_sweep();
    this->_complete_sweep();
    this->m_visitor->after_sweep();
  }

  /*! runs the sweep-line algorithm on a range of x-monotone curves and a range
   * of action event points (if a curve passed through an action point, it will
   * be split).
   * \param curves_begin  An iterator for the first x-monotone curve in the
   *                      range.
   * \param curves_end A past-the-end iterator for this range.
   * \param action_points_begin An iterator for the first point in the range.
   * \param action_points_end A past-the-end iterator for this range.
   * \pre The value-type of XCurveInputIterator is the traits-class
   *      X_monotone_curve_2, and the value-type of PointInputIterator is the
   *      traits-class Point_2.
   */
  template <typename CurveInputIterator, class PointInputIterator>
  void do_intersect_sweep(CurveInputIterator curves_begin,
                          CurveInputIterator curves_end,
                          PointInputIterator action_points_begin,
                          PointInputIterator action_points_end) {
    this->m_visitor->before_sweep();
    bool overlap = _init_do_intersect_sweep(curves_begin, curves_end);
    if (! overlap) overlap = _init_do_intersect_points(action_points_begin, action_points_end, Event::ACTION);
    if (! overlap) this->_sweep();
    this->_complete_sweep();
    this->m_visitor->after_sweep();
  }

protected:
  using Subcurve_vector = typename std::vector<Subcurve*>;

  /*! creates an event object for each input point.
   * \return (true) if an overlap occurred; (false) otherwise.
   */
  template <typename PointInputIterator>
  bool _init_do_intersect_points(PointInputIterator points_begin, PointInputIterator points_end, Attribute type) {
    for (auto pit = points_begin; pit != points_end; ++pit) {
      auto overlap = this->_init_point(*pit, type);
      if (overlap) return true;
    }
    return false;
  }

  /*! creates a Subcurve object and two Event objects for each curve.
   * \return (true) if an overlap occurred; (false) otherwise.
   */
  template <typename CurveInputIterator>
  bool _init_do_intersect_curves(CurveInputIterator curves_begin, CurveInputIterator curves_end) {
    std::size_t index = 0;
    for (auto cit = curves_begin; cit != curves_end; ++cit, ++index) {
      auto overlap = this->_init_curve(*cit);
      if (overlap) return true;
    }
    return false;
  }

  /*! initializes the sweep algorithm.
   */
  template <typename CurveInputIterator>
  bool _init_do_intersect_sweep(CurveInputIterator curves_begin, CurveInputIterator curves_end) {
    this->m_num_subcurves = std::distance(curves_begin, curves_end);
    this->_init_structures();
    bool overlap = _init_do_intersect_curves(curves_begin, curves_end);     // initialize the curves
    return overlap;
  }

  /*! completes the sweep (complete data structures). */
  virtual void _complete_sweep();

  /*! Handle the subcurves to the left of the current event point. */
  virtual void _handle_left_curves();

  /*! Handle the subcurves to the right of the current event point. */
  virtual void _handle_right_curves();

  /*! Add a subcurve to the right of an event point.
   * \param event The event point.
   * \param curve The subcurve to add.
   */
  virtual bool _add_curve_to_right(Event* event, Subcurve* curve);

  /*! Compute intersections between the two given curves.
   * If the two curves intersect, create a new event (or use the event that
   * already exits in the intersection point) and insert the curves to the
   * event.
   * \param curve1 The first curve.
   * \param curve2 The second curve.
   */
  bool _do_intersect(Subcurve* c1, Subcurve* c2);

  /*! When a curve is removed from the status line for good, its top and
   * bottom neighbors become neighbors. This method finds these cases and
   * looks for the intersection point, if one exists.
   * \param leftCurve A pointer to the curve that is about to be deleted.
   * \param remove_for_good Whether the aubcurve is removed for good.
   */
  bool _remove_curve_from_status_line(Subcurve* leftCurve);

private:
  bool m_closed;
};

} // namespace Surface_sweep_2
} // namespace CGAL

// The member-function definitions can be found in:
#include <CGAL/Surface_sweep_2/Do_intersect_surface_sweep_2_impl.h>

#endif
