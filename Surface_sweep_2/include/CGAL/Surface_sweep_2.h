// Copyright (c) 2005,2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Baruch Zukerman <baruchzu@post.tau.ac.il>
//             Efi Fogel <efifogel@gmail.com>
//               (based on old version by Tali Zvi)

#ifndef CGAL_SURFACE_SWEEP_2_H
#define CGAL_SURFACE_SWEEP_2_H

#include <CGAL/license/Surface_sweep_2.h>

/*! \file
 *
 * Definition of the Surface_sweep_2 class.
 */

#include <list>
#include <vector>

#include <CGAL/Object.h>
#include <CGAL/No_intersection_surface_sweep_2.h>
#include <CGAL/Surface_sweep_2/Random_access_output_iterator.h>
#include <CGAL/algorithm.h>

namespace CGAL {
namespace Surface_sweep_2 {

/*! \class
 *
 * Surface_sweep_2 is a class that implements the sweep line algorithm based
 * on the algorithm of Bentley and Ottmann.
 * It extends the algorithm to support not only segments but general x-monotone
 * curves as well and isolated points.
 * The X_monotone_curve_2 type and Point_2 are defined by the traits class that
 * is one of the template arguments.
 *
 * The algorithm is also extended to support the following degenerate cases:
 * - vertical segments
 * - multiple (more then two) curves intersecting at one point
 * - curves beginning and ending on other curves.
 * - overlapping curves
 *
 * General flow:
 * After the initialization stage, the events are handled from left to right.
 *
 * For each event
 *
 *  First pass - handles special cases in which curves start or end
 *               at the interior of another curve
 *  Handle left curves - iterate over the curves that intersect
 *               at the event point and defined lexicographically to the left
 *               of the event.
 *  Handle right curves - iterate over the curves that intersect
 *               the event point and defined lexicographically to the right
 *               of the event point. This is where new intersection points
 *               are calculated.
 * End
 *
 * Convensions through out the code:
 * In order to make the code as readable as possible, some convensions were
 * made in regards to variable naming:
 *
 * xp - is the intersection point between two curves
 * slIter - an iterator to the status line, always points to a curve.
 *
 */

template <typename Visitor_>
class Surface_sweep_2 : public No_intersection_surface_sweep_2<Visitor_> {
public:
  typedef Visitor_                                      Visitor;

private:
  typedef No_intersection_surface_sweep_2<Visitor>      Base;

public:
  typedef typename Base::Geometry_traits_2              Geometry_traits_2;
  typedef typename Base::Event                          Event;
  typedef typename Base::Subcurve                       Subcurve;
  typedef typename Base::Allocator                      Allocator;

  typedef typename Base::Traits_adaptor_2               Traits_adaptor_2;
  typedef typename Traits_adaptor_2::Point_2            Point_2;
  typedef typename Traits_adaptor_2::X_monotone_curve_2 X_monotone_curve_2;
  typedef typename Traits_adaptor_2::Multiplicity       Multiplicity;

  typedef typename Base::Event_queue_iterator           Event_queue_iterator;
  typedef typename Event::Subcurve_iterator             Event_subcurve_iterator;

  typedef typename Event::Attribute                     Attribute;

  typedef std::list<Subcurve*>                          Subcurve_container;
  typedef typename Subcurve_container::iterator         Subcurve_iterator;

  typedef typename Base::Status_line_iterator           Status_line_iterator;

  typedef std::pair<Point_2, Multiplicity>              Intersection_point;
  typedef boost::variant<Intersection_point, X_monotone_curve_2>
                                                        Intersection_result;
  typedef std::vector<Intersection_result>              Intersection_vector;
  typedef Random_access_output_iterator<Intersection_vector>
                                                        vector_inserter;
  typedef typename Base::Subcurve_alloc                 Subcurve_alloc;
protected:
  // Data members:
  Subcurve_container m_overlap_subCurves;
                                     // Contains all of the new sub-curves
                                     // creaed by an overlap.

  Intersection_vector m_x_objects;   // Auxiliary vector for storing the
                                     // intersection objects.

  X_monotone_curve_2 sub_cv1;        // Auxiliary varibales
  X_monotone_curve_2 sub_cv2;        // (for splitting curves).

public:
  /*! Constructor.
   * \param visitor A pointer to a sweep-line visitor object.
   */
  Surface_sweep_2(Visitor* visitor) : Base(visitor) {}

  /*!
   * Construct.
   * \param traits A pointer to a sweep-line traits object.
   * \param visitor A pointer to a sweep-line visitor object.
   */
  Surface_sweep_2(const Geometry_traits_2* traits, Visitor* visitor) :
    Base(traits, visitor)
  {}

  /*! Destrcut. */
  virtual ~Surface_sweep_2() {}

protected:
  typedef typename std::vector<Subcurve*>               Subcurve_vector;

  /*! Initialize the data structures for the sweep-line algorithm. */
  virtual void _init_structures();

  /*! Complete the sweep process (complete the data structures). */
  virtual void _complete_sweep();

  /*! Handle the subcurves to the left of the current event point. */
  virtual void _handle_left_curves();

  /*! Handle the overlap on the right curves of the current event point. */
  void _handle_overlaps_in_right_curves();

  /*! clip the last curve of a subcurve if it is not in the status line and with a left end not being the current event*/
  void _clip_non_active_curve_at_current_event(Subcurve*);

  /*! Handle the subcurves to the right of the current event point. */
  virtual void _handle_right_curves();

  /*! Add a subcurve to the right of an event point.
   * \param event The event point.
   * \param curve The subcurve to add.
   */
  virtual bool _add_curve_to_right(Event* event, Subcurve* curve);

  /*! Add a curve as a right curve or left curve when the event is created
   * or updated.
   */
  void _add_curve(Event* e, Subcurve* sc, Attribute type);

  /*! create an overlap subcurve from overlap_cv between c1 and c2.
   * \param overlap_cv the overlapping curve.
   * \param c1 first subcurve contributing to the overlap.
   * \param c2 second subcurve contributing to the overlap.
   * \param all_leaves_diff not empty in case c1 and c2 have common ancesters.
   *                        It contains the set of curves  not contained in first_parent
   *                        that are in the other subcurve
   * \param first_parent only used when c1 and c2 have common ancesters.
   *                     It is either c1 or c2 (the one having the more leaves)
   *
   */
  void _create_overlapping_curve(const X_monotone_curve_2& overlap_cv,
                                 Subcurve*& c1 , Subcurve*& c2,
                                 const Subcurve_vector& all_leaves_diff,
                                 Subcurve* first_parent,
                                 Event* event_on_overlap);

  /*! Compute intersections between the two given curves.
   * If the two curves intersect, create a new event (or use the event that
   * already exits in the intersection point) and insert the curves to the
   * event.
   * \param curve1 The first curve.
   * \param curve2 The second curve.
   */
  void _intersect(Subcurve* c1, Subcurve* c2, Event* event_for_overlap = nullptr);

  /*! When a curve is removed from the status line for good, its top and
   * bottom neighbors become neighbors. This method finds these cases and
   * looks for the intersection point, if one exists.
   * \param leftCurve A pointer to the curve that is about to be deleted.
   * \param remove_for_good Whether the aubcurve is removed for good.
   */
  void _remove_curve_from_status_line(Subcurve* leftCurve,
                                      bool remove_for_good);

  /*! Create an intersection-point event between two curves.
   * \param xp The intersection point.
   * \param mult Its multiplicity.
   * \param curve1 The first curve.
   * \param curve2 The second curve.
   */
  void _create_intersection_point(const Point_2& xp,
                                  unsigned int mult,
                                  Subcurve*& c1,
                                  Subcurve*& c2);
};

} // namespace Surface_sweep_2
} // namespace CGAL

// The member-function definitions can be found in:
#include <CGAL/Surface_sweep_2/Surface_sweep_2_impl.h>

#endif
