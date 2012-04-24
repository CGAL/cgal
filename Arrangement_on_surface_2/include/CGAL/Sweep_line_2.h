// Copyright (c) 2005,2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 (based on old version by Tali Zvi)

#ifndef CGAL_SWEEP_LINE_2_H
#define CGAL_SWEEP_LINE_2_H

/*! \file
 * Definition of the Sweep_line_2 class.
 */

#include <list>
#include <CGAL/Object.h>
#include <CGAL/Basic_sweep_line_2.h>
#include <CGAL/Sweep_line_2/Sweep_line_curve_pair.h>
#include <CGAL/Arrangement_2/Open_hash.h>

namespace CGAL {

/*! \class
 * Sweep_line_2 is a class that implements the sweep line algorithm based
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

template < class Traits_,
           class Visitor_,
           class Subcurve_ = Sweep_line_subcurve<Traits_>,
           class Event_ = Sweep_line_event<Traits_, Subcurve_>,
           typename Allocator_ = CGAL_ALLOCATOR(int) >
class Sweep_line_2 : public Basic_sweep_line_2<Traits_,
                                               Visitor_,
                                               Subcurve_,
                                               Event_,
                                               Allocator_>
{
public:

  typedef Traits_                                        Traits_2;
  typedef Visitor_                                       Visitor;
  typedef Event_                                         Event;
  typedef Subcurve_                                      Subcurve;
  typedef Allocator_                                     Allocator;

  typedef Basic_sweep_line_2<Traits_2,
                             Visitor,
                             Subcurve,
                             Event,
                             Allocator>                  Base;

  typedef typename Base::Traits_adaptor_2                Traits_adaptor_2;
  typedef typename Traits_adaptor_2::Point_2             Point_2;
  typedef typename Traits_adaptor_2::X_monotone_curve_2  X_monotone_curve_2;
 
  typedef typename Base::Event_queue_iterator         Event_queue_iterator;
  typedef typename Event::Subcurve_iterator           Event_subcurve_iterator;

  typedef typename Base::Base_event                   Base_event;
  typedef typename Base_event::Attribute              Attribute;

  typedef typename Base::Base_subcurve                Base_subcurve;
  
  typedef std::list<Subcurve*>                        Subcurve_container;
  typedef typename Subcurve_container::iterator       Subcurve_iterator; 

  typedef typename Base::Status_line_iterator         Status_line_iterator;
  
  typedef CGAL::Curve_pair<Subcurve>                       Curve_pair;
  typedef CGAL::Curve_pair_hasher<Subcurve>                Curve_pair_hasher;
  typedef CGAL::Equal_curve_pair<Subcurve>                 Equal_curve_pair;
  typedef Open_hash<Curve_pair,
                    Curve_pair_hasher,
                    Equal_curve_pair>                 Curve_pair_set;

  typedef random_access_input_iterator<std::vector<Object> >
                                                      vector_inserter;

protected:

  // Data members:
  Subcurve_container  m_overlap_subCurves;
                                     // Contains all of the new sub-curves
                                     // creaed by an overlap.

  Curve_pair_set m_curves_pair_set;  // A lookup table of pairs of Subcurves
                                     // that have been intersected.

  std::vector<Object> m_x_objects;   // Auxiliary vector for storing the
                                     // intersection objects.

  X_monotone_curve_2  sub_cv1;       // Auxiliary varibales
  X_monotone_curve_2  sub_cv2;       // (for splitting curves).

public:

  /*!
   * Constructor.
   * \param visitor A pointer to a sweep-line visitor object.
   */
  Sweep_line_2 (Visitor * visitor) :
    Base(visitor),
    m_curves_pair_set(0)
  {}


  /*!
   * Constructor.
   * \param traits A pointer to a sweep-line traits object.
   * \param visitor A pointer to a sweep-line visitor object.
   */
  Sweep_line_2 (const Traits_2 * traits, Visitor * visitor) :
    Base(traits, visitor),
    m_curves_pair_set(0)
  {}

  /*! Destrcutor. */
  virtual ~Sweep_line_2()
  {}

protected:

  /*! Initialize the data structures for the sweep-line algorithm. */
  virtual void _init_structures ();

  /*! Complete the sweep process (complete the data structures). */
  virtual void _complete_sweep ();

  /*! Handle the subcurves to the left of the current event point. */
  virtual void _handle_left_curves ();

  /*! Handle the subcurves to the right of the current event point. */
  virtual void _handle_right_curves ();

  /*!
   * Add a subcurve to the right of an event point.
   * \param event The event point.
   * \param curve The subcurve to add.
   * \return (true) if an overlap occured; (false) otherwise.
   */
  virtual bool _add_curve_to_right (Event * event, Subcurve * curve,
                                    bool overlap_exist = false);

  /*! Fix overlapping subcurves before handling the current event. */
  void _fix_overlap_subcurves();

  /*!
   * Handle overlap at right insertion to event.
   * \param event The event point.
   * \param curve The subcurve representing the overlap.
   * \param iter An iterator for the curves.
   * \param overlap_exist
   */
  void _handle_overlap (Event * event, Subcurve * curve,
                        Event_subcurve_iterator iter, bool overlap_exist);
  
  /*! 
   * Compute intersections between the two given curves.
   * If the two curves intersect, create a new event (or use the event that 
   * already exits in the intersection point) and insert the curves to the
   * event.
   * \param curve1 The first curve.
   * \param curve2 The second curve.
   */ 
  void _intersect (Subcurve * c1, Subcurve * c2);

  /*!
   * When a curve is removed from the status line for good, its top and
   * bottom neighbors become neighbors. This method finds these cases and
   * looks for the intersection point, if one exists.
   * \param leftCurve A pointer to the curve that is about to be deleted.
   * \param remove_for_good Whether the aubcurve is removed for good.
   */
  void _remove_curve_from_status_line (Subcurve * leftCurve,
                                       bool remove_for_good);

  /*!
   * Create an intersection-point event between two curves.
   * \param xp The intersection point.
   * \param mult Its multiplicity.
   * \param curve1 The first curve.
   * \param curve2 The second curve.
   * \param is_overlap Whether the two curves overlap at xp.
   */
  void _create_intersection_point (const Point_2 & xp,
                                   unsigned int mult,
                                   Subcurve *& c1,
                                   Subcurve *& c2,
                                   bool is_overlap = false);

  /*!
   * Fix a subcurve that represents an overlap.
   * \param sc The subcurve.
   */
  void _fix_finished_overlap_subcurve (Subcurve * sc);

};

} //namespace CGAL

// The member-function definitions can be found in:
#include <CGAL/Sweep_line_2/Sweep_line_2_impl.h>

#endif
