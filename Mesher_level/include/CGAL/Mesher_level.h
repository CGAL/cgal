// Copyright (c) 2004-2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_MESHER_LEVEL_H
#define CGAL_MESHER_LEVEL_H

#include <string>

namespace CGAL {

enum Mesher_level_conflict_status {
  NO_CONFLICT = 0,
  CONFLICT_BUT_ELEMENT_CAN_BE_RECONSIDERED,
  CONFLICT_AND_ELEMENT_SHOULD_BE_DROPPED
};

struct Null_mesher_level {

  template <typename Visitor>
  void refine(Visitor) {}

  template <typename P, typename Z>
  Mesher_level_conflict_status test_point_conflict_from_superior(P, Z)
  {
    return NO_CONFLICT;
  }

  bool is_algorithm_done() const
  {
    return true;
  }

  template <typename Visitor>
  bool try_to_insert_one_point(Visitor)
  {
    return false;
  }

  template <typename Visitor>
  bool one_step(Visitor)
  {
    return false;
  }

  std::string debug_info() const
  {
    return "";
  }
  std::string debug_info_header() const
  {
    return "";
  }

}; // end Null_mesher_level

template <
  class Tr, /**< The triangulation type. */
  class Derived, /**< Derived class, that implements methods. */
  class Element, /**< Type of elements that this level refines. */
  class Previous, /* = Null_mesher_level, */
  /**< Previous level type, defaults to
     \c Null_mesher_level. */
  class Triangulation_traits /** Traits class that defines types for the
         triangulation. */
  >
class Mesher_level
{
public:
  /** Type of triangulation that is meshed. */
  typedef Tr Triangulation;
  /** Type of point that are inserted into the triangulation. */
  typedef typename Triangulation::Point Point;
  /** Type of vertex handles that are returns by insertions into the
      triangulation. */
  typedef typename Triangulation::Vertex_handle Vertex_handle;
  /** Type of the conflict zone for a point that can be inserted. */
  typedef typename Triangulation_traits::Zone Zone;

  typedef Element Element_type;
  typedef Previous Previous_level;

private:
  /** \name Private member functions */

  /** Curiously recurring template pattern. */
  //@{
  Derived& derived()
  {
    return static_cast<Derived&>(*this);
  }

  const Derived& derived() const
  {
    return static_cast<const Derived&>(*this);
  }
  //@}

  /** \name Private member datas */

  Previous& previous_level; /**< The previous level of the refinement
                                    process. */
public:
  typedef Mesher_level<Tr,
                       Derived,
                       Element,
                       Previous_level,
                       Triangulation_traits> Self;

  /** \name CONSTRUCTORS */
  Mesher_level(Previous_level& previous)
    : previous_level(previous)
  {
  }

  /** \name FUNCTIONS IMPLEMENTED IN THE CLASS \c Derived */

  /**  Access to the triangulation */
  Triangulation& triangulation()
  {
    return derived().triangulation_ref_impl();
  }

  /**  Access to the triangulation */
  const Triangulation& triangulation() const
  {
    return derived().triangulation_ref_impl();
  }

  const Previous_level& previous() const
  {
    return previous_level;
  }

  Vertex_handle insert(Point p, Zone& z)
  {
    return derived().insert_impl(p, z);
  }

  Zone conflicts_zone(const Point& p, Element e)
  {
    return derived().conflicts_zone_impl(p, e);
  }

  /** Called before the first refinement, to initialized the queue of
      elements that should be refined. */
  void scan_triangulation()
  {
    derived().scan_triangulation_impl();
  }

  /** Tells if, as regards the elements of type \c Element, the refinement is
      done. */
  bool no_longer_element_to_refine()
  {
    return derived().no_longer_element_to_refine_impl();
  }

  /** Retrieves the next element that could be refined. */
  Element get_next_element()
  {
    return derived().get_next_element_impl();
  }

  /** Remove from the list the next element that could be refined. */
  void pop_next_element()
  {
    derived().pop_next_element_impl();
  }

  /** Gives the point that should be inserted to refine the element \c e */
  Point refinement_point(const Element& e)
  {
    return derived().refinement_point_impl(e);
  }

  /** Actions before testing conflicts for point \c p and element \c e */
  template <typename Mesh_visitor>
  void before_conflicts(const Element& e, const Point& p,
                        Mesh_visitor visitor)
  {
    visitor.before_conflicts(e, p);
    derived().before_conflicts_impl(e, p);
  }

  /** Tells if, as regards this level of the refinement process, if the
      point conflicts with something, and do what is needed. The return
      type is made of two booleans:
        - the first one tells if the point can be inserted,
        - in case of, the first one is \c false, the second one tells if
        the tested element should be reconsidered latter.
  */
  Mesher_level_conflict_status private_test_point_conflict(const Point& p,
                                                           Zone& zone)
  {
    return derived().private_test_point_conflict_impl(p, zone);
  }

  /** Tells if, as regards this level of the refinement process, if the
      point conflicts with something, and do what is needed. The return
      type is made of two booleans:
        - the first one tells if the point can be inserted,
        - in case of, the first one is \c false, the second one tells if
        the tested element should be reconsidered latter.
      This function is called by the superior level, if any.
  */
  Mesher_level_conflict_status
  test_point_conflict_from_superior(const Point& p,
                                    Zone& zone)
  {
    return derived().test_point_conflict_from_superior_impl(p, zone);
  }

  /**
   * Actions before inserting the point \c p in order to refine the
   * element \c e. The zone of conflicts is \c zone.
   */
  template <class Mesh_visitor>
  void before_insertion(Element& e, const Point& p, Zone& zone,
                        Mesh_visitor visitor)
  {
    visitor.before_insertion(e, p, zone);
    derived().before_insertion_impl(e, p, zone);
  }

  /** Actions after having inserted the point.
   *  \param vh is the vertex handle of the inserted point,
   *  \param visitor is the visitor.
   */
  template <class Mesh_visitor>
  void after_insertion(Vertex_handle vh, Mesh_visitor visitor)
  {
    derived().after_insertion_impl(vh);
    visitor.after_insertion(vh);
  }

  /** Actions after testing conflicts for point \c p and element \c e
   *  if no point is inserted. */
  template <class Mesh_visitor>
  void after_no_insertion(const Element& e, const Point& p, Zone& zone,
                          Mesh_visitor visitor)
  {
    derived().after_no_insertion_impl(e, p, zone);
    visitor.after_no_insertion(e, p, zone);
  }

  /** \name MESHING PROCESS
   *
   * The following functions use the functions that are implemented in the
   * derived classes.
   *
   */

  /**
   * Tells it the algorithm is done, regarding elements of type \c Element
   * or elements of previous levels.
   */
  bool is_algorithm_done()
  {
    return ( previous_level.is_algorithm_done() &&
             no_longer_element_to_refine() );
  }

  /** Refines elements of this level and previous levels. */
  template <class Mesh_visitor>
  void refine(Mesh_visitor visitor)
  {
    while(! is_algorithm_done() )
    {
      previous_level.refine(visitor.previous_level());
      if(! no_longer_element_to_refine() )
        process_one_element(visitor);
    }
  }

  /**
   * This function takes one element from the queue, and try to refine
   * it. It returns \c true if one point has been inserted.
   * @todo Merge with try_to_refine_element().
   */
  template <class Mesh_visitor>
  bool process_one_element(Mesh_visitor visitor)
  {
    Element e = get_next_element();

    const Mesher_level_conflict_status result
      = try_to_refine_element(e, visitor);

    if(result == CONFLICT_AND_ELEMENT_SHOULD_BE_DROPPED)
      pop_next_element();
    return result == NO_CONFLICT;
  }

  template <class Mesh_visitor>
  Mesher_level_conflict_status
  try_to_refine_element(Element e, Mesh_visitor visitor)
  {
    const Point& p = refinement_point(e);

    before_conflicts(e, p, visitor);

    Zone zone = conflicts_zone(p, e);

    const Mesher_level_conflict_status result = test_point_conflict(p, zone);
#ifdef CGAL_MESHES_DEBUG_REFINEMENT_POINTS
    std::cerr << "(" << p << ") ";
    switch( result )
    {
    case NO_CONFLICT:
      std::cerr << "accepted\n";
      break;
    case CONFLICT_BUT_ELEMENT_CAN_BE_RECONSIDERED:
      std::cerr << "rejected (temporarily)\n";
      break;
    case CONFLICT_AND_ELEMENT_SHOULD_BE_DROPPED:
      std::cerr << "rejected (permanent)\n";
      break;
    }
#endif

    if(result == NO_CONFLICT)
    {
      before_insertion(e, p, zone, visitor);

      Vertex_handle v = insert(p, zone);

      after_insertion(v, visitor);

      return NO_CONFLICT;
    }
    else
      after_no_insertion(e, p, zone, visitor);
    return result;
  }

  /** Return (can_split_the_element, drop_element). */
  Mesher_level_conflict_status
  test_point_conflict(const Point& p, Zone& zone)
  {
    const Mesher_level_conflict_status result =
      previous_level.test_point_conflict_from_superior(p, zone);

    if( result != NO_CONFLICT )
      return result;
    return private_test_point_conflict(p, zone);
  }

  /** \name STEP BY STEP FUNCTIONS */

  /**
   * Inserts exactly one point, if possible, and returns \c false if no
   * point has been inserted because the algorithm is done.
   */
  template <class Mesh_visitor>
  bool try_to_insert_one_point(Mesh_visitor visitor)
  {
    while(! is_algorithm_done() )
    {
      if( previous_level.try_to_insert_one_point(visitor.previous_level()) )
        return true;
      if(! no_longer_element_to_refine() )
        if( process_one_element(visitor) )
          return true;
    }
    return false;
  }

  /**
   * Applies one step of the algorithm: tries to refine one element of
   * previous level or one element of this level. Return \c false iff
   * <tt> is_algorithm_done()==true </tt>.
   */
  template <class Mesh_visitor>
  bool one_step(Mesh_visitor visitor)
  {
    if( ! previous_level.is_algorithm_done() )
      previous_level.one_step(visitor.previous_level());
    else
      if( ! no_longer_element_to_refine() )
        process_one_element(visitor);
    return ! is_algorithm_done();
  }

}; // end Mesher_level

}  // end namespace CGAL

#include <CGAL/Mesher_level_visitors.h>
#include <CGAL/Mesher_level_default_implementations.h>

#endif // CGAL_MESHER_LEVEL_H
