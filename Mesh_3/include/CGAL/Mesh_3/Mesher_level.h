// Copyright (c) 2004-2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent RINEAU, Clement JAMIN

#ifndef CGAL_MESH_3_MESHER_LEVEL_H
#define CGAL_MESH_3_MESHER_LEVEL_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/disable_warnings.h>

#include <string>

#ifdef CGAL_MESH_3_PROFILING
  #include <CGAL/Mesh_3/Profiling_tools.h>
#endif

#include <CGAL/atomic.h>
#include <CGAL/Mesh_3/Worksharing_data_structures.h>

#ifdef CGAL_CONCURRENT_MESH_3_PROFILING
# define CGAL_PROFILE
# include <CGAL/Profile_counter.h>
#endif

#include <algorithm>

#ifdef CGAL_LINKED_WITH_TBB
# include <tbb/task_group.h>
#endif

#include <string>

namespace CGAL { namespace Mesh_3 {

enum Mesher_level_conflict_status {
  NO_CONFLICT = 0
  , CONFLICT_BUT_ELEMENT_CAN_BE_RECONSIDERED
  , CONFLICT_AND_ELEMENT_SHOULD_BE_DROPPED
  , THE_FACET_TO_REFINE_IS_NOT_IN_ITS_CONFLICT_ZONE
  , ELEMENT_WAS_A_ZOMBIE
  , COULD_NOT_LOCK_ZONE
  , COULD_NOT_LOCK_ELEMENT
};

/************************************************
 *
 * Null_mesher_level class
 *
 ************************************************/

struct Null_mesher_level {

  template <typename Visitor>
  void refine(Visitor) {}

  // For sequential version
  template <typename P, typename Z>
  Mesher_level_conflict_status test_point_conflict_from_superior(P, Z)
  {
    return NO_CONFLICT;
  }
  // For parallel version
  template <typename P, typename Z, typename MV>
  Mesher_level_conflict_status test_point_conflict_from_superior(P, Z, MV &)
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

  //==============================================
  // For parallel version
  void add_to_TLS_lists(bool) {}
  void splice_local_lists() {}
  template <typename Mesh_visitor>
  void before_next_element_refinement_in_superior(Mesh_visitor) {}
  void before_next_element_refinement() {}
  //==============================================

  std::string debug_info_class_name_impl() const
  {
    return "Null_mesher_level";
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


/************************************************
 *
 * Mesher_level_base class
 *
 ************************************************/

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
class Mesher_level_base
{
public:
  /** Type of triangulation that is meshed. */
  typedef Tr Triangulation;
  /** Type of point that are inserted into the triangulation. */
  typedef typename Triangulation::Point Point;
  /** Type of point with no weight */
  typedef typename Triangulation::Bare_point Bare_point;
  /** Type of vertex handles that are returns by insertions into the
      triangulation. */
  typedef typename Triangulation::Vertex_handle Vertex_handle;
  /** Type of lock data structure for concurrency */
  typedef typename Triangulation::Lock_data_structure Lock_data_structure;
  /** Type of facet & cell handles */
  typedef typename Triangulation::Cell_handle Cell_handle;
  typedef typename Cell_handle::value_type Cell;
  typedef typename Triangulation::Facet Facet;
  /** Type of the conflict zone for a point that can be inserted. */
  typedef typename Triangulation_traits::Zone Zone;

  typedef Element Element_type;
  typedef Previous Previous_level;

protected:
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

  /// debug info: class name
  std::string debug_info_class_name() const
  {
    return derived().debug_info_class_name_impl();
  }

  std::string debug_info_element(const Element &e) const
  {
    return derived().debug_info_element_impl(e);
  }

  /** \name Private member datas */

  Previous_level& previous_level; /**< The previous level of the refinement
                                    process. */

#ifdef CGAL_MESH_3_PROFILING
protected:
  WallClockTimer m_timer;
#endif

public:
  typedef Mesher_level_base<Tr,
                       Derived,
                       Element,
                       Previous_level,
                       Triangulation_traits> Self;

  /** \name CONSTRUCTORS */
  Mesher_level_base(Previous_level& previous)
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

  void clear_refinement_queue()
  {
    derived().clear(); // Clear refinement queue
  }

  /** Called before the first refinement, to initialized the queue of
      elements that should be refined. */
  void scan_triangulation()
  {
    //derived().clear(); // Clear refinement queue
    derived().scan_triangulation_impl();

#if defined(CGAL_MESH_3_USE_LAZY_UNSORTED_REFINEMENT_QUEUE)\
 && defined(CGAL_MESH_3_IF_UNSORTED_QUEUE_JUST_SORT_AFTER_SCAN)
    std::cerr << "Sorting...";
    derived().sort();
    std::cerr << " done." << std::endl;
#endif
  }

  /** For diagnostics. */
  int number_of_bad_elements()
  {
    return derived().number_of_bad_elements_impl();
  }

  /** Tells if, as regards the elements of type \c Element, the refinement is
      done. */
  bool no_longer_element_to_refine()
  {
    return derived().no_longer_element_to_refine_impl();
  }

  /** It includes zombie elements */
  int number_of_elements_in_queue()
  {
    return derived().size();
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

  Bare_point circumcenter_of_element(const Element& e)
  {
    return derived().circumcenter_impl(e);
  }

  template <typename Mesh_visitor>
  void before_next_element_refinement_in_superior(Mesh_visitor visitor)
  {
    derived().before_next_element_refinement_in_superior_impl(visitor);
  }

  template <typename Mesh_visitor>
  void before_next_element_refinement(Mesh_visitor visitor)
  {
    derived().before_next_element_refinement_impl();
    previous_level.before_next_element_refinement_in_superior(
                                        visitor.previous_level());
  }

  /** Gives the point that should be inserted to refine the element \c e */
  Bare_point refinement_point(const Element& e)
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
#ifdef CGAL_MESH_3_PROFILING
    bool done = ( previous_level.is_algorithm_done() &&
                  no_longer_element_to_refine() );
    /*if (done)
    {
      std::cerr << "done in " << m_timer.elapsed() << " seconds." << std::endl;
      m_timer.reset();
    }*/
    return done;
#else
    return ( previous_level.is_algorithm_done() &&
             no_longer_element_to_refine() );
#endif
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
      = derived().try_to_refine_element(e, visitor);

    if(result == CONFLICT_AND_ELEMENT_SHOULD_BE_DROPPED)
    {
      pop_next_element();
    }

    return result == NO_CONFLICT;
  }

  void get_valid_vertices_of_element(const Cell_handle &e, Vertex_handle vertices[4]) const
  {
    for (int i = 0 ; i < 4 ; ++i)
      vertices[i] = e->vertex(i);
  }
  // Among the 4 values, one of them will be Vertex_handle() (~= nullptr)
  void get_valid_vertices_of_element(const Facet &e, Vertex_handle vertices[4]) const
  {
    for (int i = 0 ; i < 4 ; ++i)
      vertices[i] = (i != e.second ? e.first->vertex(i) : Vertex_handle());
  }

  Cell_handle get_cell_from_element(const Cell_handle &e) const
  {
    return e;
  }
  Cell_handle get_cell_from_element(const Facet &e) const
  {
    return e.first;
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

}; // end Mesher_level_base


/************************************************
// Class Mesher_level
// Two versions: sequential / parallel
************************************************/
// Sequential
template <
  class Tr, /**< The triangulation type. */
  class Derived, /**< Derived class, that implements methods. */
  class Element, /**< Type of elements that this level refines. */
  class Previous, /* = Null_mesher_level, */
  /**< Previous level type, defaults to
     \c Null_mesher_level. */
  class Triangulation_traits, /** Traits class that defines types for the
         triangulation. */
  typename Concurrency_tag>
class Mesher_level
  : public Mesher_level_base<Tr, Derived, Element, Previous,
                             Triangulation_traits>
{
public:

  typedef Mesher_level<Tr,
                       Derived,
                       Element,
                       Previous,
                       Triangulation_traits,
                       Concurrency_tag> Self;

  typedef Mesher_level_base<Tr,
                            Derived,
                            Element,
                            Previous,
                            Triangulation_traits> Base;

  typedef typename Base::Lock_data_structure Lock_data_structure;
  typedef typename Base::Zone                Zone;
  typedef typename Base::Point               Point;
  typedef typename Base::Vertex_handle       Vertex_handle;
  using Base::derived;
  using Base::is_algorithm_done;
  using Base::triangulation;
  using Base::insert;
  using Base::before_conflicts;
  using Base::previous_level;
  using Base::no_longer_element_to_refine;
  using Base::pop_next_element;
  using Base::debug_info_class_name;
  using Base::debug_info_element;

  /** \name CONSTRUCTORS */

  Mesher_level(Previous& previous)
    : Base(previous)
  {
  }

  void add_to_TLS_lists(bool) {}
  void splice_local_lists() {}
  bool no_longer_local_element_to_refine() { return true; }
  Element get_next_local_element() {return Element(); }
  void pop_next_local_element() {}

  Zone conflicts_zone(const Point& p
                      , Element e
                      , bool &facet_is_in_its_cz)
  {
    return derived().conflicts_zone_impl(p, e, facet_is_in_its_cz);
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
  test_point_conflict_from_superior(const Point& p, Zone& zone)
  {
    return derived().test_point_conflict_from_superior_impl(p, zone);
  }

  /** Refines elements of this level and previous levels (SEQUENTIAL VERSION). */
  template <class Mesh_visitor>
  void refine(Mesh_visitor visitor)
  {
    while(! is_algorithm_done() )
    {
      previous_level.refine(visitor.previous_level());
      if(! no_longer_element_to_refine() )
      {
        this->process_one_element(visitor);
      }
    }
  }

  template <class Mesh_visitor>
  Mesher_level_conflict_status
  try_to_refine_element(Element e, Mesh_visitor visitor)
  {
    const Tr& tr = derived().triangulation_ref_impl();
    const Point& p = tr.geom_traits().construct_weighted_point_3_object()(
                       this->refinement_point(e));

#ifdef CGAL_MESH_3_VERY_VERBOSE
    std::cerr << "Trying to insert point: " << p <<
      " inside element " << debug_info_element(e) << std::endl;
#endif

    Mesher_level_conflict_status result;
    Zone zone;

    before_conflicts(e, p, visitor);

    bool facet_is_in_its_cz = true;
    zone = conflicts_zone(p, e, facet_is_in_its_cz);
    if (!facet_is_in_its_cz)
      result = THE_FACET_TO_REFINE_IS_NOT_IN_ITS_CONFLICT_ZONE;
    else
      result = test_point_conflict(p, zone);

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
    case THE_FACET_TO_REFINE_IS_NOT_IN_ITS_CONFLICT_ZONE:
      std::cerr << "the facet to refine was not in the conflict zone "
        "(switching to exact)\n";
      break;
    case ELEMENT_WAS_A_ZOMBIE:
      std::cerr << "element was a zombie\n";
      break;
    case   COULD_NOT_LOCK_ELEMENT:
      std::cerr << "could not lock element\n";
      break;
    case COULD_NOT_LOCK_ZONE:
      std::cerr << "could not lock zone\n";
      break;
    }
#endif

    if(result == NO_CONFLICT)
    {
      this->before_insertion(e, p, zone, visitor);

      Vertex_handle vh = insert(p, zone);

      this->after_insertion(vh, visitor);
    }
    else
    {
      this->after_no_insertion(e, p, zone, visitor);
    }

    return result;
  }

  /** Refines elements of this level and previous levels.
  *   Stops when algorithm is done
  *   or when num vertices > approx_max_num_mesh_vertices
  */
  template <class Mesh_visitor>
  void
  refine_sequentially_up_to_N_vertices(Mesh_visitor visitor,
                                            int approx_max_num_mesh_vertices)
  {
    while(! is_algorithm_done()
      && triangulation().number_of_vertices() < approx_max_num_mesh_vertices)
    {
      previous_level.refine(visitor.previous_level());
      if(! no_longer_element_to_refine() )
      {
        this->process_one_element(visitor);
      }
    }
  }

  /** Return (can_split_the_element, drop_element). */
  Mesher_level_conflict_status
  test_point_conflict(const Point& p, Zone& zone)
  {
    const Mesher_level_conflict_status result =
      previous_level.test_point_conflict_from_superior(p, zone);

    if( result != NO_CONFLICT )
      return result;

    return this->private_test_point_conflict(p, zone);
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
    else if( ! no_longer_element_to_refine() )
    {
      this->process_one_element(visitor);
    }
    return ! is_algorithm_done();
  }

  // Useless here
  void set_lock_ds(Lock_data_structure *) {}
  void set_worksharing_ds(WorksharingDataStructureType *) {}
#ifndef CGAL_NO_ATOMIC
  void set_stop_pointer(CGAL::cpp11::atomic<bool>*) {}
#endif

protected:
};

// Parallel
#ifdef CGAL_LINKED_WITH_TBB
template <
  class Tr, /**< The triangulation type. */
  class Derived, /**< Derived class, that implements methods. */
  class Element, /**< Type of elements that this level refines. */
  class Previous, /* = Null_mesher_level, */
  /**< Previous level type, defaults to
     \c Null_mesher_level. */
  class Triangulation_traits> /** Traits class that defines types for the
                                  triangulation. */
class Mesher_level<Tr, Derived, Element, Previous,
                   Triangulation_traits, Parallel_tag>
  : public Mesher_level_base<Tr, Derived, Element, Previous,
                             Triangulation_traits>
{
private:
  typedef Derived Derived_;
  template <typename ML, typename Container_element,
            typename Quality, typename Mesh_visitor> class Enqueue_element;
public:

  typedef Mesher_level<Tr,
                       Derived,
                       Element,
                       Previous,
                       Triangulation_traits,
                       Parallel_tag> Self;

  typedef Mesher_level_base<Tr,
                            Derived,
                            Element,
                            Previous,
                            Triangulation_traits> Base;


  typedef typename Base::Lock_data_structure  Lock_data_structure;
  typedef typename Base::Zone                 Zone;
  typedef typename Base::Point                Point;
  typedef typename Base::Vertex_handle        Vertex_handle;
  using Base::derived;
  using Base::is_algorithm_done;
  using Base::triangulation;
  using Base::insert;
  using Base::before_conflicts;
  using Base::before_insertion;
  using Base::after_insertion;
  using Base::previous_level;
  using Base::no_longer_element_to_refine;
  using Base::pop_next_element;
  using Base::debug_info_class_name;
  using Base::debug_info_element;

  /** \name CONSTRUCTORS */

  Mesher_level(Previous& previous)
  : Base(previous),
    FIRST_GRID_LOCK_RADIUS(
    Concurrent_mesher_config::get().first_grid_lock_radius)
    , MESH_3_REFINEMENT_GRAINSIZE(
    Concurrent_mesher_config::get().first_grid_lock_radius)
    , REFINEMENT_BATCH_SIZE(
    Concurrent_mesher_config::get().refinement_batch_size)
    , m_lock_ds(0)
    , m_worksharing_ds(0)
    , m_task_group(0)
#ifndef CGAL_NO_ATOMIC
    , m_stop_ptr(0)
#endif
  {
  }

  void add_to_TLS_lists(bool add)
  {
    derived().add_to_TLS_lists_impl(add);
  }
  void splice_local_lists()
  {
    derived().splice_local_lists_impl();
  }

  bool no_longer_local_element_to_refine()
  {
    return derived().no_longer_local_element_to_refine_impl();
  }

  Element get_next_local_element()
  {
    return derived().get_next_local_element_impl();
  }

  void pop_next_local_element()
  {
    derived().pop_next_local_element_impl();
  }

  Zone conflicts_zone(const Point& p
                      , Element e
                      , bool &facet_is_in_its_cz
                      , bool &could_lock_zone)
  {
    return derived().conflicts_zone_impl(p, e, facet_is_in_its_cz,
                                         could_lock_zone);
  }

  template <typename Mesh_visitor>
  void treat_local_refinement_queue(Mesh_visitor visitor)
  {
    // We treat the elements of the local (TLS) refinement queue
    while (no_longer_local_element_to_refine() == false)
    {
      typedef typename Derived::Container::Element Container_element;
      Container_element ce = derived().get_next_local_raw_element_impl().second;

      Mesher_level_conflict_status status;
      do
      {
        status = try_lock_and_refine_element(ce, visitor);
      }
      while (status != NO_CONFLICT
        && status != CONFLICT_AND_ELEMENT_SHOULD_BE_DROPPED
        && status != ELEMENT_WAS_A_ZOMBIE);

      pop_next_local_element();
    }
  }

  /** Tells if, as regards this level of the refinement process, if the
      point conflicts with something, and do what is needed. The return
      type is made of two booleans:
        - the first one tells if the point can be inserted,
        - in case of, the first one is \c false, the second one tells if
        the tested element should be reconsidered latter.
      This function is called by the superior level, if any.
  */
  template <class Mesh_visitor>
  Mesher_level_conflict_status test_point_conflict_from_superior(
    const Point& p, Zone& zone, Mesh_visitor &visitor)
  {
    return derived().test_point_conflict_from_superior_impl(p, zone, visitor);
  }

  /** Refines elements of this level and previous levels (PARALLEL VERSION). */
  template <class Mesh_visitor>
  void refine(Mesh_visitor visitor)
  {
    while(! is_algorithm_done() )
    {
      previous_level.refine(visitor.previous_level());
      if(! no_longer_element_to_refine() )
      {
        refine_mesh_in_parallel(visitor);
      }
    }
  }

  /** Refines elements of this level and previous levels.
  *   Stops when algorithm is done
  *   or when num vertices > approx_max_num_mesh_vertices
  */
  template <class Mesh_visitor>
  void
  refine_sequentially_up_to_N_vertices(Mesh_visitor visitor,
                                            int approx_max_num_mesh_vertices)
  {

    CGAL_assertion_msg(triangulation().get_lock_data_structure() == 0,
      "In refine_sequentially_up_to_N_vertices, the triangulation's locking data structure should be nullptr");

    while(! is_algorithm_done()
      && triangulation().number_of_vertices() < approx_max_num_mesh_vertices)
    {
      previous_level.refine(visitor.previous_level());
      if(! no_longer_element_to_refine() )
      {
        this->process_one_element(visitor);
      }
    }
  }

  void unlock_all_thread_local_elements()
  {
    if (m_lock_ds)
    {
      m_lock_ds->unlock_all_points_locked_by_this_thread();
    }
  }

  template <typename Container_element, typename Quality, typename Mesh_visitor>
  void enqueue_task(
    const Container_element &ce, const Quality &quality, Mesh_visitor visitor)
  {
    CGAL_assertion(m_task_group != 0);

    m_worksharing_ds->enqueue_work(
      Enqueue_element<Self, Container_element, Quality, Mesh_visitor>(
        *this, ce, quality, visitor),
      quality,
      *m_task_group
      // NOTE: if you uncomment this line (Load_based_worksharing_ds), the element may
      // be a zombie at this point => thus, it may be "infinite" and cause an assertion error
      // in debug mode when computing the circumcenter
      //, circumcenter_of_element(derived().extract_element_from_container_value(ce))
    );
  }

  /**
    * This function takes N elements from the queue, and try to refine
    * it in parallel.
    */
  template <class Mesh_visitor>
  void refine_mesh_in_parallel(Mesh_visitor visitor)
  {
    typedef typename Derived::Container::value_type Container_quality_and_element;

#ifdef CGAL_CONCURRENT_MESH_3_VERBOSE
    std::cerr << "Refining elements...";
#endif

    previous_level.add_to_TLS_lists(true);
    add_to_TLS_lists(true);

    m_task_group = new tbb::task_group;

    while (!no_longer_element_to_refine())
    {
      Container_quality_and_element qe = derived().get_next_raw_element_impl();
      pop_next_element();
      enqueue_task(qe.second, qe.first, visitor);
    }

    m_task_group->wait();

#if defined(CGAL_MESH_3_VERBOSE) || defined(CGAL_MESH_3_PROFILING)
    std::cerr << " Flushing";
#endif
    bool keep_flushing = true;
    while (keep_flushing)
    {
      keep_flushing = m_worksharing_ds->flush_work_buffers(*m_task_group);
      m_task_group->wait();
#if defined(CGAL_MESH_3_VERBOSE) || defined(CGAL_MESH_3_PROFILING)
      std::cerr << ".";
#endif
    }

    delete m_task_group;
    m_task_group = 0;

    splice_local_lists();
    //previous_level.splice_local_lists(); // useless
    previous_level.add_to_TLS_lists(false);
    add_to_TLS_lists(false);

#ifdef CGAL_CONCURRENT_MESH_3_VERBOSE
    std::cerr << " done." << std::endl;
#endif

  }

  template <class Mesh_visitor>
  Mesher_level_conflict_status
  try_to_refine_element(Element e, Mesh_visitor visitor)
  {
    const Tr& tr = derived().triangulation_ref_impl();
    const Point& p = tr.geom_traits().construct_weighted_point_3_object()(
      this->refinement_point(e));

#ifdef CGAL_MESH_3_VERY_VERBOSE
    std::cerr << "Trying to insert point: " << p <<
      " inside element " << debug_info_element(e) << std::endl;
#endif

    before_conflicts(e, p, visitor);

    bool could_lock_zone;
    bool facet_is_in_its_cz = true;
    Zone zone = conflicts_zone(p, e, facet_is_in_its_cz, could_lock_zone);
    Mesher_level_conflict_status result;
    if (!could_lock_zone)
      result = COULD_NOT_LOCK_ZONE;
    else if (!facet_is_in_its_cz)
      result = THE_FACET_TO_REFINE_IS_NOT_IN_ITS_CONFLICT_ZONE;
    else
      result = test_point_conflict(p, zone, visitor);

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
    case THE_FACET_TO_REFINE_IS_NOT_IN_ITS_CONFLICT_ZONE:
      std::cerr << "the facet to refine was not in the conflict zone "
        "(switching to exact)\n";
      break;
    case ELEMENT_WAS_A_ZOMBIE:
      std::cerr << "element was a zombie\n";
      break;
    case   COULD_NOT_LOCK_ELEMENT:
      std::cerr << "could not lock element\n";
      break;
    case COULD_NOT_LOCK_ZONE:
      std::cerr << "could not lock zone\n";
      break;
    }
#endif

    if(result == NO_CONFLICT)
    {
      this->before_insertion(e, p, zone, visitor);

      Vertex_handle vh = insert(p, zone);

      if (vh == Vertex_handle())
      {
        this->after_no_insertion(e, p, zone, visitor);
        result = COULD_NOT_LOCK_ZONE;
      }
      else
      {
        this->after_insertion(vh, visitor);
      }
    }
    else
    {
      this->after_no_insertion(e, p, zone, visitor);
    }

    return result;
  }

  template <typename Container_element, typename Mesh_visitor>
  Mesher_level_conflict_status
  try_lock_and_refine_element(const Container_element &ce, Mesh_visitor visitor)
  {
#ifdef CGAL_CONCURRENT_MESH_3_PROFILING
    static Profile_branch_counter_3 bcounter(
      std::string("early withdrawals / late withdrawals / successes [") + debug_info_class_name() + "]");
#endif

    Mesher_level_conflict_status result;
    Derived &derivd = derived();
    if( !derivd.is_zombie(ce) )
    {
      // Lock the element area on the grid
      Element element = derivd.extract_element_from_container_value(ce);

      // This is safe to do with the concurrent compact container because even if the element `ce`
      // gets removed from the TDS at this point, it is not actually deleted in the cells container and
      // `ce->vertex(0-3)` still points to a vertex of the vertices container whose `.point()`
      // can be safely accessed (even if that vertex has itself also been removed from the TDS).
      bool locked = derivd.try_lock_element(element, FIRST_GRID_LOCK_RADIUS);

      if( locked )
      {
        // Test it again as it may have changed in the meantime
        if( !derivd.is_zombie(ce) )
        {
          result = try_to_refine_element(element, visitor);

          //lock.release();

          if (result == CONFLICT_BUT_ELEMENT_CAN_BE_RECONSIDERED
            || result == THE_FACET_TO_REFINE_IS_NOT_IN_ITS_CONFLICT_ZONE)
          {
#ifdef CGAL_CONCURRENT_MESH_3_PROFILING
            ++bcounter; // It's not a withdrawal
#endif
          }
          else if (result == COULD_NOT_LOCK_ZONE)
          {
#ifdef CGAL_CONCURRENT_MESH_3_PROFILING
            bcounter.increment_branch_1(); // THIS is a late withdrawal!
#endif
          }
          else
          {
#ifdef CGAL_CONCURRENT_MESH_3_PROFILING
            ++bcounter;
#endif
          }
        }
        else
        {
          result = ELEMENT_WAS_A_ZOMBIE;
        }

        // Unlock
        unlock_all_thread_local_elements();
      }
      // else, we will try it again
      else
      {
#ifdef CGAL_CONCURRENT_MESH_3_PROFILING
        bcounter.increment_branch_2(); // THIS is an early withdrawal!
#endif
        // Unlock
        unlock_all_thread_local_elements();

        std::this_thread::yield();
        result = COULD_NOT_LOCK_ELEMENT;
      }
    }
    else
    {
      result = ELEMENT_WAS_A_ZOMBIE;
    }

    return result;
  }

  /** Return (can_split_the_element, drop_element). */
  template <class Mesh_visitor>
  Mesher_level_conflict_status
  test_point_conflict(const Point& p, Zone& zone, Mesh_visitor &visitor)
  {
    const Mesher_level_conflict_status result =
      previous_level.test_point_conflict_from_superior(
        p, zone, visitor.previous_level());

    if( result != NO_CONFLICT )
      return result;

    return this->private_test_point_conflict(p, zone);
  }

  /**
   * Applies one step of the algorithm: tries to refine one element of
   * previous level or one element of this level. Return \c false iff
   * <tt> is_algorithm_done()==true </tt>.
   * Note that when parallelism is activated, this is not "one step"
   * but the full refinement.
   */
  template <class Mesh_visitor>
  bool one_step(Mesh_visitor visitor)
  {
    if( ! previous_level.is_algorithm_done() )
      previous_level.one_step(visitor.previous_level());
    else if( ! no_longer_element_to_refine() )
    {
      refine_mesh_in_parallel(visitor);
    }
    return ! is_algorithm_done();
  }

  void set_lock_ds(Lock_data_structure *p)
  {
    m_lock_ds = p;
  }

  void set_worksharing_ds(WorksharingDataStructureType *p)
  {
    m_worksharing_ds = p;
  }

#ifndef CGAL_NO_ATOMIC
  void set_stop_pointer(CGAL::cpp11::atomic<bool>* stop_ptr)
  {
    m_stop_ptr = stop_ptr;
  }
#endif

  bool forced_stop() const {
#ifndef CGAL_NO_ATOMIC
    if(m_stop_ptr != 0 &&
       m_stop_ptr->load(CGAL::cpp11::memory_order_acquire) == true)
    {
      CGAL_assertion(m_task_group != 0);
      m_task_group->cancel();
      return true;
    }
#endif // not defined CGAL_NO_ATOMIC
    return false;
  }

protected:

  // Member variables
  const int FIRST_GRID_LOCK_RADIUS;
  const int MESH_3_REFINEMENT_GRAINSIZE;
  const int REFINEMENT_BATCH_SIZE;
  Lock_data_structure *m_lock_ds;
  WorksharingDataStructureType *m_worksharing_ds;

  tbb::task_group *m_task_group;
#ifndef CGAL_NO_ATOMIC
  CGAL::cpp11::atomic<bool>* m_stop_ptr;
#endif

private:

  // Functor for enqueue_task function
  template <typename ML, typename Container_element, typename Quality, typename Mesh_visitor>
  class Enqueue_element
  {
    ML                & m_mesher_level;
    Container_element   m_container_element;
    Quality             m_quality;
    Mesh_visitor        m_visitor;

  public:
    // Constructor
    Enqueue_element(ML &ml,
                    const Container_element &ce,
                    const Quality &quality,
                    Mesh_visitor visitor)
    : m_mesher_level(ml),
      m_container_element(ce),
      m_quality(quality),
      m_visitor(visitor)
    {
    }

    // operator()
    void operator()() const
    {
      typedef typename ML::Derived_::Container::value_type
                                                 Container_quality_and_element;

      Mesher_level_conflict_status status;
      do
      {
        if(m_mesher_level.forced_stop()) {
          return;
        }

        status = m_mesher_level.try_lock_and_refine_element(m_container_element,
                                                            m_visitor);
      }
      while (status != NO_CONFLICT
        && status != CONFLICT_AND_ELEMENT_SHOULD_BE_DROPPED
        && status != CONFLICT_BUT_ELEMENT_CAN_BE_RECONSIDERED
        && status != ELEMENT_WAS_A_ZOMBIE);

      // Refine the new bad facets
      m_mesher_level.before_next_element_refinement(m_visitor);

      // We can now reconsider the element if requested
      if (status == CONFLICT_BUT_ELEMENT_CAN_BE_RECONSIDERED)
        m_mesher_level.enqueue_task(m_container_element, m_quality, m_visitor);

      // Finally we add the new local bad_elements to the feeder
      while (m_mesher_level.no_longer_local_element_to_refine() == false)
      {
        Container_quality_and_element qe =
          m_mesher_level.derived().get_next_local_raw_element_impl();
        m_mesher_level.pop_next_local_element();
        m_mesher_level.enqueue_task(qe.second, qe.first, m_visitor);
      }
    }
  };

};
#endif // CGAL_LINKED_WITH_TBB

} }  // end namespace CGAL::Mesh_3

#include <CGAL/Mesher_level_visitors.h>
#include <CGAL/Mesh_3/Mesher_level_default_implementations.h>

#include <CGAL/enable_warnings.h>

#endif // CGAL_MESH_3_MESHER_LEVEL_H
