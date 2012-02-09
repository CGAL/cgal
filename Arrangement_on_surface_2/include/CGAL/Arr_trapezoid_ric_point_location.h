// Copyright (c) 2005,2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
// $URL$
// $Id$
// 
//
// Author(s)     : Idit Haran   <haranidi@post.tau.ac.il>
//                 (based on old version by Oren Nechushtan and Iddo Hanniel)

#ifndef CGAL_ARR_TRAPEZOID_RIC_POINT_LOCATION_H
#define CGAL_ARR_TRAPEZOID_RIC_POINT_LOCATION_H

/*! \file
 * Definition of the Arr_trapezoid_ric_point_location<Arrangement> template.
 */

#include <CGAL/Arr_point_location_result.h>
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>
#include <CGAL/Arr_point_location/Trapezoidal_decomposition_2.h>
#include <CGAL/Arr_point_location/Td_traits.h>
#include <CGAL/Arr_observer.h>

namespace CGAL {

/*!
 * \class
 * Mapping of an x-monotone curve to the halfedge associated with it.
 */
template <typename Arrangement_>
class PL_X_curve_plus: public Arrangement_::X_monotone_curve_2 {
public:
  typedef Arrangement_                                  Arrangement_2;
  typedef typename Arrangement_2::Geometry_traits_2     Geometry_traits_2;
  typedef typename Arrangement_2::Traits_adaptor_2      Traits_adaptor_2;
  typedef typename Arrangement_2::Halfedge_handle       Halfedge_handle;
  typedef typename Arrangement_2::X_monotone_curve_2    X_monotone_curve_2;

protected:
  //Data members
  Halfedge_handle       m_parent;    // The halfedge associated with the curve.

public:
  /*! Default constructor. */
  PL_X_curve_plus() : 
    X_monotone_curve_2(),
    m_parent() 
  {}

  /*! Constructor from a curve and a halfedge. */
  PL_X_curve_plus(const X_monotone_curve_2& cv, const Halfedge_handle& p) : 
    X_monotone_curve_2(cv), 
    m_parent(p) 
  {}

  /*! Constrtuctor from a halfedge only. */
  PL_X_curve_plus(const Halfedge_handle& p) : 
    X_monotone_curve_2(p->curve()),
    m_parent(p)
  {}

  /*! Constrtuctor from a curve only. */
  PL_X_curve_plus(const X_monotone_curve_2 &cv) : 
    X_monotone_curve_2(cv),
    m_parent() 
  {}

  /*! Get the parent halfedge. */
  Halfedge_handle get_parent() const
  { return m_parent; }
};


/*! \class
 * A class that answers point-location and queries
 * on a planar arrangement using the trapezoid_ric algorithm.
 * The Arrangement parameter corresponds to an arrangement instantiation.
 */
template <typename Arrangement_>
class Arr_trapezoid_ric_point_location : public Arr_observer <Arrangement_> {
public:
  typedef Arrangement_                                    Arrangement_2;
  typedef typename Arrangement_2::Geometry_traits_2       Geometry_traits_2;
  typedef typename Arrangement_2::Traits_adaptor_2        Traits_adaptor_2;

  typedef typename Arrangement_2::Vertex_const_handle     Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_const_handle   Halfedge_const_handle;
  typedef typename Arrangement_2::Face_const_handle       Face_const_handle;
  typedef typename Arrangement_2::Vertex_handle		  Vertex_handle;
  typedef typename Arrangement_2::Halfedge_handle	  Halfedge_handle;
  typedef typename Arrangement_2::Face_handle		  Face_handle;
  typedef typename Arrangement_2::Halfedge_iterator	  Halfedge_iterator;

  typedef typename Arrangement_2::Vertex_const_iterator   Vertex_const_iterator;
  typedef typename Arrangement_2::Edge_const_iterator     Edge_const_iterator;
  typedef typename Arrangement_2::Hole_const_iterator     Hole_const_iterator;
  typedef typename Arrangement_2::Halfedge_const_iterator
    Halfedge_const_iterator;
  typedef typename Arrangement_2::Halfedge_around_vertex_const_circulator 
    Halfedge_around_vertex_const_circulator;
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator 
    Ccb_halfedge_const_circulator;
  typedef typename Arrangement_2::Ccb_halfedge_circulator 
    Ccb_halfedge_circulator;
  typedef typename Arrangement_2::Isolated_vertex_const_iterator
    Isolated_vertex_const_iterator;

  typedef typename Geometry_traits_2::Point_2             Point_2;
  typedef typename Geometry_traits_2::X_monotone_curve_2  X_monotone_curve_2;

  typedef std::list<Halfedge_const_handle>                Edge_list;
  typedef typename Edge_list::iterator                    Std_edge_iterator;

  typedef PL_X_curve_plus<Arrangement_2>                  X_curve_plus;

  typedef CGAL::Td_traits<Traits_adaptor_2, X_curve_plus> Td_traits;
  typedef Trapezoidal_decomposition_2<Td_traits>
    Trapezoidal_decomposition;
  typedef std::vector<Halfedge_const_handle>
    Halfedge_handle_container;
  typedef typename Halfedge_handle_container::iterator
    Halfedge_handle_iterator;
 
  typedef Arr_point_location_result<Arrangement_2>        Result;
  typedef typename Result::Type                           Result_type;

  // Support boost::result_of
  typedef Result_type                                     result_type;

protected:
  typedef Trapezoidal_decomposition                       TD;

  // Data members:
  const Traits_adaptor_2*   m_traits;  // Its associated traits object.

  TD                        td;        // instance of trapezoidal decomposition
  const Td_traits*          td_traits; // instance of the TD traits
                                       //for the notification functions
  X_monotone_curve_2        m_curve_before_split; 
  X_monotone_curve_2        m_curve_before_merge1;
  X_monotone_curve_2        m_curve_before_merge2;

  template<typename T>
  Result_type result_return(T t) const { return Result()(t); }
  inline Result_type result_return() const { return Result()(); }

public:
  /*! Default constructor. */
  Arr_trapezoid_ric_point_location(bool rebuild = true) : 
    m_traits(NULL),
    td_traits(NULL)
  {
    td.set_needs_update(rebuild);
  }

  /*! Constructor given an arrangement. */
  Arr_trapezoid_ric_point_location(const Arrangement_2& arr) :
    Arr_observer<Arrangement_2>(const_cast<Arrangement_2 &>(arr))
  {
    m_traits = static_cast<const Traits_adaptor_2*>(arr.traits());
    td_traits = new Td_traits(*m_traits);
    td.init_traits(td_traits);

    build_trapezoid_ric();
  }

  /*! Destructor. */
  ~Arr_trapezoid_ric_point_location() 
  {
    if (td_traits)
      delete (td_traits);
  }
 
  /*!
   * Locate the arrangement feature containing the given point.
   * \param p The query point.
   * \return An object representing the arrangement feature containing the
   *         query point. This object is either a Face_const_handle or a
   *         Halfedge_const_handle or a Vertex_const_handle.
   */
  result_type locate(const Point_2& p) const;

  /*!
   * Locate the arrangement feature which a upward vertical ray emanating from
   * the given point hits.
   * \param p The query point.
   * \return An object representing the arrangement feature the ray hits.
   *         This object is either an empty object or a
   *         Halfedge_const_handle or a Vertex_const_handle.
   */
  result_type ray_shoot_up(const Point_2& p) const
  { return (_vertical_ray_shoot(p, true)); }

  /*!
   * Locate the arrangement feature which a downward vertical ray emanating
   * from the given point hits.
   * \param p The query point.
   * \return An object representing the arrangement feature the ray hits.
   *         This object is either an empty object or a
   *         Halfedge_const_handle or a Vertex_const_handle.
   */
  result_type ray_shoot_down(const Point_2& p) const
  { return (_vertical_ray_shoot(p, false)); }

  /// \name Notification functions, inherited and overloaded from the
  //        base observer.
  //@{

  virtual void before_assign(const Arrangement_2& arr)
  {
    clear_trapezoid_ric();
    m_traits = static_cast<const Traits_adaptor_2*>(arr.traits());
  }

  virtual void after_assign()
  { build_trapezoid_ric(); }

  virtual void before_clear()
  { clear_trapezoid_ric(); }

  virtual void after_clear()
  { build_trapezoid_ric(); }

  virtual void before_attach(const Arrangement_2& arr)
  {
    clear_trapezoid_ric();
    m_traits = static_cast<const Traits_adaptor_2*>(arr.traits());
    td_traits = new Td_traits(*m_traits);
    td.init_traits(td_traits);
  }

  virtual void after_attach()
  { build_trapezoid_ric(); }

  virtual void before_detach()
  { clear_trapezoid_ric(); }

  virtual void after_create_edge(Halfedge_handle e)
  {
    // Postcondition: h->curve() with a reference back to h
    // is inserted into TD.
    td.insert(X_curve_plus(e));
  }

  //TODO IDIT OREN: what can be done in order to avoid the need 
  //to save the original curve is to find the common endpoint of the 
  //two new halfedges, locate it in the trapezoid in order to find the 
  //curve it lies on, which is the curve that was split, and then remove 
  //this curve.
  virtual void before_split_edge(Halfedge_handle e,
                                 Vertex_handle /* v */,
                                 const X_monotone_curve_2& /* c1 */,
                                 const X_monotone_curve_2& /* c2 */)
  {
    //save this curve for the "after" function.
    m_curve_before_split = e->curve();
  }

  virtual void after_split_edge(Halfedge_handle e1,
                                Halfedge_handle e2)
  {
    td.split_edge(X_curve_plus(m_curve_before_split),
                  X_curve_plus(e1),
                  X_curve_plus(e2));
  }

  //TODO IDIT OREN: create a merged X_curve_plus withput a halfedge,
  // and in the "after" function update the halfedge.
  // think ...
  virtual void before_merge_edge(Halfedge_handle e1,
                                 Halfedge_handle e2,
                                 const X_monotone_curve_2& /* c */)
  {
    //save the curves for the "after" function.
    m_curve_before_merge1 = e1->curve();
    m_curve_before_merge2 = e2->curve();
  }

  virtual void after_merge_edge(Halfedge_handle e)
  {
    td.merge_edge(X_curve_plus(m_curve_before_merge1),
                  X_curve_plus(m_curve_before_merge2),
                  X_curve_plus(e));
  }

  virtual void before_remove_edge(Halfedge_handle e)
  {
    //called before combinatoric deletion
    td.remove(X_curve_plus(e));
  }
  //@}

public:
#ifdef CGAL_TD_DEBUG
  void debug()
  { td.debug(); }
#endif

protected:
  /*! Clear the trapezoidal decomposition. */
  inline void clear_trapezoid_ric()
  { td.clear(); }

  /*! Construct the trapezoidal decomposition. */
  void build_trapezoid_ric()
  {
    td.clear();

    Halfedge_handle_container c; 
    Edge_const_iterator eit;
    Halfedge_const_handle hh;
    Arrangement_2* arr = this->arrangement();

    for (eit = arr->edges_begin(); eit != arr->edges_end(); ++eit) {
      hh = eit;
      c.push_back(hh);
    }

    // Random shuffle of the halfedges.
    std::random_shuffle(c.begin(), c.end());

    Halfedge_handle_iterator cit;
    Halfedge_handle          he;

    for (cit = c.begin(); cit < c.end(); cit++) {
      hh = *cit;
      he = arr->non_const_handle(hh);
      td.insert(X_curve_plus(he));
    }
  }

  /*!
   * Locate the arrangement feature which a vertical ray emanating from the
   * given point hits, considering isolated vertices.
   * \param p The query point.
   * \param shoot_up Indicates whether the ray is directed upward or downward.
   * \return An object representing the arrangement feature the ray hits.
   *         This object is either a Halfedge_const_handle,
   *         a Vertex_const_handle or an empty object.
   */
  result_type _vertical_ray_shoot(const Point_2& p, bool shoot_up) const;

  /*! In vertical ray shoot, when the closest halfedge is found
   * (or unbounded face)
   * we check the isolated vertices inside the face to check whether there
   * is an isolated vertex right above/below the query point.
   */ 
  result_type
  _check_isolated_for_vertical_ray_shoot(Halfedge_const_handle halfedge_found, 
                                         const Point_2& p, bool shoot_up) const;
};

} //namespace CGAL

// The member-function definitions can be found under:
#include <CGAL/Arr_point_location/Arr_trapezoid_ric_pl_impl.h>

#endif
