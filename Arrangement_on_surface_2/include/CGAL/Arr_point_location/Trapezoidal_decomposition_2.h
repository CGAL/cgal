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
// Author(s)     : Oren Nechushtan <theoren@math.tau.ac.il>
//                 Iddo Hanniel <hanniel@math.tau.ac.il>

#ifndef CGAL_TRAPEZOIDAL_DECOMPOSITION_2_H
#define CGAL_TRAPEZOIDAL_DECOMPOSITION_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>


//#define CGAL_NO_TRAPEZOIDAL_DECOMPOSITION_2_OPTIMIZATION

#include <CGAL/Arr_tags.h>
#include <CGAL/basic.h>
#include <CGAL/Arr_point_location/Td_predicates.h>
#include <CGAL/Arr_point_location/Trapezoidal_decomposition_2_misc.h>

#include <boost/optional.hpp>
#include <boost/variant.hpp>
#include <boost/shared_ptr.hpp>

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <cmath>
#include <ctime>
#include <list>
#include <vector>
#include <map>

namespace CGAL {

/*! \class Trapezoidal_decomposition_2
 * parameters    Traits
 * Description   Implementation for a planar trapezoidal map also known as
 *   trapezoidal decomposition and vertical decomposition.
 *
 * For requirements on Traits and X_curve classes see
 *   Trapezoidal_decomposition_2 documentation.
 */
template < class Td_traits>
class Trapezoidal_decomposition_2
{
public:
  enum Locate_type {
    POINT=0,
    CURVE,
    TRAPEZOID,
    UNBOUNDED_TRAPEZOID=8
  };

  //forward declarations & friend classes declarations
  class Base_map_item_iterator;
  class In_face_iterator;
  friend class In_face_iterator;

  //type of trapezoidal decomposition traits
  typedef Td_traits Traits;

  //type of the class itself
  typedef Trapezoidal_decomposition_2<Traits> Self;

  //typedef of arrangement on surface
  typedef typename Traits::Arrangement_on_surface_2 Arrangement_on_surface_2;

  //type of traits adaptor
  typedef typename Traits::Arrangement_on_surface_2::Traits_adaptor_2
    Traits_adaptor_2;

  //type of point
  typedef typename Traits::Point Point;

  //!type of Halfedge_handle
  typedef typename Traits::Halfedge_handle
    Halfedge_handle;

  //!type of Halfedge_const_handle
  typedef typename Traits::Halfedge_const_handle
    Halfedge_const_handle;

  //!type of Vertex_const_handle
  typedef typename Traits::Vertex_const_handle
    Vertex_const_handle;

  //type of X_monotone_curve
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;

  //type of Curve end: X_monotone_curve_2 ref & Arr_curve_end
  typedef typename Traits::Curve_end    Curve_end;

  //type of td_map_item
  typedef typename Traits::Td_map_item Td_map_item;

  //type of Td_nothing
  typedef typename Traits::Td_nothing Td_nothing;

  //type of Td_active_trapezoid
  typedef typename Traits::Td_active_trapezoid Td_active_trapezoid;

  //type of Td_inactive_trapezoid
  typedef typename Traits::Td_inactive_trapezoid Td_inactive_trapezoid;

  //type of Td_active_edge
  typedef typename Traits::Td_active_edge Td_active_edge;

  //type of Td_inactive_edge
  typedef typename Traits::Td_inactive_edge Td_inactive_edge;

  //type of Td_active_vertex
  typedef typename Traits::Td_active_vertex Td_active_vertex;

  //type of Td_active_fictitious_vertex
  typedef typename Traits::Td_active_fictitious_vertex
    Td_active_fictitious_vertex;

  //type of Td_inactive_vertex
  typedef typename Traits::Td_inactive_vertex Td_inactive_vertex;

  //type of Td_inactive_fictitious_vertex
  typedef typename Traits::Td_inactive_fictitious_vertex
    Td_inactive_fictitious_vertex;

  //type of Curve end pair
  typedef typename Traits::Curve_end_pair Curve_end_pair;

  //type of Halfedge_const_handle-s' vector
  typedef std::vector<Halfedge_const_handle> Halfedge_container;

  //predicates
  //typedef CGAL::Td_active_map_item<Td_map_item> Td_active_map_item;
  typedef CGAL::Td_active_edge_item<Td_map_item, Traits> Td_active_edge_item;

  //type of search structure DAG node
  typedef Td_dag_node< Traits > Dag_node;

  //type of map of DAG nodes
  typedef std::map< int,Dag_node > Nodes_map;

  //type of trapezoids comparison function - for the map
  typedef Td_map_item_handle_less<const Td_map_item* const>
    Td_map_item_ptr_less;

  //type of trapezoids ptr map
  typedef std::map<const Td_map_item*, Td_map_item*, Td_map_item_ptr_less>
    Td_map_item_ptr_map;

public:

  /*! \class Base_map_item_iterator
   * member of Trapezoidal_decomposition_2<Traits>
   * Description Implements a basic Trapezoid iterator
   */
  class Base_map_item_iterator
  {
  public:
    //constructors
    Base_map_item_iterator() : traits(0), m_cur_item(Td_map_item(0)){ }

    Base_map_item_iterator(const Traits* traits_,
                           boost::optional<Td_map_item&> curr = boost::none)
      :traits(traits_), m_cur_item((curr) ? *curr : Td_map_item(0) ) { }

    Base_map_item_iterator(const Base_map_item_iterator &it)
          :traits(it.traits), m_cur_item(it.m_cur_item) { }

    //operator overloading
    Base_map_item_iterator  & operator=(const Base_map_item_iterator &it)
    {
      traits = it.traits;
      m_cur_item = it.m_cur_item;
      return *this;
    }

    bool operator==(const Base_map_item_iterator &it) const
    {
      return (m_cur_item == it.m_cur_item);
    }

    bool operator!=(const Base_map_item_iterator &it) const
    {
      return !operator==(it);
    }

    Td_map_item& operator*() //const
    {
      CGAL_precondition(!traits->is_empty_item(m_cur_item));
      return m_cur_item;
    }

    bool operator!() const
    {
      return traits->is_empty_item(m_cur_item);//!m_cur_item;
    }

  protected:
    const Traits *traits; //pointer to the traits
    Td_map_item m_cur_item; //the current map item (or none)
  };

/*! \class In_face_iterator
   * member of Trapezoidal_decomposition_2<Traits>
   * Derived from Base_map_item_iterator class
   * Description Implements a Trapezoid iterator along a Halfedge
   */
  class In_face_iterator : public Base_map_item_iterator
  {

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
    using Base_map_item_iterator::m_cur_item;
    using Base_map_item_iterator::traits;
#endif

  protected:
    //reference to the seperating X_monotone_curve_2
    const X_monotone_curve_2& m_sep;

  public:
    //constructors
    In_face_iterator(const Traits* traits_, Halfedge_const_handle sep,
                     boost::optional<Td_map_item&> curr = boost::none)
            :Base_map_item_iterator(traits_,curr), m_sep(sep->curve())
    { }

    In_face_iterator(const Traits* traits_, const X_monotone_curve_2& sep,
                     boost::optional<Td_map_item&> curr = boost::none)
            :Base_map_item_iterator(traits_,curr), m_sep(sep)
    { }


    In_face_iterator(const In_face_iterator &it)
            :Base_map_item_iterator((Base_map_item_iterator&)it),
              m_sep(it.m_sep)
    { }

    //operatoror overloading
    bool operator==(const In_face_iterator &it) const
    {
      return ( Base_map_item_iterator::operator==(it) &&
               traits->equal_2_object()(m_sep,it.m_sep));
    }

    /*
      destription:
      advances m_cur_item to one of the right neighbours according to the relation
      between the seperating Halfedge (m_sep) and the right() trapezoid point.
      precoditions:
      m_sep doesn't intersect any existing edges except possibly on common end
      points.
      postconditions:
      if the rightmost trapezoid was traversed m_cur_item is set to NULL.
      remark:
      if the seperator is vertical, using the precondition assumptions it
      follows that there is exactly one trapezoid to travel.
    */
    In_face_iterator& operator++()
    {
      if (traits->is_empty_item(m_cur_item))
        return *this;// end reached, do nothing!

#ifndef CGAL_TD_DEBUG
      CGAL_warning(traits != NULL);
#else
      CGAL_assertion(traits != NULL);
      CGAL_assertion(traits->is_active(m_cur_item));
      //m_cur_item should be a trapezoid or an edge
      CGAL_assertion(!traits->is_td_vertex(m_cur_item));
#endif

      if (traits->is_td_trapezoid(m_cur_item))
      {
        //if the map item is a trapezoid
        Td_active_trapezoid tr (boost::get<Td_active_trapezoid>(m_cur_item));

#ifndef CGAL_TD_DEBUG
        CGAL_warning_code(Dag_node* tt = tr.dag_node();)
        CGAL_warning(!tt->is_inner_node());
#else
        CGAL_assertion_code(Dag_node* tt = tr.dag_node();)
        CGAL_assertion(tt);
        CGAL_assertion(!tt->is_inner_node());
#endif

        // handle degeneracies
        typename Traits::Compare_curve_end_xy_2 compare_xy =
          traits->compare_curve_end_xy_2_object();
        if (compare_xy (traits->vtx_to_ce(tr.left()),
                        Curve_end(m_sep,ARR_MAX_END)) != SMALLER)
        {
          //if the trapezoid's left end point is equal to or larger from the
          //  max end of sep, we reached the end of the iterator
          m_cur_item = Td_map_item(0);
        }
        else
        {
          //if the trapezoid's left end point is smaller from the sep's max end

          //comparing the y value of the trapezoid's right end point and sep
          //   (at the trapezoid's right x value), in order to select the
          //    next trapezoid in the iterator
          typename Traits::Compare_curve_end_y_at_x_2 compare_y_at_x =
            traits->compare_curve_end_y_at_x_2_object();
          switch (compare_y_at_x (traits->vtx_to_ce(tr.right()), m_sep))
          {
           case SMALLER:
              m_cur_item = tr.rt();
            break;
           case LARGER:
              m_cur_item = tr.rb();
            break;
           case EQUAL:
            // end reached
             m_cur_item = Td_map_item(0);
            break;
           default:
             m_cur_item = Td_map_item(0);
            break;
          }
        }
      }
      else
      {
        //if the map item is an edge

        Td_active_edge e (boost::get<Td_active_edge>(m_cur_item));
        CGAL_assertion_code(Dag_node* tt = e.dag_node();)
        CGAL_assertion(tt != NULL);
        CGAL_assertion(tt->is_inner_node());

        //go to next() of the current edge.
        // as long as there is an edge fragment of the same
        //  edge - next() exists.
        // If next() does not exist we reached the last fragment of the edge
        m_cur_item = e.next();
        if (!traits->is_empty_item(m_cur_item))
        {
          //if next() exists, find the next real edge fragment trapezoid
          //    (skip points)
          while(traits->is_td_vertex(m_cur_item))
          {
            Dag_node* node =
              boost::apply_visitor(dag_node_visitor(),m_cur_item);
            m_cur_item = node->left_child().get_data();
          }

          //make sure we stopped in an edge
          CGAL_warning(traits->is_td_edge(m_cur_item));
        }
      }
      return *this;
    }

    In_face_iterator operator++(int)
    {
      In_face_iterator tmp = *this;
      ++*this;
      return tmp;
    }

    const X_monotone_curve_2& seperator()
    {
      return m_sep;
    }

    Td_active_trapezoid& trp()
    {
      CGAL_precondition (!traits->is_empty_item(m_cur_item));
      CGAL_precondition (traits->is_active(m_cur_item) &&
                         traits->is_td_trapezoid(m_cur_item));
      return boost::get<Td_active_trapezoid>(m_cur_item);
    }

    Td_active_edge& e()
    {
      CGAL_precondition (!traits->is_empty_item(m_cur_item));
      CGAL_precondition (traits->is_active(m_cur_item) &&
                         traits->is_td_edge(m_cur_item));
      return boost::get<Td_active_edge>(m_cur_item);
    }

  };

  /*!  Visitors for accessing td map items methods */
  class rb_visitor : public boost::static_visitor<Td_map_item>
  {
  public:
    Td_map_item operator()(Td_active_trapezoid& t) const
    {
      return t.rb();
    }

    template < typename T >
    Td_map_item operator()(T& /* t */) const
    {
      CGAL_assertion(false);
      return Td_map_item(0);
    }
  };

  class set_rb_visitor : public boost::static_visitor<void>
  {
  public:
    set_rb_visitor (const Td_map_item& rb) : m_rb(rb) {}


    void operator()(Td_active_trapezoid& t) const
    {
      t.set_rb(m_rb);
    }

    template < typename T >
    void operator()(T& /* t */) const
    {
      CGAL_assertion(false);
    }

  private:
    const Td_map_item& m_rb;
  };

  class rt_visitor : public boost::static_visitor<Td_map_item>
  {
  public:
    Td_map_item operator()(Td_active_trapezoid& t) const
    {
      return t.rt();
    }

    template < typename T >
    Td_map_item operator()(T& /* t */) const
    {
      CGAL_assertion(false);
      return Td_map_item(0);
    }
  };

  class set_rt_visitor : public boost::static_visitor<void>
  {
  public:
    set_rt_visitor (const Td_map_item& rt) : m_rt(rt) {}

    void operator()(Td_active_trapezoid& t) const
    {
      t.set_rt(m_rt);
    }

    template < typename T >
    void operator()(T& /*t*/) const
    {
      CGAL_assertion(false);
    }

  private:
    const Td_map_item& m_rt;
  };

  class lb_visitor : public boost::static_visitor<Td_map_item>
  {
  public:
    Td_map_item operator()(Td_active_trapezoid& t) const
    {
      return t.lb();
    }

    template < typename T >
    Td_map_item operator()(T& /* t */) const
    {
      CGAL_assertion(false);
      return Td_map_item(0);
    }
  };

  class set_lb_visitor : public boost::static_visitor<void>
  {
  public:
    set_lb_visitor (const Td_map_item& lb) : m_lb(lb) {}

    void operator()(Td_active_trapezoid& t) const
    {
      return t.set_lb(m_lb);
    }

    template < typename T >
    void operator()(T& /* t */) const
    {
      CGAL_assertion(false);
    }

  private:
    const Td_map_item& m_lb;
  };

  class set_lt_visitor : public boost::static_visitor<void>
  {
  public:
    set_lt_visitor (const Td_map_item& lt) : m_lt(lt) {}

    void operator()(Td_active_trapezoid& t) const
    {
      t.set_lt(m_lt);
    }

    template < typename T >
    void operator()(T& /*t*/) const
    {
      CGAL_assertion(false);
    }

  private:
    const Td_map_item& m_lt;
  };

  class bottom_he_visitor : public boost::static_visitor<Halfedge_const_handle>
  {
  public:
    Halfedge_const_handle operator()(Td_active_trapezoid& t) const
    {
      return t.bottom();
    }

    template < typename T >
    Halfedge_const_handle operator()(T& /*t*/) const
    {
      CGAL_assertion(false);
      return Halfedge_const_handle();
    }
  };

  class set_bottom_he_visitor : public boost::static_visitor< void  >
  {
  public:
    set_bottom_he_visitor (Halfedge_const_handle he) : m_bottom_he(he) {}

    void operator()(Td_active_trapezoid& t) const
    {
      t.set_bottom(m_bottom_he);
    }

    template < typename T >
    void operator()(T& /*t*/) const
    {
      CGAL_assertion(false);
    }
  private:
    Halfedge_const_handle m_bottom_he;
  };

  class top_he_visitor : public boost::static_visitor<Halfedge_const_handle>
  {
  public:
     Halfedge_const_handle operator()(Td_active_trapezoid& t) const
    {
      return t.top();
    }

    template < typename T >
    Halfedge_const_handle operator()(T& /* t */) const
    {
      CGAL_assertion(false);
      return Halfedge_const_handle();
    }
  };

  class set_top_he_visitor : public boost::static_visitor<void>
  {
  public:
    set_top_he_visitor (Halfedge_const_handle he) : m_top_he(he) {}

    void operator()(Td_active_trapezoid& t) const
    {
      t.set_top(m_top_he);
    }

    template < typename T >
    void operator()(T& /* t */) const
    {
      CGAL_assertion(false);
    }
  private:
    Halfedge_const_handle m_top_he;
  };

  class cw_he_visitor : public boost::static_visitor< Halfedge_const_handle  >
  {
  public:
    Halfedge_const_handle operator()(Td_active_vertex& t) const
    {
      return t.cw_he();
    }
    Halfedge_const_handle operator()(Td_active_fictitious_vertex& t) const
    {
      return t.cw_he();
    }

    template < typename T >
    Halfedge_const_handle operator()(T& /*t*/) const
    {
      CGAL_assertion(false);
      return Halfedge_const_handle();
    }
  };

  class set_cw_he_visitor : public boost::static_visitor<void>
  {
  public:
    set_cw_he_visitor (Halfedge_const_handle he) : m_cw_he(he) {}

    void operator()(Td_active_vertex& t) const
    {
      t.set_cw_he(m_cw_he);
    }
    void operator()(Td_active_fictitious_vertex& t) const
    {
      t.set_cw_he(m_cw_he);
    }

    template < typename T >
    void operator()(T& /*t*/) const
    {
      CGAL_assertion(false);
    }
  private:
    Halfedge_const_handle m_cw_he;
  };

  class dag_node_visitor : public boost::static_visitor<Dag_node*>
  {
  public:
    Dag_node* operator()(Td_nothing& /* t */) const
    {
      CGAL_assertion(false);
      return NULL;
    }
    Dag_node* operator()(Td_inactive_trapezoid& /* t */) const
    {
      CGAL_assertion(false);
      return NULL;
    }

    template < typename T >
    Dag_node* operator()(T& t) const
    {
      return t.dag_node();
    }
  };

  class set_dag_node_visitor : public boost::static_visitor<void>
  {
  public:
    set_dag_node_visitor(Dag_node* node):m_node(node) {}

    void operator()(Td_nothing& /*t*/) const
    {
      CGAL_assertion(false);
    }
    void operator()(Td_inactive_trapezoid& /*t*/) const
    {
      CGAL_assertion(false);
    }

    template < typename T >
    void operator()(T& t) const
    {
      t.set_dag_node(m_node);
    }

  private:
    Dag_node* m_node;
  };

  class curve_end_for_fict_vertex_visitor :
    public boost::static_visitor<boost::optional<Curve_end> >
  {
  public:
    boost::optional<Curve_end> operator()(Td_active_fictitious_vertex& t) const
    {
      return t.curve_end();
    }

    boost::optional<Curve_end>
    operator()(Td_inactive_fictitious_vertex& t) const
    {
      return t.curve_end();
    }

    template < typename T >
    boost::optional<Curve_end> operator()(T& /* t */) const
    {
      CGAL_assertion(false);
      return boost::none;
    }
  };

  class point_for_vertex_visitor : public boost::static_visitor< Point  >
  {
  public:
    Point operator()(Td_active_vertex& t) const
    {
      return t.point();
    }

    Point operator()(Td_inactive_vertex& t) const
    {
      return t.point();
    }

    template < typename T >
    Point operator()(T& /* t */) const
    {
      CGAL_assertion(false);
      return Point();
    }
  };

  class curve_end_for_active_vertex_visitor :
    public boost::static_visitor<boost::optional<Curve_end> >
  {
    public:
    boost::optional<Curve_end> operator()(Td_active_vertex& t) const
    {
      return t.curve_end();
    }

    boost::optional<Curve_end> operator()(Td_active_fictitious_vertex& t) const
    {
      return t.curve_end();
    }

    template < typename T >
    boost::optional<Curve_end> operator()(T& /* t */) const
    {
      CGAL_assertion(false);
      return boost::none;
    }
  };

  class vertex_for_active_vertex_visitor :
    public boost::static_visitor<Vertex_const_handle>
  {
    public:
    Vertex_const_handle operator()(Td_active_vertex& t) const
    {
      return t.vertex();
    }

    Vertex_const_handle operator()(Td_active_fictitious_vertex& t) const
    {
      return t.vertex();
    }

    template < typename T >
    Vertex_const_handle operator()(T& /*t*/) const
    {
      CGAL_assertion(false);
      return Vertex_const_handle();
    }
  };

  class cv_for_edge_visitor :
    public boost::static_visitor<boost::optional<const X_monotone_curve_2&> >
  {
  public:
    boost::optional<const X_monotone_curve_2&>
    operator()(Td_active_edge& t) const
    {
      return t.halfedge()->curve();
    }

    boost::optional<const X_monotone_curve_2&>
    operator()(Td_inactive_edge& t) const
    {
      return t.curve();
    }

    template <typename T>
    boost::optional<const X_monotone_curve_2&> operator()(T& /* t */) const
    {
      CGAL_assertion(false);
      return boost::none;
    }
  };

  ////MICHAL: currently not in use since split is implemented as removed and insert two
  //struct Before_split_data
  //{
  //  X_monotone_curve_2 m_cv_before_split;
  //  Td_map_item* m_p_old_t;
  //  Td_map_item* m_p_t1;
  //  Td_map_item* m_p_t2;
  //  In_face_iterator* m_p_btm_it;
  //  In_face_iterator* m_p_mid_it;
  //  In_face_iterator* m_p_top_it;
  //
  //};

  //////////////////////////////////////////////
  //Trapezoidal_decomposition_2 member functions:
  //////////////////////////////////////////////


protected:

  /*!  is_edge_to_right variants:
      returning true if the given edge is on the right side
      of the given point / curve-end */

  bool is_edge_to_right(Halfedge_const_handle he, const Point& p) const
  {
    typename Traits::Equal_curve_end_2 equal =
      traits->equal_curve_end_2_object();
    //p is either min or max end of he
    CGAL_precondition(equal(Curve_end(he,ARR_MIN_END), p) ||
                      equal(Curve_end(he,ARR_MAX_END), p));

    return equal(Curve_end(he,ARR_MIN_END), p);
  }

  bool is_edge_to_right(Halfedge_const_handle he, const Curve_end& ce) const
  {
    typename Traits::Equal_curve_end_2 equal =
      traits->equal_curve_end_2_object();
    //p is either min or max end of he
    CGAL_precondition(equal(Curve_end(he,ARR_MIN_END), ce) ||
                      equal(Curve_end(he,ARR_MAX_END), ce));

    //if the curve end ce is on the right boundary - return false;
    if (traits->parameter_space_in_x_2_object()
          (ce.cv(), ce.ce()) == ARR_RIGHT_BOUNDARY)
    {
      return false;
    }

    return equal(Curve_end(he,ARR_MIN_END), ce);
  }

  //returns true if the given curve is on the right side of the given point
  bool is_curve_to_right(const X_monotone_curve_2& cv, const Point& p) const
  {
    typename Traits::Equal_curve_end_2 equal =
      traits->equal_curve_end_2_object();
    //p is either min or max end of he
    CGAL_precondition(equal(Curve_end(cv,ARR_MIN_END), p) ||
                      equal(Curve_end(cv,ARR_MAX_END), p));

    return equal(Curve_end(cv,ARR_MIN_END), p);
  }

  /*!  is_end_point_left_low variants:
      returning true if the first curve-end is left-low
      of the second curve-end */
  bool is_end_point_left_low(const Point& p1, const Point& p2) const
  {
    return (traits->compare_xy_2_object()(p1, p2) == SMALLER);
  }

  bool is_end_point_left_low(const Point& p, const Curve_end& ce) const
  {
    return (traits->compare_curve_end_xy_2_object()(p, ce) == SMALLER);
  }

  bool is_end_point_left_low(const Curve_end& ce, const Point& p) const
  {
    return (traits->compare_curve_end_xy_2_object()(p, ce) == LARGER);
  }

  bool is_end_point_left_low(const Curve_end& ce1, const Curve_end& ce2) const
  {
    return (traits->compare_curve_end_xy_2_object()(ce1, ce2) == SMALLER);
  }

  template <typename T>
  bool is_end_point_left_low(const T& t, const Dag_node& node) const
  {
    typename Traits::Compare_curve_end_xy_2 compare =
      traits->compare_curve_end_xy_2_object();
    Td_map_item vtx_item (node.get_data());
    bool is_fict_vtx = traits->is_fictitious_vertex(vtx_item);
    if (is_fict_vtx) {
      return (compare(t,
                      *(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),
                                             vtx_item))) == SMALLER);
    }
    else {
      return (compare(t, boost::apply_visitor(point_for_vertex_visitor(),
                                              vtx_item)) == SMALLER);
    }
  }

  /*!  is_end_point_right_top variants:
      returning true if the first curve-end is right-top
      of the second curve-end */
  bool is_end_point_right_top(const Point& p1, const Point& p2) const
  {
    return (traits->compare_xy_2_object()(p1, p2) == LARGER);
  }

  bool is_end_point_right_top(const Point& p, const Curve_end& ce) const
  {
    return (traits->compare_curve_end_xy_2_object()(p, ce) == LARGER);
  }

  bool is_end_point_right_top(const Curve_end& ce, const Point& p) const
  {
    return (traits->compare_curve_end_xy_2_object()(p, ce) == SMALLER);
  }

  bool is_end_point_right_top(const Curve_end& ce1, const Curve_end& ce2) const
  {
    return (traits->compare_curve_end_xy_2_object()(ce1, ce2) == LARGER);
  }

  template <typename T>
  bool is_end_point_right_top(const T& t, const Dag_node& node) const
  {
    typename Traits::Compare_curve_end_xy_2 compare =
      traits->compare_curve_end_xy_2_object();
    Td_map_item vtx_item (node.get_data());
    bool is_fict_vtx = traits->is_fictitious_vertex(vtx_item);
    if (is_fict_vtx) {
      return (compare(t,
                 *(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),
                                        vtx_item)))  == LARGER);
    }
    else {
      return (compare(t,
                  boost::apply_visitor(point_for_vertex_visitor(),
                                       vtx_item)) == LARGER);
    }
  }

  template <typename T>
  bool are_equal_end_points(const T& t, const Dag_node& node) const
  {
    typename Traits::Equal_curve_end_2 equal =
      traits->equal_curve_end_2_object();
    Td_map_item vtx_item (node.get_data());
    bool is_fict_vtx = traits->is_fictitious_vertex(vtx_item);
    if (is_fict_vtx) {
      return equal(t,
                   *(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),
                                          vtx_item)));
    }
    else {
      return equal(t,
                   boost::apply_visitor(point_for_vertex_visitor(),
                                        vtx_item));
    }
  }

  /*!
   * finds the node of the leftmost trapezoid with respect to a curve.
   * \param left_cv_end_node The dag node representing the left endpoint of
   *        the cv
   * \param cv The curve
   * \param cres SMALLER/EQUAL/LARGER (searching for the leftmost trapezoid
   *        which is below/on/above cv)
   * \return The required DAG node
   */
  Dag_node find_leftmost_dag_node_of_curve(const Dag_node& left_cv_end_node,
                                           const X_monotone_curve_2& cv,
                                           Comparison_result cres) const
  {
    CGAL_assertion(traits != NULL);
    Td_map_item& item = left_cv_end_node.get_data();
    CGAL_precondition(traits->is_td_vertex(item));
    CGAL_precondition (are_equal_end_points(Curve_end(cv,ARR_MIN_END),
                                            left_cv_end_node));

    //if ( traits->is_fictitious_vertex(item) )
    //{
    //  CGAL_precondition(traits->equal_curve_end_2_object()
    //    (Curve_end(cv,ARR_MIN_END), *(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),item))));
    //}
    //else
    //{
    //  CGAL_precondition(traits->equal_curve_end_2_object()
    //     (Curve_end(cv,ARR_MIN_END), boost::apply_visitor(point_for_vertex_visitor(), item)));
    //}
    //find the node of the curve's leftmost trapezoid
    Dag_node cv_leftmost_node(left_cv_end_node.right_child());
    if (traits->is_fictitious_vertex(item) )
    {
      Curve_end ce( *(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),
                                           item)));
      search_using_dag_with_cv(cv_leftmost_node, traits, ce, &cv, cres);
    }
    else
    {
      Point p( boost::apply_visitor(point_for_vertex_visitor(), item));
      search_using_dag_with_cv(cv_leftmost_node, traits, p, &cv, cres);
    }
    return cv_leftmost_node;
  }

   /*!
   * follow_curve variants:
   * follows trapezoids along a curve (below/on/above it)
   * \param left_cv_end_node The dag node representing the left endpoint of
   *        the cv
   * \param he/cv The halfedge / The curve
   * \param cres SMALLER/EQUAL/LARGER (indicating the position with respect to the curve)
   * \return An iterator for td map items along a curve
   */
  In_face_iterator follow_curve(const Dag_node& left_cv_end_node,
                                Halfedge_const_handle he,
                                Comparison_result up) const
  {
    return follow_curve(left_cv_end_node, he->curve(), up);
  }

  In_face_iterator follow_curve(const Dag_node& left_cv_end_node,
                                const X_monotone_curve_2& cv,
                                Comparison_result up) const
  {
    Dag_node cv_leftmost_node(find_leftmost_dag_node_of_curve(left_cv_end_node,cv,up));
    //return a trapezoid iterator that starts from this trapezoid
    //  and continues according to the curve cv
    return In_face_iterator(traits, cv, cv_leftmost_node.get_data());
  }

  //-----------------------------------------------------------------------------
  // Description:
  //  Input: pointer to left trapezoid, pointer to right trapezoid
  //  Output: true iff the merging took place
  //  If the two input trapezoids can be merged they are ,
  //  with one copy destroyed(the right one).
  // Preconditions:
  //  the right trapezoid is to the right of the left one
  bool merge_if_possible(Td_map_item& left_item, Td_map_item& right_item)
  {
    CGAL_precondition(traits->is_empty_item(left_item) ||
                      (traits->is_active(left_item) &&
                       traits->is_td_trapezoid(left_item)));
    CGAL_precondition(traits->is_empty_item(right_item) ||
                      (traits->is_active(right_item) &&
                       traits->is_td_trapezoid(right_item)));

    if (traits->is_empty_item(left_item) || traits->is_empty_item(right_item))
      return false;

    Td_active_trapezoid& left  (boost::get<Td_active_trapezoid>(left_item));
    Td_active_trapezoid& right (boost::get<Td_active_trapezoid>(right_item));

    if (traits->is_trapezoids_top_equal(left,right) &&
        traits->is_trapezoids_bottom_equal(left,right) &&
        traits->equal_curve_end_2_object()
         (traits->vtx_to_ce(left.right()), traits->vtx_to_ce(right.left())))
    {
      left.merge_trapezoid(right);
      //set the depth to be the max of the two merged nodes
      left.dag_node()->depth() = (std::max)(left.dag_node()->depth(),
                                            right.dag_node()->depth());
      CGAL_postcondition(
        left.is_on_right_boundary() == right.is_on_right_boundary());

      return true;
    }
    return false;
  }

  //---------------------------------------------------------------------------
  // Description:
  //  splits the trapezoid with vertical line through v
  //  assuming that he (the first cw halfedge starting at 12 o'clock) is in the
  //  desired direction, such that v is their source
  // Precondition:
  //  The trapezoid is active and contains v in its closure
  //
  Dag_node& split_trapezoid_by_vertex(Dag_node& tt,
                                      Vertex_const_handle v,
                                      Halfedge_const_handle he);

  Td_map_item build_vertex_map_item(Vertex_const_handle v,
                                    Halfedge_const_handle he,
                                    Dag_node* node);
  //---------------------------------------------------------------------------
  // Description:
  //  the opposite operation for spliting the trapezoid with
  //  vertical line through ce
  // Precondition:
  //  The root trapezoid is degenerate point (ce) and is active
  void undo_split_trapezoid_by_vertex(Dag_node& tr_node, const Curve_end& ce);

  void deactivate_trapezoid (Dag_node& trpz_node, Dag_node* active_node) const;

  void deactivate_vertex (Dag_node& vtx_node) const;

  void deactivate_edge (boost::shared_ptr<X_monotone_curve_2>& cv, Dag_node& edge_node) const;

  //-----------------------------------------------------------------------------
  // Description:
  //  splits the trapezoid that corresponds to the root of the
  //  trapezoidal tree with an input halfedge he
  // Precondition:
  //  The root trapezoid is active
  //  The root trapezoid is devided by he or is equal to it and is vertical.
  Dag_node& split_trapezoid_by_halfedge(Dag_node& split_node,
                                        Td_map_item& prev_e,
                                        Td_map_item& prev_bottom_tr,
                                        Td_map_item& prev_top_tr,
                                        Halfedge_const_handle he);


  //---------------------------------------------------------------------------
  // Description:
  //  update
  //   tr.bottom()
  //   vertical_ray_shoot downward from tr
  //   tr.top()
  //    vertical_ray_shoot upward from tr
  //  update all the curves incident to the vertex that there's a new curve
  //  starting from this vertex
  //  this point must be an interior point and not a point on the boundaries,
  //  since a point on the boundaries is related to one curve only
  Td_map_item&
  update_vtx_with_new_edge(Halfedge_const_handle he,
                           const Curve_end& ce,
                           Td_map_item& vtx_item,
                           const Locate_type& CGAL_precondition_code(lt));

  Td_map_item& insert_curve_at_vtx_using_dag(Halfedge_const_handle he,
                                             Vertex_const_handle v,
                                             Td_map_item& tr,
                                             const Locate_type&
                                             CGAL_precondition_code(lt));


  //void set_trp_params_after_halfedge_update(Halfedge_const_handle old_he,
  //                                          Halfedge_const_handle new_he,
  //                                          Td_map_item& vtx_item); //MICHAL: not in use


  void update_vtx_cw_he_after_merge(const X_monotone_curve_2& old_cv,
                                    Halfedge_const_handle new_he,
                                    Td_map_item& vtx_item);

  ////MICHAL: currently not in use since split is implemented as: remove and insert two
  //void set_trp_params_after_split_halfedge_update(Halfedge_const_handle new_he,
  //                                                Td_map_item& vtx_item,
  //                                                Halfedge_const_handle he1,
  //                                                Halfedge_const_handle he2);

  //-----------------------------------------------------------------------------
  // Description:
  //  update map items traveled along an iterator till end reached
  //   with the new halfedge
  // precondition:
  //  end==0 or end is on the path of the iterator
  // postcondition:
  //  end is pointer to the last trapezoid encountered,if any
  void update_map_items_after_merge(In_face_iterator& it,
                                    Halfedge_const_handle old_he,
                                    Halfedge_const_handle new_he,
                                    Vertex_const_handle min_v,
                                    Vertex_const_handle max_v,
                                    Td_map_item& end);


  //---------------------------------------------------------------------------
  // Description:
  //  advances input Data structure using data structure,input point p and
  //  possibly Halfedge p_he till
  //  p is found(if p_he hadn't been given)
  //  p_he is found(if p_he was given)
  //  or
  //  leaf node reached
  // postcondition:
  //  output is the closest active trapezoid to ce/p_he
  // remark:
  //  use this function with care!
  Locate_type search_using_dag(Dag_node& curr_node,
                               const Traits* traits,
                               const Point& p,
                               Halfedge_const_handle he,
                               Comparison_result up = EQUAL) const;

  ////-------------------------------------------------------------------------
  //// Description:
  ////  advances input Data structure using data structure,input point p and
  ////  possibly Halfedge p_he till
  ////  p is found(if p_he hadn't been given)
  ////  p_he is found(if p_he was given)
  ////  or
  ////  leaf node reached
  //// postcondition:
  ////  output is the closest active trapezoid to ce/p_he
  //// remark:
  ////  use this function with care!
  //void search_and_print_using_dag (std::ostream& out,
  //                                Dag_node& curr_node,
  //                                const Traits* traits,
  //                                const Point& p,
  //                                Halfedge_const_handle he,
  //                                Comparison_result up = EQUAL) const;

  //---------------------------------------------------------------------------
  // Description:
  //  advances input Data structure using data structure,input point ce and
  //  possibly Halfedge p_he till
  //  ce is found(if p_he hadn't been given)
  //  p_he is found(if p_he was given)
  //  or
  //  leaf node reached
  // postcondition:
  //  output is the closest active trapezoid to ce/p_he
  // remark:
  //  use this function with care!
  Locate_type search_using_dag (Dag_node& curr_node,
                                const Traits* traits,
                                const Curve_end& ce,
                                Halfedge_const_handle he,
                                Comparison_result up = EQUAL) const;

  //---------------------------------------------------------------------------
  // Description:
  //  advances input Data structure using data structure,input point ce and
  //  possibly X_monotone_curve_2 p_cv till
  //  ce is found(if p_cv hadn't been given)
  //  p_cv is found(if p_cv was given)
  //  or
  //  leaf node reached
  // postcondition:
  //  output is the closest active trapezoid to ce/p_cv
  // remark:
  //  use this function with care!
  Locate_type search_using_dag_with_cv(Dag_node& curr_node,
                                       const Traits* traits,
                                       const Curve_end& ce,
                                       const X_monotone_curve_2* p_cv,
                                       Comparison_result up = EQUAL) const;

  //---------------------------------------------------------------------------
  // Description:
  //  advances input Data structure using data structure,input point ce and
  //  possibly X_monotone_curve_2 p_cv till
  //  p is found(if p_cv hadn't been given)
  //  p_cv is found(if p_cv was given)
  //  or
  //  leaf node reached
  // postcondition:
  //  output is the closest active trapezoid to ce/p_cv
  // remark:
  //  use this function with care!
  Locate_type search_using_dag_with_cv (Dag_node& curr_node,
                                        const Traits* traits,
                                        const Point& p,
                                        const X_monotone_curve_2* p_cv,
                                        Comparison_result up = EQUAL) const;


  Dag_node container2dag(Nodes_map& ar, int left, int right,
                         int& num_of_new_nodes) const;

  bool is_last_edge(Halfedge_const_handle he, Td_map_item& vtx_item);

  /*==============================================
    Trapezoidal_decomposition_2 public member functions
    ==============================================*/
public:

  Trapezoidal_decomposition_2(bool with_guarantees = true) :
    m_largest_leaf_depth(0),
    m_number_of_dag_nodes(1),
    m_number_of_curves(0),
    traits(0),
    m_arr(0),
    m_depth_threshold(CGAL_TD_DEFAULT_DEPTH_THRESHOLD),
    m_size_threshold(CGAL_TD_DEFAULT_SIZE_THRESHOLD)
  {
    init();
    set_with_guarantees(with_guarantees);
  }

  Trapezoidal_decomposition_2(const double& depth_th, const double& size_th,
                              bool with_guarantees = true) :
    m_largest_leaf_depth(0),
    m_number_of_curves(0),
    m_number_of_dag_nodes(1),
    traits(0),
    m_arr(0),
    m_depth_threshold(depth_th),
    m_size_threshold(size_th)
  {
    init();
    set_with_guarantees(with_guarantees);
  }

  //MICHAL: problematic, should not be used
  //Trapezoidal_decomposition_2(const Self& td)
  //  : m_with_guarantees(td.m_with_guarantees),
  //    m_number_of_curves(td.m_number_of_curves),
  //    m_largest_leaf_depth(td.m_largest_leaf_depth),
  //    m_number_of_dag_nodes(td.m_number_of_dag_nodes),
  //    traits(td.traits),
  //    m_arr(td.m_arr),
  //    last_cv(Td_map_item(0)), prev_cv(Td_map_item(0)),
  //    m_depth_threshold(td.m_depth_threshold),
  //    m_size_threshold(td.m_size_threshold)
  //{
  //  Td_map_item_ptr_map htr;
  //  /*! \todo allocate hash_map size according to content.
  //   * \todo change vector<> to in_place_list and pointer hash to trapezoidal
  //   * hash..
  //   */
  //  //X_trapezoid_vector vtr;
  //  std::vector<Td_map_item> vtr;
  //  Td_active_map_item pr;
  //  int sz = Td_map_item_filter(vtr, &td.dag_root());
  //  //! \todo Reduce the 3 iterations to 1 (or 2) iterator.
  //  // First iteration: filter out the active trapezoids.
  //  typename std::vector<Td_map_item>::const_iterator it;
  //  for (it = vtr.begin(); it != vtr.end(); ++it)
  //  {
  //    Dag_node* ds_copy = new Dag_node(*it);
  //    const Td_map_item* cur = &*it;
  //    Td_map_item* tr_copy = &*(*ds_copy);
  //    tr_copy->set_dag_node(ds_copy);
  //    CGAL_assertion(&*(*tr_copy->dag_node()) == tr_copy);
  //    ds_copy->depth() = cur->dag_node()->depth();
  //    // We cheat a little with the depth.
  //    htr.insert(typename Td_map_item_ptr_map::value_type(cur, tr_copy));
  //    // Second iteration: generate new copies of trapezoids and nodes.
  //  }
  //
  //  for (it = vtr.begin(); it!=vtr.end(); ++it)
  //  {
  //    const Td_map_item* cur = &*it;
  //    Td_map_item* tr_copy = htr.find(cur)->second;
  //    const Dag_node* child;
  //    CGAL_assertion(tr_copy);
  //    tr_copy->set_rt(cur->rt() ?
  //                    htr.find(cur->rt())->second : NULL);
  //    tr_copy->set_rb(cur->rb() ?
  //                    htr.find(cur->rb())->second : NULL);
  //    tr_copy->set_lt(cur->lt() ?
  //                    htr.find(cur->lt())->second : NULL);
  //    tr_copy->set_lb(cur->lb() ?
  //                    htr.find(cur->lb())->second : NULL);

  //    if (cur->dag_node()->is_inner_node())
  //    {
  //      child = &cur->dag_node()->right_child();
  //      while (child && child->is_inner_node() && !pr(*(*child)))
  //        child = &child->left_child();
  //      tr_copy->dag_node()->set_right_child(*child);
  //      child = &cur->dag_node()->left_child();
  //      while (child && child->is_inner_node() && !pr(*(*child)))
  //        child = &child->left_child();
  //      tr_copy->dag_node()->set_left_child(*child);
  //    }
  //    // Third iteration: generate links in-between trapezoids
  //    //  and in-between nodes .
  //  }
  //  m_dag_root = htr.find(&*(*td.m_dag_root))->second->dag_node();
  //}
  //

  /*
    TODO: Should we add another constructor with non const argument that
    rebuild the trapezoidal decomposition prior to copy construction?
  */
  virtual ~Trapezoidal_decomposition_2()
  {
    CGAL_warning(m_dag_root != NULL);
    if (!m_dag_root) return;

    delete m_dag_root;

    if (traits)
      delete traits;
  }

  //---------------------------------------------------------------------------
  // Description:
  //  if Halfedge or twin already inserted the latter is returned.
  //  otherwise the left-low most edge-degenerate trapezoid that represents the
  //  input Halfedge is returned
  // Remark:
  //  Given an edge-degenerate trapezoid representing a Halfedge,
  //  all the other trapezoids representing the Halfedge can be extracted
  //  via moving continously to the left and right neighbours.
  Td_map_item insert(Halfedge_const_handle he);


  //---------------------------------------------------------------------------
  // Description:
  // inserts a range of halfedges into the Search structure.
  // First it randomly shuffles the container and then it inserts the Halfedges
  //  according to the new order
  // Precondition: the data structure is empty
  template <class Halfedge_iterator>
  void insert(Halfedge_iterator begin, Halfedge_iterator end)
  {
    //Precondition: the data structure is empty
    CGAL_precondition(m_number_of_curves == 0);

    if (begin == end)
      return;

    //insert the shuffled halfedges into the search structure

    //disable the rebuild check from within the halfedge insert and check here
    //  for rebuild
    bool do_rebuild = set_with_guarantees(false);

    bool start_over = true;
    while (start_over)
    {
      start_over = false;

      //random_shuffle the range
      std::random_shuffle(begin,end);

      Halfedge_const_handle he_cst;
      Halfedge_iterator it = begin;
      for (; it < end ; ++it)
      {
        if (do_rebuild && not_within_limits())
        {
          std::cout << "starting over after " << number_of_curves() << std::flush;
          start_over = true;
          clear();
          break;
        }

        he_cst = *it;
        insert(he_cst);
      }
      if (it != end)
        continue;

      //after inserting the last halfedge in the range
      //  perform another rebuild check
      if (do_rebuild && not_within_limits()) //MICHAL: should I use needs_update() instead (with the random check)?
      {
        start_over = true;
        clear();
      }
    }

    //enable the rebuild from within the halfedge insert
    set_with_guarantees(do_rebuild);
  }


  // removal functions

  //---------------------------------------------------------------------------
  // Description:
  //
  void remove(Halfedge_const_handle he);

  ////-------------------------------------------------------------------------
  //// Description:
  ////
  //template <class curve_iterator>
  //void remove(curve_iterator begin, curve_iterator end)
  //{
  //  if(begin == end)
  //    return;
  //
  //  std::random_shuffle(begin,end);
  //
  //  curve_iterator it=begin,next=it;
  //  while(it!=end)
  //  {
  //    ++next;
  //    remove(*it);
  //    it=next;
  //  }
  //}

  void clear()
  {
    delete m_dag_root;
    init();
  }


//  //-------------------------------------------------------------------------
//  // Description:
//  //  returns the active trapezoid representing the input point.
//  // Precondition:
//  //  The trapezoidal tree is not empty
//  // Postcondition:
//  //  the input locate type is set to the type of the output trapezoid.
//  // Remark:
//  //  locate call may change the class
//  Td_map_item& locate_and_print(std::ostream& out, const Point& p) const
//  {
//
//#ifdef CGAL_TD_DEBUG
//
//    CGAL_assertion(traits);
//    CGAL_assertion(m_dag_root);
//
//#endif
//
//    Dag_node curr = *m_dag_root; //MICHAL: is it ok to add &?
//
//#ifdef CGAL_TD_DEBUG
//
//    CGAL_precondition(!!curr);
//
//#endif
//    //the actual locate. curr is the DAG root, the traits,
//    //the point to location, and 0 - indicates point location
//    search_and_print_using_dag(out, curr,traits,p,Halfedge_const_handle());//m_empty_he_handle);
//
//
//#ifndef CGAL_NO_TRAPEZOIDAL_DECOMPOSITION_2_OPTIMIZATION
//
//    locate_opt_push(curr.get_data());
//
//#endif
//
//    return *curr;
//  }

  //---------------------------------------------------------------------------
  // Description:
  //  returns the active trapezoid representing the input point.
  // Precondition:
  //  The trapezoidal tree is not empty
  // Postcondition:
  //  the input locate type is set to the type of the output trapezoid.
  // Remark:
  //  locate call may change the class
  Td_map_item& locate(const Point& p,Locate_type &t) const
  {
    //print_dag_addresses(*m_dag_root);
#ifdef CGAL_TD_DEBUG

    CGAL_assertion(traits);
    CGAL_assertion(m_dag_root);

#endif

    Dag_node curr = *m_dag_root; //MICHAL: is it ok to add &?

#ifdef CGAL_TD_DEBUG

    CGAL_precondition(!!curr);

#endif
    //the actual locate. curr is the DAG root, the traits,
    //the point to location, and 0 - indicates point location
    t = search_using_dag(curr,traits,p,Halfedge_const_handle());

#ifdef CGAL_TD_DEBUG

    CGAL_postcondition(t == POINT || t == CURVE || t == TRAPEZOID ||
                       t == UNBOUNDED_TRAPEZOID);

#endif

#ifndef CGAL_NO_TRAPEZOIDAL_DECOMPOSITION_2_OPTIMIZATION

    locate_opt_push(curr.get_data());

#endif

    return *curr;
  }

  //---------------------------------------------------------------------------
  // Description:
  //  returns the active trapezoid representing the input point.
  // Precondition:
  //  The trapezoidal tree is not empty
  // Postcondition:
  //  the input locate type is set to the type of the output trapezoid.
  // Remark:
  //  locate call may change the class
  Td_map_item& locate(const Curve_end& ce, Locate_type& lt) const
  {

#ifdef CGAL_TD_DEBUG

    CGAL_assertion(traits);
    CGAL_assertion(m_dag_root);

#endif

    Dag_node curr = *m_dag_root; //MICHAL: is it ok to add &?

    //the actual locate. curr is the DAG root, the traits,
    //the end point to locate,
    //and NULL as cv ptr - indicates point location
    lt = search_using_dag (curr, traits, ce, Halfedge_const_handle());

#ifdef CGAL_TD_DEBUG

    CGAL_postcondition(lt == POINT || lt == CURVE || lt == TRAPEZOID ||
                       lt == UNBOUNDED_TRAPEZOID);

#endif

#ifndef CGAL_NO_TRAPEZOIDAL_DECOMPOSITION_2_OPTIMIZATION

    locate_opt_push(curr.get_data());

#endif

    return curr.get_data();
  }

  //---------------------------------------------------------------------------
  // Description:
  //  returns the active trapezoid containing the point represented by vertex.
  // Precondition:
  //  The trapezoidal tree is not empty
  // Postcondition:
  //  the input locate type is set to the type of the output trapezoid.
  // Remark:
  //  locate call may change the class
  Td_map_item& locate( Vertex_const_handle v, Locate_type& lt) const
  {
    CGAL_precondition(traits != NULL);
    return locate(traits->vtx_to_ce(v), lt);
  }

  //---------------------------------------------------------------------------
  // Description:
  //
  // preconditions:
  //  p is not on an edge or a vertex.
  Td_map_item& vertical_ray_shoot(const Point& p,Locate_type& t,
                                  const bool up_direction = true) const;


  ////MICHAL: commented due to inefficient depth update, remove and insert instead
  //void before_split_edge(const X_monotone_curve_2& cv,
  //                       const X_monotone_curve_2& cv1,
  //                       const X_monotone_curve_2& cv2);

  ////MICHAL: commented due to inefficient depth update, remove and insert instead
  ////-------------------------------------------------------------------------
  //// Description:
  //// Input:
  ////  1 whole curves
  ////  2 partial halfedge_handle-s
  //// precondition:
  ////  The two halfedges are valid
  ////  The first input curve is the union of the two halfedges.
  ////  The intersection of the latter is a point inside the
  ////  interior of the former.
  ////  The latter are ordered from left-down to right-up
  //// postcondition:
  ////  The first input curve is broken into two halfedges
  ////  corresponding to the input.
  //void split_edge(const X_monotone_curve_2& cv, Halfedge_const_handle he1,
  //                Halfedge_const_handle he2);


  void merge_edge(Halfedge_const_handle he1, Halfedge_const_handle he2,
                  const X_monotone_curve_2& cv);


  void
  after_merge_edge(Halfedge_const_handle CGAL_precondition_code(merged_he),
                   Halfedge_const_handle CGAL_precondition_code(before_mrg_he))
  {
    //Precondition:
    // the merge uses the suspected halfedge before the arrangement merge
    CGAL_precondition(merged_he == before_mrg_he ||
                      merged_he == before_mrg_he->twin());
  }


  unsigned long size() const
  {
    return m_dag_root->size();
  }

  unsigned long number_of_curves() const
  {
    return m_number_of_curves;
  }

  void init_arrangement_and_traits(const Arrangement_on_surface_2* arr,
                                   bool allocate_traits = true)
  {
    m_arr = arr;
    m_trts_adaptor =
      static_cast<const Traits_adaptor_2*> (arr->geometry_traits());

    if (allocate_traits)
      traits = new Td_traits(*m_trts_adaptor);
  }



#ifdef CGAL_TD_DEBUG
  /*------------------------------------------------------------------
    description:
    returns whether the Trapezoidal Dag is valid
  */
  bool is_valid(const Dag_node& ds) const
  {
    if ( !ds ) return true;
    if (ds->is_valid(traits) && ds->dag_node() &&
        is_valid(ds.left_child()) && is_valid(ds.right_child()))
      return true;
    CGAL_warning(ds->is_valid(traits));
    CGAL_warning(ds->dag_node());
    CGAL_warning(is_valid(ds.left_child()));
    CGAL_warning(is_valid(ds.right_child()));
    return false;
  }
  /*------------------------------------------------------------------
    description:
    returns whether the member Trapezoidal data structure is valid
  */
  bool is_valid() const
  {
    return is_valid(*m_dag_root);
  }
  //void debug() const
  //{
  //  std::cout << "\nTrapezoidal_decomposition_2<Traits>::debug()\n" << *this
  //            << std::endl;
  //  Td_map_item x;
  //  x.debug(); //MICHAL: will not work!
  //}
#endif

  /*------------------------------------------------------------------
    description:
    Rebuilds the trapezoid data structure by reinserting the curves
    in a random order in an empty data structure.

    postcondition:
    The old and new data structures agree on their member curves.
    ------------------------------------------------------------------*/
  Self& rebuild()
  {
#ifdef CGAL_TD_DEBUG
    std::cout << "\nrebuild!  " << m_number_of_curves << std::endl
              << std::flush;
#endif

    Halfedge_container container;

#ifdef CGAL_TD_DEBUG
    unsigned long rep = Halfedge_filter(container, &dag_root());
#endif

    clear();

    //// initialize container to point to curves in Td_map_item Tree
    //if (rep>0)
    //{
    //  bool o = set_with_guarantees(false);
    //  typename std::vector<Halfedge_const_handle>::iterator
    //    it = container.begin(),
    //    it_end = container.end();
    //  while(it!=it_end)
    //  {
    //    insert(*it);
    //    ++it;
    //  }
    //  set_with_guarantees(o);
    //}

    //insert the already inserted curves from scratch in order to build a
    // search structure guaranteeing logarithmic query time and linear size
    insert(container.begin(), container.end());

#ifdef CGAL_TD_DEBUG
    CGAL_assertion(is_valid());
    unsigned long sz = number_of_curves();
    if (sz != rep)
    {
      std::cerr << "\nnumber_of_curves()=" << sz;
      std::cerr << "\nrepresentatives.size()=" << rep;
      CGAL_assertion(number_of_curves() == rep);
    }
#endif

    container.clear();
    return *this;
  }

  /*
     Input:
     a list of pointers to Td_map_items and a Td_map_item boolean predicate.
     Output:
     void
     Postcondition:
     the list pointers correspond to all the Td_map_items in the data
     structure for which the predicate value is true.
  */

  template <typename Container, typename Predicate>
  void filter(Container& c, const Predicate& pr, const Dag_node * ds) const
  {
    CGAL_assertion(ds);
    ds->filter(c,pr);
  }

  template <typename Container, typename Predicate>
  void filter(Container& c, const Predicate& pr) const
  {
    filter(c, pr, &dag_root());
  }

  ////MICHAL: not in use
  //template <class Container>
  //unsigned long Td_map_item_filter(Container& container,
  //                                 const Dag_node* ds) const
  ///* Return a container for all active map items */
  //{
  //  ds->filter(container, Td_active_map_item());
  //  return container.size();
  //}

  template <typename Halfedge_container>
  unsigned long Halfedge_filter(Halfedge_container& container,
                                const Dag_node* ds) const
  /* Return a container for all active curves */
  {
    unsigned long sz = number_of_curves();
    std::list<Td_map_item> representatives;
    //X_trapezoid_list representatives;
    ds->filter(representatives, Td_active_edge_item(*traits));

#ifndef CGAL_TD_DEBUG

    CGAL_warning(sz==representatives.size());

#else

    unsigned long rep=representatives.size();
    if (sz != rep)
    {
      std::cerr << "\nnumber_of_curves()=" << sz;
      std::cerr << "\nrepresentatives.size()=" << rep;
      CGAL_assertion(number_of_curves()==representatives.size());
    }

#endif

    if (sz > 0)
    {
      typename std::list<Td_map_item>::iterator it = representatives.begin(),
        it_end = representatives.end();
      //typename X_trapezoid_list::iterator it = representatives.begin(),
      //  it_end = representatives.end();
      while(!(it==it_end))
      {
        Td_active_edge e (boost::get<Td_active_edge>(*it));
        container.push_back(e.halfedge()); //it represents an active trapezoid
        ++it;
      }
    }
    if(! container.empty()) {
      std::random_shuffle(container.begin(),container.end());
    }
    return sz;
  }



  /*------------------------------------------------------------------
    Input: None
    Output: bool
    Description:
    determines according to pre defined conditions whether the
    current Trapezoidal_decomposition_2<Traits> needs update
    Postconditions:
    The output is true iff the depth of the Trapezoidal Tree is more then
    DepthThreshold times log of the X_curve count or the Trapezoidal Tree's
    size
    is more then SizeThreshold times the log of the last count.
  */
  bool set_with_guarantees(bool u)
  {
    bool old = m_with_guarantees;
    m_with_guarantees = u;
    return old;
  }

  //This method occasionaly(!) checks the guarantees
  // It is currently not in use, since the guarantees are constantly checked in O(1) time
  bool needs_update()
  {
    unsigned long num_of_cv = number_of_curves();
    //to avoid signed / unsigned conversion warnings
    // rand() returns an int but a non negative one.
    if (static_cast<unsigned long>(std::rand()) >
        RAND_MAX / ( num_of_cv + 1))
      return false;
    /*       INTERNAL COMPILER ERROR overide
             #ifndef __GNUC__
    */
#ifdef CGAL_TD_REBUILD_DEBUG
    std::cout << "\n|heavy!" << std::flush;
#endif

    return not_within_limits();
  }

  /*------------------------------------------------------------------
    input: None
    output: bool
    Description:
    uses needs_update to determine whether the
    Trapezoidal_decomposition_2<Traits> needs update
    and calls rebuild accordingly
    Postcondition:
    the return value is true iff rebuilding took place.
  */
  bool update()
  {
    //if the structure violates the guarantees - rebuild
    if (not_within_limits()) //needs_update())
    {
      rebuild();
      return true;
    }

    return false;
  }

  bool not_within_limits()
  {
    unsigned long num_of_cv = number_of_curves();

    //Cond 1: Depth is greater than threshold*log(number of curves)
    bool cond1 = largest_leaf_depth() >
                   (depth_threshold()*(std::log(double(num_of_cv+1))));
    //Cond 2: Number of nodes is greater than threshold*number of curves
    bool cond2 = number_of_dag_nodes() > (size_threshold()*(num_of_cv + 1));

    //return true if at least one of the conditions is true
    return cond1 || cond2;
  }

  /* returns a reference to the internal data structure */
  const Dag_node& dag_root() const {return *m_dag_root;}

  /* returns a reference to the internal data structure */
  const Traits& get_traits() const {return *traits;}

  /* returns a reference to the internal depth threshold constant */
  const double& depth_threshold() const
  {
    return m_depth_threshold;
  }

  /* sets the internal depth threshold constant to the parameter and
     returns its reference */
  void depth_threshold(const double& depth_th)
  {
    m_depth_threshold = depth_th;
  }

  /* returns a reference to the internal size threshold constant */
  const double& size_threshold() const
  {
    return m_size_threshold;
  }

  /* sets the internal size threshold constant to the parameter and
     returns its reference */
  void size_threshold(const double& size_th)
  {
    m_size_threshold = size_th;
  }

  void update_largest_leaf_depth( unsigned long depth )
  {
    if(m_largest_leaf_depth < depth )
        m_largest_leaf_depth = depth;
  }

  unsigned long largest_leaf_depth()
  {
    //CGAL_assertion((m_largest_leaf_depth + 1) == m_dag_root->rec_depth());
    return m_largest_leaf_depth;
  }

  unsigned long number_of_dag_nodes()
  {
    //CGAL_assertion(m_number_of_dag_nodes == m_dag_root->size());
    return m_number_of_dag_nodes;
  }

  unsigned long longest_query_path_length()
  {
    return longest_query_path_length_rec(true, *m_dag_root,
                                         true, *m_dag_root, *m_dag_root);
  }


protected:

  //Trapezoidal Decomposition data members
  Dag_node* m_dag_root;
  unsigned long m_largest_leaf_depth; //holds the leargest depth of a leaf in the DAG
  unsigned long m_number_of_dag_nodes; //holds the number of nodes in the DAG
  bool m_with_guarantees; //whether the structure holds logarithmic query time and linear size guarantees //m_needs_update;
  unsigned long m_number_of_curves;
  const Traits* traits;
  //Before_split_data m_before_split;
  const Arrangement_on_surface_2* m_arr;
  const Traits_adaptor_2* m_trts_adaptor;

  Halfedge_const_handle m_empty_he_handle;

private:

#ifndef CGAL_NO_TRAPEZOIDAL_DECOMPOSITION_2_OPTIMIZATION

  mutable Td_map_item last_cv;
  mutable Td_map_item prev_cv;

#endif

  unsigned long longest_query_path_length_rec(bool minus_inf,
                                              Dag_node& min_node,
                                              bool plus_inf,
                                              Dag_node& max_node,
                                              Dag_node& node);

  void init()
  {
    // traits may be initialized later
    m_dag_root = new Dag_node(Td_active_trapezoid());
    //(*m_dag_root)->set_dag_node(m_dag_root);
    boost::apply_visitor(set_dag_node_visitor(m_dag_root),
                         m_dag_root->get_data());

    m_number_of_curves = 0;
    m_largest_leaf_depth = 0;
    m_number_of_dag_nodes = 1; //the root is the only node in the DAG

#ifndef CGAL_NO_TRAPEZOIDAL_DECOMPOSITION_2_OPTIMIZATION

    locate_opt_empty();

#endif

  }

#ifndef CGAL_NO_TRAPEZOIDAL_DECOMPOSITION_2_OPTIMIZATION

  void locate_opt_push(Td_map_item& cv_tr) const
  {
    prev_cv = last_cv;
    last_cv = cv_tr;
  }
  void locate_opt_empty() const
  {
    last_cv = Td_map_item(0);
    prev_cv = Td_map_item(0);
  }
  bool locate_opt_swap(Td_map_item& item) const
  {
    item = last_cv;
    last_cv = prev_cv;
    prev_cv = item;
    return (!traits->is_empty_item(item));
  }
  void locate_optimization(const Curve_end& ce, Td_map_item& item,
                            Locate_type& lt) const
  {
    bool res = false;
    // optimization
    if (locate_opt_swap(item) && traits->is_active(item) )
    {
      if (traits->is_td_vertex(item))
        res = are_equal_end_points(ce, item);

      //if (traits->is_td_vertex(item))
      //{
      //  if ( traits->is_fictitious_vertex(item) )
      //  {
      //    res = traits->equal_curve_end_2_object()(ce, *(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),item)));
      //  }
      //  else
      //  {
      //    res = traits->equal_curve_end_2_object()(ce, boost::apply_visitor(point_for_vertex_visitor(), item));
      //  }
      //}
      if (traits->is_td_trapezoid(item))
        res = traits->is_inside(item,ce);
    }
    if (!res && locate_opt_swap(item) && traits->is_active(item) )
    {
      if (traits->is_td_vertex(item))
        res = are_equal_end_points(ce, item);
      //if (traits->is_td_vertex(item))
      //{
      //  if ( traits->is_fictitious_vertex(item) )
      //  {
      //    res = traits->equal_curve_end_2_object()(ce, *(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),item)));
      //  }
      //  else
      //  {
      //    res = traits->equal_curve_end_2_object()(ce, boost::apply_visitor(point_for_vertex_visitor(), item));
      //  }
      //}
      if (traits->is_td_trapezoid(item))
        res = traits->is_inside(item,ce);
    }
    if (res)
    {
      if (traits->is_td_vertex(item))
        lt=POINT;
      else
      {
        Td_active_trapezoid tr (boost::get<Td_active_trapezoid>(item));
        lt = tr.is_on_boundaries()? UNBOUNDED_TRAPEZOID : TRAPEZOID;
      }
    }
    else
      item = locate(ce,lt);
  }

#endif


  void print_cv_data(const X_monotone_curve_2& cv,
                     std::ostream& out = std::cout) const
  {
    out << "min end: " << std::endl;
    print_ce_data(cv, ARR_MIN_END, out);
    out << std::endl << "max end: " << std::endl;
    print_ce_data(cv, ARR_MAX_END, out);
    out << std::endl << std::endl ;
  }

  void print_point_data(const Point& p, std::ostream& out = std::cout) const
  {
    out << "x: " << CGAL::to_double(p.x()) << ", y: "
        << CGAL::to_double(p.y()) << std::endl;
  }

  void print_ce_data(const X_monotone_curve_2& cv, Arr_curve_end ce,
                     std::ostream& out = std::cout) const
  {
    Arr_parameter_space ps_x = traits->parameter_space_in_x_2_object()(cv, ce);
    Arr_parameter_space ps_y = traits->parameter_space_in_y_2_object()(cv, ce);

    if (ps_x == ARR_INTERIOR && ps_y == ARR_INTERIOR)
    {
      if (ce == ARR_MIN_END)
        out << "x: "
            << CGAL::to_double(traits->construct_min_vertex_2_object()(cv).x())
            << ", y: "
            << CGAL::to_double(traits->construct_min_vertex_2_object()(cv).y())
            << std::endl;
      else
        out << "x: "
            << CGAL::to_double(traits->construct_max_vertex_2_object()(cv).x())
            << ", y: "
            << CGAL::to_double(traits->construct_max_vertex_2_object()(cv).y())
            << std::endl;
    }
    else if (ps_x == ARR_INTERIOR && ps_y != ARR_INTERIOR)
    {
      out << " vertical asymptote, " ;
      if (ps_y == ARR_TOP_BOUNDARY)
        out << " y -> +oo " << std::endl;
      else
        out << " y -> -oo " << std::endl;
    }
    else if (ps_x != ARR_INTERIOR && ps_y == ARR_INTERIOR)
    {
      out << " horizontal asymptote, " ;
      if (ps_x == ARR_RIGHT_BOUNDARY)
        out << " x -> +oo " << std::endl;
      else
        out << " x -> -oo " << std::endl;
    }
    else //both are not interior
    {
      if (ps_x == ARR_RIGHT_BOUNDARY)
        out << " x -> +oo " ;
      else
        out << " x -> -oo " ;
      if (ps_y == ARR_TOP_BOUNDARY)
        out << " , y -> +oo " << std::endl;
      else
        out << " , y -> -oo " << std::endl;

    }
  }


  void print_dag_addresses(const Dag_node& curr) const
  {

    std::cout << "----------------- DAG ----------------" <<std::endl
              << "--------------------------------------" <<std::endl;

    print_dag_addresses_rec(curr, 0);
    std::cout << "----------------- END OF DAG ----------------" <<std::endl
              << "---------------------------------------------" <<std::endl;

  }

  void print_dag_addresses_rec(const Dag_node& curr ,int level,
                               std::ostream& out = std::cout) const
  {

    out << "------ level " << level << ", depth " << curr.depth()
        << " ------\n";
    out << " (void *)curr : " << (void *)(&curr) << std::endl;
    out << "      (void *)curr->TRPZ : " << (void *)(curr.operator->())
        << std::endl;

    //curr is the current pointer to node in the data structure
    //item holds the trapezoidal map item connected to curr.
    Td_map_item item = curr.get_data();
    if (traits->is_td_vertex(item))
    {
      out << " POINT : " ;

      // if the map item represents a fictitious vertex
      if (traits->is_fictitious_vertex(item))
      {
        const Curve_end left_ce(*(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),item)));
        print_ce_data(left_ce.cv(), left_ce.ce(), out);
      }
      else // if the map item represents a vertex
      {
        Point p = boost::apply_visitor(point_for_vertex_visitor(),item);
        print_point_data(p, out);
      }
      out << "          (void *)left_child: " << (void*)(&(curr.left_child()))
          << std::endl;
      out << "          (void *)right_child: " << (void*)(&(curr.right_child()))
          << std::endl;
      print_dag_addresses_rec(curr.left_child(), level+1, out);
      print_dag_addresses_rec(curr.right_child(), level+1, out);
      return;
    }
    if (traits->is_td_edge(item))
    {
      // bool is_active = traits->is_active(item);
      // if the map item represents an edge
      const X_monotone_curve_2& he_cv = *(boost::apply_visitor(cv_for_edge_visitor(), item));

      //   so top() is a real Halfedge with a curve() if curr is active
      //   or curr holds the curve if curr is not active
      out << " CURVE : " ;
      print_cv_data(he_cv, out);
      out << "          (void *)left_child: " << (void*)(&(curr.left_child()))
          << std::endl;
      out << "          (void *)right_child: " << (void*)(&(curr.right_child()))
          << std::endl;
      print_dag_addresses_rec(curr.left_child(), level+1, out);
      print_dag_addresses_rec(curr.right_child(), level+1, out);
      return;
    }
    else
    {
      // if ithe map item represents a trapezoid
      if (traits->is_active(item))
        out << " TRAPEZOID \n";
      else //trapezoid is removed - may have a left child
      {
        out << " REMOVED TRAPEZOID \n";
        if (curr.left_child().is_null())
          return;
        out << "          (void *)left_child: " << (void*)(&(curr.left_child()))
            << std::endl;
        print_dag_addresses_rec(curr.left_child(), level+1, out);
      }
    }
  }

public:
  void print_dag(std::ostream& out) const
  {

    out << "----------------- DAG ----------------" << std::endl
        << "--------------------------------------" << std::endl;

    print_dag_addresses_rec(*m_dag_root , 0, out);
    out << "----------------- END OF DAG ----------------" << std::endl
        << "---------------------------------------------" << std::endl;

  }

protected:
  double m_depth_threshold;
  double m_size_threshold;
};

} //namespace CGAL

#include <CGAL/Arr_point_location/Td_active_trapezoid.h>
#include <CGAL/Arr_point_location/Td_inactive_trapezoid.h>
#include <CGAL/Arr_point_location/Td_active_edge.h>
#include <CGAL/Arr_point_location/Td_inactive_edge.h>
#include <CGAL/Arr_point_location/Td_active_vertex.h>
#include <CGAL/Arr_point_location/Td_active_fictitious_vertex.h>
#include <CGAL/Arr_point_location/Td_inactive_vertex.h>
#include <CGAL/Arr_point_location/Td_inactive_fictitious_vertex.h>

// The member-function definitions can be found under:
#include <CGAL/Arr_point_location/Trapezoidal_decomposition_2_impl.h>

#endif
