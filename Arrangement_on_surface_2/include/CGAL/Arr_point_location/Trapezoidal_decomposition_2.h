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
// Author(s)     : Oren Nechushtan <theoren@math.tau.ac.il>
//                 Iddo Hanniel <hanniel@math.tau.ac.il>

#ifndef CGAL_TRAPEZOIDAL_DECOMPOSITION_2_H
#define CGAL_TRAPEZOIDAL_DECOMPOSITION_2_H



#include <CGAL/Arr_tags.h>
#include <CGAL/basic.h>
#include <CGAL/Arr_point_location/Td_predicates.h>
#include <CGAL/Arr_point_location/Trapezoidal_decomposition_2_misc.h>

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <cmath>
#include <ctime>
#include <list>
#include <vector>
#include <map>

namespace CGAL {

#define CGAL_CURVE_COMPARE_Y_AT_X(p,cv)         \
  (traits->compare_y_at_x_2_object()((p),(cv)))   //MICHAL: should be removed

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
  class Base_trapezoid_iterator;
  class In_face_iterator;
  friend class In_face_iterator;
  class Around_point_circulator;

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

  //type of Curve end: a reference to X_monotone_curve_2 & Arr_curve_end
  typedef typename Traits::Curve_end    Curve_end;

  //type of trapezoid
  typedef typename Traits::X_trapezoid X_trapezoid;

  //type of Curve end pair
  typedef typename Traits::Curve_end_pair  Curve_end_pair;

  //type of trapezoids' list
  typedef std::list<X_trapezoid> X_trapezoid_list;
  
  //type of trapezoids' vector
  typedef std::vector<X_trapezoid> X_trapezoid_vector;
  
  //type of Halfedge_const_handle-s' vector
  typedef std::vector<Halfedge_const_handle> Halfedge_container;

  //type of base trapezoid circulator
  typedef class Base_trapezoid_iterator Base_trapezoid_circulator;
  

  //predicates
  typedef CGAL::Td_active_trapezoid<X_trapezoid> Td_active_trapezoid;
  typedef CGAL::Td_active_non_degenerate_trapezoid<X_trapezoid,Traits> 
  Td_active_non_degenerate_trapezoid; //MICHAL: not in use
  typedef CGAL::Td_active_right_degenerate_curve_trapezoid<X_trapezoid,Traits> 
  Td_active_right_degenerate_curve_trapezoid;

  //type of search structure DAG node
  typedef Td_dag_node< X_trapezoid> Dag_node; 

#if 0
  //type of vertex Trapezoidal map item 
  typedef typename Traits::Td_vertex Td_vertex;

  //type of halfedge Trapezoidal map item
  typedef typename Traits::Td_halfedge Td_halfedge;

  //type of trapezoid Trapezoidal map item
  typedef typename Traits::Td_trapezoid Td_trapezoid;

  //type of Td_map item (Td_halfedge, Td_vertex , or Td_trapezoid)
  typedef typename Traits::Td_map_item Td_map_item;
#endif

  //type of map of DAG nodes
  typedef std::map<int,Dag_node> Nodes_map;

  //type of trapezoids comparison function - for the map
  typedef Trapezoid_handle_less<const X_trapezoid* const> Trapezoid_ptr_less;
  
  //type of trapezoids ptr map
  typedef std::map<const X_trapezoid*, X_trapezoid*, Trapezoid_ptr_less> 
  X_trapezoid_ptr_map;
  
  typedef typename Traits::Left_side_category   Left_side_category;
  typedef typename Traits::Bottom_side_category Bottom_side_category;
  typedef typename Traits::Top_side_category    Top_side_category;   
  typedef typename Traits::Right_side_category  Right_side_category; 

protected:

  typedef typename Arr_are_all_sides_oblivious_tag< 
                     Left_side_category, Bottom_side_category, 
                     Top_side_category, Right_side_category >::result
  Are_all_sides_oblivious_tag;

public:
  
  /*! \class Base_trapezoid_iterator
   * member of Trapezoidal_decomposition_2<Traits>
   * Description Implements a basic Trapezoid iterator
   */
  class Base_trapezoid_iterator
  {
  public:
    //constructors
    Base_trapezoid_iterator() : traits(0), m_cur_trp(0){ }

    Base_trapezoid_iterator(const Traits* traits_, 
                            X_trapezoid* currt=0)
          :traits(traits_), m_cur_trp(currt) { }

    Base_trapezoid_iterator(const Base_trapezoid_iterator &it)
          :traits(it.traits), m_cur_trp(it.m_cur_trp) { }

    //operator overloading
    Base_trapezoid_iterator  & operator=(const Base_trapezoid_iterator &it)
    {
      traits = it.traits;
      m_cur_trp = it.m_cur_trp;
      return *this;
    }

    bool operator==(const Base_trapezoid_iterator &it) const
    {
      return (m_cur_trp == it.m_cur_trp);
    }

    bool operator!=(const Base_trapezoid_iterator &it) const
    {
      return !operator==(it);
    }
      
    X_trapezoid& operator*() const
    {
      CGAL_precondition(m_cur_trp != NULL);
      return *m_cur_trp;
    }

    X_trapezoid* operator->() const
    {
      return m_cur_trp;
    }

    bool operator!() const
    {
      return m_cur_trp == 0;
    }
    
  protected:
    const Traits* traits; //pointer to the traits
    X_trapezoid* m_cur_trp; //pointer to the current trapezoid
  };

/*! \class In_face_iterator
   * member of Trapezoidal_decomposition_2<Traits>
   * Derived from Base_trapezoid_iterator class
   * Description Implements a Trapezoid iterator along a Halfedge
   */
  class In_face_iterator : public Base_trapezoid_iterator
  {

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
    using Base_trapezoid_iterator::m_cur_trp;
    using Base_trapezoid_iterator::traits;
#endif

  protected:
    //reference to the seperating X_monotone_curve_2
    const X_monotone_curve_2& m_sep;

  public:
    //constructors
    In_face_iterator(const Traits* traits_,
                     Halfedge_const_handle sep, X_trapezoid* currt=0) 
            :Base_trapezoid_iterator(traits_,currt), m_sep(sep->curve()) 
    { }

    In_face_iterator(const Traits* traits_,
                     const X_monotone_curve_2& sep, X_trapezoid* currt=0) 
            :Base_trapezoid_iterator(traits_,currt), m_sep(sep) 
    { }


    In_face_iterator(const In_face_iterator &it) 
            :Base_trapezoid_iterator((Base_trapezoid_iterator&)it),
              m_sep(it.m_sep)                                   
    { 
    }

    //operatoror overloading
    bool operator==(const In_face_iterator &it) const
    {
      return ( Base_trapezoid_iterator::operator==(it) && 
               traits->equal_2_object()(m_sep,it.m_sep));
    }
    
    /*
      destription:
      advances m_cur_trp to one of the right neighbours according to the relation
      between the seperating Halfedge and the right() trapezoid point.
      precoditions:
      sep doesn't intersect no existing edges except possibly on common end
      points.
      postconditions:
      if the rightest trapezoid was traversed m_cur_trp is set to NULL.
      remark:
      if the seperator is vertical, using the precondition assumptions it
      follows that there is exactly one trapezoid to travel.
    */
    In_face_iterator& operator++()
    {
      if (!m_cur_trp) 
        return *this;// end reached, do nothing!
      
#ifndef CGAL_TD_DEBUG
      CGAL_warning(traits != NULL);
#else
      CGAL_assertion(traits != NULL);
      CGAL_assertion(m_cur_trp->is_active());
      //m_cur_trp should be a regular trapezoid or an edge
      CGAL_assertion(!traits->is_degenerate_point(*m_cur_trp));
#endif
      
      if (!traits->is_degenerate(*m_cur_trp))
      {
        //if the trapezoid is a regular trapezoid (not an edge)

#ifndef NDEBUG
#ifndef CGAL_TD_DEBUG
        CGAL_warning_code(Dag_node* tt = m_cur_trp->dag_node();)
          CGAL_warning(!tt->is_inner_node());
#else
        CGAL_assertion_code(Dag_node* tt = m_cur_trp->dag_node();)
          CGAL_assertion(tt);
        CGAL_assertion(!tt->is_inner_node());
#endif
#endif
          
        // handle degeneracies
        if (traits->compare_curve_end_xy_2_object()
                     (m_cur_trp->left()->curve_end(),
                      Curve_end(m_sep,ARR_MAX_END)) != SMALLER)
        {
          //if the trapezoid's left end point is equal to or larger from the 
          //  max end of sep, we reached the end of the iterator
          m_cur_trp = 0;
        }
        else
        {
          //if the trapezoid's left end point is smaller from the sep's max end

          //comparing the y value of the trapezoid's right end point and sep 
          //   (at the trapezoid's right x value), in order to select the
          //    next trapezoid in the iterator
          switch(traits->compare_curve_end_y_at_x_2_object()
                                     (m_cur_trp->right()->curve_end(), m_sep))
          {
           case SMALLER:
              m_cur_trp = m_cur_trp->rt();
            break;
           case LARGER:
              m_cur_trp = m_cur_trp->rb();
            break;
           case EQUAL:
            // end reached
              m_cur_trp = 0;
            break;
           default:       
              m_cur_trp = 0;
            break;
          }
        }
      }
      else 
      {
        //if the trapezoid is an edge (not a regular trapezoid)

#ifndef NDEBUG          
#ifndef CGAL_TD_DEBUG
          
        CGAL_warning_code(Dag_node* tt = m_cur_trp->dag_node();)
        CGAL_warning(tt != NULL);
        CGAL_warning(tt->is_inner_node());
          
#else

        CGAL_assertion_code(Dag_node* tt = m_cur_trp->dag_node();)
        CGAL_assertion(tt != NULL);
        CGAL_assertion(tt->is_inner_node());
#endif
#endif
          
        //go to rb() neighbour of the current edge.
        // as long as there is an edge fragment of the same edge - rb() 
        //  is not null. If it's null we reached the last fragment of the edge
        m_cur_trp = m_cur_trp->rb();
        if (m_cur_trp)
        {
          //if rb() is not null, find the next real edge fragment trapezoid
          //    (skip points)
          while(traits->is_degenerate_point(*m_cur_trp))
            m_cur_trp = m_cur_trp->dag_node()->left_child().operator->();
              
          //make sure we stopped in an edge
#ifndef CGAL_TD_DEBUG
              
          CGAL_warning(traits->is_degenerate_curve(*m_cur_trp));
              
#else
              
          CGAL_precondition(traits->is_degenerate_curve(*m_cur_trp));
              
#endif
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
  };
  

  /*
   * class Around_point_circulator
   * member of Trapezoidal_decomposition_2<Traits>
   * Description Implements a Trapezoid circulator around a point
   */
  class Around_point_circulator : public Base_trapezoid_circulator
  {

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
    using Base_trapezoid_circulator::m_cur_trp;
    using Base_trapezoid_circulator::traits;
#endif

  protected:
    const Curve_end& m_fixed;
  public:
    
#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
    using Base_trapezoid_circulator::operator!;
#endif

    //constructor
    Around_point_circulator(const Traits* traits_, 
                            const Curve_end& fixed, X_trapezoid* currt) 
                    :Base_trapezoid_iterator(traits_,currt),m_fixed(fixed) 
    { }
    
    Around_point_circulator(const Around_point_circulator &it) 
                    :Base_trapezoid_iterator(it),m_fixed(it.m_fixed)
    { }
    
    Around_point_circulator& operator++()
    {
      if (operator!()) return *this;
      
#ifndef CGAL_TD_DEBUG
      
      CGAL_warning(
          (!m_cur_trp->is_on_left_boundary() &&
          traits->equal_curve_end_2_object()
                    (m_fixed, m_cur_trp->left()->curve_end())) ||
          (!m_cur_trp->is_on_right_boundary() &&
          traits->equal_curve_end_2_object()
                    (m_fixed, m_cur_trp->right()->curve_end())));
      
#else
      
      CGAL_precondition(
          (!m_cur_trp->is_on_left_boundary() &&
          traits->equal_curve_end_2_object()
                    (m_fixed, m_cur_trp->left()->curve_end())) ||
          (!m_cur_trp->is_on_right_boundary() &&
          traits->equal_curve_end_2_object()
                    (m_fixed, m_cur_trp->right()->curve_end())));
      
#endif
      
      m_cur_trp = operator->();
      return *this;
    }

    Around_point_circulator operator++(int)
    {
      Around_point_circulator tmp = *this;
      ++*this;
      return tmp;
    }

    X_trapezoid* operator[](int i) const
    {
      Around_point_circulator c = *this;
      while(i-->0) c++;
      return c.m_cur_trp;
    }

    /* returns reference to the next trapezoid
       on a clockwise orientation rotation with centre
       taken as the fixed point (edge_end)
       preconditions:
       ciruclator is not empty*/
    X_trapezoid& operator*() const
    {  
      CGAL_precondition(!operator!());
      return *operator->();
    }

    /* returns pointer to the next trapezoid
       on a clockwise orientation rotation with centre
       taken as the fixed point */
    X_trapezoid* operator->() const
    {
      X_trapezoid* cand;
      if (operator!()) return m_cur_trp;
      
      cand = is_right_rotation() ? m_cur_trp->rt() : m_cur_trp->lb();

      if (traits->is_degenerate_curve(*cand)) return cand;
      
      // cand was split by a point
      while(traits->is_degenerate_point(*cand))
      {
        if ((cand->is_active() && traits->compare_curve_end_xy_2_object()
                            (cand->left()->curve_end(),m_fixed) == SMALLER ) ||
            (!cand->is_active() && 
              ((!cand->is_on_boundaries() && 
                traits->compare_curve_end_xy_2_object()
                   (cand->point_for_inner_rem_vtx(),m_fixed) == SMALLER ) ||
               (cand->is_on_boundaries() && 
                traits->compare_curve_end_xy_2_object()
                   (cand->curve_end_for_boundary_rem_vtx(),m_fixed) == SMALLER ))))
        {
          // move right using data structure
          cand = cand->dag_node()->right_child().operator->();
        }
        else
        {
          // move left using data structure
          cand = cand->dag_node()->left_child().operator->();
        }
      }
      return cand;
    }

    /* returns true if the current trapezoid is NULL 
       or if the fixed point is on its left
       or if the fixed point is on its right */
    bool is_valid() const
    {
      if ((!m_cur_trp)||
          (!m_cur_trp->is_on_left_boundary() &&
            traits->equal_curve_end_2_object()(m_fixed,m_cur_trp->left())) ||
          (!m_cur_trp->is_on_right_boundary() &&
            traits->equal_curve_end_2_object()(m_fixed,m_cur_trp->right()))) 
      {
        return true;
      }
      else 
      {
#ifdef CGAL_TD_DEBUG
        std::cerr << "\nthis=";
        write(std::cerr,*m_cur_trp,*traits,false) << std::flush;
        std::cerr << "\nfixed=" << m_fixed << std::flush;
        CGAL_warning(!(m_cur_trp && m_cur_trp->is_on_left_boundary() &&
                       m_cur_trp->is_on_right_boundary()));
#endif
        return false;
      }
    }
    
    /* description:
       inserts the input trapezoid between the
       current trapezoid and the next trapezoid
       on a clockwise orientation rotation with
       centre taken as the fixed point.
       preconditions:
       current trapezoid exist
       input trapezoid is adjacent to fixed point
    */
    void insert(X_trapezoid& tr)
    {
#ifndef CGAL_TD_DEBUG
      
      CGAL_precondition(m_cur_trp != NULL);
      CGAL_precondition(tr.is_active());
      CGAL_warning(
            (!tr.is_on_left_boundary() &&
              traits->equal_curve_end_2_object()
                    (tr.left()->curve_end(), m_fixed)) ||
            (!tr.is_on_right_boundary() &&
              traits->equal_curve_end_2_object()
                    (tr.right()->curve_end(), m_fixed)));
      
#else
      
      CGAL_precondition(m_cur_trp != NULL);
      CGAL_precondition(tr.is_active());
      CGAL_precondition(
            (!tr.is_on_left_boundary() &&
              traits->equal_curve_end_2_object()
                    (tr.left()->curve_end(), m_fixed)) ||
            (!tr.is_on_right_boundary() &&
              traits->equal_curve_end_2_object()
                    (tr.right()->curve_end(), m_fixed)));
      
#endif

      if (!tr.is_on_left_boundary() &&
          traits->equal_curve_end_2_object()
                    (tr.left()->curve_end(),m_fixed))
        tr.set_lb(operator->());
      else
        tr.set_rt(operator->());

      if (is_right_rotation())
        m_cur_trp->set_rt(&tr);
      else
        m_cur_trp->set_lb(&tr);
    }

    /* precondition:
       m_cur_trp!=NULL
    */
    void remove()
    {
      
#ifndef CGAL_TD_DEBUG
      
      CGAL_precondition(m_cur_trp != NULL);
      CGAL_precondition(m_cur_trp->is_active());
      CGAL_warning( 
          ( !m_cur_trp->is_on_left_boundary() &&
            traits->equal_curve_end_2_object()
                    (m_cur_trp->left()->curve_end(), m_fixed) ) ||
          ( !m_cur_trp->is_on_right_boundary() &&
            traits->equal_curve_end_2_object()
                    (m_cur_trp->right()->curve_end(), m_fixed) ) );
      
#else
      
      CGAL_precondition(m_cur_trp != NULL);
      CGAL_precondition(m_cur_trp->is_active());
      CGAL_warning( 
          ( !m_cur_trp->is_on_left_boundary() &&
            traits->equal_curve_end_2_object()
                    (m_cur_trp->left()->curve_end(), m_fixed) ) ||
          ( !m_cur_trp->is_on_right_boundary() &&
            traits->equal_curve_end_2_object()
                    (m_cur_trp->right()->curve_end(), m_fixed) ) );
      
#endif
      
      Around_point_circulator old = *this;
      old++;
      X_trapezoid* next = old.operator->();
      // handle 1-cycle and 2-cycles seperately
      if (m_cur_trp != next)
      {
      }
      // 2-cycle
      else if (*this != old)
      {
        next = m_cur_trp;
      }
      // 1-cycle
      else
      {
        if (is_right_rotation())
          m_cur_trp->set_rt(0);
        else
          m_cur_trp->set_lb(0);
        m_cur_trp = 0;
        return;
      }
      if (is_right_rotation())
        m_cur_trp->set_rt(next);
      else
        m_cur_trp->set_lb(next);
      if (old.is_right_rotation())
        old[0]->set_rt(0);
      else
        old[0]->set_lb(0);
    }

    void replace(X_trapezoid& tr)
    {
      
#ifndef CGAL_TD_DEBUG
      
      CGAL_precondition(m_cur_trp != NULL);
      CGAL_precondition(m_cur_trp->is_active());
      CGAL_warning(
          (!m_cur_trp->is_on_left_boundary() &&
            traits->equal_curve_end_2_object()
                  (m_cur_trp->left()->curve_end(),m_fixed)) ||
          (!m_cur_trp->is_on_right_boundary() &&
            traits->equal_curve_end_2_object()
                  (m_cur_trp->right()->curve_end(),m_fixed)));
      
#else
      
      CGAL_precondition(m_cur_trp != NULL);
      CGAL_precondition(m_cur_trp->is_active());
      CGAL_precondition(
          (!m_cur_trp->is_on_left_boundary() &&
            traits->equal_curve_end_2_object()
                  (m_cur_trp->left()->curve_end(),m_fixed)) ||
          (!m_cur_trp->is_on_right_boundary() &&
            traits->equal_curve_end_2_object()
                  (m_cur_trp->right()->curve_end(),m_fixed)));
      
#endif
      
      Around_point_circulator old = *this;
      old++;
      X_trapezoid* next = old.operator->();
      // handle 1-cycle and 2-cycles seperately
      if (m_cur_trp != next)
      {
      }
      // 2-cycle
      else if (*this != old)
      {
        next = m_cur_trp;
      }
      // 1-cycle
      else
      {
        m_cur_trp = &tr;
        if (is_right_rotation())
          m_cur_trp->set_rt(m_cur_trp);
        else
          m_cur_trp->set_lb(m_cur_trp);
        return;
      }
      if (!tr.is_on_right_boundary() && 
          traits->equal_curve_end_2_object()
                      (tr.right()->curve_end(),m_fixed))
      {
        tr.set_rt(next);
      }
      else
        tr.set_lb(next);
    
      if (is_right_rotation())
        m_cur_trp->set_rt(&tr);
    else
        m_cur_trp->set_lb(&tr);
  }
  
    bool is_right_rotation() const
  {
      return !m_cur_trp->is_on_right_boundary() &&
              traits->equal_curve_end_2_object()
                        (m_cur_trp->right()->curve_end(),m_fixed);
    }
    
    const Curve_end& fixed() const
  {
      return m_fixed;
  }
  };
  

  struct Before_split_data
  {
    X_monotone_curve_2 m_cv_before_split;
    X_trapezoid* m_p_old_t;
    X_trapezoid* m_p_t1;
    X_trapezoid* m_p_t2;
    In_face_iterator* m_p_btm_it;
    In_face_iterator* m_p_mid_it;
    In_face_iterator* m_p_top_it;
    
  };
    
        
        
        
  //////////////////////////////////////////////
  //Trapezoidal_decomposition_2 member functions:
  //////////////////////////////////////////////
    
#ifndef CGAL_TD_DEBUG
    
protected:
    
#else
    
public:
    
#endif
    
    //returns true if the given edge is on the right side of the given point
    bool is_edge_to_right(Halfedge_const_handle he, const Point& p) const
    {
      //MICHAL: add precondition that the min or max of he are equal to p
      return traits->equal_curve_end_2_object()(Curve_end(he,ARR_MIN_END), p);
    }
    
    //returns true if the given edge is on the right side of the given edge end
    bool is_edge_to_right(Halfedge_const_handle he, const Curve_end& ce) const
    {
      //MICHAL: add precondition that the min or max of he are equal to p
    
      //if the curve end ce is on the right boundary - return false;
      if (traits->parameter_space_in_x_2_object()
           (ce.cv(), ce.ce()) == ARR_RIGHT_BOUNDARY)
      {
        return false;
      }
    
      //else
      return traits->equal_curve_end_2_object()(Curve_end(he,ARR_MIN_END), ce);
    }
    
    //returns true if the given curve is on the right side of the given point
    bool is_curve_to_right(const X_monotone_curve_2& cv, const Point& p) const
    {
      //MICHAL: add precondition that the min or max of cv are equal to p
      return traits->equal_curve_end_2_object()(Curve_end(cv,ARR_MIN_END), p);
    }
    
    bool is_end_point_left_low(const Point& p1, const Point& p2) const
    {
      return (traits->compare_xy_2_object()(p1, p2) == SMALLER);
    }
    
    bool is_end_point_left_low(const Point& p, const Curve_end& ce) const
    {
      return (traits->compare_curve_end_xy_2_object()(p, ce) == SMALLER);
    }
    
    bool is_end_point_left_low(const Curve_end& ce1, const Curve_end& ce2) const
    {
      return (traits->compare_curve_end_xy_2_object()(ce1, ce2) == SMALLER);
    }
    
    bool is_end_point_right_top(const Point& p1, const Point& p2) const
    {
      return (traits->compare_xy_2_object()(p1, p2) == LARGER);
    }
    
    bool is_end_point_right_top(const Point& p, const Curve_end& ce) const
    {
      return (traits->compare_curve_end_xy_2_object()(p, ce) == LARGER);
    }
    
    bool is_end_point_right_top(const Curve_end& ce1, const Curve_end& ce2) const
    {
      return (traits->compare_curve_end_xy_2_object()(ce1, ce2) == LARGER);
    }
      
  /* input: Halfedge_const_handle,
     a dag node that coresponds the
     Halfedge's source degenerate trapezoid 
     output: trapezoid iterator
      
     Description:
     the output (trapezoid iterator) is initialized with 
     the leftmost (non degenerate) trapezoids in the trapezoid 
     interval that corresponds to the input from either the top 
     side or the bottom side, depending on the up flag.
     preconditions:
     There exist non degenerate trapezoids between the roots of the input m_dag_root
     and the node corresponding to the target degenerate trapezoid
  */
  In_face_iterator follow_curve(const Dag_node& left_cv_end,
                                Halfedge_const_handle he,
                                Comparison_result up) const
  {
    return follow_curve(left_cv_end, he->curve(), up);
  }
      
  In_face_iterator follow_curve(const Dag_node& left_cv_end,
                                const X_monotone_curve_2& cv,
                                Comparison_result up) const
  {
    CGAL_assertion(traits != NULL);
    CGAL_precondition(traits->is_degenerate_point(*left_cv_end));
    if ( left_cv_end->is_active() )
    {
      CGAL_precondition(traits->equal_curve_end_2_object()
             (Curve_end(cv,ARR_MIN_END), left_cv_end->left()->curve_end()));
    }
    else 
    { //if not active
      if ( left_cv_end->is_on_boundaries())
        CGAL_precondition(traits->equal_curve_end_2_object()
             (Curve_end(cv,ARR_MIN_END), left_cv_end->curve_end_for_boundary_rem_vtx()));
      else //interior point
        CGAL_precondition(traits->equal_curve_end_2_object()
             (Curve_end(cv,ARR_MIN_END), left_cv_end->point_for_inner_rem_vtx()));
    }
    
    //find the node of the curve's leftmost trapezoid 
    Dag_node cv_leftmost_node(left_cv_end.right_child());
    if (left_cv_end->is_active())
    {
      Curve_end ce( left_cv_end->left()->curve_end());
      search_using_dag_with_cv(cv_leftmost_node, traits, ce, &cv, up);
    }
    else
    {
      Curve_end ce( left_cv_end->curve_end_for_rem_vtx());
      search_using_dag_with_cv(cv_leftmost_node, traits, ce, &cv, up);
    }
    //return a trapezoid iterator that starts from this trapezoid
    //  and continues according to the curve cv
    return In_face_iterator(traits, cv, cv_leftmost_node.operator->());
  }
    
  //-----------------------------------------------------------------------------
  // Description:
  //  Input: pointer to left trapezoid, pointer to right trapezoid
  //  Output: true iff the merging took place
  //  If the two input trapezoids can be merged they are ,
  //  with one copy destroyed(the right one).
  // Preconditions:
  //  the right trapezoid is to the right of the left one
  bool merge_if_possible(X_trapezoid* left, X_trapezoid* right)
  {
    if (left && right &&
        traits->is_trpz_top_equal(*left,*right) &&
        traits->is_trpz_bottom_equal(*left,*right) &&
        traits->equal_curve_end_2_object()
            (left->right()->curve_end(), right->left()->curve_end()))
    {
      left->merge_trapezoid(*right);
      //set the depth to be the max of the two merged nodes
      left->dag_node()->depth() = std::max ( left->dag_node()->depth(),
                                             right->dag_node()->depth());
#ifdef CGAL_TD_DEBUG
      CGAL_assertion(
            left->is_on_right_boundary() == right->is_on_right_boundary());
#endif
    
      return true;
    }
    return false;
  }
    
  //-----------------------------------------------------------------------------
  // Description:
  //  splits the trapezoid with vertical line through v 
  //  assuming that he_bottom_ray_shoot & he_top_ray_shoot are in the
  //  desired direction, such that v is their source
  // Precondition:
  //  The trapezoid is active and contains v in its closure
  //
  Dag_node& split_trapezoid_by_vertex(Dag_node& tt,
                                     Vertex_const_handle v,
                                     Halfedge_const_handle he_bottom_ray_shoot,
                                     Halfedge_const_handle he_top_ray_shoot);

  //-----------------------------------------------------------------------------
  // Description:
  //  the opposite operation for spliting the trapezoid with 
  //  vertical line through ce 
  // Precondition:
  //  The root trapezoid is degenerate point (ce) and is active 
  void undo_split_trapezoid_by_vertex(Dag_node& tr_node, const Curve_end& ce);
 
  //-----------------------------------------------------------------------------
  // Description:
  //  splits the trapezoid that corresponds to the root of the
  //  trapezoidal tree with an input halfedge he
  // Precondition:
  //  The root trapezoid is active
  //  The root trapezoid is devided by he or is equal to it and is vertical.
  Dag_node&
  split_trapezoid_by_halfedge(Dag_node& tt, X_trapezoid*& prev,
                           X_trapezoid*& prev_bottom,
                           X_trapezoid*& prev_top, Halfedge_const_handle he);
    
#if 0 
  /* replace X_curve-point adjacency in the data structure with
     a new one
     precondition:
     the X_curve represented by t is top-right
     relative to the point represented by sep
     if and only if top=true
  */
  //MICHAL: I am not sure when is this method called - maybe should be removed?
  void replace_curve_at_point_using_geometry(X_trapezoid& cv_tr, 
                                             const X_trapezoid& sep,
                                             bool cv_top_right=true);
#endif //if 0

  //-----------------------------------------------------------------------------
  // Description:
  //  replace halfedge-vertex adjacency in the data structure with a new one
  // precondition:
  //  the halfedge represented by he_tr is top-right
  //  relative to the vertex represented by sep if and only if top=true
  void set_neighbours_after_merge_halfedge_update (X_trapezoid& he_tr, 
                                              const X_trapezoid& sep,
                                              const X_monotone_curve_2& mrg_cv,
                                              bool he_top_right=true);
  


  //-----------------------------------------------------------------------------
  // Description:
  //  replace halfedge-vertex adjacency in the data structure with a new one
  // precondition:
  //  the halfedge represented by he_tr is top-right
  //  relative to the vertex represented by sep if and only if top=true
  void set_neighbours_after_split_halfedge_update(X_trapezoid& he_tr, 
                                                const X_trapezoid& sep,
                                                Halfedge_const_handle he1, 
                                                Halfedge_const_handle he2, 
                                                bool he_top_right=true);
  

  
  void set_neighbours_after_halfedge_insertion (X_trapezoid& he_tr,
                                                X_trapezoid& v_tr);
  

  //-----------------------------------------------------------------------------
  // Description:
  //  Update top(),bottom() for trapezoid
  //  Update rt,lb
  // remarks:
  //  The point degenerate trapezoid representing a point (edge_end) ee holds as its top and
  //  bottom curves
  //  the output for a vertical ray shoot queries immidiately below the point
  //  toward up and
  //  immediately above the point toward down respectively.
  //optimization:
  //  Each degenerate X_curve trapezoid emanating from the point p holds a pointer
  //  to the next
  //  trapezoid in a clockwise sweep around ee(possibly to itself).
  //  This pointer is stored in rt or lb depending on the trapezoid is top right
  //  or bottom left of ee.
  //  For the trapezoid representing ee, rt and lb hold the previous X_curve
  //  degenerate trapezoid
  //  in a clockwise sweep to the first top right and bottom left respectively.
  void remove_halfedge_at_vertex_using_geometry(const X_trapezoid& he_tr,
                                                X_trapezoid& v_tr);
  
  //-----------------------------------------------------------------------------
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
  X_trapezoid&
  set_trp_params_after_halfedge_insertion (Halfedge_const_handle he,
                                            const Curve_end& ce,
                                            X_trapezoid*& v_tr,
                                            const Locate_type&
                                              CGAL_precondition_code(lt));
  
  
  X_trapezoid&
  insert_curve_at_point_using_dag(Halfedge_const_handle he,
                                             Vertex_const_handle v,
                                             X_trapezoid*& tr,
                                             const Locate_type&
                                             CGAL_precondition_code(lt));
        

  void set_trp_params_after_halfedge_update (Halfedge_const_handle old_he,
                                             Halfedge_const_handle new_he,
                                             X_trapezoid&   sep);
  

  void set_trp_params_after_halfedge_update(const X_monotone_curve_2& old_cv,
                                             Halfedge_const_handle new_he,
                                             X_trapezoid&   sep);
  
  
  void set_trp_params_after_split_halfedge_update(Halfedge_const_handle new_he,
                                                X_trapezoid& v_tr,
                                                Halfedge_const_handle he1, 
                                                Halfedge_const_handle he2);
  
  //-----------------------------------------------------------------------------
  // Description:
  //  update geometric boundary(top and bottom) for trapezoids
  //  traveled along an iterator till end reached
  //  precondition:
  //  end==0 or end is on the path of the iterator
  // postcondition:
  //  end is pointer to the last trapezoid encountered,if any
  void set_trp_params_after_halfedge_update(In_face_iterator& it,
                                            Halfedge_const_handle old_he,
                                            Halfedge_const_handle new_he,
                                            Vertex_const_handle min_v,
                                            Vertex_const_handle max_v,
                                            X_trapezoid*& end);
  

  //-----------------------------------------------------------------------------
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
  Locate_type search_using_dag (Dag_node& curr,
                                const Traits* traits,
                                const Point& p,
                                Halfedge_const_handle* p_he,
                                Comparison_result up = EQUAL) const;
    
      
      
  //-----------------------------------------------------------------------------
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
  Locate_type search_using_dag (Dag_node& curr,
                                const Traits* traits,
                                const Curve_end& ce,
                                Halfedge_const_handle* p_he,
                                Comparison_result up = EQUAL) const;
      
      
      
  //-----------------------------------------------------------------------------
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
  Locate_type search_using_dag_with_cv (Dag_node& curr,
                                        const Traits* traits,
                                        const Curve_end& ce,
                                        const X_monotone_curve_2* p_cv,
                                        Comparison_result up = EQUAL) const;
      
      
      
    
  Dag_node container2dag (Nodes_map& ar, int left, int right,
                          int& num_of_new_nodes) const;
    
    
    
    
    
  /*==============================================
    Trapezoidal_decomposition_2 public member functions
    ==============================================*/
public:
    
  Trapezoidal_decomposition_2 (bool rebuild = true) 
    : depth_threshold(CGAL_TD_DEFAULT_DEPTH_THRESHOLD),
      size_threshold(CGAL_TD_DEFAULT_SIZE_THRESHOLD),
      m_arr(0), traits(0), m_largest_leaf_depth(0),
      m_number_of_curves(0), m_number_of_dag_nodes(1) 
  {
    init();
    set_needs_update(rebuild);
  }
    
  Trapezoidal_decomposition_2(const double& depth_th,const double& size_th)
    : depth_threshold(depth_th),size_threshold(size_th), m_arr(0), traits(0),
      m_largest_leaf_depth(0), m_number_of_curves(0), m_number_of_dag_nodes(1)
    {
    init();
    set_needs_update(rebuild);
  }
      
  Trapezoidal_decomposition_2(const Self& td) 
    : m_needs_update(td.m_needs_update),
      m_number_of_curves(td.m_number_of_curves),
      m_largest_leaf_depth(td.m_largest_leaf_depth),
      m_number_of_dag_nodes(td.m_number_of_dag_nodes),
      traits(td.traits),
      m_arr(td.m_arr),
      last_cv(NULL), prev_cv(NULL), 
      depth_threshold(td.depth_threshold),
      size_threshold(td.size_threshold)
  {
    X_trapezoid_ptr_map htr;
    /*! \todo allocate hash_map size according to content.
     * \todo change vector<> to in_place_list and pointer hash to trapezoidal
     * hash..
     */
    X_trapezoid_vector vtr;
    Td_active_trapezoid pr;
    int sz = X_trapezoid_filter(vtr, &td.dag_root());
    //! \todo Reduce the 3 iterations to 1 (or 2) iterator.
    // First iteration: filter out the active trapezoids.
    typename X_trapezoid_vector::const_iterator it;
    for (it = vtr.begin(); it != vtr.end(); ++it) 
    {
      Dag_node* ds_copy = new Dag_node(*it);
      const X_trapezoid* cur = &*it;
      X_trapezoid* tr_copy = &*(*ds_copy);
      tr_copy->set_dag_node(ds_copy);
      CGAL_assertion(&*(*tr_copy->dag_node()) == tr_copy);
      ds_copy->depth() = cur->dag_node()->depth();
      // We cheat a little with the depth.
      htr.insert(typename X_trapezoid_ptr_map::value_type(cur, tr_copy));
      // Second iteration: generate new copies of trapezoids and nodes.
    }
      
    for (it = vtr.begin(); it!=vtr.end(); ++it) 
    {
      const X_trapezoid* cur = &*it;
      X_trapezoid* tr_copy = htr.find(cur)->second;
      const Dag_node* child;
      CGAL_assertion(tr_copy);
      tr_copy->set_rt(cur->rt() ? 
                      htr.find(cur->rt())->second : NULL);
      tr_copy->set_rb(cur->rb() ?
                      htr.find(cur->rb())->second : NULL);
      tr_copy->set_lt(cur->lt() ? 
                      htr.find(cur->lt())->second : NULL);
      tr_copy->set_lb(cur->lb() ? 
                      htr.find(cur->lb())->second : NULL);

      if (cur->dag_node()->is_inner_node()) 
      {
        child = &cur->dag_node()->right_child();
        while (child && child->is_inner_node() && !pr(*(*child)))
          child = &child->left_child();
        tr_copy->dag_node()->set_right_child(*child);
        child = &cur->dag_node()->left_child();
        while (child && child->is_inner_node() && !pr(*(*child))) 
          child = &child->left_child();
        tr_copy->dag_node()->set_left_child(*child);
      }
      // Third iteration: generate links in-between trapezoids 
      //  and in-between nodes .
    }
    m_dag_root = htr.find(&*(*td.m_dag_root))->second->dag_node();
  }
      
  /*
    TODO: Should we add another constructor with non const argument that 
    rebuild the trapezoidal decomposition prior to copy construction?
  */
  virtual ~Trapezoidal_decomposition_2()
  {
      
#ifndef CGAL_TD_DEBUG
      
    CGAL_warning(m_dag_root != NULL);
    if (!m_dag_root) return;
      
#else
    
    CGAL_assertion(m_dag_root);
    
#endif
    
    delete m_dag_root;
    
    if (traits)
      delete traits;
  }
  
  //-----------------------------------------------------------------------------
  // Description:
  //  if Halfedge or twin already inserted the latter is returned.
  //  otherwise the left-low most edge-degenerate trapezoid that represents the
  //  input Halfedge is returned
  // Remark:
  //  Given an edge-degenerate trapezoid representing a Halfedge,
  //  all the other trapezoids representing the Halfedge can be extracted
  //  via moving continously to the left and right neighbours.
  X_trapezoid insert(Halfedge_const_handle he);
    
  
  //-----------------------------------------------------------------------------
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
    bool do_rebuild = set_needs_update(false);
    
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
        if (do_rebuild && needs_update()) 
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
    set_needs_update(do_rebuild);
  }
    

  // removal functions
    
  //-----------------------------------------------------------------------------
  // Description:
  //
  void remove(Halfedge_const_handle he);
    
#if 0
  //-----------------------------------------------------------------------------
  // Description:
  //
  template <class curve_iterator>
  void remove(curve_iterator begin, curve_iterator end)
  {
    if(begin == end)
      return;
    
    std::random_shuffle(begin,end);
    
    curve_iterator it=begin,next=it;
    while(it!=end) 
    {
      ++next;
      remove(*it);
      it=next;
    }
  }
#endif //if 0
    
  void clear()
    {
    delete m_dag_root;
    init();
    }
 
  void rebuild_if_necessary(Arr_all_sides_oblivious_tag)
    {
    }
    
  void rebuild_if_necessary(Arr_not_all_sides_oblivious_tag)
    {
    rebuild(); //MICHAL: added this to avoid location problem. SHould find a better solution since it is not efficient
    }    
    

  //-----------------------------------------------------------------------------
  // Description:
  //  returns the active trapezoid representing the input point.
  // Precondition:
  //  The trapezoidal tree is not empty
  // Postcondition:
  //  the input locate type is set to the type of the output trapezoid.
  // Remark:
  //  locate call may change the class
  X_trapezoid& locate(const Point& p,Locate_type &t) const
  {
    
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
    t = search_using_dag(curr,traits,p,0);
    
#ifdef CGAL_TD_DEBUG
    
    CGAL_postcondition(t == POINT || t == CURVE || t == TRAPEZOID ||
                       t == UNBOUNDED_TRAPEZOID);
    
#endif
    
#ifndef CGAL_NO_TRAPEZOIDAL_DECOMPOSITION_2_OPTIMIZATION
    
    locate_opt_push(curr.operator->());
    
#endif
    
    return *curr;
  }
    
  //-----------------------------------------------------------------------------
  // Description:
  //  returns the active trapezoid representing the input point.
  // Precondition:
  //  The trapezoidal tree is not empty
  // Postcondition:
  //  the input locate type is set to the type of the output trapezoid.
  // Remark:
  //  locate call may change the class
  X_trapezoid& locate(const Curve_end& ce, Locate_type& lt) const
  {
    
#ifdef CGAL_TD_DEBUG
    
    CGAL_assertion(traits);
    CGAL_assertion(m_dag_root);
    
#endif
    
    Dag_node curr = *m_dag_root; //MICHAL: is it ok to add &?
#ifdef CGAL_TD_DEBUG
    
    CGAL_precondition(!!curr);
    
#endif
    
    //the actual locate. curr is the DAG root, the traits, 
    //the end point to locate, 
    //and NULL as cv ptr - indicates point location 
    lt = search_using_dag (curr, traits, ce, NULL);
    
#ifdef CGAL_TD_DEBUG
    
    CGAL_postcondition(lt == POINT || lt == CURVE || lt == TRAPEZOID ||
                       lt == UNBOUNDED_TRAPEZOID);
    
#endif
    
#ifndef CGAL_NO_TRAPEZOIDAL_DECOMPOSITION_2_OPTIMIZATION
    
    locate_opt_push(curr.operator->());
    
#endif
    
    return *curr;
  }
  

  //-----------------------------------------------------------------------------
  // Description:
  // 
  // preconditions:
  //  p is not on an edge or a vertex.
  X_trapezoid& vertical_ray_shoot(const Point & p,Locate_type & t,
                                  const bool up_direction = true) const;
  
  
  void before_split_edge(const X_monotone_curve_2& cv,
                         const X_monotone_curve_2& cv1, 
                         const X_monotone_curve_2& cv2);


  //-----------------------------------------------------------------------------
  // Description:
  // Input:
  //  1 whole curves
  //  2 partial halfedge_handle-s
  // precondition:
  //  The two halfedges are valid
  //  The first input curve is the union of the two halfedges.
  //  The intersection of the latter is a point inside the 
  //  interior of the former.
  //  The latter are ordered from left-down to right-up
  // postcondition:
  //  The first input curve is broken into two halfedges 
  //  corresponding to the input.
  void split_edge(const X_monotone_curve_2& cv, Halfedge_const_handle he1, 
                  Halfedge_const_handle he2);

 
  void merge_edge(Halfedge_const_handle he1, Halfedge_const_handle he2,
                  const X_monotone_curve_2& cv);

  
  void after_merge_edge(Halfedge_const_handle merged_he, 
                        Halfedge_const_handle before_mrg_he)
  {
    //Precondition:
    // the merge uses the suspected halfedge before the arrangement merge
    CGAL_precondition(merged_he == before_mrg_he || 
                      merged_he == before_mrg_he->twin());
    //print_dag_addresses(*m_dag_root);
    //rebuild();//MICHAL: added this to avoid point-nodes in the DAG that hold no longer existing vertex
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
  
  void init_traits(const Traits* t)
  {
    traits = t;
    
#ifdef CGAL_TD_DEBUG
    
    CGAL_assertion(!!*m_dag_root);
    
#endif
    
  }
  
  /* update geometric boundary(top and bottom) for trapezoids
     traveled along an iterator till end reached
     precondition:
     end==0 or end is on the path of the iterator
     postcondition:
     end is pointer to the last trapezoid encountered,if any
  */
  /*------------------------------------------------------------------
    description:
    returns whether the Trapezoidal Dag is valid
  */
  
#ifdef CGAL_TD_DEBUG
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
  void debug() const
  {
    std::cout << "\nTrapezoidal_decomposition_2<Traits>::debug()\n" << *this
              << std::endl;
    X_trapezoid x;
    x.debug();
  }
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
    //std::cout << "\nrebuild!  " << m_number_of_curves << "\n\n";
#ifdef CGAL_TD_DEBUG
    std::cout << "\nrebuild()" << std::flush;
#endif
    
    Halfedge_container container;
    unsigned long rep = Halfedge_filter(container, &dag_root());
    clear();
    
    // initialize container to point to curves in X_trapezoid Tree
    if (rep>0)
    {
      bool o = set_needs_update(false);
      typename std::vector<Halfedge_const_handle>::iterator 
        it = container.begin(),
        it_end = container.end();
      while(it!=it_end) 
      {
        insert(*it);
        ++it;
      }
      set_needs_update(o);
    }
    
#ifdef CGAL_TD_DEBUG
    
    CGAL_assertion(is_valid());
    unsigned long sz = number_of_curves();
    if(sz!=rep)
    {
      std::cerr << "\nnumber_of_curves()=" << sz;
      std::cerr << "\nrepresentatives.size()=" << rep;
      CGAL_assertion(number_of_curves()==rep);
    }
    
#endif
    
    container.clear();
    return *this;
  }
  
  /* 
     Input:  
     a list of pointers to X_trapezoids and a X_trapezoid boolean predicate.
     Output: 
     void
     Postcondition:
     the list pointers correspond to all the X_trapezoids in the data
     structure for which the predicate value is true. 
  */
  
  template <class Container, class Predicate>
  void filter(Container& c, const Predicate& pr, 
              const Dag_node * ds) const
  {
    CGAL_assertion(ds);
    ds->filter(c,pr);
  }

  template <class Container, class Predicate>
  void filter(Container& c, const Predicate& pr) const
  {
    filter(c, pr, &dag_root());
  }

  template <class Container>
  unsigned long X_trapezoid_filter(Container& container, 
                                   const Dag_node* ds) const
  /* Return a container for all active trapeozoids */
  {
    ds->filter(container, Td_active_trapezoid());
    return container.size();
  }

  template <class Halfedge_container>
  unsigned long Halfedge_filter(Halfedge_container& container, 
                               const Dag_node* ds) const
  /* Return a container for all active curves */
  {
    unsigned long sz=number_of_curves();
    X_trapezoid_list representatives;
    ds->filter(representatives,
               Td_active_right_degenerate_curve_trapezoid(*traits));
    
#ifndef CGAL_TD_DEBUG
    
    CGAL_warning(sz==representatives.size());
    
#else
    
    unsigned long rep=representatives.size();
    if (sz!=rep)
    {
      std::cerr << "\nnumber_of_curves()=" << sz;
      std::cerr << "\nrepresentatives.size()=" << rep;
      CGAL_assertion(number_of_curves()==representatives.size());
    }
    
#endif
    
    if (sz>0)
    {
      typename X_trapezoid_list::iterator it = representatives.begin(),
        it_end = representatives.end();
      while(it!=it_end)
      {
        container.push_back(it->top()); //it represents an active trapezoid
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
  bool set_needs_update(bool u)
  {
    bool old = m_needs_update;
    m_needs_update = u;
    return old;
  }

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
#ifdef CGAL_TD_REBUILD_DEBUG
    std::cout << "\n|" << needs_update() << std::flush;
#endif
    if (needs_update()) 
    {
      rebuild(); 
      return true;
    }
    
    return false;
  }
  
  bool not_within_limits()
  {
    unsigned long num_of_cv = number_of_curves();
    bool cond1 = largest_leaf_depth() > 
                   (get_depth_threshold()*(std::log(double(num_of_cv+1)))); //MICHAL: should we add +1 to the largest depth?
    //bool cond1 = (m_dag_root->rec_depth()-1) >
    //              (get_depth_threshold()*(std::log(double(num_of_cv+1)))); //MICHAL: should we add +1 to the largest depth?
    bool cond2 = number_of_dag_nodes() > (get_size_threshold()*(num_of_cv + 1));
    //bool cond2 = m_dag_root->size() > (get_size_threshold()*(num_of_cv + 1));
    char c1 = cond1 ? 't' : 'f';
    char c2 = cond2 ? 't' : 'f';
    //std::cout << "\n" << c1 <<"," << c2 << " --> #" << num_of_cv << "\n"; 
    
    return cond1 || cond2;
  }

  /* returns a reference to the internal data structure */
  const Dag_node& dag_root() const {return *m_dag_root;}
  
  /* returns a reference to the internal data structure */
  const Traits& get_traits() const {return *traits;}
  
  /* returns a reference to the internal depth threshold constant */
  const double& get_depth_threshold() const
  {
    return depth_threshold;
  }
  /* returns a reference to the internal size threshold constant */
  const double& get_size_threshold() const
  {
    return size_threshold;
  }
  /* sets the internal depth threshold constant to the parameter and 
     returns its reference */
  const double& set_depth_threshold(const double& depth_th)
  {
    return depth_threshold=depth_th;
  }
  
  /* sets the internal size threshold constant to the parameter and 
     returns its reference */
  const double& set_size_threshold(const double& size_th)
  {
    return size_threshold=size_th;
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
#if 0
  unsigned long rec_check()
  {
    return m_dag_root->rec_check(largest_leaf_depth()+1);
  }
#endif //0

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
  bool m_needs_update;
  unsigned long m_number_of_curves;
  const Traits* traits;
  Before_split_data m_before_split;
  const Arrangement_on_surface_2* m_arr;
  const Traits_adaptor_2* m_trts_adaptor;
  
private:
  
#ifndef CGAL_NO_TRAPEZOIDAL_DECOMPOSITION_2_OPTIMIZATION
  
  mutable X_trapezoid* last_cv;
  mutable X_trapezoid* prev_cv;
  
#endif
 
  unsigned long longest_query_path_length_rec( 
                      bool minus_inf, Dag_node& min_node, 
                      bool plus_inf, Dag_node& max_node,
                      Dag_node& node);

  unsigned char build_boundaries_flag(const Curve_end& ce)
  {
    unsigned char bndry_flag = CGAL_TD_INTERIOR;

    Arr_parameter_space x_prm_spc = 
         traits->parameter_space_in_x_2_object()(ce.cv(), ce.ce());
    Arr_parameter_space y_prm_spc = 
         traits->parameter_space_in_y_2_object()(ce.cv(), ce.ce());
    
    if (x_prm_spc != ARR_INTERIOR)
    {
      bndry_flag |= (x_prm_spc == ARR_LEFT_BOUNDARY) 
        ? CGAL_TD_ON_LEFT_BOUNDARY : CGAL_TD_ON_RIGHT_BOUNDARY;
    }
    else //if x_prm_spc == ARR_INTERIOR
    {
      if (y_prm_spc != ARR_INTERIOR)
      {
        bndry_flag |= (y_prm_spc == ARR_BOTTOM_BOUNDARY) 
          ? CGAL_TD_ON_BOTTOM_BOUNDARY : CGAL_TD_ON_TOP_BOUNDARY;
      }
    }
    return bndry_flag;
  }
    
    
  void init()
  {
    // traits may be initialized later
    m_dag_root = new Dag_node(X_trapezoid());
    (*m_dag_root)->set_dag_node(m_dag_root);
    
#ifdef CGAL_TD_DEBUG
    
    CGAL_warning(!!*m_dag_root);
    
#endif
    
    m_number_of_curves = 0;
    m_largest_leaf_depth = 0; 
    m_number_of_dag_nodes = 1; //the root is the only node in the DAG
    
#ifndef CGAL_NO_TRAPEZOIDAL_DECOMPOSITION_2_OPTIMIZATION
    
    locate_opt_empty();
    
#endif
    
  }
  
#ifndef CGAL_NO_TRAPEZOIDAL_DECOMPOSITION_2_OPTIMIZATION
  
  void locate_opt_push(X_trapezoid* cv_tr) const
  {
    prev_cv=last_cv;
    last_cv=cv_tr;
  }
  void locate_opt_empty() const
  {
    last_cv=prev_cv=0;
  }
  bool locate_opt_swap(X_trapezoid*& cv_tr) const
  {
    cv_tr=last_cv;
    last_cv=prev_cv;
    prev_cv=cv_tr;
    return (cv_tr!=0);
  }
  void locate_optimization(const Curve_end& ce,X_trapezoid*& tr,
                            Locate_type& lt) const
  {
    // optimization
    if ((locate_opt_swap(tr) && tr->is_active() &&
         ((traits->is_degenerate_point(*tr) &&
           traits->equal_curve_end_2_object()(tr->left()->curve_end(),ce)) ||
          (!traits->is_degenerate(*tr) && traits->is_inside(*tr,ce)))) ||
        (locate_opt_swap(tr) && tr->is_active() &&
         ((traits->is_degenerate_point(*tr) &&
           traits->equal_curve_end_2_object()(tr->left()->curve_end(),ce)) ||
          (!traits->is_degenerate(*tr) && traits->is_inside(*tr,ce)))))
    {
      if (traits->is_degenerate_point(*tr)) 
        lt=POINT;
    else
        lt=tr->is_on_boundaries()? UNBOUNDED_TRAPEZOID : TRAPEZOID;
    }
    else
      tr=&locate(ce,lt);
  }
  
#endif

 
  void print_cv_data(const X_monotone_curve_2& cv) const
  {
    std::cout << "min end: " << std::endl;
    
    print_ce_data(cv, ARR_MIN_END);
    
    std::cout << std::endl << "max end: " << std::endl;
    
    print_ce_data(cv, ARR_MAX_END);
    
    std::cout << std::endl << std::endl ;
  }

  void print_ce_data(const X_monotone_curve_2& cv, Arr_curve_end ce) const
  {
    Arr_parameter_space ps_x = traits->parameter_space_in_x_2_object()(cv, ce);
    Arr_parameter_space ps_y = traits->parameter_space_in_y_2_object()(cv, ce);
    
    if (ps_x == ARR_INTERIOR && ps_y == ARR_INTERIOR)
    {
      if (ce == ARR_MIN_END)
        std::cout << "x: " << CGAL::to_double(traits->construct_min_vertex_2_object()(cv).x())
         << ", y: " << CGAL::to_double(traits->construct_min_vertex_2_object()(cv).y()) << std::endl;
      else
        std::cout  << "x: " << CGAL::to_double(traits->construct_max_vertex_2_object()(cv).x())
         << ", y: " << CGAL::to_double(traits->construct_max_vertex_2_object()(cv).y()) << std::endl;
    }
    else if (ps_x == ARR_INTERIOR && ps_y != ARR_INTERIOR)
    {
      std::cout << " vertical asymptote, " ;
      if (ps_y == ARR_TOP_BOUNDARY)
        std::cout << " y -> +oo " << std::endl;
      else
        std::cout << " y -> -oo " << std::endl;
    }
    else if (ps_x != ARR_INTERIOR && ps_y == ARR_INTERIOR)
    {
      std::cout << " horizontal asymptote, " ;
      if (ps_x == ARR_RIGHT_BOUNDARY)
        std::cout << " x -> +oo " << std::endl;
      else
        std::cout << " x -> -oo " << std::endl;
    }
    else //both are not interior
    {
      if (ps_x == ARR_RIGHT_BOUNDARY)
        std::cout << " x -> +oo " ;
      else
        std::cout << " x -> -oo " ;
      if (ps_y == ARR_TOP_BOUNDARY)
        std::cout << " , y -> +oo " << std::endl;
      else
        std::cout << " , y -> -oo " << std::endl;
    
    }
  }

  void print_dag_addresses(const Dag_node& curr)
  {
    
    std::cout << "----------------- DAG ----------------" <<std::endl
              << "--------------------------------------" <<std::endl;
    
    print_dag_addresses_rec(curr, 0);
    std::cout << "----------------- END OF DAG ----------------" <<std::endl
              << "---------------------------------------------" <<std::endl;
    
  }
  void print_dag_addresses_rec(const Dag_node& curr ,int level)
  {
    std::cout << "------ level " << level << ", depth " << curr.depth() << " ------\n";
    std::cout << " (void *)curr : " << (void *)(&curr) << std::endl;
    std::cout << "      (void *)curr->TRPZ : " << (void *)(curr.operator->()) << std::endl;
    //curr is the current pointer to node in the data structure
    if (traits->is_degenerate_point(*curr))
    { // if the trapezoid (curr) represents a point
      const Curve_end left_ce(curr->is_active()? 
        curr->left()->curve_end() : curr->curve_end_for_rem_vtx());
      std::cout << " POINT : " ;
      print_ce_data(left_ce.cv(), left_ce.ce());
      std::cout << "          (void *)left_child: " << (void *)(&(curr.left_child())) << std::endl;
      std::cout << "          (void *)right_child: " << (void *)(&(curr.right_child())) << std::endl;
      print_dag_addresses_rec(curr.left_child(), level+1);
      print_dag_addresses_rec(curr.right_child(), level+1);
      return;
    }
    if (traits->is_degenerate_curve(*curr))
    { 
      const X_monotone_curve_2* p_he_cv = 
      (curr->is_active()) ? &curr->top()->curve() : &curr->curve_for_rem_he();
    
      // if the trapezoid (curr) represents a curve, 
      //   so top() is a real Halfedge with a curve() if curr is active
      //   or curr holds the curve if curr is not active 
      std::cout << " CURVE : " ;
      print_cv_data(*p_he_cv);
      std::cout << "          (void *)left_child: " << (void *)(&(curr.left_child())) << std::endl;
      std::cout << "          (void *)right_child: " << (void *)(&(curr.right_child())) << std::endl;
      print_dag_addresses_rec(curr.left_child(), level+1);
      print_dag_addresses_rec(curr.right_child(), level+1);
      return;
    }
    else
    {
      // if is_degenerate() == 0, meaning: the trapezoid (curr)
      // is neither a point nor a curve , but a real trapezoid
      if (curr->is_active())
        std::cout << " TRAPEZOID \n";
      else //trapezoid is removed - may have a left child
      {
        std::cout << " REMOVED TRAPEZOID \n";
        if (!curr.left_child())
          return;
        std::cout << "          (void *)left_child: " << (void *)(&(curr.left_child())) << std::endl;
        print_dag_addresses_rec(curr.left_child(), level+1);
      }
    }
  }

protected:
  double depth_threshold,size_threshold;
};

} //namespace CGAL

#ifndef CGAL_TD_X_TRAPEZOID_H
#include <CGAL/Arr_point_location/Td_X_trapezoid.h>
#endif

#ifdef CGAL_TD_DEBUG
#ifndef CGAL_TRAPEZOIDAL_DECOMPOSITION_2_IOSTREAM_H
#include <CGAL/Arr_point_location/Trapezoidal_decomposition_2_iostream.h>
#endif
#endif

// The member-function definitions can be found under:
#include <CGAL/Arr_point_location/Trapezoidal_decomposition_2_impl.h>

#endif
