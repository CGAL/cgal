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
//trapezoid
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://balasmic@scm.gforge.inria.fr/svn/cgal/branches/features/Arrangement_on_surface_2-RIC_pl_for_unbounded-balasmic/Arrangement_on_surface_2/include/CGAL/Arr_point_location/Td_active_trapezoid.h $
// $Id: Td_active_trapezoid.h 65793 2011-10-10 17:04:38Z balasmic $
// 
//
// Author(s)	 : Oren Nechushtan <theoren@math.tau.ac.il>
//               updated by: Michal Balas <balasmic@post.tau.ac.il>

#ifndef CGAL_TD_ACTIVE_TRAPEZOID_H
#define CGAL_TD_ACTIVE_TRAPEZOID_H

/*! \file
 * Defintion of the Td_active_trapezoid<Td_traits> class.
 */

#include <CGAL/Arr_point_location/Trapezoidal_decomposition_2.h>
#include <CGAL/tuple.h>

#include <boost/variant.hpp>
#include <boost/optional.hpp>


#ifdef CGAL_TD_DEBUG
#define CGAL_TD_INLINE
#else
#define CGAL_TD_INLINE inline
#endif

namespace CGAL {

/*! \class
 * Implementation of an active trapezoid as two halfedges(top,bottom)
 * and two vertices(left,right).
 * Trapezoids are represented as two vertices called right and left and
 * two halfedges called top and bottom. The curve-ends (points) lie on the 
 * right and left boundaries of the trapezoid respectively and the halfedges 
 * bound the trapezoid from above and below.
 * There exist degenerate trapezoids called infinite trapezoid; this happens 
 * when one of the four sides is on the parameter space boundary.
 * Each trapezoid has at most four neighbouring trapezoids.
 */
template <typename Td_traits_>
class Td_active_trapezoid : public Handle
{
public:
  
  //type of traits class
  typedef Td_traits_                                    Traits;
  
  //type of point (Point_2)
  typedef typename Traits::Point                        Point;

  //type of X_monotone_curve_2
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;

  //type of Curve_end
  typedef typename Traits::Curve_end                    Curve_end;

  //type of Halfedge_const_handle (trapezoid edge)
  typedef typename Traits::Halfedge_const_handle        Halfedge_const_handle;
  
  //type of Vertex_const_handle (trapezoid vertex)
  typedef typename Traits::Vertex_const_handle          Vertex_const_handle;

  //type of Td_active_trapezoid (Self)
  typedef typename Traits::Td_active_trapezoid          Self;

  typedef typename Traits::Td_map_item                  Td_map_item;
  
  //type of Trapezoidal decomposition
  typedef Trapezoidal_decomposition_2<Traits>           TD;
  
  //type of In face iterator
  typedef typename TD::In_face_iterator                 In_face_iterator;

  //type of Trapezoidal map search structure
  typedef typename TD::Dag_node                         Dag_node;

  //friend class declarations:

  friend class Trapezoidal_decomposition_2<Traits>;
  
#ifdef CGAL_PM_FRIEND_CLASS
#if defined(__SUNPRO_CC) || defined(__PGI) || defined(__INTEL_COMPILER)
  friend class Trapezoidal_decomposition_2<Traits>::In_face_iterator;
#elif defined(__GNUC__)

#if ((__GNUC__ < 3) || ((__GNUC__ == 3) && (__GNUC_MINOR__ <= 2)))
  friend typename Trapezoidal_decomposition_2<Traits>::In_face_iterator;
#else
  friend class Trapezoidal_decomposition_2<Traits>::In_face_iterator;
#endif
  
#else
  friend class In_face_iterator;
#endif
#endif
  
  protected:

   /*! \class
   * Inner class Data derived from Rep class
   */
  class Data : public Rep
  {
    friend class Td_active_trapezoid<Td_traits_>;

  public:
    //c'tors
    Data (Vertex_const_handle _left_v,   
          Vertex_const_handle _right_v,
          Halfedge_const_handle _bottom_he,
          Halfedge_const_handle _top_he,
          const Td_map_item& _lb,
          const Td_map_item& _lt,
          const Td_map_item& _rb,
          const Td_map_item& _rt,
          Dag_node* _p_node)
          : left_v(_left_v), right_v(_right_v),
            bottom_he(_bottom_he), top_he(_top_he),
            lb(_lb), lt(_lt), rb(_rb), rt(_rt), p_node(_p_node)
    { }
    
    ~Data() { }

  protected:
    Vertex_const_handle left_v; 
    Vertex_const_handle right_v;
    Halfedge_const_handle bottom_he;
    Halfedge_const_handle top_he;
    Td_map_item lb;
    Td_map_item lt;
    Td_map_item rb; 
    Td_map_item rt;
    Dag_node* p_node;
  };
  
 private:
  
  Data* ptr() const { return (Data*)(PTR); }
	
#ifndef CGAL_TD_DEBUG
#ifdef CGAL_PM_FRIEND_CLASS
protected:
#else
public: // workaround
#endif
#else //CGAL_TD_DEBUG
public:
#endif //CGAL_TD_DEBUG
	
  //Dag_node* m_dag_node; //pointer to the search structure (DAG) node
	
  /*! Initialize the trapezoid's neighbours. */
  inline void init_neighbours(boost::optional<Td_map_item&> lb,
                              boost::optional<Td_map_item&> lt,
                              boost::optional<Td_map_item&> rb,
                              boost::optional<Td_map_item&> rt)
  {
    set_lb((lb) ? *lb : Td_map_item(0));
    set_lt((lt) ? *lt : Td_map_item(0));
    set_rb((rb) ? *rb : Td_map_item(0));
    set_rt((rt) ? *rt : Td_map_item(0));
  }

  /*! Set the DAG node. */
  inline void set_dag_node(Dag_node* p) 
  {
    ptr()->p_node = p;
//    m_dag_node = p;
//  
//#ifdef CGAL_TD_DEBUG
//  
//    CGAL_assertion(!p || **p == *this);
//  
//#endif		
  }
  
  /*! Set the trapezoid's left (Vertex_const_handle). */
  inline void set_left(Vertex_const_handle v) { ptr()->left_v = v; }
  
  /*! Set the trapezoid's right (Vertex_const_handle). */
  inline void set_right(Vertex_const_handle v) { ptr()->right_v = v; }
  
  /*! Set the trapezoid's bottom (Halfedge_const_handle). */
  inline void set_bottom(Halfedge_const_handle he) 
  {
    if (!is_on_bottom_boundary() && (bottom()->direction() != he->direction()))
    {
      ptr()->bottom_he = he->twin(); 
    }
    else
    {
      ptr()->bottom_he = he; 
    }
  }
  
  /*! Set the trapezoid's top (Halfedge_const_handle). */
  inline void set_top(Halfedge_const_handle he) 
  {
    if (!is_on_top_boundary() && (top()->direction() != he->direction()))
    {
      ptr()->top_he = he->twin(); 
    }
    else
    {
      ptr()->top_he = he; 
    }
  }
 
  
  /*! Set left bottom neighbour. */
  inline void set_lb(const Td_map_item& lb) { ptr()->lb = lb; }
  
  /*! Set left top neighbour. */
  inline void set_lt(const Td_map_item& lt) { ptr()->lt = lt; }
  
  /*! Set right bottom neighbour. */
  inline void set_rb(const Td_map_item& rb) { ptr()->rb = rb; }
  
  /*! Set right top neighbour. */
  inline void set_rt(const Td_map_item& rt) { ptr()->rt = rt; }

 public:
  
  /// \name Constructors.
  //@{

  /*! Default constructor. */
  Td_active_trapezoid()
  {
    //define the initial trapezoid: left, right, btm, top are at infinity.
    // has no neighbours
    PTR = new Data(Traits::empty_vtx_handle(),
                   Traits::empty_vtx_handle(),
                   Traits::empty_he_handle(),
                   Traits::empty_he_handle(),
                   Td_map_item(0), Td_map_item(0) ,
                   Td_map_item(0), Td_map_item(0), NULL);
   
    //m_dag_node = 0;
  }

  /*! Constructor given Vertex & Halfedge handles. */ 
  Td_active_trapezoid (Vertex_const_handle l, Vertex_const_handle r,
                       Halfedge_const_handle b, Halfedge_const_handle t,
                       boost::optional<Td_map_item&> lb = boost::none, 
                       boost::optional<Td_map_item&> lt = boost::none,
                       boost::optional<Td_map_item&> rb = boost::none, 
                       boost::optional<Td_map_item&> rt = boost::none,
                       Dag_node* node = 0)
  {
    PTR = new Data(l, r, b, t, (lb) ? *lb : Td_map_item(0),
                   (lt) ? *lt : Td_map_item(0), (rb) ? *rb : Td_map_item(0),
                   (rt) ? *rt : Td_map_item(0), node);
    //m_dag_node = node;
  }
 
  
  /*! Copy constructor. */
  Td_active_trapezoid (const Self& tr) : Handle(tr)
  {
    //m_dag_node = tr.m_dag_node;
  }
  
  //@}
  
  /// \name Operator overloading.
  //@{

  /*! Assignment operator. 
  *   operator= should not copy m_dag_node (or otherwise update 
  *     Dag_node::replace)
    */
  inline Self& operator=(const Self& t2)
  {
    Handle::operator=(t2);
    return *this;
  }

  /*! Operator==. */
  inline bool operator==(const Self& t2) const
  {
    return (ptr() == t2.ptr());
  }

  /*! Operator!=. */
  inline bool operator!=(const Self& t2) const
  {
    return !(operator==(t2));
  }

  //@}


  /// \name Access methods.
  //@{

  inline Self& self() { return *this; }

  inline const Self& self() const { return *this; }

  /*! Access the trapezoid id (PTR). */
  inline unsigned long id() const { return (unsigned long) PTR; }

  /*! Access trapezoid left. 
  *   filters out the infinite case which returns predefined dummy values
  */
  inline Vertex_const_handle left() const { return ptr()->left_v; }

  /*! Access trapezoid right. 
  *   filters out the infinite case which returns predefined dummy values
  */
  inline Vertex_const_handle right() const { return ptr()->right_v; }
  
  /*! Access trapezoid bottom. 
  *   filters out the infinite case which returns predefined dummy values
  */
  inline Halfedge_const_handle bottom() const { return  ptr()->bottom_he; }
 
  /*! Access trapezoid top. 
  *   filters out the infinite case which returns predefined dummy values
  */
  inline Halfedge_const_handle top() const { return ptr()->top_he; }

  /*! Access is on left boundary. */
  inline bool is_on_left_boundary() const 
  { return (left() == Traits::empty_vtx_handle()); }

  /*! Access is on right boundary. */
  inline bool is_on_right_boundary() const 
  { return (right() == Traits::empty_vtx_handle()); }

  /*! Access is on bottom boundary. */
  inline bool is_on_bottom_boundary() const 
  { return (bottom() == Traits::empty_he_handle()); }

  /*! Access is on top boundary. */
  inline bool is_on_top_boundary() const 
  { return (top() == Traits::empty_he_handle()); }
  
  /*! Access is on at least one boundary. */
  inline bool is_on_boundaries() const
  {			
    return (is_on_left_boundary() || is_on_right_boundary() || 
            is_on_bottom_boundary() || is_on_top_boundary());
  }
  
  /*! Access left bottom neighbour. */
  Td_map_item& lb() const { return ptr()->lb; }
  
  /*! Access left top neighbour. */
  Td_map_item& lt() const { return ptr()->lt; }
  
  /*! Access right bottom neighbour. */
  Td_map_item& rb() const { return ptr()->rb; }
  
  /*! Access right top neighbour. */
  Td_map_item& rt() const { return ptr()->rt; }
  
  /*! Access DAG node. */
  Dag_node* dag_node() const {return ptr()->p_node; }
  
  void clear_neighbors()
  {
    set_lb(Td_map_item(0));
    set_lt(Td_map_item(0));
    set_rb(Td_map_item(0));
    set_rt(Td_map_item(0));
  }
  
  //@}
  
  
  // Merge this trapezoid with the input trapezoid.
  // Precondition:
  //   Both trapezoids are active and have the same bounding edges from
  //   above and below and the trapezoids are adjacent to one another
  //   with the first to the left.
  // Postcondition:
  //   This trapezoid is the union of the old this trapezoid and the input
  //   trapezoid.
  inline void merge_trapezoid(Self& right)
  {
    //precondition: the left trapezoid is not on the right boundary
    CGAL_assertion(!is_on_right_boundary());

    // bool on_right_boundary = right.is_on_right_boundary();

    ptr()->left_v = left();
    ptr()->right_v = right.right();
    ptr()->bottom_he = bottom();
    ptr()->top_he = top();
    ptr()->lb = lb();
    ptr()->lt = lt();
    ptr()->rb = right.rb(); 
    ptr()->rt = right.rt();
    
    Td_map_item item (*this);

    if (ptr()->rb.which() != 0)
    {
      Self tr(boost::get<Self>(rb()));     
      tr.set_lb(item); 
    }
    if (ptr()->rt.which() != 0)
    {
      Self tr(boost::get<Self>(rt()));     
      tr.set_lt(item);
    }
    CGAL_assertion(is_on_right_boundary() == right.is_on_right_boundary());
  }
};

} //namespace CGAL

#endif
