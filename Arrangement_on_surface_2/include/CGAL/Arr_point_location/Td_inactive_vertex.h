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
// $URL: svn+ssh://balasmic@scm.gforge.inria.fr/svn/cgal/branches/features/Arrangement_on_surface_2-RIC_pl_for_unbounded-balasmic/Arrangement_on_surface_2/include/CGAL/Arr_point_location/Td_inactive_vertex.h $
// $Id: Td_inactive_vertex.h 65793 2011-10-10 17:04:38Z balasmic $
// 
//
// Author(s)	 : Oren Nechushtan <theoren@math.tau.ac.il>
//               updated by: Michal Balas <balasmic@post.tau.ac.il>

#ifndef CGAL_TD_INACTIVE_VERTEX_H
#define CGAL_TD_INACTIVE_VERTEX_H

/*! \file
 * Defintion of the Td_inactive_vertex<Td_traits> class.
 */

#include <CGAL/Arr_point_location/Trapezoidal_decomposition_2.h>
#include <boost/variant.hpp>
#include <boost/shared_ptr.hpp>


#ifdef CGAL_TD_DEBUG
#define CGAL_TD_INLINE
#else
#define CGAL_TD_INLINE inline
#endif

namespace CGAL {

/*! \class
 * Implementation of a pseudo-trapezoid as two halfedges(top,bottom)
 * and two curve-ends(left,right).
 * Trapezoids are represented as two curve-ends called right and left and
 * two halfedges called top and bottom. The curve-ends (points) lie on the 
 * right and left boundaries of the trapezoid respectively and the halfedges 
 * bound the trapezoid from above and below.
 * There exist degenerate trapezoids called infinite trapezoid; this happens 
 * when one of the four sides is on the parameter space boundary.
 * Trapezoids are created as active and become inactive when Remove() member
 * function called.
 * Each trapezoid has at most four neighbouring trapezoids.
 * X_trapezoid structure can represent a real trapezoid, a Td-edge or an 
 * edge-end (end point).
 */
template <class Td_traits_>
class Td_inactive_vertex : public Handle
{
public:
  
  //type of trapezoid type
  enum Type 
  {
      TD_TRAPEZOID,
      TD_EDGE,
      TD_VERTEX
  };

  //type of traits class
  typedef Td_traits_                                   Traits;
  
  //type of point (Point_2)
  typedef typename Traits::Point                       Point;

  //type of X_monotone_curve_2
  typedef typename Traits::X_monotone_curve_2     X_monotone_curve_2;

  //type of Curve_end
  typedef typename Traits::Curve_end              Curve_end;

  //type of Halfedge_const_handle (trapezoid edge)
  typedef typename Traits::Halfedge_const_handle  Halfedge_const_handle;
  
  //type of Vertex_const_handle (trapezoid vertex)
  typedef typename Traits::Vertex_const_handle    Vertex_const_handle;

  //type of Td_inactive_vertex (Self)
  typedef typename Traits::Td_inactive_vertex            Self;
  
  typedef typename Traits::Td_map_item            Td_map_item;

  //type of Trapezoid parameter space
  // Ninetuple which represents the Trapezoid:
  //  - for regular & edge trapezoids or active point trapezoids:
  //      left vertex, right vertex, bottom halfedge, top halfedge
  //  - for removed point trapezoids:
  //      point or X_monotone_curve_2+ cv end
  //  type flag + on boundaries flags,
  //  left-bottom neighbor trapezoid, left-top neighbor trapezoid,
  //  right-bottom neighbor trapezoid, right-top neighbor trapezoid
  typedef Td_ninetuple<Point, 
                       Vertex_const_handle,
                       Halfedge_const_handle,
                       Halfedge_const_handle, 
                       unsigned char,
                       boost::optional<Td_map_item>, boost::optional<Td_map_item>,
                       boost::optional<Td_map_item>, boost::optional<Td_map_item> >            Trpz_parameter_space;
  
  //type of Trapezoidal decomposition
  typedef Trapezoidal_decomposition_2<Traits>          TD;
  
  //type of Around point circulator
  typedef typename TD::Around_point_circulator         Around_point_circulator;
  
  //type of In face iterator
  typedef typename TD::In_face_iterator                In_face_iterator;

  //type of Trapezoidal map search structure
  typedef typename TD::Dag_node                 Dag_node;


  //friend class declarations:

  friend class Trapezoidal_decomposition_2<Traits>;
  
#ifdef CGAL_PM_FRIEND_CLASS
#if defined(__SUNPRO_CC) || defined(__PGI) || defined(__INTEL_COMPILER)
  friend class Trapezoidal_decomposition_2<Traits>::Around_point_circulator;
  friend class Trapezoidal_decomposition_2<Traits>::In_face_iterator;
#elif defined(__GNUC__)

#if ((__GNUC__ < 3) || ((__GNUC__ == 3) && (__GNUC_MINOR__ <= 2)))
  friend typename Trapezoidal_decomposition_2<Traits>::Around_point_circulator;
  friend typename Trapezoidal_decomposition_2<Traits>::In_face_iterator;
#else
  friend class Trapezoidal_decomposition_2<Traits>::Around_point_circulator;
  friend class Trapezoidal_decomposition_2<Traits>::In_face_iterator;
#endif
  
#else
  friend class Around_point_circulator;
  friend class In_face_iterator;
#endif
#endif
  
  
  
 private:
  
  Trpz_parameter_space* ptr() const { return (Trpz_parameter_space*)(PTR);  }
	
	
#ifndef CGAL_TD_DEBUG
#ifdef CGAL_PM_FRIEND_CLASS
 protected:
#else
 public: // workaround
#endif
#else //CGAL_TD_DEBUG
 public:
#endif //CGAL_TD_DEBUG
	
  Dag_node* m_dag_node; //pointer to the search structure (DAG) node
	
  /*! Initialize the trapezoid's neighbours. */
  inline void init_neighbours(boost::optional<Td_map_item> lb = boost::none, boost::optional<Td_map_item> lt = boost::none,
                              boost::optional<Td_map_item> rb = boost::none, boost::optional<Td_map_item> rt = boost::none)
  {
    set_lb(lb);
    set_lt(lt);
    set_rb(rb);
    set_rt(rt);
  }

  /*! Set the DAG node. */
  inline void set_dag_node(Dag_node* p) 
  {
    m_dag_node = p;
  
#ifdef CGAL_TD_DEBUG
  
    CGAL_assertion(!p || **p == *this);
  
#endif	
	
  }
  
  ///*! Set the trapezoid's left (Vertex_const_handle). */
  //CGAL_TD_INLINE void set_left(Vertex_const_handle v) 
  //{
  //  CGAL_precondition(is_active());
  //  ptr()->e0 = v;
  //}
  
  inline void set_point(Vertex_const_handle v_before_rem)
  {
    ptr()->e0 = v_before_rem->point();
  }

  ///*! Set the trapezoid's right (Vertex_const_handle). */
  //CGAL_TD_INLINE void set_right(Vertex_const_handle v) 
  //{
  //  CGAL_precondition(is_active());
  //  ptr()->e1 = v;
  //}
  //
  ///*! Set the trapezoid's bottom (Halfedge_const_handle). */
  //CGAL_TD_INLINE void set_bottom(Halfedge_const_handle he) 
  //{
  //  CGAL_precondition(is_active());
  //  if (!is_on_bottom_boundary() &&
  //      bottom_unsafe()->direction() != he->direction())
  //  {
  //    ptr()->e2 = he->twin();
  //  }
  //  else
  //  {
  //    ptr()->e2 = he;
  //  }
  //}
  //
  ///*! Set the trapezoid's top (Halfedge_const_handle). */
  //CGAL_TD_INLINE void set_top(Halfedge_const_handle he) 
  //{
  //  CGAL_precondition(is_active());
  //  if (!is_on_top_boundary() &&
  //      top_unsafe()->direction() != he->direction())
  //  {
  //    ptr()->e3 = he->twin();
  //  }
  //  else
  //  {
  //    ptr()->e3 = he;
  //  }
  //}
  
  
  
  /*! Set left bottom neighbour. */
  inline void set_lb(boost::optional<Td_map_item> lb) { ptr()->e5 = lb; }
  
  /*! Set left top neighbour. */
  inline void set_lt(boost::optional<Td_map_item> lt) { ptr()->e6 = lt; }
  
  /*! Set right bottom neighbour. */
  inline void set_rb(boost::optional<Td_map_item> rb) { ptr()->e7 = rb; }
  
  /*! Set right top neighbour. */
  inline void set_rt(boost::optional<Td_map_item> rt) { ptr()->e8 = rt; }

 public:
  
  /// \name Constructors.
  //@{

  ///*! Default constructor. */
  //Td_inactive_vertex()
  //{
  //  //define the initial trapezoid: left, right, btm, top are at infinity.
  //  // its type is TD_TRAPEZOID ,it is on all boundaries, and has no neighbours
  //  PTR = new Trpz_parameter_space
  //    (Traits::vtx_at_left_infinity(),
  //     Traits::vtx_at_right_infinity(),
  //     Traits::he_at_bottom_infinity(),
  //     Traits::he_at_top_infinity(),
  //     CGAL_TD_TRAPEZOID | CGAL_TD_ON_ALL_BOUNDARIES ,
  //     0, 0, 0, 0);

  //  m_dag_node = 0;
  //}
  
  /*! Constructor given Vertex & Halfedge handles. */
  Td_inactive_vertex (Vertex_const_handle v_before_rem,
                      boost::optional<Td_map_item> lb = boost::none, boost::optional<Td_map_item> lt = boost::none,
                      boost::optional<Td_map_item> rb = boost::none, boost::optional<Td_map_item> rt = boost::none,
                      Dag_node* node = 0)
    //Vertex_const_handle l, Vertex_const_handle r,
    //              Halfedge_const_handle b, Halfedge_const_handle t,
    //              Type tp = TD_TRAPEZOID,
    //              unsigned char boundness_flag = CGAL_TD_INTERIOR,
    //              Self* lb = 0, Self* lt = 0,
    //              Self* rb = 0, Self* rt = 0,
    //              Dag_node* node = 0)
  {
    
    //build the type flag
    unsigned char type_flag = 0;
    type_flag |= CGAL_TD_VERTEX;

    PTR = new Trpz_parameter_space
      (v_before_rem->point(),
        Traits::empty_vtx_handle(), 
       Traits::empty_he_handle(), 
       Traits::empty_he_handle(), 
       type_flag | CGAL_TD_INTERIOR, lb, lt, rb, rt);
    m_dag_node = node;
  }
  
  
  /*! Copy constructor. */
  Td_inactive_vertex (const Self& tr) : Handle(tr)
    {
    m_dag_node = tr.m_dag_node;
    }
  
  //@}
  
  /// \name Operator overloading.
  //@{

  /*! Assignment operator. 
  *   operator= should not copy m_dag_node (or otherwise update 
  *     Dag_node::replace)
    */
  CGAL_TD_INLINE Self& operator= (const Self& t2)
      {
	Handle::operator=(t2);
	return *this;
      }

  /*! Operator==. */
  CGAL_TD_INLINE bool operator== (const Self& t2) const
  {
      return CGAL::identical(*this,t2);
  }

  /*! Operator!=. */
  CGAL_TD_INLINE bool operator!= (const Self& t2) const
  {
    return !(operator==(t2));
  }

  //@}


  /// \name Access methods.
  //@{

  CGAL_TD_INLINE Self& self() 
    {
      return *this;
    }
  
  CGAL_TD_INLINE const Self& self() const 
    {
      return *this;
    }

  /*! Access the trapezoid id (PTR). */
  CGAL_TD_INLINE unsigned long id() const
    {
      return (unsigned long) PTR;
    }

  ///*! Access trapezoid left. */
  //CGAL_TD_INLINE Vertex_const_handle left_unsafe() const
  //  {
  //  CGAL_precondition(is_active()); 
  //  CGAL_assertion(boost::get<Vertex_const_handle>(&(ptr()->e0)) != NULL);
  //  return boost::get<Vertex_const_handle>(ptr()->e0);
  //  }

  ///*! Access trapezoid left. 
  //*   filters out the infinite case which returns predefined dummy values
  //*/
  //CGAL_TD_INLINE Vertex_const_handle left() const
  //  {
  //  CGAL_precondition(is_active()); 
  //  if (is_on_left_boundary() && is_on_bottom_boundary()
  //      && is_on_top_boundary())
  //  {
  //    return Traits::vtx_at_left_infinity();
  //  }
  //  //else
  //  return left_unsafe();
  //  }

  ///*! Access trapezoid right. */
  //CGAL_TD_INLINE Vertex_const_handle right_unsafe() const
  //  {
  //  CGAL_precondition(is_active()); 
  //  CGAL_assertion(boost::get<Vertex_const_handle>(&(ptr()->e1)) != NULL);
  //  return boost::get<Vertex_const_handle>(ptr()->e1);
  //  }

  ///*! Access trapezoid right. 
  //*   filters out the infinite case which returns predefined dummy values
  //*/
  //CGAL_TD_INLINE Vertex_const_handle right () const
  //{
  //  CGAL_precondition(is_active()); 
  //  if (is_on_right_boundary() && is_on_bottom_boundary()
  //      && is_on_top_boundary())
  //  {
  //    return Traits::vtx_at_right_infinity();
  //  }
  //  //else
  //  return right_unsafe();
  //}
  //
  ///*! Access trapezoid bottom. */
  //CGAL_TD_INLINE Halfedge_const_handle bottom_unsafe () const
  //{
  //  CGAL_precondition(is_active()); 
  //  CGAL_assertion(boost::get<Halfedge_const_handle>(&(ptr()->e2)) != NULL);
  //  return boost::get<Halfedge_const_handle>(ptr()->e2);
  //  }
  //
  ///*! Access trapezoid bottom. 
  //*   filters out the infinite case which returns predefined dummy values
  //*/
  //CGAL_TD_INLINE Halfedge_const_handle bottom () const
  //  {
  //  CGAL_precondition(is_active()); 
  //  return !is_on_bottom_boundary() ?  
  //          bottom_unsafe() : Traits::he_at_bottom_infinity();
  //  }

  ///*! Access trapezoid top. */
  //CGAL_TD_INLINE Halfedge_const_handle top_unsafe () const
  //  {
  //  CGAL_precondition(is_active()); 
  //  return ptr()->e3;
  //  }
  //
  ///*! Access trapezoid top. 
  //*   filters out the infinite case which returns predefined dummy values
  //*/
  //CGAL_TD_INLINE Halfedge_const_handle top () const
  //  {
  //  CGAL_precondition(is_active()); 
  //  return !is_on_top_boundary() ?	
  //          top_unsafe() : Traits::he_at_top_infinity();
  //  }

  inline Point& point() const  
  {
    return ptr()->e0;
    //CGAL_assertion(boost::get<Point>( &(ptr()->e0)) != NULL);
    //return boost::get<Point>( ptr()->e0 );
  }
  
  //CGAL_TD_INLINE std::pair<X_monotone_curve_2*,Arr_curve_end> curve_end_pair_for_boundary_rem_vtx() const  
  //  {
  //  CGAL_precondition(!is_active()); 
  //  CGAL_precondition(type() == TD_VERTEX);
  //  CGAL_precondition(is_on_boundaries());

  //  CGAL_assertion(boost::get<unsigned char>( &(ptr()->e1)) != NULL);
  //  CGAL_assertion(boost::get<boost::shared_ptr<X_monotone_curve_2> >(&(ptr()->e2)) != NULL);
  //  X_monotone_curve_2* cv_ptr = (boost::get<boost::shared_ptr<X_monotone_curve_2> >(ptr()->e2)).get();
  //  CGAL_assertion(cv_ptr != NULL);
  // 
  //  Arr_curve_end ce = 
  //    (boost::get<unsigned char>(ptr()->e1) == CGAL_TD_CV_MIN_END) ?
  //      ARR_MIN_END : ARR_MAX_END;
  //
  //  return std::make_pair(cv_ptr, ce);
  //}

  //CGAL_TD_INLINE Curve_end curve_end_for_boundary_rem_vtx() const  
  //{
  //  CGAL_precondition(!is_active()); 
  //  CGAL_precondition(type() == TD_VERTEX);
  //  CGAL_precondition(is_on_boundaries());

  //  CGAL_assertion(boost::get<unsigned char>( &(ptr()->e1)) != NULL);
  //  CGAL_assertion(boost::get<boost::shared_ptr<X_monotone_curve_2> >(&(ptr()->e2)) != NULL);
  //  X_monotone_curve_2* cv_ptr = (boost::get<boost::shared_ptr<X_monotone_curve_2> >(ptr()->e2)).get();
  //  CGAL_assertion(cv_ptr != NULL);
  // 
  //  Arr_curve_end ce = 
  //    (boost::get<unsigned char>(ptr()->e1) == CGAL_TD_CV_MIN_END) ?
  //      ARR_MIN_END : ARR_MAX_END;
  //
  //  return Curve_end(*cv_ptr, ce);
  //}

  //CGAL_TD_INLINE Curve_end curve_end_for_rem_vtx() const  
  //{
  //  CGAL_precondition(!is_active()); 
  //  CGAL_precondition(type() == TD_VERTEX);

  //  CGAL_assertion(boost::get<unsigned char>( &(ptr()->e1)) != NULL);
  //  CGAL_assertion(boost::get<boost::shared_ptr<X_monotone_curve_2> >(&(ptr()->e2)) != NULL);
  //  X_monotone_curve_2* cv_ptr = (boost::get<boost::shared_ptr<X_monotone_curve_2> >(ptr()->e2)).get();
  //  CGAL_assertion(cv_ptr != NULL);
  //  
  //  Arr_curve_end ce = 
  //    (boost::get<unsigned char>(ptr()->e1) == CGAL_TD_CV_MIN_END) ?
  //      ARR_MIN_END : ARR_MAX_END;
  //
  //  return Curve_end(*cv_ptr, ce);
  //}

  //CGAL_TD_INLINE X_monotone_curve_2& curve_for_rem_he() const  
  //{
  //  CGAL_precondition(!is_active() && type() == TD_EDGE); 
  //  
  //  CGAL_assertion(boost::get<boost::shared_ptr<X_monotone_curve_2> >(&(ptr()->e2)) != NULL);
  //  X_monotone_curve_2* cv_ptr = (boost::get<boost::shared_ptr<X_monotone_curve_2> >(ptr()->e2)).get();
  //  CGAL_assertion(cv_ptr != NULL);
  //  return *cv_ptr;
  //}

  ///*! Access trapezoid type. */
  //CGAL_TD_INLINE Type type() const 
  //{
  //  switch(ptr()->e4 & CGAL_TD_TYPE_MASK)
  //  {
  //  case CGAL_TD_TRAPEZOID:
  //    return TD_TRAPEZOID;
  //  case CGAL_TD_EDGE:
  //    return TD_EDGE;
  //  case CGAL_TD_VERTEX:
  //    return TD_VERTEX;
  //  default:
  //     CGAL_assertion(false);
  //     return TD_TRAPEZOID;
  //  } 
  //  }

  ///*! Access trapezoid type flag. */
  //CGAL_TD_INLINE unsigned char type_flag() const 
  //  {
  //  return (ptr()->e4 & CGAL_TD_TYPE_MASK);
  //  }
  //
  ///*! Access on boundaries flag. */
  //CGAL_TD_INLINE unsigned char on_boundaries_flag() const 
  //{
  //  return (ptr()->e4 & CGAL_TD_ON_ALL_BOUNDARIES);
  //}
  //
  ///*! Access is on left boundary. */
  //CGAL_TD_INLINE bool is_on_left_boundary() const 
  //  {			
  //  return (ptr()->e4 & CGAL_TD_ON_LEFT_BOUNDARY) != 0;
  //  }

  ///*! Access is on right boundary. */
  //CGAL_TD_INLINE bool is_on_right_boundary() const 
  //{
  //  return (ptr()->e4 & CGAL_TD_ON_RIGHT_BOUNDARY) != 0;
  //  }

  ///*! Access is on bottom boundary. */
  //CGAL_TD_INLINE bool is_on_bottom_boundary() const 
  //{
  //  return (ptr()->e4 & CGAL_TD_ON_BOTTOM_BOUNDARY) != 0;
  //  }

  ///*! Access is on top boundary. */
  //CGAL_TD_INLINE bool is_on_top_boundary() const 
  //{
  //  return (ptr()->e4 & CGAL_TD_ON_TOP_BOUNDARY) != 0;
  //}
  //
  ///*! Access is on at least one boundary. */
  //CGAL_TD_INLINE bool is_on_boundaries() const
  //{			
  //  return (ptr()->e4 & CGAL_TD_ON_ALL_BOUNDARIES) != 0;
  //}
  
  /*! Access left bottom neighbour. */
  boost::optional<Td_map_item> lb() const    { return ptr()->e5; }
  
  /*! Access left top neighbour. */
  boost::optional<Td_map_item> lt() const    { return ptr()->e6; }
  
  /*! Access right bottom neighbour. */
  boost::optional<Td_map_item> rb() const    { return ptr()->e7; }
  
  /*! Access right top neighbour. */
  boost::optional<Td_map_item> rt() const    { return ptr()->e8; }
  
  /*! Access DAG node. */
  Dag_node* dag_node() const            {return m_dag_node;}
  
  
  //@}
  
	

};

} //namespace CGAL

#endif
