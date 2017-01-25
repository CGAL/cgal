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
// $URL: svn+ssh://balasmic@scm.gforge.inria.fr/svn/cgal/branches/features/Arrangement_on_surface_2-RIC_pl_for_unbounded-balasmic/Arrangement_on_surface_2/include/CGAL/Arr_point_location/Td_inactive_fictitious_vertex.h $
// $Id: Td_inactive_fictitious_vertex.h 65793 2011-10-10 17:04:38Z balasmic $
// 
//
// Author(s)	 : Oren Nechushtan <theoren@math.tau.ac.il>
//               updated by: Michal Balas <balasmic@post.tau.ac.il>

#ifndef CGAL_TD_INACTIVE_FICTITIOUS_VERTEX_H
#define CGAL_TD_INACTIVE_FICTITIOUS_VERTEX_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
 * Defintion of the Td_inactive_fictitious_vertex<Td_traits> class.
 */

#include <CGAL/Arr_point_location/Trapezoidal_decomposition_2.h>
#include <boost/variant.hpp>


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
class Td_inactive_fictitious_vertex : public Handle
{
public:
  
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

  //type of Halfedge_around_vertex_const_circulator
  typedef typename Traits::Halfedge_around_vertex_const_circulator  
    Halfedge_around_vertex_const_circulator;

  //type of Td_inactive_fictitious_vertex (Self)
  typedef typename Traits::Td_inactive_fictitious_vertex            Self;

  typedef typename Traits::Td_map_item            Td_map_item;

  //type of Trapezoidal decomposition
  typedef Trapezoidal_decomposition_2<Traits>          TD;
  
  //type of In face iterator
  typedef typename TD::In_face_iterator                In_face_iterator;

  //type of Trapezoidal map search structure
  typedef typename TD::Dag_node                 Dag_node;


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
  
     /*! \class
   * Inner class Data derived from Rep class
   */
  class Data : public Rep
  {
    friend class Td_inactive_fictitious_vertex<Td_traits_>;

  public:
    //c'tors
    Data (const X_monotone_curve_2& _cv,   
          Arr_curve_end _ce, 
          Dag_node* _p_node): cv(_cv),ce(_ce),p_node(_p_node)  
    { }
    
    ~Data() { }

  protected:
    X_monotone_curve_2 cv; 
    Arr_curve_end ce;
    Dag_node* p_node;
  };
  
 private:
  
  Data* ptr() const { return (Data*)(PTR);  }
	
	Curve_end vtx_to_ce(Vertex_const_handle v) const
  {
    //the circulator is of incoming halfedges
    Halfedge_around_vertex_const_circulator he = v->incident_halfedges(); 
    //if the vertex is associated with a point on the bounded coords,
    // we can take any incident halfedge. o/w if the vertex lies at infinity,
    //  it has 2 fictitious incident halfedges
    if (v->is_at_open_boundary() && he->source()->is_at_open_boundary()) ++he;
    if (v->is_at_open_boundary() && he->source()->is_at_open_boundary()) ++he;

    return Curve_end(he->curve(),
                     (he->direction() == ARR_RIGHT_TO_LEFT)? 
                      ARR_MIN_END : ARR_MAX_END);
  }

#ifndef CGAL_TD_DEBUG
#ifdef CGAL_PM_FRIEND_CLASS
 protected:
#else
 public: // workaround
#endif
#else //CGAL_TD_DEBUG
 public:
#endif //CGAL_TD_DEBUG
	
  /*! Set the DAG node. */
  inline void set_dag_node(Dag_node* p) 
  {
    ptr()->p_node = p;
  }
  
  inline void set_curve_end(Vertex_const_handle v_before_rem)
  {
    Curve_end v_ce(vtx_to_ce(v_before_rem));
    ptr()->cv(v_ce.cv());
    ptr()->ce(v_ce.ce());
  }
  
 

 public:
  
  /// \name Constructors.
  //@{

  /*! Constructor given Vertex & Halfedge handles. */
  Td_inactive_fictitious_vertex (Vertex_const_handle v_before_rem, Dag_node* node = NULL)
  {
    Curve_end v_ce(vtx_to_ce(v_before_rem));
   
    PTR = new Data( v_ce.cv(), v_ce.ce(), node);
    
  }
  
  /*! Copy constructor. */
  Td_inactive_fictitious_vertex (const Self& tr) : Handle(tr)
  {
   
  }
  
  //@}
  
  /// \name Operator overloading.
  //@{

  /*! Assignment operator. 
  *   operator= should not copy m_dag_node (or otherwise update 
  *     Dag_node::replace)
    */
  inline Self& operator= (const Self& t2)
  {
	  Handle::operator=(t2);
	  return *this;
  }

  /*! Operator==. */
  inline bool operator== (const Self& t2) const
  {
    return (ptr() == t2.ptr());
  }

  /*! Operator!=. */
  inline bool operator!= (const Self& t2) const
  {
    return !(operator==(t2));
  }

  //@}


  /// \name Access methods.
  //@{

  inline Self& self() 
  {
    return *this;
  }
  
  inline const Self& self() const 
  {
    return *this;
  }

  /*! Access the trapezoid id (PTR). */
  inline unsigned long id() const
  {
    return (unsigned long) PTR;
  }

  
  inline Curve_end curve_end() const  
  {
    return Curve_end(ptr()->cv, ptr()->ce);
  }

  
  /*! Access DAG node. */
  Dag_node* dag_node() const  
  {
    return ptr()->p_node; 
  } 
  
  
  //@}
  
  
	

};

} //namespace CGAL

#endif
