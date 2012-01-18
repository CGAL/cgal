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
// $URL: svn+ssh://balasmic@scm.gforge.inria.fr/svn/cgal/branches/features/Arrangement_on_surface_2-RIC_pl_for_unbounded-balasmic/Arrangement_on_surface_2/include/CGAL/Arr_point_location/Td_inactive_fictitious_vertex.h $
// $Id: Td_inactive_fictitious_vertex.h 65793 2011-10-10 17:04:38Z balasmic $
// 
//
// Author(s)	 : Oren Nechushtan <theoren@math.tau.ac.il>
//               updated by: Michal Balas <balasmic@post.tau.ac.il>

#ifndef CGAL_TD_INACTIVE_FICTITIOUS_VERTEX_H
#define CGAL_TD_INACTIVE_FICTITIOUS_VERTEX_H

/*! \file
 * Defintion of the Td_inactive_fictitious_vertex<Td_traits> class.
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

  //type of Curve_end_pair
  typedef typename Traits::Curve_end_pair         Curve_end_pair;

  //type of Halfedge_const_handle (trapezoid edge)
  typedef typename Traits::Halfedge_const_handle  Halfedge_const_handle;
  
  //type of Vertex_const_handle (trapezoid vertex)
  typedef typename Traits::Vertex_const_handle    Vertex_const_handle;

  //type of Td_inactive_fictitious_vertex (Self)
  typedef typename Traits::Td_inactive_fictitious_vertex            Self;

  typedef typename Traits::Td_map_item            Td_map_item;

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
  
     /*! \class
   * Inner class Data derived from Rep class
   */
  class Data : public Rep
  {
    friend class Td_inactive_fictitious_vertex<Td_traits_>;

  public:
    //c'tors
    Data (boost::shared_ptr<X_monotone_curve_2> _cv,   
          unsigned char _chr, 
          Dag_node* _p_node): cv(_cv),chr(_chr),p_node(_p_node)  //MICHAL: Do we need neighbours for inactive fict vertex?
    { }
    
    ~Data() { }

  protected:
    boost::shared_ptr<X_monotone_curve_2> cv; 
    unsigned char chr;
    Dag_node* p_node;
  };
  
 private:
  
  Data* ptr() const { return (Data*)(PTR);  }
	
	
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
	
  ///*! Initialize the trapezoid's neighbours. */
  //inline void init_neighbours(boost::optional<Td_map_item> lb = boost::none, boost::optional<Td_map_item> lt = boost::none,
  //                            boost::optional<Td_map_item> rb = boost::none, boost::optional<Td_map_item> rt = boost::none)
  //{
  //  set_lb(lb);
  //  set_lt(lt);
  //  set_rb(rb);
  //  set_rt(rt);
  //}

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
  
  inline void set_curve_end(Vertex_const_handle v_before_rem)
  {
    Curve_end v_ce(v_before_rem->curve_end());
    ptr()->cv = (boost::shared_ptr<X_monotone_curve_2>)(new X_monotone_curve_2(v_ce.cv()));
    //CGAL_assertion(boost::get<boost::shared_ptr<X_monotone_curve_2>>( &(ptr()->e2)) != NULL);
    ptr()->chr = (v_ce.ce() == ARR_MIN_END ) ? CGAL_TD_CV_MIN_END : CGAL_TD_CV_MAX_END;
  }
  
 ///*! Set left bottom neighbour. */
 // inline void set_lb(boost::optional<Td_map_item> lb) {  }
 // 
 // /*! Set left top neighbour. */
 // inline void set_lt(boost::optional<Td_map_item> lt) {  }
 // 
 // /*! Set right bottom neighbour. */
 // inline void set_rb(boost::optional<Td_map_item> rb) {  }
 // 
 // /*! Set right top neighbour. */
 // inline void set_rt(boost::optional<Td_map_item> rt) {  }

 public:
  
  /// \name Constructors.
  //@{

  /*! Constructor given Vertex & Halfedge handles. */
  Td_inactive_fictitious_vertex (Vertex_const_handle v_before_rem, Dag_node* node = NULL)
  {
    Curve_end v_ce(v_before_rem->curve_end());
   
    PTR = new Data((boost::shared_ptr<X_monotone_curve_2>)(new X_monotone_curve_2(v_ce.cv())),
                   (v_ce.ce() == ARR_MIN_END ) ? CGAL_TD_CV_MIN_END : CGAL_TD_CV_MAX_END, node);
    //m_dag_node = node;
  }
  
  /*! Copy constructor. */
  Td_inactive_fictitious_vertex (const Self& tr) : Handle(tr)
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
    Curve_end_pair pair(curve_end_pair());
    return Curve_end(pair);
  }

  inline Curve_end_pair curve_end_pair() const  
  {
    X_monotone_curve_2* cv_ptr = (ptr()->cv).get();
    CGAL_assertion(cv_ptr != NULL);
   
    Arr_curve_end ce = 
      (ptr()->chr == CGAL_TD_CV_MIN_END) ?
        ARR_MIN_END : ARR_MAX_END;
  
    return std::make_pair(cv_ptr, ce);
  }

  ///*! Access left bottom neighbour. */
  //boost::optional<Td_map_item> lb() const    { return boost::none; }
  //
  ///*! Access left top neighbour. */
  //boost::optional<Td_map_item> lt() const    { return boost::none; }
  //
  ///*! Access right bottom neighbour. */
  //boost::optional<Td_map_item> rb() const    { return boost::none; }
  //
  ///*! Access right top neighbour. */
  //boost::optional<Td_map_item> rt() const    { return boost::none; }
  
  /*! Access DAG node. */
  Dag_node* dag_node() const            {return ptr()->p_node; } //m_dag_node;}
  
  
  //@}
  
  
	

};

} //namespace CGAL

#endif
