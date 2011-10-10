// Copyright (c) 1997  Tel-Aviv University (Israel).
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
// $URL: svn+ssh://balasmic@scm.gforge.inria.fr/svn/cgal/trunk/Arrangement_on_surface_2/include/CGAL/Arr_point_location/Td_vertex.h $
// $Id: Td_vertex.h 56667 2010-06-09 07:37:13Z sloriot $
// 
//
// Author(s)	 : Michal Balas <balasmic@post.tau.ac.il>

#ifndef CGAL_TD_VERTEX_H
#define CGAL_TD_VERTEX_H

/*! \file
 * Defintion of the Td_vertex<Td_traits> class.
 */

#include <CGAL/Arr_point_location/Trapezoidal_decomposition_2.h>

#ifdef CGAL_TD_DEBUG
#define CGAL_TD_INLINE
#else
#define CGAL_TD_INLINE inline
#endif

namespace CGAL {

/*! \class
 * Implementation of a vertex trapezoidal map item
 * Td_vertex holds the vertex and bottm & top halfedges (used for locating
 * unbounded trapezoids).
 * Td_vertex-s are created as active and become inactive when Remove() member
 * function called.
 * Each Td_vertex has at most four neighbouring Td_halfedges.
 */
template <class Td_traits_>
class Td_vertex : public Handle
{
public:
  //forward declarations
 // class Td_halfedge;


  //type of traits class
  typedef Td_traits_                              Traits;
  
  //type of Halfedge_const_handle 
  typedef typename Traits::Halfedge_const_handle  Halfedge_const_handle;

  //type of Vertex_const_handle 
  typedef typename Traits::Vertex_const_handle    Vertex_const_handle;
  
  //type of Td_vertex (Self)
  typedef Td_vertex                               Self;

  typedef typename Traits::Td_halfedge            Td_halfedge;
  
  //type of Td_vertex parameter space
  // seventuple which represents the Td_vertex map item:
  //  vertex, 
  //  bottom halfedge, top halfedge, (needed for locating unbounded trapezoids)
  //  left-bottom neighbor Td_halfedge, left-top neighbor Td_halfedge,
  //  right-bottom neighbor Td_halfedge, right-top neighbor Td_halfedge
  typedef Td_seventuple<Vertex_const_handle,
                        Halfedge_const_handle,
                        Halfedge_const_handle,
                        Td_halfedge*, Td_halfedge*,
                        Td_halfedge*, Td_halfedge*>   Vtx_param_space;

  //type of Trapezoidal decomposition
  typedef Trapezoidal_decomposition_2<Traits>     TD;
  
  //type of Around point circulator
  typedef typename TD::Around_point_circulator    Around_point_circulator;
  
  //type of In face iterator
  typedef typename TD::In_face_iterator           In_face_iterator;

  //type of Trapezoidal map search structure
  typedef typename TD::Dag_node                   Dag_node;


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
  
  Vtx_param_space* ptr() const { return (Vtx_param_space*)(PTR);  }
	
	
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
	
  /*! Initialize the Td_vertex neighbours. */
  CGAL_TD_INLINE void init_neighbours(Td_halfedge* lb = 0, Td_halfedge* lt = 0,
                                      Td_halfedge* rb = 0, Td_halfedge* rt = 0)
  {
    set_lb(lb);
    set_lt(lt);
    set_rb(rb);
    set_rt(rt);
  }

  /*! Set the DAG node. */
  CGAL_TD_INLINE void set_dag_node(Dag_node* p_node) 
  {
    m_dag_node = p_node;
  
#ifdef CGAL_TD_DEBUG
  
    CGAL_assertion(!p_node || **p_node == *this);
  
#endif	
	
  }

  /*! Set the Td_vertex vertex itself (Vertex_const_handle). */
  CGAL_TD_INLINE void set_v(Vertex_const_handle v) 
  {
    ptr()->e0 = v;
  }
  
  /*! Set the Td_vertex bottom halfedge (Halfedge_const_handle). */
  CGAL_TD_INLINE void set_btm_he(Halfedge_const_handle he) 
  {
    if (btm_he()->direction() != he->direction())
    {
      ptr()->e1 = he->twin();
    }
    else
    {
      ptr()->e1 = he;
    }
  }
  
  /*! Set the Td_vertex top halfedge (Halfedge_const_handle). */
  CGAL_TD_INLINE void set_top_he(Halfedge_const_handle he) 
  {
    if (top_he()->direction() != he->direction())
    {
      ptr()->e2 = he->twin();
    }
    else
    {
      ptr()->e2 = he;
    }
  }
  
  /*! Set left bottom neighbour. */
  CGAL_TD_INLINE void set_lb(Td_halfedge* lb) { ptr()->e3 = lb; }

  /*! Set left top neighbour. */
  CGAL_TD_INLINE void set_lt(Td_halfedge* lt) { ptr()->e4 = lt; }
  
  /*! Set right bottom neighbour. */
  CGAL_TD_INLINE void set_rb(Td_halfedge* rb) { ptr()->e5 = rb; }
  
  /*! Set right top neighbour. */
  CGAL_TD_INLINE void set_rt(Td_halfedge* rt) { ptr()->e6 = rt; }

 public:
  
  /// \name Constructors.
  //@{

  /*! Default constructor. */
  Td_vertex()
  {
    PTR = new Vtx_param_space(Traits::vtx_at_left_infinity(),
                              Traits::he_at_bottom_infinity(),
                              Traits::he_at_top_infinity(), 0, 0, 0, 0);
    m_dag_node = 0;
  }
  
  /*! Constructor given Vertex & Halfedge handles. */
  Td_vertex (Vertex_const_handle v,
             Halfedge_const_handle btm_he,
             Halfedge_const_handle top_he,
             Td_halfedge* lb = 0, Td_halfedge* lt = 0,
             Td_halfedge* rb = 0, Td_halfedge* rt = 0,
             Dag_node* node = 0)
  {
    PTR = new Vtx_param_space(v, btm_he, top_he, lb, lt, rb, rt);
    m_dag_node = node;
  }
  
  /*! Constructor given indicator , Vertex & Halfedge handles. */
  Td_vertex (bool is_v_given,
             Vertex_const_handle v,
             Halfedge_const_handle btm_he,
             Halfedge_const_handle top_he,
             Td_halfedge* lb = 0, Td_halfedge* lt = 0,
             Td_halfedge* rb = 0, Td_halfedge* rt = 0,
             Dag_node* node = 0)
  {
    PTR = new Vtx_param_space
                   (is_v_given? v : Traits::vtx_at_left_infinity(),
                    is_v_given? btm_he : Traits::he_at_bottom_infinity(),
                    is_v_given? top_he : Traits::he_at_top_infinity(),
                    lb, lt, rb, rt);
    m_dag_node = node;
  }
  
  /*! Copy constructor. */
  Td_vertex (const Self& td_v) : Handle(td_v)
  {
    m_dag_node = td_v.m_dag_node;
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

  /*! Access the Td_vertex id (PTR). */
  CGAL_TD_INLINE unsigned long id() const
  {
    return (unsigned long) PTR;
  }

  /*! Access the vertex itself. */
  CGAL_TD_INLINE Vertex_const_handle v() const
  {
    return ptr()->e0;
  }

  /*! Access Td_vertex bottom halfedge. */
  CGAL_TD_INLINE Halfedge_const_handle btm_he() const
  {
    return ptr()->e1;
  }

  /*! Access Td_vertex top halfedge. */
  CGAL_TD_INLINE Halfedge_const_handle top_he() const
  {
    return ptr()->e2;
  }

  /*! Access left bottom neighbour. */
  CGAL_TD_INLINE Td_halfedge* lb() const    { return ptr()->e3; }
  
  /*! Access left top neighbour. */
  CGAL_TD_INLINE Td_halfedge* lt() const    { return ptr()->e4; }
  
  /*! Access right bottom neighbour. */
  CGAL_TD_INLINE Td_halfedge* rb() const    { return ptr()->e5; }
  
  /*! Access right top neighbour. */
  CGAL_TD_INLINE Td_halfedge* rt() const    { return ptr()->e6; }
  
  /*! Access DAG node. */
  CGAL_TD_INLINE Dag_node* dag_node() const  { return m_dag_node;  }

  
  //@}

  /*! is Td_vertex active */
  bool is_active() const 
  {
    return (rb() != (Self*)CGAL_TD_DELETE_SIGNATURE);
  }

  /*! Removing this Td_vertex (defining it as in-active) */
  CGAL_TD_INLINE void remove(Dag_node* left=0)
  {
    CGAL_precondition(is_active());
 
    // mark Td_vertex as deleted,
    set_rb((Self*)CGAL_TD_DELETE_SIGNATURE);
  		
    // resets left son in data structure depending on input.
    if (left)
      m_dag_node->set_left_child(*left);
  }								
	
#ifdef CGAL_TD_DEBUG
  
  bool is_valid(const Traits* traits) const
  {
    if (!is_active()) return true;

    if (get_node() && **get_node()!=*this)
    {
      std::cerr << "\nthis=";
      write(std::cerr,*this,*traits,false);
      std::cerr << "\nget_node= ";
      write(std::cerr,**get_node(),*traits,false) << std::flush;
      CGAL_warning(**get_node()==*this);
      return false;
    }
    
    //precondition: v is in the x-range of btm_he
    Comparison_result res = traits->compare_curve_end_y_at_x_2_object()
                                      (v->curve_end(), btm_he()); 
    if (res != EQUAL)
    {
      std::cerr << "\nthis=";
      write(std::cerr, *this, *traits, false) << std::flush;
      std::cerr << "\nres==" << res << std::flush;
      CGAL_warning(res == EQUAL);
      return false;
    }
      
    //precondition: v is in the x-range of top_he
    res = traits->compare_curve_end_y_at_x_2_object()
                    (v->curve_end(), top_he()); 
    if (res != EQUAL)
    {
      std::cerr << "\nthis=";
      write(std::cerr, *this, *traits, false) << std::flush;
      std::cerr << "\nres==" << res << std::flush;
      CGAL_warning(res == EQUAL);
      return false;
    }

    //rt & lb should be Td_halfedge if they exist
    if ((rt() && !rt()->is_active()) ||
        (lb() && !lb()->is_active()) )
    {
      CGAL_warning(!rt() || rt()->is_active());
      CGAL_warning(!lb() || lb()->is_active());
      return false;
    }
    
    return true;
  }

  CGAL_TD_INLINE void debug() const // instantiate ptr functions.
  {
    ptr();
    v();
    btm_he();
    top_he();
  }

#endif //CGAL_TD_DEBUG

};

} //namespace CGAL

#endif
