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
// $URL: svn+ssh://balasmic@scm.gforge.inria.fr/svn/cgal/trunk/Arrangement_on_surface_2/include/CGAL/Arr_point_location/Td_halfedge.h $
// $Id: Td_halfedge.h 56667 2010-06-09 07:37:13Z sloriot $
// 
//
// Author(s)	 : Michal Balas <balasmic@post.tau.ac.il>

#ifndef CGAL_TD_HALFEDGE_H
#define CGAL_TD_HALFEDGE_H

/*! \file
 * Defintion of the Td_halfedge<Td_traits> class.
 */

#include <CGAL/Arr_point_location/Trapezoidal_decomposition_2.h>

#ifdef CGAL_TD_DEBUG
#define CGAL_TD_INLINE
#else
#define CGAL_TD_INLINE inline
#endif

namespace CGAL {

/*! \class
 * Implementation of a halfedge trapezoidal map item
 * Td_halfedges are created as active and become inactive when Remove() member
 * function called.
 * Each Td_halfedge has at most four neighbouring Td_halfedges.
 */
template <class Td_traits_>
class Td_halfedge : public Handle
{
public:

  //type of traits class
  typedef Td_traits_                              Traits;
  
  //type of Halfedge_const_handle 
  typedef typename Traits::Halfedge_const_handle  Halfedge_const_handle;
  
  //type of Td_halfedge (Self)
  typedef Td_halfedge                             Self;
  
  //type of Td_halfedge parameter space
  // fivetuple which represents the Td_halfedge map item:
  //  halfedge,
  //  left-bottom neighbor Td_halfedge, left-top neighbor Td_halfedge,
  //  right-bottom neighbor Td_halfedge, right-top neighbor Td_halfedge
  typedef Td_fivetuple<Halfedge_const_handle,
                       Self*, Self*,
                       Self*, Self*>              He_param_space;

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
  
  He_param_space* ptr() const { return (He_param_space*)(PTR);  }
	
	
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
	
  /*! Initialize the Td_halfedge's neighbours. */
  CGAL_TD_INLINE void init_neighbours(Self* lb = 0, Self* lt = 0,
                                      Self* rb = 0, Self* rt = 0)
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
  
  /*! Set the Td_halfedge's  halfedge itself (Halfedge_const_handle). */
  CGAL_TD_INLINE void set_he(Halfedge_const_handle new_he) 
  {
    if (he()->direction() != new_he->direction())
    {
      ptr()->e0 = new_he->twin();
    }
    else
    {
      ptr()->e0 = new_he;
    }
  }
  
  /*! Set left bottom neighbour. */
  CGAL_TD_INLINE void set_lb(Self* lb) { ptr()->e1 = lb; }

  /*! Set left top neighbour. */
  CGAL_TD_INLINE void set_lt(Self* lt) { ptr()->e2 = lt; }
  
  /*! Set right bottom neighbour. */
  CGAL_TD_INLINE void set_rb(Self* rb) { ptr()->e3 = rb; }
  
  /*! Set right top neighbour. */
  CGAL_TD_INLINE void set_rt(Self* rt) { ptr()->e4 = rt; }

 public:
  
  /// \name Constructors.
  //@{

  /*! Default constructor. */
  Td_halfedge()
  {
    PTR = new He_param_space(Traits::he_at_bottom_infinity(),0, 0, 0, 0);
    m_dag_node = 0;
  }
  
  /*! Constructor given a Halfedge handle. */
  Td_halfedge (Halfedge_const_handle he,
               Self* lb = 0, Self* lt = 0,
               Self* rb = 0, Self* rt = 0,
               Dag_node* node = 0)
  {
    PTR = new He_param_space(he, lb, lt, rb, rt);
    m_dag_node = node;
  }
  
  /*! Constructor given indicator & Halfedge handle. */
  Td_halfedge (bool is_he_given,
               Halfedge_const_handle he,
               Self* lb = 0, Self* lt = 0,
               Self* rb = 0, Self* rt = 0,
               Dag_node* node = 0)
  {
    PTR = new He_param_space(is_he_given? he : Traits::he_at_bottom_infinity(),
                             lb, lt, rb, rt);
    m_dag_node = node;
  }
  
  /*! Copy constructor. */
  Td_halfedge (const Self& td_he) : Handle(td_he)
  {
    m_dag_node = td_he.m_dag_node;
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

  /*! Access the Td_halfedge id (PTR). */
  CGAL_TD_INLINE unsigned long id() const
  {
    return (unsigned long) PTR;
  }

  /*! Access the halfedge itself. */
  CGAL_TD_INLINE Halfedge_const_handle he() const
  {
    return ptr()->e0;
  }

  /*! Access left bottom neighbour. */
  CGAL_TD_INLINE Self* lb() const    { return ptr()->e1; }
  
  /*! Access left top neighbour. */
  CGAL_TD_INLINE Self* lt() const    { return ptr()->e2; }
  
  /*! Access right bottom neighbour. */
  CGAL_TD_INLINE Self* rb() const    { return ptr()->e3; }
  
  /*! Access right top neighbour. */
  CGAL_TD_INLINE Self* rt() const    { return ptr()->e4; }
  
  /*! Access DAG node. */
  CGAL_TD_INLINE Dag_node* dag_node() const  { return m_dag_node;  }

  
  //@}

  /*! is Td_halfedge active */
  bool is_active() const 
  {
    return (rb() != (Self*)CGAL_TD_DELETE_SIGNATURE);
  }

  /*! Removing this Td_halfedge (defining it as in-active) */
  CGAL_TD_INLINE void remove(Dag_node* left=0)
  {
    CGAL_precondition(is_active());
 
    // mark Td_halfedge as deleted,
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
    
    if ((rt() && !rt()->is_active()) ||
        (rb() && !rb()->is_active()) ||
        (lb() && !lb()->is_active()) )
    {
      CGAL_warning(!rt() || rt()->is_active());
      CGAL_warning(!rb() || rb()->is_active());
      CGAL_warning(!lb() || lb()->is_active());
      return false;
    }
    
    return true;
  }

  CGAL_TD_INLINE void debug() const // instantiate ptr functions.
  {
    ptr();
    he();
  }

#endif //CGAL_TD_DEBUG

};

} //namespace CGAL

#endif
