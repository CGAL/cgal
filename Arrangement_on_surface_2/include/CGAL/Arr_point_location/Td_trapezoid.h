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
// $URL: svn+ssh://balasmic@scm.gforge.inria.fr/svn/cgal/trunk/Arrangement_on_surface_2/include/CGAL/Arr_point_location/Td_trapezoid.h $
// $Id: Td_trapezoid.h 56667 2010-06-09 07:37:13Z sloriot $
// 
//
// Author(s)	 : Michal Balas <balasmic@post.tau.ac.il> 
//        (according to previous version of Oren Nechushtan <theoren@math.tau.ac.il>)

#ifndef CGAL_TD_TRAPEZOID_H
#define CGAL_TD_TRAPEZOID_H

/*! \file
 * Defintion of the Td_trapezoid<Td_traits> class.
 */

#include <CGAL/Arr_point_location/Trapezoidal_decomposition_2.h>

#ifdef CGAL_TD_DEBUG
#define CGAL_TD_INLINE
#else
#define CGAL_TD_INLINE inline
#endif

namespace CGAL {

/*! \class
 * Implementation of a td-trapezoid as two halfedges(top,bottom)
 * and two vertices(left,right).
 * Trapezoids are represented as two vertices called right and left and
 * two halfedges called top and bottom. The vertices lie on the 
 * right and left boundaries of the trapezoid respectively and the halfedges 
 * bound the trapezoid from above and below.
 * There exist degenerate trapezoids called infinite trapezoid; this happens 
 * when one of the four sides is on the parameter space boundary.
 * Trapezoids are created as active and become inactive when Remove() member
 * function called.
 * Each trapezoid has at most four neighbouring trapezoids.
 */
template <class Td_traits_>
class Td_trapezoid : public Handle
{
public:
  
  //type of traits class
  typedef Td_traits_                              Traits;
  
  //type of Halfedge_const_handle 
  typedef typename Traits::Halfedge_const_handle  Halfedge_const_handle;

  //type of Vertex_const_handle 
  typedef typename Traits::Vertex_const_handle    Vertex_const_handle;
  
  //type of Td_trapezoid (Self)
  typedef Td_trapezoid                            Self;

  //type of Td_trapezoid parameter space
  // ninetuple which represents the Td_trapezoid map item:
  //  left vertex, right vertex,
  //  bottom halfedge, top halfedge, 
  //  on boundaries flags,
  //  left-bottom neighbor trapezoid, left-top neighbor trapezoid,
  //  right-bottom neighbor trapezoid, right-top neighbor trapezoid
  typedef Td_ninetuple<Vertex_const_handle, Vertex_const_handle, 
                       Halfedge_const_handle, Halfedge_const_handle, 
                       unsigned char,
                       Self*, Self*,
                       Self*, Self*>              Trpz_param_space;

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
  
  Trpz_param_space* ptr() const { return (Trpz_param_space*)(PTR);  }
	
	
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
	
  /*! Initialize the Td_trapezoid neighbours. */
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

  /*! Set the Td_trapezoid left vertex (Vertex_const_handle). */
  CGAL_TD_INLINE void set_left_v(Vertex_const_handle v) 
  {
    ptr()->e0 = v;
  }
  
  /*! Set the Td_trapezoid right vertex (Vertex_const_handle). */
  CGAL_TD_INLINE void set_right_v(Vertex_const_handle v) 
  {
    ptr()->e1 = v;
  }
  
  /*! Set the Td_trapezoid bottom halfedge (Halfedge_const_handle). */
  CGAL_TD_INLINE void set_btm_he(Halfedge_const_handle he) 
  {
    if (btm_he()->direction() != he->direction())
    {
      ptr()->e2 = he->twin();
    }
    else
    {
      ptr()->e2 = he;
    }
  }
  
  /*! Set the Td_trapezoid top halfedge (Halfedge_const_handle). */
  CGAL_TD_INLINE void set_top_he(Halfedge_const_handle he) 
  {
    if (top_he()->direction() != he->direction())
    {
      ptr()->e3 = he->twin();
    }
    else
    {
      ptr()->e3 = he;
    }
  }

  /*! Set is on left boundary flag. */
  CGAL_TD_INLINE void set_is_on_left_boundary(bool b) 
  {
    if (b)
      ptr()->e4 |= CGAL_TD_ON_LEFT_BOUNDARY;
    else
      ptr()->e4 &= ~CGAL_TD_ON_LEFT_BOUNDARY;
  }
  
  /*! Set is on right boundary flag. */
  CGAL_TD_INLINE void set_is_on_right_boundary(bool b) 
  {
    if (b)
      ptr()->e4 |= CGAL_TD_ON_RIGHT_BOUNDARY;
    else
      ptr()->e4 &= ~CGAL_TD_ON_RIGHT_BOUNDARY;
  }
  
  /*! Set is on bottom boundary flag. */
  CGAL_TD_INLINE void set_is_on_bottom_boundary(bool b) 
  {
    if (b)
      ptr()->e4 |= CGAL_TD_ON_BOTTOM_BOUNDARY;
    else
      ptr()->e4 &= ~CGAL_TD_ON_BOTTOM_BOUNDARY;
  }
  
  /*! Set is on top boundary flag. */
  CGAL_TD_INLINE void set_is_on_top_boundary(bool b) 
  {
    if (b)
      ptr()->e4 |= CGAL_TD_ON_TOP_BOUNDARY;
    else
      ptr()->e4 &= ~CGAL_TD_ON_TOP_BOUNDARY;
  }
  
  /*! Set left bottom neighbour. */
  CGAL_TD_INLINE void set_lb(Self* lb) { ptr()->e5 = lb; }

  /*! Set left top neighbour. */
  CGAL_TD_INLINE void set_lt(Self* lt) { ptr()->e6 = lt; }
  
  /*! Set right bottom neighbour. */
  CGAL_TD_INLINE void set_rb(Self* rb) { ptr()->e7 = rb; }
  
  /*! Set right top neighbour. */
  CGAL_TD_INLINE void set_rt(Self* rt) { ptr()->e8 = rt; }

 public:
  
  /// \name Constructors.
  //@{

  /*! Default constructor. */
  Td_trapezoid()
  {
    //define the initial trapezoid: left, right, btm, top are at infinity.
    // it is on all boundaries, and has no neighbours
    PTR = new Trpz_param_space(Traits::vtx_at_left_infinity(),
                               Traits::vtx_at_right_infinity(),
                               Traits::he_at_bottom_infinity(),
                               Traits::he_at_top_infinity(),
                               CGAL_TD_ON_ALL_BOUNDARIES ,
                               0, 0, 0, 0);
    m_dag_node = 0;
  }
  
  /*! Constructor given Vertex & Halfedge handles. */
  Td_trapezoid (Vertex_const_handle l, Vertex_const_handle r,
                Halfedge_const_handle b, Halfedge_const_handle t,
                unsigned char boundness_flag = CGAL_TD_INTERIOR,
                Self* lb = 0, Self* lt = 0,
                Self* rb = 0, Self* rt = 0,
                Dag_node* node = 0)
  {
    PTR = new Trpz_param_space(l, r, b, t, boundness_flag, lb, lt, rb, rt);
    m_dag_node = node;
  }
  
  /*! Constructor given on boundary flags , Vertex & Halfedge handles. */
  Td_trapezoid (Vertex_const_handle l, Vertex_const_handle r,
                Halfedge_const_handle b, Halfedge_const_handle t,
                bool  on_left_bndry, bool  on_right_bndry,
                bool  on_btm_bndry, bool  on_top_bndry,
                Self* lb = 0, Self* lt = 0,
                Self* rb = 0, Self* rt = 0,
                Dag_node* node = 0)
  {
    PTR = new Trpz_param_space
               (on_left_bndry ? Traits::vtx_at_left_infinity()  : l,
                on_right_bndry? Traits::vtx_at_right_infinity() : r, 
                on_btm_bndry  ? Traits::he_at_bottom_infinity() : b,
                on_top_bndry  ? Traits::he_at_top_infinity()    : t,
                ((on_left_bndry  ? CGAL_TD_ON_LEFT_BOUNDARY   : 0) | 
                 (on_right_bndry ? CGAL_TD_ON_RIGHT_BOUNDARY  : 0) | 
                 (on_btm_bndry   ? CGAL_TD_ON_BOTTOM_BOUNDARY : 0) | 
                 (on_top_bndry   ? CGAL_TD_ON_TOP_BOUNDARY    : 0) ),
                lb, lt, rb, rt);
    m_dag_node = node;
  }
  
  /*! Copy constructor. */
  Td_trapezoid (const Self& td_trpz) : Handle(td_trpz)
  {
    m_dag_node = td_trpz.m_dag_node;
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

  /*! Access the Td_trapezoid id (PTR). */
  CGAL_TD_INLINE unsigned long id() const
  {
    return (unsigned long) PTR;
  }

  /*! Access trapezoid left. */
  CGAL_TD_INLINE Vertex_const_handle left_v_unsafe() const
  {
    return ptr()->e0;
  }

  /*! Access trapezoid left. 
  *   filters out the infinite case which returns predefined dummy values
  */
  CGAL_TD_INLINE Vertex_const_handle left_v() const
  {
    if (is_on_left_boundary() && is_on_bottom_boundary()
        && is_on_top_boundary())
    {
      return Traits::vtx_at_left_infinity();
    }
    //else
    return left_v_unsafe();
  }

  /*! Access trapezoid right. */
  CGAL_TD_INLINE Vertex_const_handle right_v_unsafe() const
  {
    return ptr()->e1;
  }

  /*! Access trapezoid right. 
  *   filters out the infinite case which returns predefined dummy values
  */
  CGAL_TD_INLINE Vertex_const_handle right_v () const
  {
    if (is_on_right_boundary() && is_on_bottom_boundary()
        && is_on_top_boundary())
    {
      return Traits::vtx_at_right_infinity();
    }
    //else
    return right_v_unsafe();
  }

  /*! Access trapezoid bottom. */
  CGAL_TD_INLINE Halfedge_const_handle btm_he_unsafe () const
  {
    return ptr()->e2;
  }

  /*! Access trapezoid bottom. 
  *   filters out the infinite case which returns predefined dummy values
  */
  CGAL_TD_INLINE Halfedge_const_handle btm_he() const
  {
    return !is_on_bottom_boundary() ?  
            btm_he_unsafe() : Traits::he_at_bottom_infinity();
  }

  /*! Access trapezoid top. */
  CGAL_TD_INLINE Halfedge_const_handle top_he_unsafe () const
  {
    return ptr()->e3;
  }

  /*! Access trapezoid top. 
  *   filters out the infinite case which returns predefined dummy values
  */
  CGAL_TD_INLINE Halfedge_const_handle top_he() const
  {
    return !is_on_top_boundary() ?	
            top_he_unsafe() : Traits::he_at_top_infinity();
  }

  /*! Access on boundaries flag. */
  CGAL_TD_INLINE unsigned char on_boundaries_flag() const 
  {
    return (ptr()->e4 & CGAL_TD_ON_ALL_BOUNDARIES);
  }
  
  /*! Access is on left boundary. */
  CGAL_TD_INLINE bool is_on_left_boundary() const 
  {
    return (ptr()->e4 & CGAL_TD_ON_LEFT_BOUNDARY) != 0;
  }

  /*! Access is on right boundary. */
  CGAL_TD_INLINE bool is_on_right_boundary() const 
  {
    return (ptr()->e4 & CGAL_TD_ON_RIGHT_BOUNDARY) != 0;
  }

  /*! Access is on bottom boundary. */
  CGAL_TD_INLINE bool is_on_bottom_boundary() const 
  {
    return (ptr()->e4 & CGAL_TD_ON_BOTTOM_BOUNDARY) != 0;
  }

  /*! Access is on top boundary. */
  CGAL_TD_INLINE bool is_on_top_boundary() const 
  {
    return (ptr()->e4 & CGAL_TD_ON_TOP_BOUNDARY) != 0;
  }
  
  /*! Access is on at least one boundary. */
  CGAL_TD_INLINE bool is_on_boundaries() const
  {			
    CGAL_assertion(is_active());
    return (ptr()->e4 & CGAL_TD_ON_ALL_BOUNDARIES) != 0;
  }

  /*! Access left bottom neighbour. */
  CGAL_TD_INLINE Self* lb() const    { return ptr()->e3; }
  
  /*! Access left top neighbour. */
  CGAL_TD_INLINE Self* lt() const    { return ptr()->e4; }
  
  /*! Access right bottom neighbour. */
  CGAL_TD_INLINE Self* rb() const    { return ptr()->e5; }
  
  /*! Access right top neighbour. */
  CGAL_TD_INLINE Self* rt() const    { return ptr()->e6; }
  
  /*! Access DAG node. */
  CGAL_TD_INLINE Dag_node* dag_node() const  { return m_dag_node;  }

  
  //@}

  /*! is Td_trapezoid active */
  bool is_active() const 
  {
    return (rb() != (Self*)CGAL_TD_DELETE_SIGNATURE);
  }

  /*! Removing this Td_trapezoid (defining it as in-active) */
  CGAL_TD_INLINE void remove(Dag_node* left=0)
  {
    CGAL_precondition(is_active());
 
    // mark Td_trapezoid as deleted,
    set_rb((Self*)CGAL_TD_DELETE_SIGNATURE);
  		
    // resets left son in data structure depending on input.
    if (left)
      m_dag_node->set_left_child(*left);
  }		

  /* Merge this trapezoid with the input trapezoid.
     Precondition:
      both trapezoids are active and have the same
      bounding edges from above and below and the trapezoids are adjacent to
      one another with the first to the left
     Postcondition:
      this trapezoid is the union of the old this trapezoid
      and the input trapezoid
  */
  CGAL_TD_INLINE void merge_trapezoid( Self& right)
  {
    //precondition: the left trapezoid is not on the right boundary
    CGAL_assertion(!is_on_right_boundary());

    bool on_right_boundary = right.is_on_right_boundary();

    *this = Self
      (!is_on_left_boundary() ? left_v() : Traits::vtx_at_left_infinity(),
		   !on_right_boundary ? right.right_v() : Traits::vtx_at_right_infinity(),
		   !is_on_bottom_boundary() ? btm_he() : Traits::he_at_bottom_infinity(),
		   !is_on_top_boundary() ? top_he() : Traits::he_at_top_infinity(),
       is_on_left_boundary(), on_right_boundary,
       is_on_bottom_boundary(), is_on_top_boundary(),
		   lb(),lt(), right.rb(), right.rt());
		              
    if (rb())  rb()->set_lb(this);
  
    if (rt())  rt()->set_lt(this);

    CGAL_assertion(is_on_right_boundary() == right.is_on_right_boundary());
  }
	
	
#ifdef CGAL_TD_DEBUG
  
  bool is_valid(const Traits* traits) const
  {
    if (!is_active()) return true;

    Comparison_result res;

    if (get_node() && **get_node()!=*this)
    {
      std::cerr << "\nthis=";
      write(std::cerr,*this,*traits,false);
      std::cerr << "\nget_node= ";
      write(std::cerr, **get_node(), *traits, false) << std::flush;
      CGAL_warning(**get_node()==*this);
      return false;
    }

    //check that right_v is not lexicographically smaller than left_v
    if (traits->compare_curve_end_xy_2_object()
                  (right_v()->curve_end(), left_v()->curve_end()) == SMALLER)
    {
      std::cerr << "\nthis=";
      write(std::cerr, *this, *traits, false) << std::flush;
      CGAL_warning(traits->compare_curve_end_xy_2_object()
                    (right_v()->curve_end(), left_v()->curve_end()) != SMALLER);
      return false;
    }
    
    if (!is_on_bottom_boundary())
    {
      //check that left_v is not below btm_he
      res = traits->compare_curve_end_y_at_x_2_object()
                        (left_v()->curve_end(), btm_he()); 
      if (res == SMALLER)
      {
        std::cerr << "\nthis=";
        write(std::cerr, *this, *traits, false) << std::flush;
        std::cerr << "\nres==" << res << std::flush;
        CGAL_warning(res != SMALLER);
        return false;
      }
        
      //check that right_v is not below btm_he
      res = traits->compare_curve_end_y_at_x_2_object()
                         (right_v()->curve_end(), btm_he()); 
      if (res == SMALLER)
      {
        std::cerr << "\nthis=";
        write(std::cerr, *this, *traits, false) << std::flush;
        std::cerr << "\nres==" << res << std::flush;
        CGAL_warning(res != SMALLER);
        return false;
      }
    }

    if (!is_on_top_boundary())
    {
      //check that left_v is not above top_he
      res = traits->compare_curve_end_y_at_x_2_object()
                        (left_v()->curve_end(), top_he()); 
      if (res == LARGER)
      {
        std::cerr << "\nthis=";
        write(std::cerr, *this, *traits, false) << std::flush;
        std::cerr << "\nres==" << res << std::flush;
        CGAL_warning(res != LARGER);
        return false;
      }
        
      //check that right_v is not above top_he
      res = traits->compare_curve_end_y_at_x_2_object()
                         (right_v()->curve_end(), top_he()); 
      if (res == LARGER)
      {
        std::cerr << "\nthis=";
        write(std::cerr, *this, *traits, false) << std::flush;
        std::cerr << "\nres==" << res << std::flush;
        CGAL_warning(res != LARGER);
        return false;
      }
    }

    if ((rt() && 
         !( (is_on_top_boundary() && rt()->is_on_top_boundary()) ||
            (traits->equal_2_object()(top_he()->curve(), 
                                      rt()->top_he()->curve()))   ) )   ||
        (lt() && 
         !( (is_on_top_boundary() && lt()->is_on_top_boundary()) ||
            (traits->equal_2_object()(top_he()->curve(), 
                                      lt()->top_he()->curve()))   ) )   ||
        (rb() && 
         !( (is_on_bottom_boundary() && rb()->is_on_bottom_boundary())||
            (traits->equal_2_object()(btm_he()->curve(), 
                                      rb()->btm_he()->curve()))   ) )   ||
        (lb() && 
         !( (is_on_bottom_boundary() && lb()->is_on_bottom_boundary())||
            (traits->equal_2_object()(btm_he()->curve(), 
                                      lb()->btm_he()->curve()))   ) )   )
    {
      std::cerr << "\nthis=";
      write(std::cerr, *this, *traits, false) << std::flush;
      CGAL_warning(!(rt() && 
                    !( (is_on_top_boundary() && rt()->is_on_top_boundary()) ||
                       (traits->equal_2_object()(top_he()->curve(), 
                                              rt()->top_he()->curve()))   ) ));
      CGAL_warning(!(lt() && 
                    !( (is_on_top_boundary() && lt()->is_on_top_boundary()) ||
                        (traits->equal_2_object()(top_he()->curve(), 
                                              lt()->top_he()->curve()))   ) ));
      CGAL_warning(!(rb() && 
                    !( (is_on_bottom_boundary() && 
                          rb()->is_on_bottom_boundary())||
                       (traits->equal_2_object()(btm_he()->curve(), 
                                              rb()->btm_he()->curve()))   ) ));
      CGAL_warning(!(lb() && 
                    !( (is_on_bottom_boundary() && 
                          lb()->is_on_bottom_boundary())||
                       (traits->equal_2_object()(btm_he()->curve(), 
                                              lb()->btm_he()->curve()))   ) ));
      return false;
    }
          
    if ( (rt() && !rt()->is_active() ) ||
         (lt() && !lt()->is_active() ) ||
         (rb() && !rb()->is_active() ) ||
         (lb() && !lb()->is_active() )  )
    {
      std::cerr << "\nleft=" << left_v() << " right=" << right)v()
                << " bottom=" << btm_he() << " top=" << top_he()
                << std::flush;
      CGAL_warning(!(rt() && !rt()->is_active()));
      CGAL_warning(!(lt() && !lt()->is_active()));
      CGAL_warning(!(rb() && !rb()->is_active()));
      CGAL_warning(!(lb() && !lb()->is_active()));
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
