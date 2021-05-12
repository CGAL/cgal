// Copyright (c) 2005,2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)         : Oren Nechushtan <theoren@math.tau.ac.il>
//               updated by: Michal Balas <balasmic@post.tau.ac.il>

#ifndef CGAL_TD_X_TRAPEZOID_H
#define CGAL_TD_X_TRAPEZOID_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
 * Defintion of the Td_X_trapezoid<Td_traits> class.
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
class Td_X_trapezoid : public Handle
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

  //type of Trapezoid (Self)
  typedef typename Traits::X_trapezoid            Self;

  //type of Trapezoid parameter space
  // Ninetuple which represents the Trapezoid:
  //  - for regular & edge trapezoids or active point trapezoids:
  //      left vertex, right vertex, bottom halfedge, top halfedge
  //  - for removed point trapezoids:
  //      point or X_monotone_curve_2+ cv end
  //  type flag + on boundaries flags,
  //  left-bottom neighbor trapezoid, left-top neighbor trapezoid,
  //  right-bottom neighbor trapezoid, right-top neighbor trapezoid
  typedef Td_ninetuple<boost::variant<Vertex_const_handle,Point>,
                       boost::variant<Vertex_const_handle,unsigned char>,
                       boost::variant<Halfedge_const_handle,
                                      boost::shared_ptr<X_monotone_curve_2> >,
                       Halfedge_const_handle,
                       unsigned char,
                       Self*, Self*,
                       Self*, Self*>            Trpz_parameter_space;

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
  CGAL_TD_INLINE void init_neighbours(Self* lb_ = 0, Self* lt_ = 0,
                                      Self* rb_ = 0, Self* rt_ = 0)
  {
    set_lb(lb_);
    set_lt(lt_);
    set_rb(rb_);
    set_rt(rt_);
  }

  /*! Set the DAG node. */
  CGAL_TD_INLINE void set_dag_node(Dag_node* p)
  {
    m_dag_node = p;

#ifdef CGAL_TD_DEBUG

    CGAL_assertion(!p || **p == *this);

#endif

  }

  /*! Set the trapezoid's left (Vertex_const_handle). */
  CGAL_TD_INLINE void set_left(Vertex_const_handle v)
  {
    CGAL_precondition(is_active());
    ptr()->e0 = v;
  }

  /*! Set the trapezoid's right (Vertex_const_handle). */
  CGAL_TD_INLINE void set_right(Vertex_const_handle v)
  {
    CGAL_precondition(is_active());
    ptr()->e1 = v;
  }

  /*! Set the trapezoid's bottom (Halfedge_const_handle). */
  CGAL_TD_INLINE void set_bottom(Halfedge_const_handle he)
  {
    CGAL_precondition(is_active());
    if (!is_on_bottom_boundary() &&
        bottom_unsafe()->direction() != he->direction())
    {
      ptr()->e2 = he->twin();
    }
    else
    {
      ptr()->e2 = he;
    }
  }

  /*! Set the trapezoid's top (Halfedge_const_handle). */
  CGAL_TD_INLINE void set_top(Halfedge_const_handle he)
  {
    CGAL_precondition(is_active());
    if (!is_on_top_boundary() &&
        top_unsafe()->direction() != he->direction())
    {
      ptr()->e3 = he->twin();
    }
    else
    {
      ptr()->e3 = he;
    }
  }

  CGAL_TD_INLINE void update_removed_trpz()
  {
    CGAL_precondition(is_active());

    if (type() == TD_EDGE)
    {
      //ptr()->e2 = (boost::shared_ptr<X_monotone_curve_2>)(new X_monotone_curve_2(top()->curve()));
      set_curve_for_rem_he(top()->curve());
      return;
    }

    //else if (type() == TD_VERTEX)

    Curve_end v_ce(left()->curve_end());
    ptr()->e2 = (boost::shared_ptr<X_monotone_curve_2>)(new X_monotone_curve_2(v_ce.cv()));
    //CGAL_assertion(boost::get<boost::shared_ptr<X_monotone_curve_2>>( &(ptr()->e2)) != nullptr);

    ptr()->e1 = (v_ce.ce() == ARR_MIN_END ) ? CGAL_TD_CV_MIN_END : CGAL_TD_CV_MAX_END;

    if (!is_on_boundaries())
    { //if the trapezoid respresents an inner vertex
      ptr()->e0 = left()->point();
    }
  }


  /*! Set the x_monotone_curve_2 for removed edge degenerate trapezoid. */
  CGAL_TD_INLINE void set_curve_for_rem_he(const X_monotone_curve_2& cv)
  {
    CGAL_precondition (type() == TD_EDGE);

    ptr()->e2 = (boost::shared_ptr<X_monotone_curve_2>)(new X_monotone_curve_2(cv));
  }

  /*! Set the trapezoid's type flag (Trapezoid/Edge/Vertex). */
  CGAL_TD_INLINE void set_type(unsigned char obj_type)
  {
    ptr()->e4 &= ~CGAL_TD_TYPE_MASK;
    ptr()->e4 |= obj_type;
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
  Td_X_trapezoid()
  {
    //define the initial trapezoid: left, right, btm, top are at infinity.
    // its type is TD_TRAPEZOID ,it is on all boundaries, and has no neighbours
    PTR = new Trpz_parameter_space
      (Traits::vtx_at_left_infinity(),
       Traits::vtx_at_right_infinity(),
       Traits::he_at_bottom_infinity(),
       Traits::he_at_top_infinity(),
       CGAL_TD_TRAPEZOID | CGAL_TD_ON_ALL_BOUNDARIES ,
       0, 0, 0, 0);

    m_dag_node = 0;
  }

  /*! Constructor given Vertex & Halfedge handles. */
  Td_X_trapezoid (Vertex_const_handle l, Vertex_const_handle r,
                  Halfedge_const_handle b, Halfedge_const_handle t,
                  Type tp = TD_TRAPEZOID,
                  unsigned char boundness_flag = CGAL_TD_INTERIOR,
                  Self* lb = 0, Self* lt = 0,
                  Self* rb = 0, Self* rt = 0,
                  Dag_node* node = 0)
  {

    //build the type flag
    unsigned char type_flag = 0;
    if (tp == TD_TRAPEZOID)
      type_flag |= CGAL_TD_TRAPEZOID;
    else if (tp == TD_EDGE)
      type_flag |= CGAL_TD_EDGE;
    else //tp == TD_VERTEX
      type_flag |= CGAL_TD_VERTEX;

    PTR = new Trpz_parameter_space
      (l, r, b, t, type_flag | boundness_flag, lb, lt, rb, rt);
    m_dag_node = node;
  }

  /*! Constructor given Pointers to Vertex & Halfedge handles. */
  Td_X_trapezoid (Vertex_const_handle* l, Vertex_const_handle* r ,
                  Halfedge_const_handle* b, Halfedge_const_handle* t,
                  unsigned char type_flag,
                  bool  on_left_bndry,
                  bool  on_right_bndry,
                  bool  on_bottom_bndry,
                  bool  on_top_bndry,
                  Self* lb = 0, Self* lt = 0,
                  Self* rb = 0, Self* rt = 0,
                  Dag_node* node = 0)
  {
    PTR = new Trpz_parameter_space
      (l ? *l : Traits::vtx_at_left_infinity(),
       r ? *r : Traits::vtx_at_right_infinity(),
       b ? *b : Traits::he_at_bottom_infinity(),
       t ? *t : Traits::he_at_top_infinity(),
       (type_flag |
        (on_left_bndry   ? CGAL_TD_ON_LEFT_BOUNDARY   : 0) |
        (on_right_bndry  ? CGAL_TD_ON_RIGHT_BOUNDARY  : 0) |
        (on_bottom_bndry ? CGAL_TD_ON_BOTTOM_BOUNDARY : 0) |
        (on_top_bndry    ? CGAL_TD_ON_TOP_BOUNDARY    : 0) ),
         lb, lt, rb, rt);
    m_dag_node = node;
    }

  /*! Copy constructor. */
  Td_X_trapezoid (const Self& tr) : Handle(tr)
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

  /*! Access trapezoid left. */
  CGAL_TD_INLINE Vertex_const_handle left_unsafe() const
    {
    CGAL_precondition(is_active());
    CGAL_assertion(boost::get<Vertex_const_handle>(&(ptr()->e0)) != nullptr);
    return boost::get<Vertex_const_handle>(ptr()->e0);
    }

  /*! Access trapezoid left.
  *   filters out the infinite case which returns predefined dummy values
  */
  CGAL_TD_INLINE Vertex_const_handle left() const
    {
    CGAL_precondition(is_active());
    if (is_on_left_boundary() && is_on_bottom_boundary()
        && is_on_top_boundary())
    {
      return Traits::vtx_at_left_infinity();
    }
    //else
    return left_unsafe();
    }

  /*! Access trapezoid right. */
  CGAL_TD_INLINE Vertex_const_handle right_unsafe() const
    {
    CGAL_precondition(is_active());
    CGAL_assertion(boost::get<Vertex_const_handle>(&(ptr()->e1)) != nullptr);
    return boost::get<Vertex_const_handle>(ptr()->e1);
    }

  /*! Access trapezoid right.
  *   filters out the infinite case which returns predefined dummy values
  */
  CGAL_TD_INLINE Vertex_const_handle right () const
  {
    CGAL_precondition(is_active());
    if (is_on_right_boundary() && is_on_bottom_boundary()
        && is_on_top_boundary())
    {
      return Traits::vtx_at_right_infinity();
    }
    //else
    return right_unsafe();
  }

  /*! Access trapezoid bottom. */
  CGAL_TD_INLINE Halfedge_const_handle bottom_unsafe () const
  {
    CGAL_precondition(is_active());
    CGAL_assertion(boost::get<Halfedge_const_handle>(&(ptr()->e2)) != nullptr);
    return boost::get<Halfedge_const_handle>(ptr()->e2);
    }

  /*! Access trapezoid bottom.
  *   filters out the infinite case which returns predefined dummy values
  */
  CGAL_TD_INLINE Halfedge_const_handle bottom () const
    {
    CGAL_precondition(is_active());
    return !is_on_bottom_boundary() ?
            bottom_unsafe() : Traits::he_at_bottom_infinity();
    }

  /*! Access trapezoid top. */
  CGAL_TD_INLINE Halfedge_const_handle top_unsafe () const
    {
    CGAL_precondition(is_active());
    return ptr()->e3;
    }

  /*! Access trapezoid top.
  *   filters out the infinite case which returns predefined dummy values
  */
  CGAL_TD_INLINE Halfedge_const_handle top () const
    {
    CGAL_precondition(is_active());
    return !is_on_top_boundary() ?
            top_unsafe() : Traits::he_at_top_infinity();
    }

  CGAL_TD_INLINE Point point_for_inner_rem_vtx() const
    {
    CGAL_precondition(!is_active());
    CGAL_precondition(type() == TD_VERTEX);
    CGAL_precondition(!is_on_boundaries());

    CGAL_assertion(boost::get<Point>( &(ptr()->e0)) != nullptr);
    return boost::get<Point>( ptr()->e0 );
    }

  CGAL_TD_INLINE std::pair<X_monotone_curve_2*,Arr_curve_end> curve_end_pair_for_boundary_rem_vtx() const
    {
    CGAL_precondition(!is_active());
    CGAL_precondition(type() == TD_VERTEX);
    CGAL_precondition(is_on_boundaries());

    CGAL_assertion(boost::get<unsigned char>( &(ptr()->e1)) != nullptr);
    CGAL_assertion(boost::get<boost::shared_ptr<X_monotone_curve_2> >(&(ptr()->e2)) != nullptr);
    X_monotone_curve_2* cv_ptr = (boost::get<boost::shared_ptr<X_monotone_curve_2> >(ptr()->e2)).get();
    CGAL_assertion(cv_ptr != nullptr);

    Arr_curve_end ce =
      (boost::get<unsigned char>(ptr()->e1) == CGAL_TD_CV_MIN_END) ?
        ARR_MIN_END : ARR_MAX_END;

    return std::make_pair(cv_ptr, ce);
  }

  CGAL_TD_INLINE Curve_end curve_end_for_boundary_rem_vtx() const
  {
    CGAL_precondition(!is_active());
    CGAL_precondition(type() == TD_VERTEX);
    CGAL_precondition(is_on_boundaries());

    CGAL_assertion(boost::get<unsigned char>( &(ptr()->e1)) != nullptr);
    CGAL_assertion(boost::get<boost::shared_ptr<X_monotone_curve_2> >(&(ptr()->e2)) != nullptr);
    X_monotone_curve_2* cv_ptr = (boost::get<boost::shared_ptr<X_monotone_curve_2> >(ptr()->e2)).get();
    CGAL_assertion(cv_ptr != nullptr);

    Arr_curve_end ce =
      (boost::get<unsigned char>(ptr()->e1) == CGAL_TD_CV_MIN_END) ?
        ARR_MIN_END : ARR_MAX_END;

    return Curve_end(*cv_ptr, ce);
  }

  CGAL_TD_INLINE Curve_end curve_end_for_rem_vtx() const
  {
    CGAL_precondition(!is_active());
    CGAL_precondition(type() == TD_VERTEX);

    CGAL_assertion(boost::get<unsigned char>( &(ptr()->e1)) != nullptr);
    CGAL_assertion(boost::get<boost::shared_ptr<X_monotone_curve_2> >(&(ptr()->e2)) != nullptr);
    X_monotone_curve_2* cv_ptr = (boost::get<boost::shared_ptr<X_monotone_curve_2> >(ptr()->e2)).get();
    CGAL_assertion(cv_ptr != nullptr);

    Arr_curve_end ce =
      (boost::get<unsigned char>(ptr()->e1) == CGAL_TD_CV_MIN_END) ?
        ARR_MIN_END : ARR_MAX_END;

    return Curve_end(*cv_ptr, ce);
  }

  CGAL_TD_INLINE X_monotone_curve_2& curve_for_rem_he() const
  {
    CGAL_precondition(!is_active() && type() == TD_EDGE);

    CGAL_assertion(boost::get<boost::shared_ptr<X_monotone_curve_2> >(&(ptr()->e2)) != nullptr);
    X_monotone_curve_2* cv_ptr = (boost::get<boost::shared_ptr<X_monotone_curve_2> >(ptr()->e2)).get();
    CGAL_assertion(cv_ptr != nullptr);
    return *cv_ptr;
  }

  /*! Access trapezoid type. */
  CGAL_TD_INLINE Type type() const
  {
    switch(ptr()->e4 & CGAL_TD_TYPE_MASK)
    {
    case CGAL_TD_TRAPEZOID:
      return TD_TRAPEZOID;
    case CGAL_TD_EDGE:
      return TD_EDGE;
    case CGAL_TD_VERTEX:
      return TD_VERTEX;
    default:
       CGAL_assertion(false);
       return TD_TRAPEZOID;
    }
    }

  /*! Access trapezoid type flag. */
  CGAL_TD_INLINE unsigned char type_flag() const
    {
    return (ptr()->e4 & CGAL_TD_TYPE_MASK);
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
    return (ptr()->e4 & CGAL_TD_ON_ALL_BOUNDARIES) != 0;
  }

  /*! Access left bottom neighbour. */
  Self* lb() const    { return ptr()->e5; }

  /*! Access left top neighbour. */
  Self* lt() const    { return ptr()->e6; }

  /*! Access right bottom neighbour. */
  Self* rb() const    { return ptr()->e7; }

  /*! Access right top neighbour. */
  Self* rt() const    { return ptr()->e8; }

  /*! Access DAG node. */
  Dag_node* dag_node() const            {return m_dag_node;}


  //@}

  /*! is trapezoid active */
    bool is_active() const
    {
    return rb()!=
            (Self*)CGAL_TD_DELETE_SIGNATURE;
    }

  /*! Removing this trapezoid (defining it as in-active) */
  CGAL_TD_INLINE void remove(Dag_node* left=0)
    {
      CGAL_precondition(is_active());

    // update vertex/edge trapezoid parameters after remove
    if (type() != TD_TRAPEZOID)
      update_removed_trpz();

      // mark trapezoid as deleted,
    set_rb((Self*)CGAL_TD_DELETE_SIGNATURE);

    if (type() == TD_VERTEX)
      curve_end_for_rem_vtx();

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
    //precondition: both are of type trapezoid
    CGAL_precondition((type() == TD_TRAPEZOID) &&
                      (right.type() == TD_TRAPEZOID));
    //precondition: both are active
    CGAL_precondition(is_active() && right.is_active());
    //precondition: the left trapezoid is not on the right boundary
    CGAL_assertion(!is_on_right_boundary());

    bool on_right_boundary = right.is_on_right_boundary();

    *this = Self (!is_on_left_boundary() ? & left() : 0,
                              !on_right_boundary ? &right.right() : 0,
                              !is_on_bottom_boundary() ? &bottom() : 0,
                              !is_on_top_boundary() ? &top() : 0,
                  CGAL_TD_TRAPEZOID,
                  is_on_left_boundary(), on_right_boundary,
                  is_on_bottom_boundary(), is_on_top_boundary(),
                              lb(),lt(),
                              right.rb(),
                              right.rt());

    if (rb())
      rb()->set_lb(this);

    if (rt())
      rt()->set_lt(this);

    CGAL_assertion(is_on_right_boundary() == right.is_on_right_boundary());
    }

#ifdef CGAL_TD_DEBUG
  //MICHAL: This method should not compile!!
    bool is_valid(const Traits* traits) const
    {
      Comparison_result t;
      bool              b;

      if (is_active())
      {
        if (get_node() && **get_node()!=*this)
        {
          std::cerr << "\nthis=";
          write(std::cerr,*this,*traits,false);
          std::cerr << "\nget_node= ";
          write(std::cerr,**get_node(),*traits,false) << std::flush;
          CGAL_warning(**get_node()==*this);
          return false;
        }
        if (!is_on_left_boundary() && !is_on_right_boundary() &&
            CGAL_POINT_IS_LEFT_LOW(right(),left()))
        {
          std::cerr << "\nthis=";
          write(std::cerr,*this,*traits,false) << std::flush;
          CGAL_warning(!CGAL_POINT_IS_LEFT_LOW(right(),left()));
          return false;
        }

        if (!is_on_bottom_boundary())
        {
          if (is_on_left_boundary() || is_on_right_boundary())
          {
            std::cerr << "\nthis=";
            write(std::cerr,*this,*traits,false) << std::flush;
            CGAL_warning(!(is_on_left_boundary() ||is_on_right_boundary()));
            return false;
          }

          b = CGAL_IS_IN_X_RANGE(bottom(),left());
          if (b) {
            t = CGAL_CURVE_COMPARE_Y_AT_X(left(), bottom());
          }
          if (!b || t == SMALLER)
          {
            std::cerr << "\nthis=";
            write(std::cerr,*this,*traits,false) << std::flush;
            std::cerr << "\nt==" << t << std::flush;
            CGAL_warning(b);
            CGAL_warning(t != SMALLER);
            return false;
          }

          b=CGAL_IS_IN_X_RANGE(bottom(),right());
          if (b) {
            t = CGAL_CURVE_COMPARE_Y_AT_X(right(), bottom());
          }
          if (!b || t == SMALLER)
          {
            std::cerr << "\nthis=";
            write(std::cerr,*this,*traits,false) << std::flush;
            std::cerr << "\nt==" << t << std::flush;
            CGAL_warning(b);
            CGAL_warning(t != SMALLER);
            return false;
          }
        }
        if (!is_on_top_boundary())
        {
          if (is_on_left_boundary() || is_on_right_boundary())
          {
            std::cerr << "\nthis=";
            write(std::cerr,*this,*traits,false) << std::flush;
            CGAL_warning(!(is_on_left_boundary() || is_on_right_boundary()));
            return false;
          }

          b=CGAL_IS_IN_X_RANGE(top(),left());
          if (b) {
            t = CGAL_CURVE_COMPARE_Y_AT_X(left(), top());
          }
          if (!b || t == LARGER)
          {
            std::cerr << "\nthis=";
            write(std::cerr,*this,*traits,false) << std::flush;
            std::cerr << "\nt==" << t << std::flush;
            CGAL_warning(b);
            CGAL_warning(t != LARGER);
            return false;
          }

          b=CGAL_IS_IN_X_RANGE(top(),right());
          if (b) {
            t = CGAL_CURVE_COMPARE_Y_AT_X(right(), top());
          }
          if (!b || t == LARGER)
          {
            std::cerr << "\nthis=";
            write(std::cerr,*this,*traits,false) << std::flush;
            std::cerr << "\nt==" << t << std::flush;
            CGAL_warning(b);
            CGAL_warning(t != LARGER);
            return false;
          }
        }
        if (!traits->is_degenerate(*this))
        {
          if (rt() &&
              (! is_top_curve_equal(*rt(), traits)) ||
              lt() &&
              (! is_top_curve_equal(*lt(), traits)) ||
              rb() &&
              (! is_bottom_curve_equal(*rb(), traits)) ||
              lb() &&
              (! is_bottom_curve_equal(*lb(), traits)) ||
              rt() &&
              traits->is_degenerate(*rt()) ||
              lt() &&
              traits->is_degenerate(*lt()) ||
              rb() &&
              traits->is_degenerate(*rb()) ||
              lb() &&
              traits->is_degenerate(*lb()))
          {
            std::cerr << "\nthis=";
            write(std::cerr,*this,*traits,false) << std::flush;
            CGAL_warning(!(rt() &&
                           (! is_top_curve_equal(*rt(), traits))));
            CGAL_warning(!(lt() &&
                           (! is_top_curve_equal(*lt(), traits))));
            CGAL_warning(!(rb() &&
                           (! is_bottom_curve_equal(*rb(), traits))));
            CGAL_warning(!(lb() &&
                           (! is_bottom_curve_equal(*lb(), traits))));
            CGAL_warning(!(rt() &&
                           traits->is_degenerate(*rt())));
            CGAL_warning(!(lt() &&
                           traits->is_degenerate(*lt())));
            CGAL_warning(!(rb() &&
                           traits->is_degenerate(*rb())));
            CGAL_warning(!(lb() &&
                           traits->is_degenerate(*lb())));
            return false;
          }
          if (rt()&&!rt()->is_active()||
              lt()&&!lt()->is_active()||
              rb()&&!rb()->is_active()||
              lb()&&!lb()->is_active())
          {
            std::cerr << "\nleft=" << left() << " right=" << right()
                      << " bottom=" << bottom() << " top=" << top()
                      << std::flush;
            CGAL_warning(!(rt() &&
                           !rt()->is_active()));
            CGAL_warning(!(lt() &&
                           !lt()->is_active()));
            CGAL_warning(!(rb() &&
                           !rb()->is_active()));
            CGAL_warning(!(lb() &&
                           !lb()->is_active()));
            return false;
          }
        }
        else
        {
          /* if the trapezoid is degenerate, the left() and right()
             points should be on the top() and bottom() curves.
             In any case none of the geometric boundaries should
             be unbounded */
          if (is_on_bottom_boundary()||
              is_on_top_boundary()||
              is_on_left_boundary()||
              is_on_right_boundary()
              )
          {
            std::cerr << "\nbottom()==" << bottom() << std::flush;
            std::cerr << "\ntop()==" << top() << std::flush;
            std::cerr << "\nleft()==" << left() << std::flush;
            std::cerr << "\nright()==" << right() << std::flush;
            CGAL_warning((!is_on_bottom_boundary()));
            CGAL_warning((!is_on_top_boundary()));
            CGAL_warning((!is_on_left_boundary()));
            CGAL_warning((!is_on_right_boundary()));
            return false;
          }
          if (!CGAL_IS_IN_X_RANGE(bottom(),left()) ||
              CGAL_CURVE_COMPARE_Y_AT_X(left(), bottom()) != EQUAL)
          {
            std::cerr << "\nbottom()==" << bottom() << std::flush;
            std::cerr << "\nleft()==" << left() << std::flush;
            CGAL_warning(CGAL_IS_IN_X_RANGE(bottom(),left()) &&
                         CGAL_CURVE_COMPARE_Y_AT_X(left(), bottom()) ==
                         EQUAL);
            return false;
          }
          if (!CGAL_IS_IN_X_RANGE(bottom(),right()) ||
              CGAL_CURVE_COMPARE_Y_AT_X(right(), bottom()) != EQUAL)
          {
            std::cerr << "\nbottom()==" << bottom() << std::flush;
            std::cerr << "\nright()==" << right() << std::flush;
            CGAL_warning(CGAL_IS_IN_X_RANGE(bottom(),right()) &&
                         CGAL_CURVE_COMPARE_Y_AT_X(right(), bottom()) ==
                         EQUAL);
            return false;
          }
          if (!CGAL_IS_IN_X_RANGE(top(),left()) ||
              CGAL_CURVE_COMPARE_Y_AT_X(left(), top()) != EQUAL)
          {
            std::cerr << "\ntop()==" << top() << std::flush;
            std::cerr << "\nleft()==" << left() << std::flush;
            CGAL_warning(!CGAL_IS_IN_X_RANGE(top(),left()) &&
                         CGAL_CURVE_COMPARE_Y_AT_X(left(), top()) == EQUAL);
            return false;
          }
          if (!CGAL_IS_IN_X_RANGE(top(),right()) ||
              CGAL_CURVE_COMPARE_Y_AT_X(right(), top()) != EQUAL)
          {
            std::cerr << "\ntop()==" << top() << std::flush;
            std::cerr << "\nright()==" << right() << std::flush;
            CGAL_warning(CGAL_IS_IN_X_RANGE(top(),right()) &&
                         CGAL_CURVE_COMPARE_Y_AT_X(right(), top()) == EQUAL);
            return false;
          }
          if (traits->is_degenerate_curve(*this))
          {
            if (rt()&&!rt()->is_active()||
                //!lt()||!lt()->is_active()||
                rb() &&
                !rb()->is_active()||
                lb() && !lb()->is_active()
                )
            {
              CGAL_warning(!rt() ||
                           rt()->is_active());
              //CGAL_warning(!lt() ||
              //lt()->is_active());
              CGAL_warning(!rb() ||
                           rb()->is_active());
              CGAL_warning(!lb() ||
                           lb()->is_active());
              return false;
            }
            if (
                /* if trapezoid is end relative to supporting X_curve, that is
                   adjacent(trapezoid's right end point,supporting X_curve right
                   end point) , rt() returns next such trapezoid
                   around right() point in clockwise oriented order
                   adjacent(trapezoid's left end point,supporting X_curve left end
                   point), lb() returns next such trapezoid
                   around left() point in clockwise oriented order */
                /* rb() points to next trapezoid on
                   supporting X_curve, if such exist */
                rt() &&
                !traits->is_degenerate_curve(*rt())||
                // !lt() ||
                // !traits->is_degenerate_curve(*lt())||
                rb() &&
                !traits->is_degenerate_curve(*rb())||
                lb() &&
                !traits->is_degenerate_curve(*lb())
                )
            {
              CGAL_warning(!rt() ||
                           traits->is_degenerate_curve(*rt()));
              //CGAL_warning(!lt() ||
              //!traits->is_degenerate_curve(*lt()));
              CGAL_warning(!rb() ||
                           traits->
                           is_degenerate_curve(*rb()));
              CGAL_warning(!lb() ||
                           traits->
                           is_degenerate_curve(*lb()));
              return false;
            }
          }
          else if (traits->is_degenerate_point(*this))
          {
            if (rt() &&
                !traits->is_degenerate_curve(*rt())||
                lb() &&
                !traits->is_degenerate_curve(*lb())
                )
            {
              CGAL_warning(!rt() ||
                           traits->is_degenerate_curve(*rt()));
              CGAL_warning(!lb() ||
                           traits->
                           is_degenerate_curve(*lb()));
              return false;
            }
            if (rt()&&!rt()->is_active()||
                lb()&&!lb()->is_active()
                )
            {
              CGAL_warning(!rt() ||
                           rt()->is_active());
              CGAL_warning(!lb() ||
                           lb()->is_active());
              return false;
            }
            if (!traits->equal_curve_end_2_object()(left(),right()))
            {
              std::cerr << "\nleft()==" << left() << std::flush;
              std::cerr << "\nright()==" << right() << std::flush;
              CGAL_warning(traits->equal_curve_end_2_object()(left(),right()));
              return false;
            }
          }
        }
      }
      return true;
    }

  CGAL_TD_INLINE void debug() const // instantiate ptr functions.
  {
    ptr();
    bottom();
    top();
    if (type() == TD_VERTEX && !is_active())
    {
      if (!is_on_boundaries())
        point_for_inner_rem_vtx();
      else
        curve_end_for_boundary_rem_vtx();
    }
    else //MICHAL: this is problematic since left, right may not exist for inactive trapezoids as well
    {
    left();
    right();
  }
  }
#endif //CGAL_TD_DEBUG

};

} //namespace CGAL

#endif
