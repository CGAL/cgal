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
// $URL$
// $Id$
// 
//
// Author(s)	 : Oren Nechushtan <theoren@math.tau.ac.il>

#ifndef CGAL_TD_X_TRAPEZOID_H
#define CGAL_TD_X_TRAPEZOID_H

/*! \file
 * Defintion of the Td_X_trapezoid<Td_traits> class.
 */

#include <CGAL/Arr_point_location/Trapezoidal_decomposition_2.h>

#ifdef CGAL_TD_DEBUG
#define CGAL_TD_INLINE
#else
#define CGAL_TD_INLINE inline
#endif

namespace CGAL {

/*! \class
 * Implementation of a pseudo-trapezoid as two curves(top,bottom)
 * and two points(left,right).
 * Trapezoids are represented as two points called right and left and
 * two curves called top and bottom. The points lie on the right and left
 * boundaries of the trapezoid respectively and the curves bound the trapezoid
 * from above and below.
 * There exist degenerate trapezoids called infinite trapezoid; this happens 
 * when one of the four sides degenerates into a point/X_curve at infinity.
 * Trapezoids are created as active and become inactive when Remove() member
 * function called.
 * Each trapezoid has at most four neighbouring trapezoids.
 */
template <class Td_traits_>
class Td_X_trapezoid : public Handle
{
 public:
  typedef Td_traits_                                   Traits;
  typedef typename Traits::Point                       Point;
  typedef typename Traits::X_curve                     X_curve;
  typedef typename Traits::X_curve_ptr                 curve_pointer;
  typedef typename Traits::X_curve_ref                 curve_ref;
  typedef typename Traits::X_curve_const_ref           curve_const_ref;
  typedef typename Traits::X_trapezoid                 X_trapezoid;
  typedef X_trapezoid                                  Self;
  typedef typename Traits::X_trapezoid_ptr             pointer;
  typedef typename Traits::X_trapezoid_ref             ref;
  typedef typename Traits::X_trapezoid_const_ref       const_ref;
  typedef Td_ninetuple<Point, Point, X_curve, X_curve, unsigned char,
                       pointer, pointer,
                       pointer, pointer>               Arr_parameter_space;
  typedef Trapezoidal_decomposition_2<Traits>          TD;
  typedef typename TD::Unbounded                       Unbounded;
  typedef typename TD::Around_point_circulator         Around_point_circulator;
  typedef typename TD::In_face_iterator                In_face_iterator;
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
  
  typedef typename TD::Data_structure Data_structure;
  
 private:
  
  Arr_parameter_space* ptr() const {return (Arr_parameter_space*)(PTR);}
	
	
#ifndef CGAL_TD_DEBUG
#ifdef CGAL_PM_FRIEND_CLASS
 protected:
#else
 public: // workaround
#endif
#else //CGAL_TD_DEBUG
 public:
#endif //CGAL_TD_DEBUG
	
  Data_structure* node;
	
  CGAL_TD_INLINE void init_neighbours(pointer lb_ = 0, pointer lt_ = 0,
                                      pointer rb_ = 0, pointer rt_ = 0)
  {set_lb(lb_);set_lt(lt_);set_rb(rb_);set_rt(rt_);}
  CGAL_TD_INLINE void set_node(Data_structure* p) {node=p;
  
#ifdef CGAL_TD_DEBUG
  
    CGAL_assertion(!p || **p==*this);
  
#endif	
	
  }
  CGAL_TD_INLINE void set_left(const Point& p) 
  {ptr()->e0=p;ptr()->e4&=~CGAL_TRAPEZOIDAL_DECOMPOSITION_2_LEFT_UNBOUNDED;}
  
  CGAL_TD_INLINE void set_right(const Point& p) 
  {ptr()->e1=p;ptr()->e4&=~CGAL_TRAPEZOIDAL_DECOMPOSITION_2_RIGHT_UNBOUNDED;}
  
  CGAL_TD_INLINE void set_bottom(const X_curve& cv) 
  {ptr()->e2=cv;ptr()->e4&=~CGAL_TRAPEZOIDAL_DECOMPOSITION_2_BOTTOM_UNBOUNDED;}
  
  CGAL_TD_INLINE void set_top(const X_curve& cv) 
  {ptr()->e3=cv;ptr()->e4&=~CGAL_TRAPEZOIDAL_DECOMPOSITION_2_TOP_UNBOUNDED;}
  
  CGAL_TD_INLINE void set_left_unbounded() 
  {ptr()->e4|=CGAL_TRAPEZOIDAL_DECOMPOSITION_2_LEFT_UNBOUNDED;}
  
  CGAL_TD_INLINE void set_right_unbounded() 
  {ptr()->e4|=CGAL_TRAPEZOIDAL_DECOMPOSITION_2_RIGHT_UNBOUNDED;}
  
  CGAL_TD_INLINE void set_bottom_unbounded() 
  {ptr()->e4|=CGAL_TRAPEZOIDAL_DECOMPOSITION_2_BOTTOM_UNBOUNDED;}
  
  CGAL_TD_INLINE void set_top_unbounded() 
  {ptr()->e4|=CGAL_TRAPEZOIDAL_DECOMPOSITION_2_TOP_UNBOUNDED;}
  
  CGAL_TD_INLINE void set_lb(X_trapezoid* lb) {ptr()->e5=lb;}
  CGAL_TD_INLINE void set_lt(X_trapezoid* lt) {ptr()->e6=lt;}
  CGAL_TD_INLINE void set_rb(X_trapezoid* rb) {ptr()->e7=rb;}
  CGAL_TD_INLINE void set_rt(X_trapezoid* rt) {ptr()->e8=rt;}

 public:
  
  Td_X_trapezoid()
  {
    PTR = new Arr_parameter_space
      (Traits::point_at_left_top_infinity(),
       Traits::point_at_right_bottom_infinity(),
       Traits::curve_at_infinity(),
       Traits::curve_at_infinity(),
       CGAL_TRAPEZOIDAL_DECOMPOSITION_2_TOTALLY_UNBOUNDED,
       0, 0, 0, 0);
    node = 0;
  }
  
  Td_X_trapezoid (const Point &l, const Point &r,
                  const X_curve &b, const X_curve &t,
                  unsigned char c = CGAL_TRAPEZOIDAL_DECOMPOSITION_2_BOUNDED,
                  X_trapezoid *lb = 0, X_trapezoid *lt = 0,
                  X_trapezoid *rb = 0, X_trapezoid *rt = 0,
                  Data_structure *p = 0)
    {
      PTR = new Arr_parameter_space(l, r, b, t, c, lb, lt, rb, rt);
      node = p;
    }
  
  Td_X_trapezoid (const Point *l, const Point *r ,
                  const X_curve *b, const X_curve *t,
                  X_trapezoid *lb = 0, X_trapezoid *lt = 0,
                  X_trapezoid *rb = 0, X_trapezoid *rt = 0,
                  Data_structure *p = 0)
    {
      PTR = new Arr_parameter_space
	(l ? *l : Traits::point_at_left_top_infinity(),
	 r ? *r : Traits::point_at_right_bottom_infinity(),
	 b ? *b : Traits::curve_at_infinity(),
	 t ? *t : Traits::curve_at_infinity(),
	 ((l ? 0 : CGAL_TRAPEZOIDAL_DECOMPOSITION_2_LEFT_UNBOUNDED) | 
	  (r ? 0 : CGAL_TRAPEZOIDAL_DECOMPOSITION_2_RIGHT_UNBOUNDED) | 
	  (b ? 0 : CGAL_TRAPEZOIDAL_DECOMPOSITION_2_BOTTOM_UNBOUNDED) | 
	  (t ? 0 : CGAL_TRAPEZOIDAL_DECOMPOSITION_2_TOP_UNBOUNDED)),
	 lb, lt, rb, rt);
      node = p;
    }
  
  Td_X_trapezoid (const X_trapezoid &tr) :
    Handle(tr)
    {
      node = tr.node;
    }
  
    /*
      remark:
      operator= should not copy node (or otherwise update 
      Data_structure::replace)
    */
    CGAL_TD_INLINE X_trapezoid& operator= (const X_trapezoid& t2)
      {
	Handle::operator=(t2);
	return *this;
      }
    CGAL_TD_INLINE ref self() 
    {
      return *this;
    }
    CGAL_TD_INLINE const_ref self() const 
    {
      return *this;
    }
    CGAL_TD_INLINE unsigned long id () const
    {
      return (unsigned long) PTR;
    }

    CGAL_TD_INLINE bool operator== (const X_trapezoid& t2) const
    {
      return CGAL::identical(*this,t2);
    }

    CGAL_TD_INLINE bool operator!= (const X_trapezoid& t2) const
    {
      return !(operator==(t2));
    }

    CGAL_TD_INLINE const Point& left_unsafe () const
    {
      return ptr()->e0;
    }

    CGAL_TD_INLINE const Point& left () const
    {
      return !is_left_unbounded() ? 
	left_unsafe() : Traits::point_at_left_top_infinity();
    }
  
    CGAL_TD_INLINE const Point& right_unsafe () const
    {
      return ptr()->e1;
    }

    CGAL_TD_INLINE const Point& right () const
    {
      return !is_right_unbounded() ? 
	right_unsafe() : Traits::point_at_right_bottom_infinity();
    }
  
    // filters out the infinite case where at returns predefined dummy values
    CGAL_TD_INLINE const X_curve& bottom_unsafe () const
    {
      return ptr()->e2;
    }

    // filters out the infinite case where at returns predefined dummy values
    CGAL_TD_INLINE const X_curve& bottom () const
    {
      return !is_bottom_unbounded() ?  
	bottom_unsafe() : Traits::curve_at_infinity();
    }
  
    CGAL_TD_INLINE const X_curve& top_unsafe () const
    {
      return ptr()->e3;
    }

    CGAL_TD_INLINE const X_curve& top () const
    {
      return !is_top_unbounded() ?	
	top_unsafe() : Traits::curve_at_infinity();
    }
  
    unsigned char boundedness() const {return ptr()->e4;}
  
    CGAL_TD_INLINE bool is_left_unbounded() const {
      return (ptr()->e4&CGAL_TRAPEZOIDAL_DECOMPOSITION_2_LEFT_UNBOUNDED)!=0;}
    CGAL_TD_INLINE bool is_right_unbounded() const {
      return (ptr()->e4&CGAL_TRAPEZOIDAL_DECOMPOSITION_2_RIGHT_UNBOUNDED)!=0;}
    CGAL_TD_INLINE bool is_bottom_unbounded() const {
      return (ptr()->e4&CGAL_TRAPEZOIDAL_DECOMPOSITION_2_BOTTOM_UNBOUNDED)!=0;}
    CGAL_TD_INLINE bool is_top_unbounded() const {
      return (ptr()->e4&CGAL_TRAPEZOIDAL_DECOMPOSITION_2_TOP_UNBOUNDED)!=0;}
    CGAL_TD_INLINE bool is_unbounded() const
    {			
      CGAL_assertion(is_active());
      return (ptr()->e4&CGAL_TRAPEZOIDAL_DECOMPOSITION_2_TOTALLY_UNBOUNDED)!=0;
    }

    CGAL_TD_INLINE bool is_top_curve_equal(const Self& t,
                                           const Traits* traits) const {
      if (is_top_unbounded()) return t.is_top_unbounded();
      else if (t.is_top_unbounded()) return false;
      return traits->equal_2_object()(top_unsafe(),t.top_unsafe());
    }
    CGAL_TD_INLINE bool is_bottom_curve_equal(const Self& t,
                                              const Traits* traits) const {
      if (is_bottom_unbounded()) return t.is_bottom_unbounded();
      else if (t.is_bottom_unbounded()) return false;
      return traits->equal_2_object()(bottom_unsafe(),t.bottom_unsafe());
    }

    pointer get_lb() const {return ptr()->e5;}
  
    pointer get_lt() const {return ptr()->e6;}
  
    pointer get_rb() const {return ptr()->e7;}
  
    pointer get_rt() const {return ptr()->e8;}
  
    pointer left_bottom_neighbour() const {return get_lb();}
  
    pointer left_top_neighbour() const {return get_lt();}
  
    pointer right_bottom_neighbour() const {return get_rb();}
  
    pointer right_top_neighbour() const {return get_rt();}
  
    Data_structure* get_node() const {return node;}
  
    bool is_active() const 
    {
      return right_bottom_neighbour()!=
	(pointer)CGAL_TRAPEZOIDAL_DECOMPOSITION_2_DELETE_SIGNATURE;
    }
  
    CGAL_TD_INLINE void remove(Data_structure* left=0)
    {
      CGAL_precondition(is_active());
 
      // mark trapezoid as deleted,
      set_rb((pointer)CGAL_TRAPEZOIDAL_DECOMPOSITION_2_DELETE_SIGNATURE);
    		
      // resets left son in data structure depending on input.
      if (left)
	node->set_left(*left);
    }								

    /* precondition:
       both trapezoidal are active and have the same
       bounding edges from above and below and the trapezoids are adjacent to
       one another with the first to the left
       postcondition:
       this trapezoid is the union of the old this trapezoid
       and the input trapezoids
    */
    CGAL_TD_INLINE void merge_trapezoid(X_trapezoid& right)
    {
      CGAL_assertion(!is_right_unbounded());

      bool right_unbounded = right.is_right_unbounded();
      *this=X_trapezoid (!is_left_unbounded() ? &left() : 0,
			 !right_unbounded ? &right.right() : 0,
			 !is_bottom_unbounded() ? &bottom() : 0,
			 !is_top_unbounded() ? &top() : 0,
			 left_bottom_neighbour(),left_top_neighbour(),
			 right.right_bottom_neighbour(),
			 right.right_top_neighbour());

      if (right_bottom_neighbour())
	right_bottom_neighbour()->set_lb(this);
    
      if (right_top_neighbour())
	right_top_neighbour()->set_lt(this);

      CGAL_assertion(is_right_unbounded()==right.is_right_unbounded());
    }
	
#ifdef CGAL_TD_DEBUG
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
        if (!is_left_unbounded() && !is_right_unbounded() &&
            CGAL_POINT_IS_LEFT_LOW(right(),left()))
        {
          std::cerr << "\nthis=";
          write(std::cerr,*this,*traits,false) << std::flush;
          CGAL_warning(!CGAL_POINT_IS_LEFT_LOW(right(),left()));
          return false;
        }
        
        if (!is_bottom_unbounded())
        {
          if (is_left_unbounded() || is_right_unbounded())
          {
            std::cerr << "\nthis=";
            write(std::cerr,*this,*traits,false) << std::flush;
            CGAL_warning(!(is_left_unbounded() ||is_right_unbounded()));
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
        if (!is_top_unbounded())
        {
          if (is_left_unbounded() || is_right_unbounded())
          {
            std::cerr << "\nthis=";
            write(std::cerr,*this,*traits,false) << std::flush;
            CGAL_warning(!(is_left_unbounded() || is_right_unbounded()));
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
          if (right_top_neighbour() && 
              (! is_top_curve_equal(*right_top_neighbour(), traits)) ||
              left_top_neighbour() && 
              (! is_top_curve_equal(*left_top_neighbour(), traits)) ||
              right_bottom_neighbour() &&
              (! is_bottom_curve_equal(*right_bottom_neighbour(), traits)) ||
              left_bottom_neighbour() &&
              (! is_bottom_curve_equal(*left_bottom_neighbour(), traits)) ||
              right_top_neighbour() &&
              traits->is_degenerate(*right_top_neighbour()) ||
              left_top_neighbour() &&
              traits->is_degenerate(*left_top_neighbour()) ||
              right_bottom_neighbour() &&
              traits->is_degenerate(*right_bottom_neighbour()) ||
              left_bottom_neighbour() &&
              traits->is_degenerate(*left_bottom_neighbour()))
          {
            std::cerr << "\nthis=";
            write(std::cerr,*this,*traits,false) << std::flush;
            CGAL_warning(!(right_top_neighbour() &&
                           (! is_top_curve_equal(*right_top_neighbour(), traits))));
            CGAL_warning(!(left_top_neighbour() &&
                           (! is_top_curve_equal(*left_top_neighbour(), traits))));
            CGAL_warning(!(right_bottom_neighbour() &&
                           (! is_bottom_curve_equal(*right_bottom_neighbour(), traits))));
            CGAL_warning(!(left_bottom_neighbour() &&
                           (! is_bottom_curve_equal(*left_bottom_neighbour(), traits))));
            CGAL_warning(!(right_top_neighbour() &&
                           traits->is_degenerate(*right_top_neighbour())));
            CGAL_warning(!(left_top_neighbour() &&
                           traits->is_degenerate(*left_top_neighbour())));
            CGAL_warning(!(right_bottom_neighbour() &&
                           traits->is_degenerate(*right_bottom_neighbour())));
            CGAL_warning(!(left_bottom_neighbour() &&
                           traits->is_degenerate(*left_bottom_neighbour())));
            return false;
          }
          if (right_top_neighbour()&&!right_top_neighbour()->is_active()||
              left_top_neighbour()&&!left_top_neighbour()->is_active()||
              right_bottom_neighbour()&&!right_bottom_neighbour()->is_active()||
              left_bottom_neighbour()&&!left_bottom_neighbour()->is_active())
          {
            std::cerr << "\nleft=" << left() << " right=" << right()
                      << " bottom=" << bottom() << " top=" << top()
                      << std::flush;
            CGAL_warning(!(right_top_neighbour() &&
                           !right_top_neighbour()->is_active()));
            CGAL_warning(!(left_top_neighbour() &&
                           !left_top_neighbour()->is_active()));
            CGAL_warning(!(right_bottom_neighbour() &&
                           !right_bottom_neighbour()->is_active()));
            CGAL_warning(!(left_bottom_neighbour() &&
                           !left_bottom_neighbour()->is_active()));
            return false;
          }
        }
        else
        {
          /* if the trapezoid is degenerate, the left() and right()
             points should be on the top() and bottom() curves.
             In any case none of the geometric boundaries should 
             be unbounded */
          if (is_bottom_unbounded()||
              is_top_unbounded()||
              is_left_unbounded()||
              is_right_unbounded()
              )
          {
            std::cerr << "\nbottom()==" << bottom() << std::flush;
            std::cerr << "\ntop()==" << top() << std::flush;
            std::cerr << "\nleft()==" << left() << std::flush;
            std::cerr << "\nright()==" << right() << std::flush;
            CGAL_warning((!is_bottom_unbounded()));
            CGAL_warning((!is_top_unbounded()));
            CGAL_warning((!is_left_unbounded()));
            CGAL_warning((!is_right_unbounded()));
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
            if (right_top_neighbour()&&!right_top_neighbour()->is_active()||
                //!left_top_neighbour()||!left_top_neighbour()->is_active()||
                right_bottom_neighbour() &&
                !right_bottom_neighbour()->is_active()||
                left_bottom_neighbour() && !left_bottom_neighbour()->is_active()
                )
            {
              CGAL_warning(!right_top_neighbour() ||
                           right_top_neighbour()->is_active());
              //CGAL_warning(!left_top_neighbour() ||
              //left_top_neighbour()->is_active());
              CGAL_warning(!right_bottom_neighbour() ||
                           right_bottom_neighbour()->is_active());
              CGAL_warning(!left_bottom_neighbour() ||
                           left_bottom_neighbour()->is_active());
              return false;
            }
            if (
                /* if trapezoid is end relative to supporting X_curve, that is
                   adjacent(trapezoid's right end point,supporting X_curve right
                   end point) , right_top_neighbour() returns next such trapezoid
                   around right() point in clockwise oriented order
                   adjacent(trapezoid's left end point,supporting X_curve left end
                   point), left_bottom_neighbour() returns next such trapezoid
                   around left() point in clockwise oriented order */
                /* right_bottom_neighbour() points to next trapezoid on
                   supporting X_curve, if such exist */
                right_top_neighbour() &&
                !traits->is_degenerate_curve(*right_top_neighbour())||
                // !left_top_neighbour() ||
                // !traits->is_degenerate_curve(*left_top_neighbour())||
                right_bottom_neighbour() &&
                !traits->is_degenerate_curve(*right_bottom_neighbour())||
                left_bottom_neighbour() &&
                !traits->is_degenerate_curve(*left_bottom_neighbour())
                )
            {
              CGAL_warning(!right_top_neighbour() ||
                           traits->is_degenerate_curve(*right_top_neighbour()));
              //CGAL_warning(!left_top_neighbour() ||
              //!traits->is_degenerate_curve(*left_top_neighbour()));
              CGAL_warning(!right_bottom_neighbour() ||
                           traits->
                           is_degenerate_curve(*right_bottom_neighbour()));
              CGAL_warning(!left_bottom_neighbour() ||
                           traits->
                           is_degenerate_curve(*left_bottom_neighbour()));
              return false;
            }
          }
          else if (traits->is_degenerate_point(*this))
          {
            if (right_top_neighbour() &&
                !traits->is_degenerate_curve(*right_top_neighbour())||
                left_bottom_neighbour() &&
                !traits->is_degenerate_curve(*left_bottom_neighbour())
                )
            {
              CGAL_warning(!right_top_neighbour() ||
                           traits->is_degenerate_curve(*right_top_neighbour()));
              CGAL_warning(!left_bottom_neighbour() ||
                           traits->
                           is_degenerate_curve(*left_bottom_neighbour()));
              return false;
            }
            if (right_top_neighbour()&&!right_top_neighbour()->is_active()||
                left_bottom_neighbour()&&!left_bottom_neighbour()->is_active()
                )
            {
              CGAL_warning(!right_top_neighbour() ||
                           right_top_neighbour()->is_active());
              CGAL_warning(!left_bottom_neighbour() ||
                           left_bottom_neighbour()->is_active());
              return false;
            }
            if (!traits->equal_2_object()(left(),right()))
            {
              std::cerr << "\nleft()==" << left() << std::flush;
              std::cerr << "\nright()==" << right() << std::flush;
              CGAL_warning(traits->equal_2_object()(left(),right()));
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
    left();
    right();
  }
#endif

};

} //namespace CGAL

#endif
