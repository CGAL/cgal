// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.2-I-26 $
// release_date  : $CGAL_Date: 2000/07/11 $
//
// file          : include/CGAL/Td_X_trapezoid.h
// package       : Trapezoidal_decomposition (1.06)
// maintainer    : Oren Nechushtan <theoren@math.tau.ac.il>
// source		 : 
// revision 	 : 
// revision_date : 
// author(s)	 : Oren Nechushtan <theoren@math.tau.ac.il>
//
//
// maintainer(s) : Oren Nechushtan <theoren@math.tau.ac.il>
//
//
// coordinator	 : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter		 : 
// ======================================================================

#ifndef CGAL_TD_X_TRAPEZOID_H
#define CGAL_TD_X_TRAPEZOID_H

/* ------------------------------------------------------------------------

  class 			CGAL::TD_X_trapezoid.
  parameters		Td_traits, implementation of a trapezoid traits.
  inherited from	CGAL::Handle.
  description		Implements a trapezoid as two curves(top,bottom) 
                          and two points(left,right).
  
------------------------------------------------------------------------ */

#ifndef CGAL_TRAPEZOIDAL_DECOMPOSITION_2_H
#include <CGAL/Trapezoidal_decomposition_2.h>
#endif

CGAL_BEGIN_NAMESPACE

template < class Td_traits_>
class Td_X_trapezoid : public Handle
{
public:
  typedef Td_traits_ Traits;
  typedef typename Traits::Point Point;
  typedef typename Traits::X_curve X_curve;
  typedef typename Traits::X_curve_ptr curve_pointer;
  typedef typename Traits::X_curve_ref curve_ref;
  typedef typename Traits::X_curve_const_ref curve_const_ref;
  typedef typename Traits::X_trapezoid X_trapezoid;
  typedef typename Traits::X_trapezoid_ptr pointer;
  typedef typename Traits::X_trapezoid_ref reference;
  typedef typename Traits::X_trapezoid_const_ref const_ref;
  typedef Td_ninetuple<Point,Point,X_curve,X_curve,unsigned char,
    pointer,pointer,pointer,pointer> Boundary_type;
  typedef Trapezoidal_decomposition_2<Traits> TD;
  typedef typename TD::Unbounded Unbounded;
  typedef typename TD::Around_point_circulator Around_point_circulator;
  typedef typename TD::In_face_iterator In_face_iterator;
  friend class Trapezoidal_decomposition_2<Traits>;
  
#ifdef CGAL_PM_FRIEND_CLASS
  
  friend class Trapezoidal_decomposition_2<Traits>::Around_point_circulator;
  friend class Trapezoidal_decomposition_2<Traits>::In_face_iterator;
  
#endif
  
  typedef typename TD::Data_structure Data_structure;
  
  // implementation for trapezoids with two sides parallel to y axis
  /* X_trapezoid abbriviates Planar trapezoid with vertical parallel 
     left and right sides.
     Trapezoids are represented as two points called right and left and
     two curves called top and bottom. The points lie on the right and left
     boundaries of the trapezoid respectively and the curves bound the 
     trapezoid
     from above and below.
     There are degenerate trapezoids called infinite trapezoid; this happens 
     when one of the four sides degenerates into a point/X_curve at infinity.
     Trapezoids are created as active and become inactive when Remove() member
     function called.
     Each X_trapezoid has at most four neighbouring trapezoids.
  */
private:
	Boundary_type* ptr() const {return (Boundary_type*)(PTR);}
	
	
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
	
  void init_neighbours(pointer lb_=0,pointer lt_=0,pointer rb_=0,pointer rt_=0)
  {set_lb(lb_);set_lt(lt_);set_rb(rb_);set_rt(rt_);}
  void set_node(Data_structure* p) {node=p;
  
#ifdef CGAL_TD_DEBUG
  
  CGAL_assertion(!p || **p==*this);
  
#endif	
	
  }
  void set_left(const Point& p) 
  {ptr()->e0=p;ptr()->e4&=~CGAL_TRAPEZOIDAL_DECOMPOSITION_2_LEFT_UNBOUNDED;}
  void set_right(const Point& p) 
  {ptr()->e1=p;ptr()->e4&=~CGAL_TRAPEZOIDAL_DECOMPOSITION_2_RIGHT_UNBOUNDED;}
  void set_bottom(const X_curve& cv) 
  {ptr()->e2=cv;ptr()->e4&=~CGAL_TRAPEZOIDAL_DECOMPOSITION_2_BOTTOM_UNBOUNDED;}
  void set_top(const X_curve& cv) 
  {ptr()->e3=cv;ptr()->e4&=~CGAL_TRAPEZOIDAL_DECOMPOSITION_2_TOP_UNBOUNDED;}
  void set_left_unbounded() 
  {ptr()->e4|=CGAL_TRAPEZOIDAL_DECOMPOSITION_2_LEFT_UNBOUNDED;}
  void set_right_unbounded() 
  {ptr()->e4|=CGAL_TRAPEZOIDAL_DECOMPOSITION_2_RIGHT_UNBOUNDED;}
  void set_bottom_unbounded() 
  {ptr()->e4|=CGAL_TRAPEZOIDAL_DECOMPOSITION_2_BOTTOM_UNBOUNDED;}
  void set_top_unbounded() 
  {ptr()->e4|=CGAL_TRAPEZOIDAL_DECOMPOSITION_2_TOP_UNBOUNDED;}
  void set_lb(X_trapezoid* lb) {ptr()->e5=lb;}
  void set_lt(X_trapezoid* lt) {ptr()->e6=lt;}
  void set_rb(X_trapezoid* rb) {ptr()->e7=rb;}
  void set_rt(X_trapezoid* rt) {ptr()->e8=rt;}
public:
  Td_X_trapezoid(){
    PTR=new Boundary_type(
                          Traits::get_point_at_left_top_infinity(),
                          Traits::get_point_at_right_bottom_infinity(),
                          Traits::get_curve_at_infinity(),
                          Traits::get_curve_at_infinity(),
                          CGAL_TRAPEZOIDAL_DECOMPOSITION_2_TOTALLY_UNBOUNDED,
                          0,0,0,0);node=0;}
  Td_X_trapezoid(
    const Point& l ,const Point&r,const X_curve& b,const X_curve &t,
    unsigned char c=CGAL_TRAPEZOIDAL_DECOMPOSITION_2_BOUNDED,
    X_trapezoid *lb=0,X_trapezoid *lt=0,X_trapezoid *rb=0,X_trapezoid *rt=0,
    Data_structure* p=0){PTR=new Boundary_type(l,r,b,t,c,lb,lt,rb,rt);node=p;}
  Td_X_trapezoid(
    const Point* l ,const Point* r , const X_curve* b ,const X_curve *t,
    X_trapezoid *lb=0,X_trapezoid *lt=0,X_trapezoid *rb=0,X_trapezoid *rt=0,
    Data_structure* p=0)
  {
    PTR=new Boundary_type(
       l ? *l : Traits::get_point_at_left_top_infinity(),
       r ? *r : Traits::get_point_at_right_bottom_infinity(),
       b ? *b : Traits::get_curve_at_infinity(),
       t ? *t : Traits::get_curve_at_infinity(),
       (l ? 0 : CGAL_TRAPEZOIDAL_DECOMPOSITION_2_LEFT_UNBOUNDED) | 
       (r ? 0 : CGAL_TRAPEZOIDAL_DECOMPOSITION_2_RIGHT_UNBOUNDED) | 
       (b ? 0 : CGAL_TRAPEZOIDAL_DECOMPOSITION_2_BOTTOM_UNBOUNDED) | 
       (t ? 0 : CGAL_TRAPEZOIDAL_DECOMPOSITION_2_TOP_UNBOUNDED),
       lb,lt,rb,rt);node=p;
  }
  Td_X_trapezoid(const X_trapezoid& tr) : Handle(tr){node=tr.node;}
  ~Td_X_trapezoid() {}
  /*
    remark:
    operator= should not copy node (or otherwise update 
    Data_structure::replace)
  */
  X_trapezoid& operator=(const X_trapezoid& t2)
  {
    Handle::operator=(t2);
    return *this;
  }
  unsigned long id() const
  {
    return (unsigned long) PTR;
  }
  bool operator==(const X_trapezoid& t2) const
  {
    return identical(*this,t2);
  }
  bool operator!=(const X_trapezoid& t2) const
  {
    return !(operator==(t2));
  }
  const Point& left() const {
    return !is_left_unbounded() ? 
      ptr()->e0 : Traits::get_point_at_left_top_infinity();}
  const Point& right() const {
    return !is_right_unbounded() ? 
      ptr()->e1 : Traits::get_point_at_right_bottom_infinity();}
  // filters out the infinite case where at returns predefined dummy values
  const X_curve& bottom() const {
    return !is_bottom_unbounded() ?  
      ptr()->e2 : Traits::get_curve_at_infinity();}
  const X_curve& top() const {
    return !is_top_unbounded() ?	
      ptr()->e3 : Traits::get_curve_at_infinity();}
  unsigned char boundedness() const {return ptr()->e4;}
  bool is_left_unbounded() const {
    return (ptr()->e4&CGAL_TRAPEZOIDAL_DECOMPOSITION_2_LEFT_UNBOUNDED)!=0;}
  bool is_right_unbounded() const {
    return (ptr()->e4&CGAL_TRAPEZOIDAL_DECOMPOSITION_2_RIGHT_UNBOUNDED)!=0;}
  bool is_bottom_unbounded() const {
    return (ptr()->e4&CGAL_TRAPEZOIDAL_DECOMPOSITION_2_BOTTOM_UNBOUNDED)!=0;}
  bool is_top_unbounded() const {
    return (ptr()->e4&CGAL_TRAPEZOIDAL_DECOMPOSITION_2_TOP_UNBOUNDED)!=0;}
  bool is_unbounded() const
  {	
#ifdef CGAL_TD_DEBUG
		
    CGAL_assertion(is_active());
    
#endif
    
    return (ptr()->e4&CGAL_TRAPEZOIDAL_DECOMPOSITION_2_TOTALLY_UNBOUNDED)!=0;
  }
  pointer left_bottom_neighbour() const {return ptr()->e5;}
  pointer left_top_neighbour() const {return ptr()->e6;}
  pointer right_bottom_neighbour() const {return ptr()->e7;}
  pointer right_top_neighbour() const {return ptr()->e8;}
  Data_structure* get_node() const {return node;}
  bool is_active() const 
  {return right_bottom_neighbour()!=
     (pointer)CGAL_TRAPEZOIDAL_DECOMPOSITION_2_DELETE_SIGNATURE;}
  
  void remove(Data_structure* left=0)
  {
    
#ifndef CGAL_TD_DEBUG
    CGAL_warning(is_active());
#else
    CGAL_precondition(is_active());
#endif
    
    // mark trapezoid as deleted,
    set_rb((pointer)CGAL_TRAPEZOIDAL_DECOMPOSITION_2_DELETE_SIGNATURE);
    
#ifdef CGAL_TD_DEBUG
    CGAL_warning(node);
#endif
		
    // resets left son in data structure depending on input.
    if (left) node->set_left(*left);
  }								

  /* precondition:
     both trapezoidal are active and have the same
     bounding edges from above and below and the trapezoids are adjacent to one
     another with the first to the left
     postcondition:
     this trapezoid is the union of the old this trapezoid
     and the input trapezoids
  */

  void merge_trapezoid(X_trapezoid& right)
  {

#ifdef CGAL_TD_DEBUG
	  CGAL_assertion(!is_right_unbounded());
#endif

	bool right_unbounded = right.is_right_unbounded();
    *this=X_trapezoid(
                      !is_left_unbounded() ? &left() : 0,
                      !right_unbounded ? &right.right() : 0,
                      !is_bottom_unbounded() ? &bottom() : 0,
                      !is_top_unbounded() ? &top() : 0,
                      left_bottom_neighbour(),left_top_neighbour(),
                      right.right_bottom_neighbour(),
                      right.right_top_neighbour());
//	if (right_unbounded) set_right_unbounded();
    if (right_bottom_neighbour()) right_bottom_neighbour()->set_lb(this);
    if (right_top_neighbour()) right_top_neighbour()->set_lt(this);

#ifdef CGAL_TD_DEBUG
	  CGAL_assertion(is_right_unbounded()==right.is_right_unbounded());
#endif

  }
	
#ifdef CGAL_TD_DEBUG
	bool is_valid(const Traits* traits) const
	{
		typename Traits::Curve_point_status t;
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
			if (!is_left_unbounded() && !is_right_unbounded() && traits->point_is_left_low(right(),left()))
			{
				std::cerr << "\nthis=";
				write(std::cerr,*this,*traits,false) << std::flush;
				CGAL_warning(!traits->point_is_left_low(right(),left()));
				return false;
			}
			
			if (!is_bottom_unbounded())
			{
				if (is_left_unbounded() ||
					is_right_unbounded())
				{
					std::cerr << "\nthis=";
					write(std::cerr,*this,*traits,false) << std::flush;
					CGAL_warning(!(is_left_unbounded() ||is_right_unbounded()));
					return false;
				}
				t=traits->curve_get_point_status(bottom(),left());
				if (!(t==traits->ABOVE_CURVE || t==traits->ON_CURVE))
				{
					std::cerr << "\nthis=";
					write(std::cerr,*this,*traits,false) << std::flush;
					std::cerr << "\nt==" << t << std::flush;
					CGAL_warning(t==Traits::ABOVE_CURVE || t==Traits::ON_CURVE);
					return false;
				}
				t=traits->curve_get_point_status(bottom(),right());
				if (!(t==traits->ABOVE_CURVE || t==traits->ON_CURVE))
				{
					std::cerr << "\nthis=";
					write(std::cerr,*this,*traits,false) << std::flush;
					std::cerr << "\nt==" << t << std::flush;
					CGAL_warning(!(!(t==Traits::ABOVE_CURVE || t==Traits::ON_CURVE)));
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
				t=traits->curve_get_point_status(top(),left());
				if (!(t==traits->UNDER_CURVE || t==traits->ON_CURVE))
				{
					std::cerr << "\nthis=";
					write(std::cerr,*this,*traits,false) << std::flush;
					std::cerr << "\nt==" << t << std::flush;
					CGAL_warning(!(!(t==Traits::UNDER_CURVE || t==Traits::ON_CURVE)));
					return false;
				}
				t=traits->curve_get_point_status(top(),right());
				if (!(t==traits->UNDER_CURVE || t==traits->ON_CURVE))
				{
					std::cerr << "\nthis=";
					write(std::cerr,*this,*traits,false) << std::flush;
					std::cerr << "\nt==" << t << std::flush;
					CGAL_warning(!(!(t==Traits::UNDER_CURVE || t==Traits::ON_CURVE)));
					return false;
				}
			}
			if (!traits->is_degenerate(*this))
			{
				if(
					right_top_neighbour() && right_top_neighbour()->top()!=top() ||
					left_top_neighbour() && left_top_neighbour()->top() != top() ||
					right_bottom_neighbour() && right_bottom_neighbour()->bottom()!=bottom() ||
					left_bottom_neighbour() && left_bottom_neighbour()->bottom() != bottom() ||
					right_top_neighbour() && traits->is_degenerate(*right_top_neighbour()) ||
					left_top_neighbour() && traits->is_degenerate(*left_top_neighbour()) ||
					right_bottom_neighbour() && traits->is_degenerate(*right_bottom_neighbour()) ||
					left_bottom_neighbour() && traits->is_degenerate(*left_bottom_neighbour())
					)
				{
					std::cerr << "\nthis=";
					write(std::cerr,*this,*traits,false) << std::flush;
					CGAL_warning(!(right_top_neighbour() && right_top_neighbour()->top()!=top()));
					CGAL_warning(!(left_top_neighbour() && left_top_neighbour()->top() != top()));
					CGAL_warning(!(right_bottom_neighbour() && right_bottom_neighbour()->bottom()!=bottom()));
					CGAL_warning(!(left_bottom_neighbour() && left_bottom_neighbour()->bottom() != bottom()));
					CGAL_warning(!(right_top_neighbour() && traits->is_degenerate(*right_top_neighbour())));
					CGAL_warning(!(left_top_neighbour() && traits->is_degenerate(*left_top_neighbour())));
					CGAL_warning(!(right_bottom_neighbour() && traits->is_degenerate(*right_bottom_neighbour())));
					CGAL_warning(!(left_bottom_neighbour() && traits->is_degenerate(*left_bottom_neighbour())));
					return false;
				}
				if (
					right_top_neighbour()&&!right_top_neighbour()->is_active()||
					left_top_neighbour()&&!left_top_neighbour()->is_active()||
					right_bottom_neighbour()&&!right_bottom_neighbour()->is_active()||
					left_bottom_neighbour()&&!left_bottom_neighbour()->is_active()
					)
				{
					std::cerr << "\nleft=" << left() << " right=" << right() << " bottom=" << bottom() << " top=" << top() << std::flush;
					CGAL_warning(!(right_top_neighbour()&&!right_top_neighbour()->is_active()));
					CGAL_warning(!(left_top_neighbour()&&!left_top_neighbour()->is_active()));
					CGAL_warning(!(right_bottom_neighbour()&&!right_bottom_neighbour()->is_active()));
					CGAL_warning(!(left_bottom_neighbour()&&!left_bottom_neighbour()->is_active()));
					return false;
				}
			}
			else
			{
			/* if the trapezoid is degenerate, the left() and right()
			points should be on the top() and bottom() curves.
			In any case none of the geometric boundaries should 
				be unbounded */
				if (
					is_bottom_unbounded()||
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
				if(!(traits->curve_get_point_status(bottom(),left())==Traits::ON_CURVE))
				{
					std::cerr << "\nbottom()==" << bottom() << std::flush;
					std::cerr << "\nleft()==" << left() << std::flush;
					CGAL_warning(traits->curve_get_point_status(bottom(),left())==Traits::ON_CURVE);
					return false;
				}
				if(!(traits->curve_get_point_status(bottom(),right())==Traits::ON_CURVE))
				{
					std::cerr << "\nbottom()==" << bottom() << std::flush;
					std::cerr << "\nright()==" << right() << std::flush;
					CGAL_warning(traits->curve_get_point_status(bottom(),right())==Traits::ON_CURVE);
					return false;
				}
				if(!(traits->curve_get_point_status(top(),left())==Traits::ON_CURVE))
				{
					std::cerr << "\ntop()==" << top() << std::flush;
					std::cerr << "\nleft()==" << left() << std::flush;
					CGAL_warning(traits->curve_get_point_status(top(),left())==Traits::ON_CURVE);
					return false;
				}
				if(!(traits->curve_get_point_status(top(),right())==Traits::ON_CURVE))
				{
					std::cerr << "\ntop()==" << top() << std::flush;
					std::cerr << "\nright()==" << right() << std::flush;
					CGAL_warning(traits->curve_get_point_status(top(),right())==Traits::ON_CURVE);
					return false;
				}
				if (traits->is_degenerate_curve(*this))
				{
					if(
						right_top_neighbour()&&!right_top_neighbour()->is_active()||
						//!left_top_neighbour()||!left_top_neighbour()->is_active()||
						right_bottom_neighbour()&&!right_bottom_neighbour()->is_active()||
						left_bottom_neighbour()&&!left_bottom_neighbour()->is_active()
						)
					{
						CGAL_warning(!right_top_neighbour()||right_top_neighbour()->is_active());
						//CGAL_warning(!left_top_neighbour()||left_top_neighbour()->is_active());
						CGAL_warning(!right_bottom_neighbour()||right_bottom_neighbour()->is_active());
						CGAL_warning(!left_bottom_neighbour()||left_bottom_neighbour()->is_active());
						return false;
					}
					if(
					/* if trapezoid is end relative to supporting X_curve, that is
					adjacent(trapezoid's right end point,supporting X_curve right end point) , right_top_neighbour() returns next such trapezoid around right() point in clockwise oriented order
					adjacent(trapezoid's left end point,supporting X_curve left end point), left_bottom_neighbour() returns next such trapezoid around left() point in clockwise oriented order */
					/* right_bottom_neighbour() points to next trapezoid on supporting X_curve, if such exist */
					right_top_neighbour()&&!traits->is_degenerate_curve(*right_top_neighbour())||
					//				!left_top_neighbour()||!traits->is_degenerate_curve(*left_top_neighbour())||
					right_bottom_neighbour()&&!traits->is_degenerate_curve(*right_bottom_neighbour())||
					left_bottom_neighbour()&&!traits->is_degenerate_curve(*left_bottom_neighbour())
					)
					{
						CGAL_warning(!right_top_neighbour()||traits->is_degenerate_curve(*right_top_neighbour()));
						//CGAL_warning(!left_top_neighbour()||!traits->is_degenerate_curve(*left_top_neighbour()));
						CGAL_warning(!right_bottom_neighbour()||traits->is_degenerate_curve(*right_bottom_neighbour()));
						CGAL_warning(!left_bottom_neighbour()||traits->is_degenerate_curve(*left_bottom_neighbour()));
						return false;
					}
				}
				else if (traits->is_degenerate_point(*this))
				{
					if(
						right_top_neighbour()&&!traits->is_degenerate_curve(*right_top_neighbour())||
						left_bottom_neighbour()&&!traits->is_degenerate_curve(*left_bottom_neighbour())
						)
					{
						CGAL_warning(!right_top_neighbour()||traits->is_degenerate_curve(*right_top_neighbour()));
						CGAL_warning(!left_bottom_neighbour()||traits->is_degenerate_curve(*left_bottom_neighbour()));
						return false;
					}
					if(
						right_top_neighbour()&&!right_top_neighbour()->is_active()||
						left_bottom_neighbour()&&!left_bottom_neighbour()->is_active()
						)
					{
						CGAL_warning(!right_top_neighbour()||right_top_neighbour()->is_active());
						CGAL_warning(!left_bottom_neighbour()||left_bottom_neighbour()->is_active());
						return false;
					}
					if (!traits->point_is_same(left(),right()))
					{
						std::cerr << "\nleft()==" << left() << std::flush;
						std::cerr << "\nright()==" << right() << std::flush;
						CGAL_warning(traits->point_is_same(left(),right()));
						return false;
					}
				}
			}
		}
		return true;
	}

	void debug() const // instantiate ptr functions.
	{
		ptr();
		bottom();
		top();
		left();
		right();
	}
#endif

};

CGAL_END_NAMESPACE

#endif // CGAL_TD_X_TRAPEZOID_H
