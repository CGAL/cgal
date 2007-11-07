// Copyright (c) 2001-2004  ENS of Paris (France).
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
// Author(s)     : Pierre Angelier, Michel Pocchiola

#ifndef CGAL_VISIBILITY_COMPLEX_2_FUNCTION_OBJECTS
#define CGAL_VISIBILITY_COMPLEX_2_FUNCTION_OBJECTS

CGAL_BEGIN_NAMESPACE
namespace Visibility_complex_2_details {
// -----------------------------------------------------------------------------

template <class Gtr_> struct Parallel_chi3 {
  typedef typename Gtr_::Bitangent_2 Bitangent_2;

  // decides whether a.target_object lies before b.source_object along the
  // supporting line of a and b, assuming that b.source_object lies after 
  // a.source_object
  bool order (const Bitangent_2&a,bool upa, const Bitangent_2&b) {
    typename Gtr_::Orientation_infinite chi2;
    typename Gtr_::Is_upward_directed upw;
    Bitangent_2 u(Bitangent_2::RL,b.source_object(),a.target_object());
    switch (chi2(a,u)) {
    case LEFT_TURN: return true;
    case COLLINEAR: 
      if ((upw(u)==upa)) return true; else return false;
    case RIGHT_TURN: return false;
    }
    CGAL_postcondition(false);
    return true;
  }

  Orientation operator () (const Bitangent_2&a,bool upa,
                           const Bitangent_2&b) {
    typename Gtr_::Orientation_object chi2;
    typename Gtr_::Is_upward_directed upw;
    CGAL_precondition(a.type()==b.type());
    if (a.is_left_xx()) {
      Bitangent_2 abll(Bitangent_2::LL,a.source_object(),b.source_object());
      Bitangent_2 ablr(abll,false,Bitangent_2::LR);
      Bitangent_2 ball;
      Bitangent_2 barl;
      Orientation ablla=chi2(abll,a);
      Orientation aablr;
      Orientation balla;
      switch (ablla) {
      case COLLINEAR:
        if (upa==upw(abll)) {
          // b is left_xx
          if (a.is_xx_left())
            if (order(a,upa,b)) return LEFT_TURN; else return RIGHT_TURN;
          else
            if (order(a,upa,b)) return RIGHT_TURN; else return LEFT_TURN;
        }
        // else -> carry on
      case LEFT_TURN:
        aablr=chi2(a,ablr);
        switch (aablr) {
        case LEFT_TURN:
        case COLLINEAR: 
          return RIGHT_TURN;
        case RIGHT_TURN:
          ball=Bitangent_2(abll,true,Bitangent_2::RR);
          balla=chi2(ball,a);
          switch (balla) {
          case LEFT_TURN: return LEFT_TURN;
          case COLLINEAR:
            return LEFT_TURN;
          case RIGHT_TURN:
            return RIGHT_TURN;
          }
        }
      case RIGHT_TURN:
        ball=Bitangent_2(abll,true,Bitangent_2::RR);
        balla=chi2(ball,a);
        switch (balla) {
        case LEFT_TURN: return LEFT_TURN;
        case COLLINEAR: return LEFT_TURN;
        case RIGHT_TURN:
          return RIGHT_TURN;
        }
      }
    } else {
      Bitangent_2 abrr(Bitangent_2::RR,a.source_object(),b.source_object());
      Bitangent_2 abrl(abrr,false,Bitangent_2::RL);
      Bitangent_2 barr;
      Bitangent_2 balr;
      Orientation abrra=chi2(abrr,a);
      Orientation aabrl;
      Orientation barra;
      switch (abrra) {
      case COLLINEAR:
        if (upa==upw(abrr)) {
          // b is right_xx
          if (a.is_xx_right())
            if (order(a,upa,b)) return RIGHT_TURN; else return LEFT_TURN;
          else
            if (order(a,upa,b)) return LEFT_TURN; else return RIGHT_TURN;
        }
        // else -> carry on
      case RIGHT_TURN:
        aabrl=chi2(a,abrl);
        switch (aabrl) {
        case RIGHT_TURN:
        case COLLINEAR: 
          return LEFT_TURN;
        case LEFT_TURN:
          barr=Bitangent_2(abrr,true,Bitangent_2::LL);
          barra=chi2(barr,a);
          switch (barra) {
          case RIGHT_TURN: 
          case COLLINEAR: 
            return RIGHT_TURN;
          case LEFT_TURN:
            return LEFT_TURN;
          }
        }
      case LEFT_TURN:
        barr=Bitangent_2(abrr,true,Bitangent_2::LL);
        barra=chi2(barr,a);
        switch (barra) {
        case RIGHT_TURN: 
        case COLLINEAR: 
          return RIGHT_TURN;
        case LEFT_TURN:
          return LEFT_TURN;
        }
      }
    }
    CGAL_postcondition(false);
    return COLLINEAR;
  }
};

// chi2 wrapper

template < class Gtr_ >
struct Less_bitangent {
    typedef typename Gtr_::Bitangent_2 Bitangent_2;
    bool operator()(const Bitangent_2& a, const Bitangent_2& b) const {

	if (a == b) return false;

	// If this and b have the same source and target object we sometimes
	// don't need to compute a determinant. This solves some degenerate
	// cases.
	if (a.source_object() == b.source_object() && 
	    a.target_object() == b.target_object()) {
	    if (a.is_right_left() || b.is_left_right()) return true;
	    if (a.is_left_right() || b.is_right_left()) return false;
	    // Case where the left-left and right-right bitangents are
	    // geometrically equal. Arbitrary choice ll < rr.
	    typename Gtr_::Equal_as_segments equal_as_segments;
	    if (equal_as_segments(a,b)) return (a.is_left_left());
	}

	// Chi2 predicate from geometric traits class
	typename Gtr_::Orientation_object orientation;
	Orientation chi2ab = orientation(a,b);

	// Non degenerate cases.
	if (chi2ab == RIGHT_TURN){ return false;}
	if (chi2ab == LEFT_TURN) { return true; }

	// Treatment of degenerate case when two bitangents are parallel.

	typename Gtr_::Is_upward_directed is_upward_directed;
	bool upa = is_upward_directed(a); 
	bool upb = is_upward_directed(b);

	// Treating the case where a and b do not have the same angle.
	// Their angle difference is \pi
	if (!upa && upb)       return false;
	else if (upa && !upb)  return true;
	// Now a and b have the same angle

	// a and b have distinct types
	if (a.type() != b.type()) {
	    if (a.is_right_left() || b.is_left_right())   return true;
	    if (a.is_left_right() || b.is_right_left())   return false;
	    return (a.is_left_left());
	}

	// Final case : a and b have the same angle and same type

        if (a.target_object() == b.source_object()) return a.is_xx_left();
        if (a.source_object() == b.target_object()) return b.is_xx_right();

        if (a.source_object()==b.source_object()) {
          Bitangent_2 u(Bitangent_2::LR,a.target_object(),b.target_object());
          Orientation uo=orientation(a,u);
          if (uo==COLLINEAR) {
            if (is_upward_directed(u)==upa) uo=LEFT_TURN; else uo=RIGHT_TURN;
          }
          if (uo==LEFT_TURN)
            return a.is_xx_left();
        }
        if (a.target_object()==b.target_object()) {
          Bitangent_2 u(Bitangent_2::LR,a.source_object(),b.source_object());
          Orientation uo=orientation(a,u);
          if (uo==COLLINEAR) {
            if (is_upward_directed(u)==upa) uo=LEFT_TURN; else uo=RIGHT_TURN;
          }
          if (uo==LEFT_TURN)
            return a.is_right_xx();
        }
        Parallel_chi3<Gtr_> chi3;
        Orientation op=chi3(a,upa,b);
        if (op == RIGHT_TURN)     return false;
        else if (op == LEFT_TURN) return true;
        CGAL_error();
        return true;
    }
};

template < class Gtr_ >
struct Greater_bitangent {
    typedef typename Gtr_::Bitangent_2 Bitangent_2;
    bool operator()(const Bitangent_2& a, const Bitangent_2& b) {
	if (a == b) return false;
	return Less_bitangent<Gtr_>()(b,a);
    }
};

// -----------------------------------------------------------------------------
// x-compare two objects|bitangents
template <class Gtr_>
class Less_object {
private:
    // -------------------------------------------------------------------------
    typedef typename Gtr_::Disk             Disk;
    typedef typename Gtr_::Bitangent_2          Bitangent_2;
  typedef Bitangent_2 BT;
    typedef typename Gtr_::Is_upward_directed   Is_upward_directed;
    typedef typename Gtr_::Compare_extreme_yx   Compare;
    typedef typename Gtr_::Orientation_infinite Orientation_infinite;
    // -------------------------------------------------------------------------
public:
    // -------------------------------------------------------------------------
    bool operator() (const Disk& o1, const Disk& o2) const {
	// ---------------------------------------------------------------------
	if (&o1 == &o2) return false;
	BT rl(BT::RL,&o1,&o2);
	return !Is_upward_directed()(rl);
	// ---------------------------------------------------------------------
    }                      
    // -------------------------------------------------------------------------
    bool operator() (const Disk& o, const Bitangent_2& b) const {
	// ---------------------------------------------------------------------
	Is_upward_directed is_upward_directed;
	CGAL_precondition(is_upward_directed(b));
	// ---------------------------------------------------------------------
	// o is the source or target object of b.
	if (b.source_object() == &o) return b.is_left_xx();
	if (b.target_object() == &o) return b.is_xx_left();
	// ---------------------------------------------------------------------
	// Otherwise o and b involve three convex objects.
	Comparison_result cmp = Compare()(false,o,false,b);
	CGAL_precondition(cmp != EQUAL);
	Orientation_infinite chi3;
	if (cmp == SMALLER) return (chi3(b,o) == LEFT_TURN);
	else                return (chi3(o,b) == RIGHT_TURN);
	// Find somtething better than the above not using chi3...
	// ---------------------------------------------------------------------
    }                      
    // -------------------------------------------------------------------------
    bool operator() (const Bitangent_2& b, const Disk& o) const 
    { return !operator()(o,b); }                      
    // -------------------------------------------------------------------------
    bool operator() (const Bitangent_2& b1, const Bitangent_2& b2) const {
	// ---------------------------------------------------------------------
	if (b1 == b2) return false;
	// ---------------------------------------------------------------------
	Is_upward_directed is_upward_directed;
	CGAL_precondition(is_upward_directed(b1));
	CGAL_precondition(is_upward_directed(b2));
	// ---------------------------------------------------------------------
	// b1 and b2 involve two objects
	if ((b1.source_object() == b2.source_object() &&
	     b1.target_object() == b2.target_object()) ||
	    (b1.source_object() == b2.target_object() &&
	     b1.target_object() == b2.source_object()))
	    return (b1.is_right_right() || b2.is_left_left());
	// ---------------------------------------------------------------------
	// Some cases when b1 and b2 involve three objects
	//if (b1.target_object() == b2.source_object()) return b2.is_left_xx();
	//if (b1.source_object() == b2.target_object()) return b1.is_right_xx();
	Less_bitangent<Gtr_> chi1;
	if (b1.target_object() == b2.target_object()) {
	    if (b1.is_xx_left() == b2.is_xx_left()) return chi1(b1,b2);
	    else                                    return b1.is_xx_right();
	}
	if (b1.source_object() == b2.source_object()) {
	    if (b1.is_left_xx() == b2.is_left_xx()) return chi1(b2,b1);
	    else                                    return b1.is_right_xx();
	}
	// ---------------------------------------------------------------------
	// The remaining cases.
	Comparison_result cmp = Compare()(false,b1,false,b2);
	CGAL_precondition(cmp != EQUAL);
	Orientation_infinite chi3;
	if (cmp == SMALLER) return (chi3(b2,b1) == LEFT_TURN);
	else                return (chi3(b1,b2) == RIGHT_TURN);
	// Find somtething better than the above not using chi3...
	// ---------------------------------------------------------------------
    }
    // -------------------------------------------------------------------------
};

// -----------------------------------------------------------------------------
// Used to compute the antichain : x-compare two faces at angle 0
template < class Vc_ >
struct Less_face_handle {
    template < class Face_handle > 
    bool operator() (const Face_handle& e1, const Face_handle& e2) const {
	if (e1 == e2 || e1->front_view() == 0) return false;
	if (e2->front_view() == 0) return true;

	typedef typename Vc_::Gt Gtr_; 
	typedef typename Vc_::Edge_handle Edge_handle;
	Edge_handle e = e1->front_view();
	Edge_handle f = e2->front_view();
	if (e->object() == 0 && f->object() == 0) 
	    return Less_object<Gtr_>()(*e->sup(),*f->sup());
	else if (e->object() == 0 && f->object() != 0) 
	    return Less_object<Gtr_>()(*e->sup(),*f->object());
	else if (e->object() != 0 && f->object() == 0)
	    return Less_object<Gtr_>()(*e->object(),*f->sup());
	return Less_object<Gtr_>()(*e->object(), *f->object());
    }                      
};

// Used to compute the antichain : y-compares the horizontal tangents
template < class Gtr_ >
struct Less_edge_handle {
    typedef bool result_type;
    template < class Edge_handle > 
    bool operator () (const Edge_handle& e1, const Edge_handle& e2) const 
    { 
        CGAL_precondition(e1->object()&&e2->object());

	// Trivial cases not needing geometry
	if (e1 == e2) {  return false; }
	if (e1->object() == e2->object()) 
	    return (e1->sign()  && !e2->sign());

	// Otherwise use predicate from geometric traits class
	// to compare lexicographically the base points
	typename Gtr_::Compare_extreme_yx compare;
	Comparison_result cmp=compare(e1->sign(),*e1->object(),
                                      e2->sign(),*e2->object());
	// distinct base points
	if (cmp != EQUAL) return (cmp == SMALLER) ? true : false;

	// equal base points and e1, e2 are regular edges
	return Less_object<Gtr_>()(*e1->object(),*e2->object());
    }
};
}
CGAL_END_NAMESPACE

#endif // VISIBILITY_COMPLEX_2_FUNCTION_OBJECTS
