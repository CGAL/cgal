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

#ifndef CGAL_VISIBILITY_COMPLEX_2_SEGMENT_TRAITS_H
#define CGAL_VISIBILITY_COMPLEX_2_SEGMENT_TRAITS_H


#include <CGAL/basic.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Visibility_complex_2/Bitangent_2.h>
#include <CGAL/Visibility_complex_2/Arc_2.h>
#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/Real_embeddable_traits.h>

CGAL_BEGIN_NAMESPACE
namespace Visibility_complex_2_details {
// -----------------------------------------------------------------------------

template <class R_,class DistanceNT,class RToDistanceNT> class Segment_traits;

template <class R_,class DistanceNT,class RToDistanceNT>
class Bitangent_2<Segment_traits<R_,DistanceNT,RToDistanceNT> >
  : public Segment_traits<R_,DistanceNT,RToDistanceNT>::Segment_2, 
    public Bitangent_base<
             typename Segment_traits<R_,DistanceNT,RToDistanceNT>::Disk>
{
public:
  typedef Segment_traits<R_,DistanceNT,RToDistanceNT> Gt;
  typedef typename Gt::Disk                  Disk;
private:
  typedef Bitangent_base<Disk> Base;
  typedef Bitangent_2<Gt> Self;
  using Base::LL;
  using Base::LR;
  using Base::RL;
  using Base::RR;
public:
  // -------------------------------------------------------------------------
  typedef R_                          R;
  typedef typename R::FT              FT;
  typedef typename Gt::Segment_2       Segment_2;
  typedef typename Gt::Point_2         Point_2;
  typedef typename Base::Disk_handle           Disk_handle;
  typedef typename Gt::Arc_2 Arc_2;
  typedef typename Base::Type                  Type;
  // -------------------------------------------------------------------------

public:
  // Constructeurs -----------------------------------------------------------
  Bitangent_2() : Base() { }
  Bitangent_2(const Point_2& v1 , const Point_2& v2 , 
              Type t, Disk_handle start, Disk_handle finish)
    : Segment_2(v1,v2) , Base(t,start,finish) { }
  Bitangent_2(Type t ,  Disk_handle o1 , Disk_handle o2) {
    if (is_bitangent(t,o1->source(),o2->source(),
                     o1->target(),o2->target()))
      *this = Bitangent_2(o1->source(),o2->source(),t,o1,o2);
    else if (is_bitangent(t,o1->source(),o2->target(),
                          o1->target(),o2->source()))
      *this = Bitangent_2(o1->source(),o2->target(),t,o1,o2);
    else if (is_bitangent(t,o1->target(),o2->source(),
                          o1->source(),o2->target()))
      *this = Bitangent_2(o1->target(),o2->source(),t,o1,o2);
    else *this = Bitangent_2(o1->target(),o2->target(),t,o1,o2);
  }
  Bitangent_2(Type t, const Arc_2& source, const Arc_2& target)
  { *this = Bitangent_2(t,source.object(),target.object()); }

  Bitangent_2(const Bitangent_2&sibling,bool reverse,Type t) {
    if (reverse) {
      *this=Bitangent_2(Visibility_complex_2_details::reverse(t),
                        sibling.target_object(),
                        sibling.source_object());
    } else {
      *this=Bitangent_2(t,sibling.source_object(),
                        sibling.target_object());
    }
  }

  //--------------------------------------------------------------------------
  bool operator==(const Bitangent_2& b) const 
  { return Base::operator==(b); }
  bool operator!=(const Bitangent_2& b) const 
  { return Base::operator!=(b); }
  // -------------------------------------------------------------------------
private:
  // b = (p1 , p2) and q1 and q2 are the two other points respectively on
  // source and target object. Returns true if the bitangent is valid.
  bool is_bitangent(Type t , const Point_2& p1, const Point_2& p2,
                    const Point_2& q1, const Point_2& q2)
  {
    return ((t == this->LL 
             && (left_turn (p1,p2,q1) ||
                 (collinear(p1,p2,q1) && are_ordered_along_line(q1,p1,p2)))
             && (left_turn (p1,p2,q2) ||
                 (collinear(p1,p2,q2) && are_ordered_along_line(p1,p2,q2))))
            ||	(t == LR 
                 && (left_turn (p1,p2,q1) ||
                     (collinear(p1,p2,q1) && are_ordered_along_line(q1,p1,p2)))
                 && (right_turn(p1,p2,q2) ||
                     (collinear(p1,p2,q2) && are_ordered_along_line(p1,p2,q2))))
            ||	(t == RR 
                 && (right_turn(p1,p2,q1) ||
                     (collinear(p1,p2,q1) && are_ordered_along_line(q1,p1,p2)))
                 && (right_turn(p1,p2,q2) ||
                     (collinear(p1,p2,q2) && are_ordered_along_line(p1,p2,q2))))
            ||	(t == RL 
                 && (right_turn(p1,p2,q1) ||
                     (collinear(p1,p2,q1) && are_ordered_along_line(q1,p1,p2)))
                 && (left_turn (p1,p2,q2) ||
                     (collinear(p1,p2,q2) && are_ordered_along_line(p1,p2,q2)))));
  }
};


  template < class R_, class DistanceNT,class RToDistanceNT >
class Segment_traits
{
  typedef Segment_traits<R_,DistanceNT,RToDistanceNT> Self;
public:
    // -------------------------------------------------------------------------
    typedef R_                            R;
    typedef typename R::FT                FT;
    typedef Point_2<R_>          Point_2;
    typedef CGAL::Segment_2<R_>           Segment_2;
    typedef Segment_2                     Disk;

    typedef Bitangent_2<Self>  Bitangent_2;
    typedef Arc_2<Self>               Arc_2;

    // -------------------------------------------------------------------------
    // The chi2 predicate
    struct Orientation_object {
	Orientation operator()(const Bitangent_2& a,const Bitangent_2& b) const{ 
	  return R().orientation_2_object()(a.source() , a.target() ,
			     a.source() + (b.target() - b.source()));
	}	
    };
    // -------------------------------------------------------------------------
    // The two follwing give the chi2 predicate with a point at infinity
    struct Compare_extreme_yx {
	const Point_2& extreme_point(bool b, const Disk& c) const { 

	    Comparison_result comp = R().compare_y_2_object()(c.source(), c.target());
	    comp = (comp != EQUAL) ? comp : R().compare_x_2_object()(c.source(), c.target());

	    return ((b == true  && comp == SMALLER) || 
		    (b == false && comp == LARGER)) ?  c.source() : c.target() ;
	}
	const Point_2& extreme_point(bool b, const Bitangent_2& c) const 
	{ return (b) ? c.source() : c.target(); }
	template < class C_ , class D_ >
	Comparison_result operator() (bool sa , const C_& a,
				      bool sb , const D_& b) const { 
	    const Point_2& ap = extreme_point(sa,a);
	    const Point_2& bp = extreme_point(sb,b);

	    Comparison_result cr = R().compare_y_2_object()(ap,bp);
	    cr = (cr != EQUAL) ? cr : R().compare_x_2_object()(ap,bp);
	    return cr;
	    
	}
    };

    struct Is_upward_directed {
	bool operator()(const Bitangent_2& b) const {
	    Comparison_result comp = R().compare_y_2_object()(b.source(), b.target());
	  comp = (comp != EQUAL) ? comp : R().compare_x_2_object()(b.source(), b.target());
	    return (comp != LARGER);
	}
    };

    // The chi3 predicate
    typedef Tag_true Supports_chi3;
    struct Orientation_infinite {
	const Point_2& top(const Disk& c) const {

	  Comparison_result comp = R().compare_y_2_object()(c.source(), c.target());
	  comp = (comp != EQUAL) ? comp : R().compare_x_2_object()(c.source(), c.target());

	    return (comp == SMALLER) ? c.target() : c.source();
	}
	const Point_2& bot(const Disk& o) const 
	{ return (top(o) == o.source()) ? o.target() : o.source(); }
	Orientation operator() (const Bitangent_2& a, 
				const Disk& o) const
	{ return R().orientation_2_object()(a.source(),a.target(),top(o)); } 
	Orientation operator() (const Disk& o, 
				const Bitangent_2& b) const
	{ return R().orientation_2_object()(bot(o),top(o),b.target()); }
	Orientation operator() (const Bitangent_2& a, 
				const Bitangent_2& b) const
	{ return R().orientation_2_object()(a.source(),a.target(),b.target()); } 
    };
    // Detection of degenerate cases
    struct Equal_as_segments {
	bool operator() (const Bitangent_2& a, const Bitangent_2& b) const {
	    return (a.source() == b.source() && a.target() == b.target());
	}
    };
    struct Is_point {
	bool operator() (const Disk& c) const 
	{ return c.source() == c.target(); }
    };
    // ----------------------------------------------------------------------
    typedef DistanceNT Distance_NT;

    struct Make_convex_from_point {
      Disk operator () (const Point_2& p) const { return Disk(p,p); }
    };

    struct Length {
      Distance_NT operator () (const Bitangent_2& b) const 
	{ return Distance()(b.source(),b.target()); }
      Distance_NT operator() (const Arc_2& a, 
                      const Bitangent_2& inf, const Bitangent_2& sup) const 
      { 
	Point_2 p = (inf.source_object() == a.object()) ? inf.source() :
          inf.target();
	Point_2 q = (sup.source_object() == a.object()) ? sup.source() :
          sup.target();
	if (p == q) return 0;
	return Distance()(p,q);
      }

    };


private:

    typedef typename Algebraic_structure_traits<Distance_NT>::Sqrt Sqrt;

    struct Distance {
      Distance_NT operator () (const Point_2& p, const Point_2& q) const {
        return
          Sqrt()(
            RToDistanceNT()(R().compute_squared_distance_2_object()(p,q)));
      }
    };
};

}

template <class R_,class DistanceNT=double,
          class RToDistanceNT=
            typename Real_embeddable_traits<typename R_::FT>::To_double>
class Visibility_complex_2_segment_traits
  :public 
  Visibility_complex_2_details::Segment_traits<R_,DistanceNT,RToDistanceNT>
{};



CGAL_END_NAMESPACE

#endif // VISIBILITY_COMPLEX_2_SEGMENT_TRAITS_H

