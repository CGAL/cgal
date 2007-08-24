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

#ifndef CGAL_VISIBILITY_COMPLEX_2_POINT_TRAITS_H
#define CGAL_VISIBILITY_COMPLEX_2_POINT_TRAITS_H

#include <CGAL/basic.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/Real_embeddable_traits.h>
#include <CGAL/Visibility_complex_2/Bitangent_2.h>
#include <CGAL/Visibility_complex_2/Arc_2.h>


CGAL_BEGIN_NAMESPACE
namespace Visibility_complex_2_details {
// -----------------------------------------------------------------------------
template<class R_,class DistanceNT, class RToDistanceNT> class Point_traits;

template<class R_,class DistanceNT, class RToDistanceNT> 
class Bitangent_2<Point_traits<R_,DistanceNT,RToDistanceNT> >
  : public Point_traits<R_,DistanceNT,RToDistanceNT>::Segment_2 ,
    public 
  Bitangent_base<typename Point_traits<R_,DistanceNT,RToDistanceNT>::Disk >
{
public:
  typedef Point_traits<R_,DistanceNT,RToDistanceNT> Gt;
  typedef typename Point_traits<R_,DistanceNT,RToDistanceNT>::Disk Disk;
private:
  typedef Bitangent_base<Disk>                Base;
public:
  // -------------------------------------------------------------------------
  typedef typename Gt::R                         R;
  typedef typename R::FT                          FT;
  typedef typename Gt::Arc_2                     Arc_2;
  typedef typename Gt::Segment_2                 Segment_2;

  typedef typename Gt::Point_2                   Point_2;
  typedef typename Base::Type                     Type;
  typedef typename Base::Disk_handle              Disk_handle;
  // Constructeurs -----------------------------------------------------------
  Bitangent_2() : Base() { }
  Bitangent_2(const Point_2& v1 , const Point_2& v2 , 
              Type t, Disk_handle start, Disk_handle finish)
    : Segment_2(v1,v2) , Base(t,start,finish) { }
  Bitangent_2(Type t, const Arc_2& source, const Arc_2& target) 
    : Segment_2(*source.object(),*target.object()) ,
      Base(t,source.object(),target.object()) { }
  Bitangent_2(Type t , Disk_handle o1 , Disk_handle o2) 
    : R_::Segment_2(*o1,*o2) , Base(t,o1,o2) { }

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
  //--------------------------------------------------------------------------
};


template <class R_,class DistanceNT, class RToDistanceNT>
class Point_traits
{
    typedef Point_traits<R_,DistanceNT,RToDistanceNT> Self;
public:
    // -------------------------------------------------------------------------
    typedef R_                            R;
    typedef typename R::FT                FT;
    typedef Point_2<R>           Point_2;
    typedef Segment_2<R>         Segment_2;
    typedef Point_2                       Disk;
    typedef Bitangent_2<Self> Bitangent_2;
    typedef Arc_2<Self>    Arc_2;

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
      Point_2 extreme_point(bool /*b*/, const Disk& c) const { return c; }
	Point_2 extreme_point(bool b, const Bitangent_2& c) const 
	{ return (b) ? c.source() : c.target(); }
	template < class C , class D >
	Comparison_result operator() (bool sa , const C& a,
				      bool sb , const D& b) const { 
	    Point_2 ap = extreme_point(sa,a);
	    Point_2 bp = extreme_point(sb,b);
	    Comparison_result res = R().compare_y_2_object()(ap,bp);
	    return (res == EQUAL) ? R().compare_x_2_object()(ap,bp) : res;
	}
    };
    // -------------------------------------------------------------------------
    struct Is_upward_directed {
	bool operator()(const Bitangent_2& b) const {
	  Comparison_result comp = R().compare_y_2_object()(b.source(), 
							    b.target());
	  comp = (comp == EQUAL) ? R().compare_x_2_object()(b.source(), 
							    b.target()) : comp;
	  return (comp != LARGER);
	}
    };
    // -------------------------------------------------------------------------
    // The chi3 predicate
    typedef Tag_true Supports_chi3;
    struct Orientation_infinite {
	Orientation operator() (const Bitangent_2& a, 
				const Disk& o) const
	{ return orientation(a.source(),a.target(),o); } 
      Orientation operator() (const Disk& /*o*/, 
			      const Bitangent_2& /*b*/) const{ return COLLINEAR; }
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
	bool operator() (const Disk& c) const { return true; }
    };

    typedef DistanceNT Distance_NT;

    struct Make_convex_from_point {
      Disk operator () (const Point_2& p) const { return p; }
    };

    struct Length {
      Distance_NT operator()(const Arc_2&, 
                         const Bitangent_2&, 
                         const Bitangent_2&) const {
        return 0;
      }
      Distance_NT operator()(const Bitangent_2& b) const 
      { return Distance()(b.source(),b.target()); }
    };
    // -------------------------------------------------------------------------
private:
    // -------------------------------------------------------------------------
    typedef typename Algebraic_structure_traits<Distance_NT>::Sqrt Sqrt;

    struct Distance {
      Distance_NT operator () (const Point_2& p, const Point_2& q) const {
        return
          Sqrt()(
            RToDistanceNT()(R().compute_squared_distance_2_object()(p,q)));
      }
    };
};
// ----------------------------------------------------------------------------- 
}

template <class R_,class DistanceNT=double,
          class RToDistanceNT=
            typename Real_embeddable_traits<typename R_::FT>::To_double>
class Visibility_complex_2_point_traits
  :public 
  Visibility_complex_2_details::Point_traits<R_,DistanceNT,RToDistanceNT> {};


CGAL_END_NAMESPACE

#endif // VISIBILITY_COMPLEX_2_POINT_TRAITS_H
