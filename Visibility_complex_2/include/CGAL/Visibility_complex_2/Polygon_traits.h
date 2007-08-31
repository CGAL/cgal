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

#ifndef CGAL_VISIBILITY_COMPLEX_2_POLYGON_TRAITS_H
#define CGAL_VISIBILITY_COMPLEX_2_POLYGON_TRAITS_H

#include <CGAL/basic.h>

#include <CGAL/Polygon_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h> 
#include <CGAL/Visibility_complex_2/Bitangent_2.h>
#include <CGAL/Visibility_complex_2/Arc_2.h>

CGAL_BEGIN_NAMESPACE
namespace Visibility_complex_2_details {
// -----------------------------------------------------------------------------

template <class R_,class Container_,class DistanceNT,class RToDistanceNT>
class Polygon_traits;

template <class R_,class Container_,class DistanceNT,class RToDistanceNT>
class Arc_2<Polygon_traits<R_,Container_,DistanceNT,RToDistanceNT> >
  : public Arc_base<
      typename Polygon_traits<R_,Container_,DistanceNT,RToDistanceNT>::Disk>
{
public:
  // -------------------------------------------------------------------------
  typedef Polygon_traits<R_,Container_,DistanceNT,RToDistanceNT> Gt;
private:
  typedef Arc_base<Polygon_2<R_,Container_> >     Base;
public:
  typedef typename Gt::Point_2 Point_2;
  typedef typename Gt::Disk                      Disk;
  typedef typename Base::Disk_handle              Disk_handle;
  typedef typename Gt::Bitangent_2               Bitangent_2;
  typedef typename Disk::Vertex_const_circulator  Vertex_const_iterator;
  typedef typename Disk::Vertex_const_circulator  Vertex_iterator;
  // -------------------------------------------------------------------------
  using Base::object;

public:
  // -------------------------------------------------------------------------
  Arc_2() : Base(0) { }
  Arc_2(Disk_handle P) 
    : Base(P) , first(Vertex_iterator()) , beyond(Vertex_iterator()) { }
  // -------------------------------------------------------------------------
  void split (Arc_2& tmp, const Bitangent_2& v)
  {
    Point_2 p=v.source_object()==object()?v.source():v.target();
    if (first == CGAL_CIRC_NULL) {
      first = object()->vertices_circulator();
      while (*first != p) ++first;
      beyond = first; //++beyond;
      tmp.first = first; tmp.beyond = beyond;
    }
    else {
      tmp.beyond = beyond;
      tmp.first  = first;
      while (*tmp.first != p) ++tmp.first;
      beyond = tmp.first; 
      ++beyond;
    }
  } 
  void split_cw(Arc_2& tmp, const Bitangent_2& v) 
  {
    Point_2 p=v.source_object()==object()?v.source():v.target();
    if (first == CGAL_CIRC_NULL) {
      first = object()->vertices_circulator();
      while (*first != p) ++first;
      beyond = first; //++beyond;
      tmp.first = first; tmp.beyond = beyond;
    }
    else {
      tmp.beyond = beyond;
      tmp.first  = first;
      while (*tmp.beyond != p) --tmp.beyond;
      first = tmp.beyond;
      ++tmp.beyond;
    }
  }
  void update_begin(const Bitangent_2& v) {
    Point_2 p=v.source_object()==object()?v.source():v.target();
    first = beyond; --first;
    while (*first != p) --first;
  }
  void update_end(const Bitangent_2& v) {
    Point_2 p=v.source_object()==object()?v.source():v.target();
    beyond = first;
    while (*beyond != p) ++beyond;
    ++beyond;
  }
  void join (Arc_2& y) { beyond  = y.beyond; }
  // -------------------------------------------------------------------------
  // these methods are specific to this specialization
  Vertex_iterator begin()          const { return first;  }
  Vertex_iterator end()            const { return beyond; }
  Vertex_iterator vertices_begin() const { return first;  }
  Vertex_iterator vertices_end()   const { return beyond; }
  // -------------------------------------------------------------------------
private:
  Vertex_iterator first;
  Vertex_iterator beyond;
};


template <class R_,class Container_,class DistanceNT,class RToDistanceNT>
class Bitangent_2<Polygon_traits<R_,Container_,DistanceNT,RToDistanceNT> >
  : public Polygon_traits<R_,Container_,DistanceNT,RToDistanceNT>::Segment_2 , 
    public Bitangent_base<
      typename Polygon_traits<R_,Container_,DistanceNT,RToDistanceNT>::Disk >
{
  // -------------------------------------------------------------------------
public:
  typedef Polygon_traits<R_,Container_,DistanceNT,RToDistanceNT> Gt;
  typedef typename Gt::R R;
  typedef typename Gt::Disk                     Disk;
private:
  typedef Bitangent_base<Disk> Base;
  typedef Bitangent_2<Gt> Self;
public:
  typedef typename R::FT                          FT;
  typedef typename Gt::Arc_2                      Arc_2;
  typedef typename Base::Disk_handle              Disk_handle;
  typedef typename Gt::Segment_2                  Segment_2;
  typedef typename Gt::Point_2                    Point_2;
  typedef typename Base::Type                     Type;
  // -------------------------------------------------------------------------
  typedef typename Disk::Vertex_const_iterator Vertex_iterator;
  typedef typename Disk::Vertex_const_iterator Vertex_const_iterator;
  typedef typename Disk::Vertex_const_circulator Vertex_circulator;
  typedef typename Disk::Vertex_const_circulator Vertex_const_circulator;

  using Base::LL;
  using Base::LR;
  using Base::RL;
  using Base::RR;
  using Base::source_object;
  using Base::target_object;

  // Constructeurs -----------------------------------------------------------
  Bitangent_2() : Base() { }
  Bitangent_2(const Point_2& v1 , const Point_2& v2 , Type t ,
              Disk_handle start, Disk_handle finish)
    : Segment_2(v1,v2) , Base(t,start,finish) { }
  Bitangent_2(Type t, const Arc_2& source, const Arc_2& target) {
    CGAL_precondition(source.object()->orientation()!=CLOCKWISE);
    CGAL_precondition(target.object()->orientation()!=CLOCKWISE);


    if (source.begin() == source.end() || target.begin() == target.end()) {
      *this = Bitangent_2(t,source.object(),target.object());
      return;
    }
    
    bool t1,t2;

    Vertex_const_iterator it_source = source.begin();
    Vertex_const_iterator it_source_succ;
    Vertex_const_iterator it_target = target.begin();
    Vertex_const_iterator it_target_succ;

    Point_2 p_source(*it_source);
    Point_2 p_target(*it_target);
    Point_2 p_target_succ,p_source_succ;

    bool finished_source = false;
    bool finished_target = false;

    do {
      if (it_source == source.end()) finished_source = true;
      if (it_target == target.end()) finished_target = true;
      if (finished_target && finished_source) break;

      p_source = *(it_source); 
      p_target = *(it_target); 
	    
      it_source_succ = it_source; ++it_source_succ; 
      if (it_source_succ != source.end()) 
        p_source_succ = *(it_source_succ);
      else p_source_succ = p_source;

      it_target_succ = it_target; ++it_target_succ; 
      if (it_target_succ != target.end()) 
        p_target_succ = *(it_target_succ);
      else p_target_succ = p_target;
	    
      t1 = ( (t == RR || t == RL) && !left_turn (p_source,p_target,p_source_succ))
        || ( (t == LL || t == LR) && !right_turn(p_source,p_target,p_source_succ));
      if ( collinear(p_source,p_target,p_source_succ) && p_source != p_source_succ)
        t1 = are_ordered_along_line(p_source_succ,p_source,p_target);
	    
      t2 = ( (t == LL || t == RL) && !right_turn(p_source,p_target,p_target_succ))
        || ( (t == RR || t == LR) && !left_turn (p_source,p_target,p_target_succ));
      if ( collinear(p_source,p_target,p_target_succ) && p_target != p_target_succ)
        t2 = are_ordered_along_line(p_source,p_target,p_target_succ);

      if (finished_target)         ++it_source;
      else if (finished_source)    ++it_target;
      else if (t1 == 1 && t2 == 0) ++it_target;      
      else if (t1 == 0 && t2 == 1) ++it_source;
      else if (t1 == 0 && t2 == 0) {
        // Le test d'orientation est different suivant que 
        // source_is_left et target_is_left sont de meme 
        // signe ou non 
        Point_2 tmp(p_target_succ + (p_source_succ - p_target));
        if ( (    (t == LR || t == RL)    
                  && right_turn(p_source_succ, p_source, tmp))
             || ( (t == LL || t == RR)
                  && left_turn (p_source_succ, p_source, tmp)))
          ++it_target;
        else ++it_source;
      }
    } while (t1 == 0 || t2 == 0);
    *this = Bitangent_2(p_source,p_target,t,
                        source_object(),target_object());
  }
  Bitangent_2(Type t ,  Disk_handle o1 , Disk_handle o2){ 
    CGAL_precondition(o1->orientation()!=CLOCKWISE);
    CGAL_precondition(o2->orientation()!=CLOCKWISE);
    bool exists = true;
    bool found = false;
    bool found_start = false;
    bool found_finish = false;
    bool found_start_succ = false;
    bool found_start_pred = false;
    bool found_finish_succ = false;
    bool found_finish_pred = false;
    int  loop_counter = o1->size() + o2->size();

    Vertex_circulator start  = o1->vertices_circulator();
    Vertex_circulator finish = o2->vertices_circulator();
    do {
      Vertex_circulator start_succ  = start;  ++start_succ;
      Vertex_circulator finish_succ = finish; ++finish_succ; 
      Vertex_circulator start_pred  = start;  --start_pred;
      Vertex_circulator finish_pred = finish; --finish_pred;

      if (t == LL || t == LR) { 
        // --------------------------------------------------------------------
        found_start_succ = 
          (left_turn(*start,*finish,*start_succ)                      || 
           *start == *start_succ                                     ||
           (collinear(*start,*finish,*start_succ) &&
            are_ordered_along_line(*start_succ,*start,*finish)));
        // --------------------------------------------------------------------
        found_start_pred = 
          (left_turn(*start,*finish,*start_pred)                      || 
           *start == *start_pred                                     || 
           (collinear(*start,*finish,*start_pred) &&
            are_ordered_along_line(*start_pred,*start,*finish)));
        // --------------------------------------------------------------------
      }
      else   {
        // --------------------------------------------------------------------
        found_start_succ = 
          (right_turn(*start,*finish,*start_succ)                     || 
           *start == *start_succ                                     ||
           (collinear(*start,*finish,*start_succ) &&
            are_ordered_along_line(*start_succ,*start,*finish)));
        // --------------------------------------------------------------------
        found_start_pred = 
          (right_turn(*start,*finish,*start_pred)                     || 
           *start == *start_pred                                     || 
           (collinear(*start,*finish,*start_pred) &&
            are_ordered_along_line(*start_pred,*start,*finish)));
        // --------------------------------------------------------------------
      }
      found_start = (found_start_pred && found_start_succ);

      if (t == LL || t == RL) {
        // --------------------------------------------------------------------
        found_finish_succ = 
          (right_turn(*finish,*start,*finish_succ)                    || 
           *finish == *finish_succ                                   ||
           (collinear(*finish,*start,*finish_succ) &&
            are_ordered_along_line(*start,*finish,*finish_succ)));
        // --------------------------------------------------------------------
        found_finish_pred = 
          (right_turn(*finish,*start,*finish_pred)                    || 
           *finish == *finish_pred                                   ||
           (collinear(*finish,*start,*finish_pred) && 
            are_ordered_along_line(*start,*finish,*finish_pred)));
        // --------------------------------------------------------------------
      }
      else  {
        // --------------------------------------------------------------------
        found_finish_succ = 
          (left_turn(*finish,*start,*finish_succ)                     || 
           *finish == *finish_succ                                   ||
           (collinear(*finish,*start,*finish_succ) &&
            are_ordered_along_line(*start,*finish,*finish_succ)));
        // --------------------------------------------------------------------
        found_finish_pred = 
          (left_turn(*finish,*start,*finish_pred)                     || 
           *finish == *finish_pred                                   ||
           (collinear(*finish,*start,*finish_pred) && 
            are_ordered_along_line(*start,*finish,*finish_pred)));
        // --------------------------------------------------------------------
      }
      found_finish = (found_finish_pred && found_finish_succ);

      /* If found == true we have found the bitangent */
      found = (found_start && found_finish);

      /* Otherwise we compute one more determinant to decide on which */
      /* object to advance */

      if (found == false) 
        {
          if (( (t == LR || t == RL) 
                && right_turn(*start_succ,*start,*finish_succ + (*start_succ - *finish)))||
              ((t == LL || t == RR) 
               && left_turn (*start_succ,*start,*finish_succ + (*start_succ - *finish))))
            {
              ++finish;
            }
          else if (((t == LR || t == RL) 
                    && left_turn (*start_succ,*start,*finish_succ + (*start_succ - *finish)))||
                   ((t == LL || t == RR) 
                    && right_turn(*start_succ,*start,*finish_succ + (*start_succ - *finish))))
            {
              ++start;
            }
          else if (found_start == true && found_finish == false) ++finish;
          else ++start;
        }
      --loop_counter; 
      if (loop_counter < -100) {exists = false; break; }
    } while (found == false);

    if (exists){ *this = Bitangent_2(*start,*finish,t,o1,o2); }
    else { *this = Bitangent_2(); }
  }
  Bitangent_2(const Bitangent_2&sibling,bool reverse,Type t) {
    if (sibling.type()==t) {
      if (reverse) {
        *this=Bitangent_2(sibling.target(),sibling.source(),
                          Base::reverse(t),
                          sibling.target_object(),sibling.source_object());
      } else {
        *this=Bitangent_2(sibling.source(),sibling.target(),
                          t,
                          sibling.source_object(),sibling.target_object());
      }
    } else {
      if (reverse) {
        *this=Bitangent_2(Base::reverse(t),
                          sibling.target_object(),
                          sibling.source_object());
      } else {
        *this=Bitangent_2(t,sibling.source_object(),
                          sibling.target_object());
      }        
    }
  }
  //--------------------------------------------------------------------------
  bool operator==(const Bitangent_2& b) const 
  { return Base::operator==(b); }
  bool operator!=(const Bitangent_2& b) const 
  { return Base::operator!=(b); }
  //--------------------------------------------------------------------------
};


template < class R_, class Container_,class DistanceNT,class RToDistanceNT>
class Polygon_traits
{
public:
    // -------------------------------------------------------------------------
    // Geometric types
    typedef R_                                R;
    typedef typename R::FT                    FT;
    typedef Container_                Polygon_container;

private:
    typedef Polygon_traits<R,Polygon_container,DistanceNT,RToDistanceNT> Self;
public:
    typedef CGAL::Point_2<R_> Point_2;
    typedef Visibility_complex_2_details::Bitangent_2<Self> Bitangent_2;
    typedef Visibility_complex_2_details::Arc_2<Self>       Arc_2;
    typedef CGAL::Segment_2<R>                              Segment_2;
    typedef Polygon_2<R,Polygon_container>    Disk;
    // -------------------------------------------------------------------------
public:
    // -------------------------------------------------------------------------
    // The chi2 predicate
    struct Orientation_object {
	Orientation operator()(const Bitangent_2& a,const Bitangent_2& b) const{ 
	    return R().orientation_2_object()(a.source() , a.target() ,
					      a.source() + (b.target() 
							    - b.source()));
	}	
    };
    // -------------------------------------------------------------------------
    // The two follwing give the chi2 predicate with a point at infinity
    struct Compare_extreme_yx {
	const Point_2& extreme_point(bool b, const Disk& c) const {
	    return (b) ? *c.bottom_vertex() : *c.top_vertex();
	}
	const Point_2& extreme_point(bool b, const Bitangent_2& c) const 
	{ return (b) ? c.source() : c.target(); }
	template < class C , class D >
	Comparison_result operator() (bool sa , const C& a,
				      bool sb , const D& b) const { 
	    const Point_2& ap = extreme_point(sa,a);
	    const Point_2& bp = extreme_point(sb,b);

	    Comparison_result cr = R().compare_y_2_object()(ap,bp);
	    cr = (cr != EQUAL) ? cr : R().compare_x_2_object()(ap,bp);
	    return cr;
	    
	}
    };
    // -------------------------------------------------------------------------
    struct Is_upward_directed {
	bool operator()(const Bitangent_2& b) const {
	  Comparison_result comp = R().compare_y_2_object()(b.source(), 
							    b.target());
	  comp = (comp != EQUAL) ? comp : R().compare_x_2_object()(b.source(), 
								   b.target());
	  return (comp != LARGER);
	}
    };
    // -------------------------------------------------------------------------
    // The chi3 predicate
    typedef Tag_true Supports_chi3;
    struct Orientation_infinite {
	Orientation operator() (const Bitangent_2& a, 
				const Disk& o) const{
	    return R().orientation_2_object()(a.source(),a.target(),
			       *top_vertex_2(o.vertices_begin(),
					     o.vertices_end()));
	} 
	Orientation operator() (const Disk& o, 
				const Bitangent_2& b) const{
	    return R().orientation_2_object()(*bottom_vertex_2(o.vertices_begin(),
					        o.vertices_end()),
			       *top_vertex_2   (o.vertices_begin(),
						o.vertices_end()),
			       b.target());
	} 
	Orientation operator() (const Bitangent_2& a, 
				const Bitangent_2& b) const
	{ return R().orientation_2_object()(a.source(),a.target(),b.target()); } 
    };
    // -------------------------------------------------------------------------
    // Detection of degenerate cases
    struct Equal_as_segments {
	bool operator() (const Bitangent_2& a, const Bitangent_2& b) const {
	    return (a.source() == b.source() && a.target() == b.target());
	}
    };
    struct Is_point {
	bool operator() (const Disk& c) const { return (c.size() == 1); }
    };

    typedef DistanceNT Distance_NT;


    struct Make_convex_from_point {
      Disk operator () (const Point_2& p) const {
	Disk c;
	c.push_back(p);
	return c;
      }
    };

    struct Length {
      Distance_NT operator () (const Bitangent_2& b) const 
	{ return Distance()(b.source(),b.target()); }
      Distance_NT operator() (const Arc_2& a, 
                      const Bitangent_2&, const Bitangent_2& sup) const {
	// Arc reduced to one point : length is 0
	if (*a.begin() == sup.source() || *a.begin() == sup.target()) 
          return 0;
	// Arc containing at least two points : length is > 0
        typedef typename Arc_2::Vertex_const_iterator Arc_const_iterator;
	Arc_const_iterator it  = a.begin();
	Arc_const_iterator its = it; ++its;
	Distance_NT l = 0;
	while (*its != sup.source() && *its != sup.target()) {
          l = l + Distance()(*it,*its);
          it = its; ++its;
	}
	l = l + Distance()(*it,*its);
	return l;
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

// ----------------------------------------------------------------------------- 
}

template <class R_,class Container_ = typename Polygon_2<R_>::Container,
  class DistanceNT=double,
  class RToDistanceNT=
  typename Real_embeddable_traits<typename R_::FT>::To_double>
class Visibility_complex_2_polygon_traits
  :public 
  Visibility_complex_2_details::
  Polygon_traits<R_,Container_,DistanceNT,RToDistanceNT>
{};


CGAL_END_NAMESPACE

#endif // CGAL_VISIBILITY_COMPLEX_2_POLYGON_TRAITS_H
