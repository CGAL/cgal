// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision:  $
// release_date  : $CGAL_Date:  $
//
// file          : include/CGAL/squared_distance_2_1.h
// source        : sqdistance_2.fw
// author(s)     : Geert-Jan Giezeman
//                 Michel Hoffmann <hoffmann@inf.ethz.ch>
//                 Andreas Fabri <Andreas.Fabri@geometryfactory.com>
//
// coordinator   : Saarbruecken
//
// ============================================================================


#ifndef CGAL_SQUARED_DISTANCE_2_1_H
#define CGAL_SQUARED_DISTANCE_2_1_H

#include <CGAL/user_classes.h>


#include <CGAL/utils.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/enum.h>
#include <CGAL/wmult.h>
#include <CGAL/squared_distance_utils.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {
  
  template <class K>
  inline typename K::FT
  squared_distance(const typename CGAL_WRAP(K)::Point_2 & pt1,
		   const typename CGAL_WRAP(K)::Point_2 & pt2,
		   const K& k)
  {
    typename K::Vector_2 vec = k.construct_vector_2_object()(pt2, pt1);
    return (typename K::FT)(vec*vec);
  }

  template <class K>
  typename K::FT
  squared_distance(const typename CGAL_WRAP(K)::Point_2 &pt,
		   const typename CGAL_WRAP(K)::Line_2 &line,
		   const K&)
  {
    typedef typename K::RT RT;
    RT a = line.a();
    RT b = line.b();
    RT w = pt.hw();
    RT n = a*pt.hx() + b*pt.hy() + wmult((K*)0, line.c(), w);
    RT d = wmult((K*)0, RT(a*a+b*b), w, w);
    return K::make_FT(n*n, d);
  }
  
  template <class K>
  inline typename K::FT
  squared_distance(const typename CGAL_WRAP(K)::Line_2 &line, 
		   const typename CGAL_WRAP(K)::Point_2 &pt,
		   const K& k)
  {
    return CGALi::squared_distance(pt, line, k);
  }
  
  template <class K>
  extern typename K::FT
  squared_distance(const typename CGAL_WRAP(K)::Point_2 &pt,
		   const typename CGAL_WRAP(K)::Ray_2 &ray,
		   const K& k)
  {
    typedef typename K::Vector_2 Vector_2;
    Vector_2 diff = pt-ray.source();
    const Vector_2 &dir = ray.direction().vector();
    if (!is_acute_angle(dir,diff, k) )
      return (typename K::FT)(diff*diff);
    return CGALi::squared_distance(pt, ray.supporting_line(), k);
  }

  template <class K>
  inline typename K::FT
  squared_distance(const typename CGAL_WRAP(K)::Ray_2 &ray,
		   const typename CGAL_WRAP(K)::Point_2 &pt,
		   const K& k)
  {
    return CGALi::squared_distance(pt, ray, k);
  }

  template <class K>
  typename K::FT
  squared_distance(const typename CGAL_WRAP(K)::Point_2 &pt,
		   const typename CGAL_WRAP(K)::Segment_2 &seg,
		   const K& k)
  {
    typedef typename K::Vector_2 Vector_2;
    typedef typename K::RT RT;
    // assert that the segment is valid (non zero length).
    Vector_2 diff = pt-seg.source();
    Vector_2 segvec = seg.target()-seg.source();
    RT d = wdot(diff,segvec, k);
    if (d <= (RT)0)
      return (typename K::FT)(diff*diff);
    RT e = wdot(segvec,segvec, k);
    if (wmult((K*)0 ,d, segvec.hw()) > wmult((K*)0, e, diff.hw()))
      return CGALi::squared_distance(pt, seg.target(), k);
    return CGALi::squared_distance(pt, seg.supporting_line(), k);
  }

  template <class K>
  inline typename K::FT
  squared_distance(const typename CGAL_WRAP(K)::Segment_2 &seg,
		   const typename CGAL_WRAP(K)::Point_2 &pt,
		   const K& k)
  {
    return CGALi::squared_distance(pt, seg, k);
  }

  template <class K>
  typename K::FT
  squared_distance_parallel(const typename CGAL_WRAP(K)::Segment_2 &seg1,
			    const typename CGAL_WRAP(K)::Segment_2 &seg2,
			    const K& k)
  {
    typedef typename K::Vector_2 Vector_2;
    bool same_direction;
    const Vector_2 &dir1 = seg1.direction().vector();
    const Vector_2 &dir2 = seg2.direction().vector();
    if (CGAL_NTS abs(dir1.hx()) > CGAL_NTS abs(dir1.hy())) {
      same_direction = (CGAL_NTS sign(dir1.hx()) == CGAL_NTS sign(dir2.hx()));
    } else {
      same_direction = (CGAL_NTS sign(dir1.hy()) == CGAL_NTS sign(dir2.hy()));
    }
    if (same_direction) {
      if (!is_acute_angle(seg1.source(), seg1.target(), seg2.source(), k))
	return CGALi::squared_distance(seg1.target(), seg2.source(), k);
      if (!is_acute_angle(seg1.target(), seg1.source(), seg2.target(), k))
	return CGALi::squared_distance(seg1.source(), seg2.target(), k);
    } else {
      if (!is_acute_angle(seg1.source(), seg1.target(), seg2.target(), k))
	return CGALi::squared_distance(seg1.target(), seg2.target(), k);
      if (!is_acute_angle(seg1.target(), seg1.source(), seg2.source(), k))
	return CGALi::squared_distance(seg1.source(), seg2.source(), k);
    }
    return CGALi::squared_distance(seg2.source(), seg1.supporting_line(), k);
  }

  template <class K>
  inline typename K::RT 
  _distance_measure_sub(const typename K::RT &startwcross, 
			const typename K::RT &endwcross,
			const typename CGAL_WRAP(K)::Point_2 &start, 
			const typename CGAL_WRAP(K)::Point_2 &end)
  {
    return  CGAL_NTS abs(wmult((K*)0, startwcross, end.hw())) -
      CGAL_NTS abs(wmult((K*)0, endwcross, start.hw()));
  }

  template <class K>
  typename K::FT
  squared_distance(const typename CGAL_WRAP(K)::Segment_2 &seg1,
		   const typename CGAL_WRAP(K)::Segment_2 &seg2,
		   const K& k)
  {
    typedef typename K::RT RT;
    typedef typename K::FT FT;
    bool crossing1, crossing2;
    RT c1s, c1e, c2s, c2e;
    if (seg1.source() == seg1.target())
      return CGALi::squared_distance(seg1.source(), seg2, k);
    if (seg2.source() == seg2.target())
      return CGALi::squared_distance(seg2.source(), seg1, k);
    c1s = wcross(seg2.source(), seg2.target(), seg1.source(), k);
    c1e = wcross(seg2.source(), seg2.target(), seg1.target(), k);
    c2s = wcross(seg1.source(), seg1.target(), seg2.source(), k);
    c2e = wcross(seg1.source(), seg1.target(), seg2.target(), k);
    if (c1s < RT(0)) {
      crossing1 = (c1e >= RT(0));
    } else {
      if (c1e <= RT(0)) {
	if (c1s == RT(0) && c1e == RT(0))
	  return CGALi::squared_distance_parallel(seg1, seg2, k);
	crossing1 = true;
      } else {
	crossing1 = (c1s == RT(0));
      }
    }
    if (c2s < RT(0)) {
      crossing2 = (c2e >= RT(0));
    } else {
      if (c2e <= RT(0)) {
	if (c2s == RT(0) && c2e == RT(0))
	  return CGALi::squared_distance_parallel(seg1, seg2, k);
	crossing2 = true;
      } else {
	crossing2 = (c2s == RT(0));
      }
    }

    if (crossing1) {
      if (crossing2)
	return (FT)0;
      RT dm;
      dm = _distance_measure_sub<K>(c2s,c2e, seg2.source(), seg2.target());
      if (dm < RT(0)) {
	return CGALi::squared_distance(seg2.source(), seg1, k);
      } else {
	if (dm > RT(0)) {
	  return CGALi::squared_distance(seg2.target(), seg1, k);
	} else {
	  // parallel, should not happen (no crossing)
	  return CGALi::squared_distance_parallel(seg1, seg2, k);
	}
      }
    } else {
      if (crossing2) {
	RT dm;
	dm =
	  _distance_measure_sub<K>(c1s, c1e,seg1.source(),seg1.target());
	if (dm < RT(0)) {
	  return CGALi::squared_distance(seg1.source(), seg2, k);
	} else {
	  if (dm > RT(0)) {
	    return CGALi::squared_distance(seg1.target(), seg2, k);
	  } else {
	    // parallel, should not happen (no crossing)
	    return CGALi::squared_distance_parallel(seg1, seg2, k);
	  }
	}
      } else {

	FT min1, min2;
	RT dm = _distance_measure_sub<K>(
				      c1s, c1e, seg1.source(), seg1.target());
	if (dm == RT(0))
	  return CGALi::squared_distance_parallel(seg1, seg2, k);
	min1 = (dm < RT(0)) ?
	  CGALi::squared_distance(seg1.source(), seg2, k):
	  CGALi::squared_distance(seg1.target(), seg2, k);
	dm = _distance_measure_sub<K>(
				   c2s, c2e, seg2.source(), seg2.target());
	if (dm == RT(0))  // should not happen.
	  return CGALi::squared_distance_parallel(seg1, seg2, k);
	min2 = (dm < RT(0)) ?
	  CGALi::squared_distance(seg2.source(), seg1, k):
	  CGALi::squared_distance(seg2.target(), seg1, k);
	return (min1 < min2) ? min1 : min2;
      }
    }
  }

  template <class K>
  inline typename K::RT 
  _distance_measure_sub(const typename K::RT &startwcross, 
			const typename K::RT &endwcross,
			const typename CGAL_WRAP(K)::Vector_2 &start, 
			const typename CGAL_WRAP(K)::Vector_2 &end)
  {
    return  CGAL_NTS abs(wmult((K*)0, startwcross, end.hw())) -
      CGAL_NTS abs(wmult((K*)0, endwcross, start.hw()));
  }

  template <class K>
  typename K::FT
  squared_distance_parallel(const typename CGAL_WRAP(K)::Segment_2 &seg,
			    const typename CGAL_WRAP(K)::Ray_2 &ray,
			    const K& k)
  {
    typedef typename K::Vector_2 Vector_2;
    bool same_direction;
    const Vector_2 &dir1 = seg.direction().vector();
    const Vector_2 &dir2 = ray.direction().vector();
    if (CGAL_NTS abs(dir1.hx()) > CGAL_NTS abs(dir1.hy())) {
      same_direction = (CGAL_NTS sign(dir1.hx()) == CGAL_NTS sign(dir2.hx()));
    } else {
      same_direction = (CGAL_NTS sign(dir1.hy()) == CGAL_NTS sign(dir2.hy()));
    }
    if (same_direction) {
      if (!is_acute_angle(seg.source(), seg.target(), ray.source(), k))
	return CGALi::squared_distance(seg.target(), ray.source(), k);
    } else {
      if (!is_acute_angle(seg.target(), seg.source(), ray.source(), k))
	return CGALi::squared_distance(seg.source(), ray.source(), k);
    }
    return CGALi::squared_distance(ray.source(), seg.supporting_line(), k);
  }

  template <class K>
  typename K::FT
  squared_distance(const typename CGAL_WRAP(K)::Segment_2 &seg,
		   const typename CGAL_WRAP(K)::Ray_2 &ray,
		   const K& k)
  {
    typedef typename K::RT RT;
    typedef typename K::FT FT;
    typedef typename K::Vector_2 Vector_2;
    const Vector_2 &raydir = ray.direction().vector();
    Vector_2 startvec(seg.source()-ray.source());
    Vector_2 endvec(seg.target()-ray.source());
    typename K::Orientation_2 orientation;

    bool crossing1, crossing2;
    RT c1s, c1e;
    Orientation ray_s_side;
    if (seg.source() == seg.target())
      return CGALi::squared_distance(seg.source(), ray, k);
    c1s = wcross(raydir, startvec, k);
    c1e = wcross(raydir, endvec, k);
    if (c1s < RT(0)) {
      crossing1 = (c1e >= RT(0));
    } else {
      if (c1e <= RT(0)) {
	if (c1s == RT(0) && c1e == RT(0))
	  return CGALi::squared_distance_parallel(seg, ray, k);
	crossing1 = true;
      } else {
	crossing1 = (c1s == RT(0));
      }
    }
    ray_s_side = orientation(seg.source(), seg.target(), ray.source());
    switch (ray_s_side) {
    case LEFT_TURN:
      crossing2 = right_turn(seg.target()-seg.source(), raydir, k);
      break;
    case RIGHT_TURN:
      crossing2 = left_turn(seg.target()-seg.source(), raydir, k);
      break;
    case COLLINEAR:
      crossing2 = true;
      break;
    }

    if (crossing1) {
      if (crossing2)
	return FT(0);
      return CGALi::squared_distance(ray.source(), seg, k);
    } else {
      if (crossing2) {
	RT dm;
	dm = _distance_measure_sub<K>(c1s, c1e, startvec, endvec);
	if (dm < RT(0)) {
	  return CGALi::squared_distance(seg.source(), ray, k);
	} else {
	  if (dm > RT(0)) {
	    return CGALi::squared_distance(seg.target(), ray, k);
	  } else {
	    // parallel, should not happen (no crossing)
	    return CGALi::squared_distance_parallel(seg, ray, k);
	  }
	}
      } else {

	FT min1, min2;
	RT dm;
	dm = _distance_measure_sub<K>(c1s, c1e, startvec, endvec);
	if (dm == RT(0))
	  return CGALi::squared_distance_parallel(seg, ray, k);
	min1 = (dm < RT(0))
	  ? CGALi::squared_distance(seg.source(), ray, k)
	  : CGALi::squared_distance(seg.target(), ray, k);
	min2 = CGALi::squared_distance(ray.source(), seg, k);
	return (min1 < min2) ? min1 : min2;
      }
    }
  }

  template <class K>
  inline typename K::FT
  squared_distance(const typename CGAL_WRAP(K)::Ray_2 &ray,
		   const typename CGAL_WRAP(K)::Segment_2 &seg,
		   const K& k)
  {
    return CGALi::squared_distance(seg, ray, k);
  }

  template <class K>
  typename K::FT
  _sqd_to_line(const typename CGAL_WRAP(K)::Vector_2 &diff,
	       const typename K::RT & wcross, 
	       const typename CGAL_WRAP(K)::Vector_2 &dir )
  {
    typedef typename K::RT RT;
    typedef typename K::FT FT;
    RT numerator = wcross*wcross;
    RT denominator = wmult((K*)0, RT(wdot(dir,dir, K())),
			   diff.hw(), diff.hw());
    FT result = K::make_FT(numerator, denominator);
    return result;
  }

  template <class K>
  typename K::FT
  squared_distance(const typename CGAL_WRAP(K)::Segment_2 &seg,
		   const typename CGAL_WRAP(K)::Line_2 &line,
		   const K& k)
  {
    typedef typename K::RT RT;
    typedef typename K::FT FT;
    typedef typename K::Vector_2 Vector_2;
    typedef typename K::Point_2  Point_2;
    const Vector_2 &linedir = line.direction().vector();
    const Point_2 &linepoint = line.point();
    Vector_2 startvec(seg.source()-linepoint);
    Vector_2 endvec(seg.target()-linepoint);

    bool crossing1;
    RT c1s, c1e;
    if (seg.source() == seg.target())
      return CGALi::squared_distance(seg.source(), line, k);
    c1s = wcross(linedir, startvec, k);
    c1e = wcross(linedir, endvec, k);
    if (c1s < RT(0)) {
      crossing1 = (c1e >= RT(0));
    } else {
      if (c1e <= RT(0)) {
	crossing1 = true;
      } else {
	crossing1 = (c1s == RT(0));
      }
    }

    if (crossing1) {
      return (FT)0;
    } else {
      RT dm;
      dm = _distance_measure_sub<K>(c1s, c1e, startvec, endvec);
      if (dm <= RT(0)) {
	return _sqd_to_line<K>(startvec, c1s, linedir);
      } else {
	return _sqd_to_line<K>(endvec, c1e, linedir);
      }
    }
  }

  template <class K>
  inline typename K::FT
  squared_distance(const typename CGAL_WRAP(K)::Line_2 &line,
		   const typename CGAL_WRAP(K)::Segment_2 &seg,
		   const K& k)
  {
    return CGALi::squared_distance(seg, line, k);
  }

  template <class K>
  typename K::FT
  ray_ray_squared_distance_parallel(
    const typename CGAL_WRAP(K)::Vector_2 &ray1dir,
    const typename CGAL_WRAP(K)::Vector_2 &ray2dir,
    const typename CGAL_WRAP(K)::Vector_2 &from1to2,
    const K& k)
  {
    typedef typename K::RT RT;
    typedef typename K::FT FT;
    if (!is_acute_angle(ray1dir, from1to2, k)) {
      bool same_direction;
      if (CGAL_NTS abs(ray1dir.hx()) > CGAL_NTS abs(ray1dir.hy())) {
	same_direction =
	  (CGAL_NTS sign(ray1dir.hx()) == CGAL_NTS sign(ray2dir.hx()));
      } else {
	same_direction =
	  (CGAL_NTS sign(ray1dir.hy()) == CGAL_NTS sign(ray2dir.hy()));
      }
      if (!same_direction)
	return (typename K::FT)(from1to2*from1to2);
    }
    RT wcr, w;
    wcr = wcross(ray1dir, from1to2, k);
    w = from1to2.hw();
    return (typename K::FT)(FT(wcr*wcr)
			    / FT(wmult((K*)0, RT(wdot(ray1dir, ray1dir, k)), w, w)));
  }

  template <class K>
  typename K::FT
  squared_distance(const typename CGAL_WRAP(K)::Ray_2 &ray1,
		   const typename CGAL_WRAP(K)::Ray_2 &ray2,
		   const K& k)
  {
    typedef typename K::Vector_2 Vector_2;
    typedef typename K::FT FT;
    const Vector_2 &ray1dir = ray1.direction().vector();
    const Vector_2 &ray2dir = ray2.direction().vector();
    Vector_2 diffvec(ray2.source()-ray1.source());

    bool crossing1, crossing2;
    Orientation dirorder;
    dirorder = orientation(ray1dir, ray2dir, k);
    switch (dirorder) {
    case COUNTERCLOCKWISE:
      crossing1 = !clockwise(diffvec, ray2dir, k);
      crossing2 = !counterclockwise(ray1dir, diffvec, k);
      break;
    case CLOCKWISE:
      crossing1 = !counterclockwise(diffvec, ray2dir, k);
      crossing2 = !clockwise(ray1dir, diffvec, k);
      break;
    case COLLINEAR:
      return ray_ray_squared_distance_parallel(ray1dir,ray2dir,diffvec,k);
    }

    if (crossing1) {
      if (crossing2)
	return (FT)0;
      return CGALi::squared_distance(ray2.source(), ray1, k);
    } else {
      if (crossing2) {
	return CGALi::squared_distance(ray1.source(), ray2, k);
      } else {

	FT min1, min2;
	min1 = CGALi::squared_distance(ray1.source(), ray2, k);
	min2 = CGALi::squared_distance(ray2.source(), ray1, k);
	return (min1 < min2) ? min1 : min2;
      }
    }
  }
  
  template <class K>
  extern typename K::FT
  squared_distance(const typename CGAL_WRAP(K)::Line_2 &line,
		   const typename CGAL_WRAP(K)::Ray_2 &ray,
		   const K& k)
  {
    typedef typename K::FT FT;
    typedef typename K::Vector_2 Vector_2;
    Vector_2 normalvec(line.a(), line.b());
    Vector_2 diff = ray.source()-line.point();
    FT sign_dist = diff*normalvec;
    if (sign_dist < FT(0)) {
      if (is_acute_angle(normalvec, ray.direction().vector(), k) )
	return (FT)0;
    } else {
      if (is_obtuse_angle(normalvec, ray.direction().vector(), k) )
	return (FT)0;
    }
    return (typename K::FT)((sign_dist*sign_dist)/(normalvec*normalvec));
  }

  template <class K>
  inline typename K::FT
  squared_distance(const typename CGAL_WRAP(K)::Ray_2 &ray,
		   const typename CGAL_WRAP(K)::Line_2 &line,
		   const K& k)
  {
    return CGALi::squared_distance(line, ray, k);
  }

  template <class K>
  inline bool
  _are_parallel(const typename CGAL_WRAP(K)::Line_2 &line1,
		const typename CGAL_WRAP(K)::Line_2 &line2)
  {
    return line1.a()*line2.b() == line2.a()*line1.b();
  }

  template <class K>
  inline typename K::FT
  squared_distance(const typename CGAL_WRAP(K)::Line_2 &line1,
		   const typename CGAL_WRAP(K)::Line_2 &line2,
		   const K& k)
  {
    typedef typename K::FT FT;
    if (_are_parallel<K>(line1,line2))
      return CGALi::squared_distance(line1.point(), line2, k);
    else
      return (FT)0;
  }

  template <class K>
  extern void
  distance_index(int &ind,
		 const typename CGAL_WRAP(K)::Point_2 &pt,
		 const typename CGAL_WRAP(K)::Ray_2 &ray,
		 const K& k)
  {
    if (!is_acute_angle(ray.direction().vector(),pt-ray.source(), k)) {
      ind = 0;
      return;
    }
    ind = -1;
  }

  template <class K>
  extern void
  distance_index(int &ind,
		 const typename CGAL_WRAP(K)::Point_2 &pt,
		 const typename CGAL_WRAP(K)::Segment_2 &seg,
		 const K& k)
  {
    if (!is_acute_angle(seg.target(),seg.source(),pt, k)) {
      ind = 0;
      return;
    }
    if (!is_acute_angle(seg.source(),seg.target(),pt, k)) {
      ind = 1;
      return;
    }
    ind = -1;
  }

  template <class K>
  inline typename K::FT
  squared_distance_indexed(const typename CGAL_WRAP(K)::Point_2 &pt,
			   const typename CGAL_WRAP(K)::Ray_2 &ray, 
			   int ind,
			   const K& k)
  {
    if (ind == 0)
      return CGALi::squared_distance(pt, ray.source(), k);
    return CGALi::squared_distance(pt, ray.supporting_line(), k);
  }

  template <class K>
  inline typename K::FT
  squared_distance_indexed(const typename CGAL_WRAP(K)::Point_2 &pt,
			   const typename CGAL_WRAP(K)::Segment_2 &seg, 
			   int ind,
			   const K& k)
  {
    if (ind == 0)
      return CGALi::squared_distance(pt, seg.source(), k);
    if (ind == 1)
      return CGALi::squared_distance(pt, seg.target(), k);
    return CGALi::squared_distance(pt, seg.supporting_line(), k);
  }
  
} // namespace CGALi

template <class K>
inline typename K::FT
squared_distance(
    const Point_2<K> & pt1,
    const Point_2<K> & pt2)
{
  return CGALi::squared_distance(pt1, pt2, K());
}

template <class K>
class Squared_distance_to_line {
    typename K::RT  a, b, c, sqnorm;
  public:
    Squared_distance_to_line(typename K::Line_2 const &line)
    : a(line.a()), b(line.b()), c(line.c())
    {
        sqnorm = a*a+b*b;
    }
    typename K::FT operator()(typename K::Point_2 const &pt) const
    {
        typedef typename K::RT RT;
        RT w = pt.hw();
        RT n = a*pt.hx() + b*pt.hy() + wmult((K*)0, c, w);
        RT d = wmult((K*)0, sqnorm, w, w);
        return K::make_FT(n*n, d);
    }
};

template <class K>
inline typename K::FT
squared_distance(
    const Point_2<K> &pt,
    const Line_2<K> &line)
{
  return CGALi::squared_distance(pt, line, K());
}


template <class K>
inline typename K::FT
squared_distance(
    const Line_2<K> & line,
    const Point_2<K> & pt)
{
    return squared_distance(pt, line);
}

template <class K>
class Squared_distance_to_ray {
    typename K::Vector_2 ray_dir;
    typename K::Point_2 ray_source;
    Squared_distance_to_line<K> supline_dist;
  public:
    Squared_distance_to_ray(typename K::Ray_2 const &ray)
    : ray_dir(ray.direction().vector()),
      ray_source(ray.source()),
      supline_dist(ray.supporting_line())
    { }
    typename K::FT operator()(typename K::Point_2 const &pt) const
    {
        typename K::Vector_2 diff = pt-ray_source;
        if (! CGALi::is_acute_angle(ray_dir,diff, K()) )
            return (typename K::FT)(diff*diff);
        return supline_dist(pt);
    }
};



template <class K>
inline typename K::FT
squared_distance(
    const Point_2<K> &pt,
    const Ray_2<K> &ray)
{
  return CGALi::squared_distance(pt, ray, K());
}


template <class K>
inline typename K::FT
squared_distance(
    const Ray_2<K> & ray,
    const Point_2<K> & pt)
{
    return squared_distance(pt, ray);
}




template <class K>
class Squared_distance_to_segment {
    typename K::Point_2 seg_source, seg_target;
    Squared_distance_to_line<K> supline_dist;
    typename K::Vector_2 segvec;
    typename K::RT e;
  public:
    Squared_distance_to_segment(typename K::Segment_2 const &seg)
    : seg_source(seg.source()), seg_target(seg.target()),
      supline_dist(seg.supporting_line())
    {
        segvec = seg_target-seg_source;
        e = CGALi::wdot(segvec,segvec, K());
    }
    typename K::FT operator()(typename K::Point_2 const &pt) const
    {
        typedef typename K::RT RT;
        // assert that the segment is valid (non zero length).
        typename K::Vector_2 diff = pt-seg_source;
        RT d = CGALi::wdot(diff,segvec, K());
        if (d <= (RT)0)
            return (typename K::FT)(diff*diff);
        if (wmult((K*)0 ,d, segvec.hw()) > wmult((K*)0, e, diff.hw()))
            return CGALi::squared_distance(pt, seg_target, K());
        return supline_dist(pt);
    }
};


template <class K>
inline typename K::FT
squared_distance(
    const Point_2<K> &pt,
    const Segment_2<K> &seg)
{
  return CGALi::squared_distance(pt, seg, K());
}


template <class K>
inline typename K::FT
squared_distance(
    const Segment_2<K> & seg,
    const Point_2<K> & pt)
{
  return CGALi::squared_distance(pt, seg, K());
}


template <class K>
inline typename K::FT
squared_distance(
    const Segment_2<K> &seg1,
    const Segment_2<K> &seg2)
{
  return CGALi::squared_distance(seg1, seg2, K());
}

template <class K>
inline typename K::FT
squared_distance(
    const Segment_2<K> &seg,
    const Ray_2<K> &ray)
{
  return CGALi::squared_distance(seg, ray, K());
}

template <class K>
inline typename K::FT
squared_distance(
    const Ray_2<K> & ray,
    const Segment_2<K> & seg)
{
  return CGALi::squared_distance(seg, ray, K());
}

template <class K>
inline typename K::FT
squared_distance(
    const Segment_2<K> &seg,
    const Line_2<K> &line)
{
  return CGALi::squared_distance(seg, line, K());
}

template <class K>
inline typename K::FT
squared_distance(
    const Line_2<K> & line,
    const Segment_2<K> & seg)
{
  return CGALi::squared_distance(seg, line, K());
}

template <class K>
inline typename K::FT
squared_distance(
    const Ray_2<K> &ray1,
    const Ray_2<K> &ray2)
{
  return CGALi::squared_distance(ray1, ray2, K());
}

template <class K>
inline typename K::FT
squared_distance(
    const Line_2<K> &line,
    const Ray_2<K> &ray)
{
  return CGALi::squared_distance(line, ray, K());
}

template <class K>
inline typename K::FT
squared_distance(
    const Ray_2<K> & ray,
    const Line_2<K> & line)
{
  return CGALi::squared_distance(line, ray, K());
}

template <class K>
inline typename K::FT
squared_distance(
    const Line_2<K> &line1,
    const Line_2<K> &line2)
{
  return CGALi::squared_distance(line1, line2, K());
}

CGAL_END_NAMESPACE

#endif
