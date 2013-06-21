// Copyright (c) 2000,2001  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Michael Seel

#ifndef CGAL_INTERSECTIONS_D_H
#define CGAL_INTERSECTIONS_D_H

#include <CGAL/basic.h>
#include <CGAL/Intersection_traits.h>

namespace CGAL {
namespace internal {

template <class R>
typename cpp11::result_of<typename R::Intersect_d(typename R::Line_d, typename R::Line_d)>::type
intersection(const typename R::Line_d& l1, const typename R::Line_d& l2, const R&)
{ 
  typedef typename R::Line_d_Line_d_pair ll_pair;
  ll_pair LL(l1, l2);
  switch (LL.intersection_type()) {
    case ll_pair::NO_INTERSECTION:
    default: 
      return intersection_return<typename R::Intersect_d, typename R::Line_d, typename R::Line_d>();
    case ll_pair::POINT: {
      typename R::Point_d pt;
      LL.intersection(pt);
      return intersection_return<typename R::Intersect_d, typename R::Line_d, typename R::Line_d>(pt);
    }
    case ll_pair::LINE:
      return intersection_return<typename R::Intersect_d, typename R::Line_d, typename R::Line_d>(l1);
  }
  return intersection_return<typename R::Intersect_d, typename R::Line_d, typename R::Line_d>(); // never reached
}

template <class R>
typename cpp11::result_of<typename R::Intersect_d(typename R::Ray_d, typename R::Ray_d)>::type
intersection(const typename R::Ray_d& l1, const typename R::Ray_d& l2, const R&)
{ 
  typedef typename R::Ray_d_Ray_d_pair ll_pair;
  ll_pair LL(l1,l2);
  switch (LL.intersection_type()) {
    case ll_pair::NO_INTERSECTION:
    default: 
      return intersection_return<typename R::Intersect_d, typename R::Ray_d, typename R::Ray_d>();
    case ll_pair::POINT: {
      typename R::Point_d p;
      LL.intersection(p);
      return intersection_return<typename R::Intersect_d, typename R::Ray_d, typename R::Ray_d>(p);
    }
    case ll_pair::RAY: {
      typename R::Ray_d r;
      LL.intersection(r);
      return intersection_return<typename R::Intersect_d, typename R::Ray_d, typename R::Ray_d>(r);
    }
    case ll_pair::SEGMENT: {
      typename R::Segment_d s;
      LL.intersection(s);
      return intersection_return<typename R::Intersect_d, typename R::Ray_d, typename R::Ray_d>(s);
    }
  }
  return intersection_return<typename R::Intersect_d, typename R::Ray_d, typename R::Ray_d>(); // never reached
}

template <class R>
typename cpp11::result_of<typename R::Intersect_d(typename R::Segment_d, typename R::Segment_d)>::type
intersection(const typename R::Segment_d& l1, const typename R::Segment_d& l2, const R&)
{
  typedef typename R::Segment_d_Segment_d_pair ll_pair;
  ll_pair LL(l1,l2);
  switch (LL.intersection_type()) {
    case ll_pair::NO_INTERSECTION:
    default: 
      return intersection_return<typename R::Intersect_d, typename R::Segment_d, typename R::Segment_d>();
    case ll_pair::POINT: {
      typename R::Point_d p;
      LL.intersection(p);
      return intersection_return<typename R::Intersect_d, typename R::Segment_d, typename R::Segment_d>(p);
    }
    case ll_pair::SEGMENT: {
      typename R::Segment_d s;
      LL.intersection(s);
      return intersection_return<typename R::Intersect_d, typename R::Segment_d, typename R::Segment_d>(s);
    }
  }
  return intersection_return<typename R::Intersect_d, typename R::Segment_d, typename R::Segment_d>(); // never reached
}

template <class R>
typename cpp11::result_of<typename R::Intersect_d(typename R::Line_d, typename R::Ray_d)>::type
intersection(const typename R::Line_d& l, const typename R::Ray_d& r, const R&)
{
  typedef typename R::Line_d_Ray_d_pair lr_pair;
  lr_pair LR(l,r);
  switch (LR.intersection_type()) {
    case lr_pair::NO_INTERSECTION:
    default:
        return intersection_return<typename R::Intersect_d, typename R::Line_d, typename R::Ray_d>();
    case lr_pair::POINT: {
        typename R::Point_d pt;
        LR.intersection(pt);
        return intersection_return<typename R::Intersect_d, typename R::Line_d, typename R::Ray_d>(pt);
    }
    case lr_pair::RAY: {
        return intersection_return<typename R::Intersect_d, typename R::Line_d, typename R::Ray_d>(r);
    }
  }
  return intersection_return<typename R::Intersect_d, typename R::Line_d, typename R::Ray_d>(); // never reached
}

template <class R>
typename cpp11::result_of<typename R::Intersect_d(typename R::Ray_d, typename R::Line_d)>::type
intersection(const typename R::Ray_d& r, const typename R::Line_d& l, const R& k)
{ return intersection(l,r,k); }

template <class R>
typename cpp11::result_of<typename R::Intersect_d(typename R::Ray_d, typename R::Segment_d)>::type
intersection(const typename R::Ray_d& r, const typename R::Segment_d& s, const R&)
{
  typedef typename R::Ray_d_Segment_d_pair rs_pair;
  rs_pair RS(r,s);
  switch (RS.intersection_type()) {
    case rs_pair::NO_INTERSECTION:
    default:
        return intersection_return<typename R::Intersect_d, typename R::Ray_d, typename R::Segment_d>();
    case rs_pair::POINT: {
        typename R::Point_d pt;
        RS.intersection(pt);
        return intersection_return<typename R::Intersect_d, typename R::Ray_d, typename R::Segment_d>(pt);
    }
    case rs_pair::SEGMENT: {
        typename R::Segment_d st;
        RS.intersection(st);
        return intersection_return<typename R::Intersect_d, typename R::Ray_d, typename R::Segment_d>(st);
    }
  }
  return intersection_return<typename R::Intersect_d, typename R::Ray_d, typename R::Segment_d>(); // never reached
}

template <class R>
typename cpp11::result_of<typename R::Intersect_d(typename R::Segment_d, typename R::Ray_d)>::type
intersection(const typename R::Segment_d& s, const typename R::Ray_d& r, const R& k)
{ return intersection(r,s, k); }

template <class R>
typename cpp11::result_of<typename R::Intersect_d(typename R::Line_d, typename R::Segment_d)>::type
intersection(const typename R::Line_d& l, const typename R::Segment_d& s, const R&)
{
  typedef typename R::Line_d_Segment_d_pair rs_pair;
  rs_pair RS(l,s);
  switch (RS.intersection_type()) {
    case rs_pair::NO_INTERSECTION:
    default:
        return intersection_return<typename R::Intersect_d, typename R::Line_d, typename R::Segment_d>();
    case rs_pair::POINT: {
        typename R::Point_d pt;
        RS.intersection(pt);
        return intersection_return<typename R::Intersect_d, typename R::Line_d, typename R::Segment_d>(pt);
    }
    case rs_pair::SEGMENT: {
        typename R::Segment_d st;
        RS.intersection(st);
        return intersection_return<typename R::Intersect_d, typename R::Line_d, typename R::Segment_d>(st);
    }
  }
  return intersection_return<typename R::Intersect_d, typename R::Line_d, typename R::Segment_d>(); // never reached
}

template <class R>
typename cpp11::result_of<typename R::Intersect_d(typename R::Segment_d, typename R::Line_d)>::type
intersection(const typename R::Segment_d& s, const typename R::Line_d& l, const R& r)
{ return intersection(l,s,r); }

template <class R>
typename cpp11::result_of<typename R::Intersect_d(typename R::Line_d, typename R::Hyperplane_d)>::type
intersection(const typename R::Line_d& l, const typename R::Hyperplane_d& h, const R&)
{
  typedef typename R::Line_d_Hyperplane_d_pair lh_pair;
  lh_pair LH(l,h);
  switch (LH.intersection_type()) {
    case lh_pair::NO_INTERSECTION:
    default:
        return intersection_return<typename R::Intersect_d, typename R::Line_d, typename R::Hyperplane_d>();
    case lh_pair::POINT: {
        typename R::Point_d pt;
        LH.intersection(pt);
        return intersection_return<typename R::Intersect_d, typename R::Line_d, typename R::Hyperplane_d>(pt);
    }
    case lh_pair::LINE:
        return intersection_return<typename R::Intersect_d, typename R::Line_d, typename R::Hyperplane_d>(l);
  }
  return intersection_return<typename R::Intersect_d, typename R::Line_d, typename R::Hyperplane_d>(); // never reached
}

template <class R>
typename cpp11::result_of<typename R::Intersect_d(typename R::Hyperplane_d, typename R::Line_d)>::type
intersection(const typename R::Hyperplane_d& h, const typename R::Line_d& l, const R& r)
{ return intersection(l,h,r); }

template <class R>
typename cpp11::result_of<typename R::Intersect_d(typename R::Ray_d, typename R::Hyperplane_d)>::type
intersection(const typename R::Ray_d& r, const typename R::Hyperplane_d& h, const R&)
{
  typedef typename R::Ray_d_Hyperplane_d_pair rh_pair;
  rh_pair RH(r,h);
  switch (RH.intersection_type()) {
    case rh_pair::NO_INTERSECTION:
    default:
        return intersection_return<typename R::Intersect_d, typename R::Ray_d, typename R::Hyperplane_d>();
    case rh_pair::POINT: {
        typename R::Point_d pt;
        RH.intersection(pt);
        return intersection_return<typename R::Intersect_d, typename R::Ray_d, typename R::Hyperplane_d>(pt);
    }
    case rh_pair::RAY:
        return intersection_return<typename R::Intersect_d, typename R::Ray_d, typename R::Hyperplane_d>(r);
  }
  return intersection_return<typename R::Intersect_d, typename R::Ray_d, typename R::Hyperplane_d>(); // never reached
}

template <class R>
typename cpp11::result_of<typename R::Intersect_d(typename R::Hyperplane_d, typename R::Ray_d)>::type
intersection(const typename R::Hyperplane_d& h, const typename R::Ray_d& r, const R& k)
{ return intersection(r,h,k); }

template <class R>
typename cpp11::result_of<typename R::Intersect_d(typename R::Segment_d, typename R::Hyperplane_d)>::type
intersection(const typename R::Segment_d& s, const typename R::Hyperplane_d& h, const R&)
{
  typedef typename R::Segment_d_Hyperplane_d_pair sh_pair;
  sh_pair SH(s,h);
  switch (SH.intersection_type()) {
    case sh_pair::NO_INTERSECTION:
    default:
        return intersection_return<typename R::Intersect_d, typename R::Segment_d, typename R::Hyperplane_d>();
    case sh_pair::POINT: {
        typename R::Point_d pt;
        SH.intersection(pt);
        return intersection_return<typename R::Intersect_d, typename R::Segment_d, typename R::Hyperplane_d>(pt);
    }
    case sh_pair::SEGMENT:
        return intersection_return<typename R::Intersect_d, typename R::Segment_d, typename R::Hyperplane_d>(s);
  }
  return intersection_return<typename R::Intersect_d, typename R::Segment_d, 
                             typename R::Hyperplane_d>(); // never reached
}

template <class R>
typename cpp11::result_of<typename R::Intersect_d(typename R::Hyperplane_d, typename R::Segment_d)>::type
intersection(const typename R::Hyperplane_d& h, const typename R::Segment_d& s, const R& r)
{ return intersection(s,h,r); }

template <class R>
inline bool do_intersect(const typename R::Line_d &l1, const typename R::Line_d &l2, const R&)
{ typedef typename R::Line_d_Line_d_pair ll_pair;
  ll_pair LL(l1,l2);
  return LL.intersection_type() != ll_pair::NO_INTERSECTION;
}

template <class R>
inline bool do_intersect(const typename R::Ray_d &l1, const typename R::Ray_d &l2, const R&)
{ typedef typename R::Ray_d_Ray_d_pair ll_pair;
  ll_pair LL(l1,l2);
  return LL.intersection_type() != ll_pair::NO_INTERSECTION;
}

template <class R>
inline bool do_intersect(const typename R::Segment_d &l1, const typename R::Segment_d &l2, const R&)
{ typedef typename R::Segment_d_Segment_d_pair ll_pair;
  ll_pair LL(l1,l2);
  return LL.intersection_type() != ll_pair::NO_INTERSECTION;
}

template <class R>
inline bool do_intersect(const typename R::Line_d& l, const typename R::Ray_d& r, const R&)
{ typedef typename R::Line_d_Ray_d_pair lr_pair;
  lr_pair LR(l,r);
  return LR.intersection_type() != lr_pair::NO_INTERSECTION;
}

template <class R>
inline bool do_intersect(const typename R::Ray_d& r, const typename R::Line_d& l, const R& k)
{ return do_intersect(l,r, k); }

template <class R>
inline bool do_intersect(const typename R::Line_d& l, const typename R::Segment_d& s, const R&)
{ typedef typename R::Line_d_Segment_d_pair ls_pair;
  ls_pair LS(l,s);
  return LS.intersection_type() != ls_pair::NO_INTERSECTION;
}

template <class R>
inline bool do_intersect(const typename R::Segment_d& s, const typename R::Line_d& l, const R& r)
{ return do_intersect(l,s, r); }

template <class R>
inline bool do_intersect(const typename R::Ray_d& r, const typename R::Segment_d& s, const R&)
{ typedef typename R::Ray_d_Segment_d_pair rs_pair;
  rs_pair RS(r,s);
  return RS.intersection_type() != rs_pair::NO_INTERSECTION;
}

template <class R>
inline bool do_intersect(const typename R::Segment_d& s, const typename R::Ray_d& r, const R& k)
{ return do_intersect(r,s,k); }

template <class R>
inline bool do_intersect(const typename R::Line_d& l, const typename R::Hyperplane_d& h, const R&)
{ typedef typename R::Line_d_Hyperplane_d_pair lh_pair;
  lh_pair LH(l,h);
  return LH.intersection_type() != lh_pair::NO_INTERSECTION;
}

template <class R>
inline bool do_intersect(const typename R::Hyperplane_d& h, const typename R::Line_d& l, const R&)
{ return do_intersect(l,h); }

template <class R>
inline bool do_intersect(const typename R::Ray_d& r, const typename R::Hyperplane_d& h, const R&)
{ typedef typename R::Ray_d_Hyperplane_d_pair rh_pair;
  rh_pair RH(r,h);
  return RH.intersection_type() != rh_pair::NO_INTERSECTION;
}

template <class R>
inline bool do_intersect(const typename R::Hyperplane_d& h, const typename R::Ray_d& r, const R& k)
{ return do_intersect(r,h,k); }

template <class R>
inline bool do_intersect(const typename R::Segment_d& s, const typename R::Hyperplane_d& h, const R&)
{ typedef typename R::Segment_d_Hyperplane_d_pair sh_pair;
  sh_pair SH(s,h);
  return SH.intersection_type() != sh_pair::NO_INTERSECTION;
}

template <class R>
inline bool do_intersect(const typename R::Hyperplane_d& h, const typename R::Segment_d& s, const R& r)
{ return do_intersect(s,h,r); }

} //namespace internal

template<typename T>
class Hyperplane_d;
template<typename T>
class Line_d;
template<typename T>
class Segment_d;
template<typename T>
class Ray_d;

// global intersection
CGAL_INTERSECTION_FUNCTION_SELF(Line_d, d)
CGAL_INTERSECTION_FUNCTION_SELF(Ray_d, d)
CGAL_INTERSECTION_FUNCTION_SELF(Segment_d, d)
CGAL_INTERSECTION_FUNCTION(Line_d, Ray_d, d)
CGAL_INTERSECTION_FUNCTION(Ray_d, Segment_d, d)
CGAL_INTERSECTION_FUNCTION(Line_d, Segment_d, d)
CGAL_INTERSECTION_FUNCTION(Line_d, Hyperplane_d, d)
CGAL_INTERSECTION_FUNCTION(Ray_d, Hyperplane_d, d)
CGAL_INTERSECTION_FUNCTION(Segment_d, Hyperplane_d, d)

// global do_intersect
CGAL_DO_INTERSECT_FUNCTION_SELF(Line_d, d)
CGAL_DO_INTERSECT_FUNCTION_SELF(Ray_d, d)
CGAL_DO_INTERSECT_FUNCTION_SELF(Segment_d, d)
CGAL_DO_INTERSECT_FUNCTION(Line_d, Ray_d, d)
CGAL_DO_INTERSECT_FUNCTION(Line_d, Segment_d, d)
CGAL_DO_INTERSECT_FUNCTION(Ray_d, Segment_d, d)
CGAL_DO_INTERSECT_FUNCTION(Line_d, Hyperplane_d, d)
CGAL_DO_INTERSECT_FUNCTION(Ray_d, Hyperplane_d, d)
CGAL_DO_INTERSECT_FUNCTION(Segment_d, Hyperplane_d, d)

} //namespace CGAL
#endif //CGAL_INTERSECTIONS_D_H
