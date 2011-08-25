// Copyright (c) 2000,2001  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
#include <CGAL/Intersection_traits_d.h>
#include <CGAL/Kernel_traits.h>

namespace CGAL {
namespace internal {
template <class R>
typename IT<R, typename R::Line_d, typename R::Line_d>::result_type
intersection(const typename R::Line_d& l1, const typename R::Line_d& l2, const R&)
{ 
  typedef typename IT<R, typename R::Line_d, typename R::Line_d>
    ::result_type result_type;
  typedef typename R::Line_d_Line_d_pair ll_pair;
  ll_pair LL(l1, l2);
  switch (LL.intersection_type()) {
    case ll_pair::NO_INTERSECTION:
    default: 
      return result_type();
    case ll_pair::POINT: {
      typename R::Point_d pt;
      LL.intersection(pt);
      return result_type(pt);
    }
    case ll_pair::LINE:
      return result_type(l1);
  }
  return result_type(); // never reached
}

template <class R>
typename IT<R, typename R::Ray_d, typename R::Ray_d>::result_type
intersection(const typename R::Ray_d& l1, const typename R::Ray_d& l2, const R&)
{ 
  typedef typename IT<R, typename R::Ray_d, typename R::Ray_d>
    ::result_type result_type;
  typedef typename R::Ray_d_Ray_d_pair ll_pair;
  ll_pair LL(l1,l2);
  switch (LL.intersection_type()) {
    case ll_pair::NO_INTERSECTION:
    default: 
      return result_type();
    case ll_pair::POINT: {
      typename R::Point_d p;
      LL.intersection(p);
      return result_type(p);
    }
    case ll_pair::RAY: {
      typename R::Ray_d r;
      LL.intersection(r);
      return result_type(r);
    }
    case ll_pair::SEGMENT: {
      typename R::Segment_d s;
      LL.intersection(s);
      return result_type(s);
    }
  }
  return result_type(); // never reached
}

template <class R>
typename IT<R, typename R::Segment_d, typename R::Segment_d>::result_type
intersection(const typename R::Segment_d& l1, const typename R::Segment_d& l2, const R&)
{
  typedef typename IT<R, typename R::Segment_d, typename R::Segment_d>
    ::result_type result_type;
  typedef typename R::Segment_d_Segment_d_pair ll_pair;
  ll_pair LL(l1,l2);
  switch (LL.intersection_type()) {
    case ll_pair::NO_INTERSECTION:
    default: 
      return result_type();
    case ll_pair::POINT: {
      typename R::Point_d p;
      LL.intersection(p);
      return result_type(p);
    }
    case ll_pair::SEGMENT: {
      typename R::Segment_d s;
      LL.intersection(s);
      return result_type(s);
    }
  }
  return result_type(); // never reached
}

template <class R>
typename IT<R, typename R::Line_d, typename R::Ray_d>::result_type
intersection(const typename R::Line_d& l, const typename R::Ray_d& r, const R&)
{
  typedef typename IT<R, typename R::Line_d, typename R::Ray_d>
    ::result_type result_type;
  typedef typename R::Line_d_Ray_d_pair lr_pair;
  lr_pair LR(l,r);
  switch (LR.intersection_type()) {
    case lr_pair::NO_INTERSECTION:
    default:
        return result_type();
    case lr_pair::POINT: {
        typename R::Point_d pt;
        LR.intersection(pt);
        return result_type(pt);
    }
    case lr_pair::RAY: {
        return result_type(r);
    }
  }
  return result_type(); // never reached
}

template <class R>
typename IT<R, typename R::Ray_d, typename R::Line_d>::result_type
intersection(const typename R::Ray_d& r, const typename R::Line_d& l, const R& k)
{ return intersection(l,r,k); }

template <class R>
typename IT<R, typename R::Ray_d, typename R::Segment_d>::result_type
intersection(const typename R::Ray_d& r, const typename R::Segment_d& s, const R&)
{
  typedef typename IT<R, typename R::Ray_d, typename R::Segment_d>
    ::result_type result_type;
  typedef typename R::Ray_d_Segment_d_pair rs_pair;
  rs_pair RS(r,s);
  switch (RS.intersection_type()) {
    case rs_pair::NO_INTERSECTION:
    default:
        return result_type();
    case rs_pair::POINT: {
        typename R::Point_d pt;
        RS.intersection(pt);
        return result_type(pt);
    }
    case rs_pair::SEGMENT: {
        typename R::Segment_d st;
        RS.intersection(st);
        return result_type(st);
    }
  }
  return result_type(); // never reached
}

template <class R>
typename IT<R, typename R::Segment_d, typename R::Ray_d>::result_type
intersection(const typename R::Segment_d& s, const typename R::Ray_d& r, const R& k)
{ return intersection(r,s, k); }

template <class R>
typename IT<R, typename R::Line_d, typename R::Segment_d>::result_type
intersection(const typename R::Line_d& l, const typename R::Segment_d& s, const R&)
{
  typedef typename IT<R, typename R::Line_d, typename R::Segment_d>
    ::result_type result_type;
  typedef typename R::Line_d_Segment_d_pair rs_pair;
  rs_pair RS(l,s);
  switch (RS.intersection_type()) {
    case rs_pair::NO_INTERSECTION:
    default:
        return result_type();
    case rs_pair::POINT: {
        typename R::Point_d pt;
        RS.intersection(pt);
        return result_type(pt);
    }
    case rs_pair::SEGMENT: {
        typename R::Segment_d st;
        RS.intersection(st);
        return result_type(st);
    }
  }
  return result_type(); // never reached
}

template <class R>
typename IT<R, typename R::Segment_d, typename R::Line_d>::result_type
intersection(const typename R::Segment_d& s, const typename R::Line_d& l, const R& r)
{ return intersection(l,s,r); }

template <class R>
typename IT<R, typename R::Line_d, typename R::Hyperplane_d>::result_type
intersection(const typename R::Line_d& l, const typename R::Hyperplane_d& h, const R&)
{
  typedef typename IT<R, typename R::Line_d, typename R::Hyperplane_d>
    ::result_type result_type;
  typedef typename R::Line_d_Hyperplane_d_pair lh_pair;
  lh_pair LH(l,h);
  switch (LH.intersection_type()) {
    case lh_pair::NO_INTERSECTION:
    default:
        return result_type();
    case lh_pair::POINT: {
        typename R::Point_d pt;
        LH.intersection(pt);
        return result_type(pt);
    }
    case lh_pair::LINE:
        return result_type(l);
  }
  return result_type(); // never reached
}

template <class R>
typename IT<R, typename R::Hyperplane_d, typename R::Line_d>::result_type
intersection(const typename R::Hyperplane_d& h, const typename R::Line_d& l, const R& r)
{ return intersection(l,h,r); }

template <class R>
typename IT<R, typename R::Ray_d, typename R::Hyperplane_d>::result_type
intersection(const typename R::Ray_d& r, const typename R::Hyperplane_d& h, const R&)
{
  typedef typename IT<R, typename R::Ray_d, typename R::Hyperplane_d>
    ::result_type result_type;
  typedef typename R::Ray_d_Hyperplane_d_pair rh_pair;
  rh_pair RH(r,h);
  switch (RH.intersection_type()) {
    case rh_pair::NO_INTERSECTION:
    default:
        return result_type();
    case rh_pair::POINT: {
        typename R::Point_d pt;
        RH.intersection(pt);
        return result_type(pt);
    }
    case rh_pair::RAY:
        return result_type(r);
  }
  return result_type(); // never reached
}

template <class R>
typename IT<R, typename R::Hyperplane_d, typename R::Ray_d>::result_type
intersection(const typename R::Hyperplane_d& h, const typename R::Ray_d& r, const R& k)
{ return intersection(r,h,k); }

template <class R>
typename IT<R, typename R::Segment_d, typename R::Hyperplane_d>::result_type
intersection(const typename R::Segment_d& s, const typename R::Hyperplane_d& h, const R&)
{
  typedef typename IT<R, typename R::Segment_d, typename R::Hyperplane_d>
    ::result_type result_type;
  typedef typename R::Segment_d_Hyperplane_d_pair sh_pair;
  sh_pair SH(s,h);
  switch (SH.intersection_type()) {
    case sh_pair::NO_INTERSECTION:
    default:
        return result_type();
    case sh_pair::POINT: {
        typename R::Point_d pt;
        SH.intersection(pt);
        return result_type(pt);
    }
    case sh_pair::SEGMENT:
        return result_type(s);
  }
  return result_type(); // never reached
}

template <class R>
typename IT<R, typename R::Hyperplane_d, typename R::Segment_d>::result_type
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
} //namespace CGAL
#endif //CGAL_INTERSECTIONS_D_H
