// Copyright (c) 1998-2004
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Geert-Jan Giezeman
//                 Michel Hoffmann <hoffmann@inf.ethz.ch>
//                 Andreas Fabri <Andreas.Fabri@geometryfactory.com>


#ifndef CGAL_SQUARED_DISTANCE_2_1_H
#define CGAL_SQUARED_DISTANCE_2_1_H

#include <CGAL/user_classes.h>


#include <CGAL/kernel_assertions.h>
#include <CGAL/enum.h>
#include <CGAL/wmult.h>
#include <CGAL/squared_distance_utils.h>
#include <CGAL/Kernel/global_functions_2.h>

namespace CGAL {

namespace internal {

  template <class K>
  inline typename K::FT
  squared_distance(const typename K::Point_2 & pt1,
                   const typename K::Point_2 & pt2,
                   const K& k)
  {
    typename K::Vector_2 vec = k.construct_vector_2_object()(pt2, pt1);
    return (typename K::FT)k.compute_squared_length_2_object()(vec);
  }

  template <class K>
  typename K::FT
  squared_distance(const typename K::Point_2 &pt,
                   const typename K::Line_2 &line,
                   const K&,
                   const Homogeneous_tag&)
  {
    typedef typename K::RT RT;
    typedef typename K::FT FT;
    const RT & a = line.a();
    const RT & b = line.b();
    const RT & w = pt.hw();
    RT n = a*pt.hx() + b*pt.hy() + w * line.c();
    RT d = (CGAL_NTS square(a) + CGAL_NTS square(b)) * CGAL_NTS square(w);
    return Rational_traits<FT>().make_rational(CGAL_NTS square(n), d);
  }

  template <class K>
  typename K::FT
  squared_distance(const typename K::Point_2 &pt,
                   const typename K::Line_2 &line,
                   const K&,
                   const Cartesian_tag&)
  {
    typedef typename K::FT FT;
    const FT & a = line.a();
    const FT & b = line.b();
    FT n = a*pt.x() + b*pt.y() + line.c();
    FT d = CGAL_NTS square(a) + CGAL_NTS square(b);
    return CGAL_NTS square(n)/d;
  }


  template <class K>
  typename K::FT
  squared_distance(const typename K::Point_2 &pt,
                   const typename K::Line_2 &line,
                   const K& k)
  {
    typedef typename K::Kernel_tag Tag;
    Tag tag;
    return squared_distance(pt, line, k, tag);
  }

  template <class K>
  inline typename K::FT
  squared_distance(const typename K::Line_2 &line,
                   const typename K::Point_2 &pt,
                   const K& k)
  {
    return internal::squared_distance(pt, line, k);
  }

  template <class K>
  typename K::FT
  squared_distance(const typename K::Point_2 &pt,
                   const typename K::Ray_2 &ray,
                   const K& k)
  {
    typedef typename K::Vector_2 Vector_2;
    typename K::Construct_vector_2 construct_vector;
    Vector_2 diff = construct_vector(ray.source(), pt);
    const Vector_2 &dir = ray.direction().vector();
    if (!is_acute_angle(dir,diff, k) )
      return (typename K::FT)k.compute_squared_length_2_object()(diff);
    return internal::squared_distance(pt, ray.supporting_line(), k);
  }

  template <class K>
  inline typename K::FT
  squared_distance(const typename K::Ray_2 &ray,
                   const typename K::Point_2 &pt,
                   const K& k)
  {
    return internal::squared_distance(pt, ray, k);
  }

  template <class K>
  typename K::FT
  squared_distance(const typename K::Point_2 &pt,
                   const typename K::Segment_2 &seg,
                   const K& k)
  {
    typename K::Construct_vector_2 construct_vector;
    typedef typename K::Vector_2 Vector_2;
    typedef typename K::RT RT;
    // assert that the segment is valid (non zero length).
    Vector_2 diff = construct_vector(seg.source(), pt);
    Vector_2 segvec = construct_vector(seg.source(), seg.target());
    RT d = wdot(diff,segvec, k);
    if (d <= (RT)0)
      return (typename K::FT)k.compute_squared_length_2_object()(diff);
    RT e = wdot(segvec,segvec, k);
    if (wmult((K*)0 ,d, segvec.hw()) > wmult((K*)0, e, diff.hw()))
      return internal::squared_distance(pt, seg.target(), k);
    return internal::squared_distance(pt, seg.supporting_line(), k);
  }

  template <class K>
  inline typename K::FT
  squared_distance(const typename K::Segment_2 &seg,
                   const typename K::Point_2 &pt,
                   const K& k)
  {
    return internal::squared_distance(pt, seg, k);
  }

  template <class K>
  typename K::FT
  squared_distance_parallel(const typename K::Segment_2 &seg1,
                            const typename K::Segment_2 &seg2,
                            const K& k)
  {
    typedef typename K::Vector_2 Vector_2;
    const Vector_2 &dir1 = seg1.direction().vector();
    const Vector_2 &dir2 = seg2.direction().vector();
    if (same_direction(dir1, dir2, k)) {
      if (!is_acute_angle(seg1.source(), seg1.target(), seg2.source(), k))
        return internal::squared_distance(seg1.target(), seg2.source(), k);
      if (!is_acute_angle(seg1.target(), seg1.source(), seg2.target(), k))
        return internal::squared_distance(seg1.source(), seg2.target(), k);
    } else {
      if (!is_acute_angle(seg1.source(), seg1.target(), seg2.target(), k))
        return internal::squared_distance(seg1.target(), seg2.target(), k);
      if (!is_acute_angle(seg1.target(), seg1.source(), seg2.source(), k))
        return internal::squared_distance(seg1.source(), seg2.source(), k);
    }
    return internal::squared_distance(seg2.source(), seg1.supporting_line(), k);
  }

  template <class K>
  inline typename K::RT
  _distance_measure_sub(const typename K::RT &startwcross,
                        const typename K::RT &endwcross,
                        const typename K::Point_2 &start,
                        const typename K::Point_2 &end)
  {
    return  CGAL_NTS abs(wmult((K*)0, startwcross, end.hw())) -
      CGAL_NTS abs(wmult((K*)0, endwcross, start.hw()));
  }

  template <class K>
  typename K::FT
  squared_distance(const typename K::Segment_2 &seg1,
                   const typename K::Segment_2 &seg2,
                   const K& k,
                   const Cartesian_tag&)
  {
    typedef typename K::RT RT;
    typedef typename K::FT FT;
    bool crossing1, crossing2;
    RT c1s, c1e, c2s, c2e;
    if (seg1.source() == seg1.target())
      return internal::squared_distance(seg1.source(), seg2, k);
    if (seg2.source() == seg2.target())
      return internal::squared_distance(seg2.source(), seg1, k);

    Orientation o1s = orientation(seg2.source(), seg2.target(), seg1.source());
    Orientation o1e = orientation(seg2.source(), seg2.target(), seg1.target());
    if (o1s == RIGHT_TURN) {
      crossing1 = (o1e != RIGHT_TURN);
    } else {
      if (o1e != LEFT_TURN) {
        if (o1s == COLLINEAR && o1e == COLLINEAR)
          return internal::squared_distance_parallel(seg1, seg2, k);
        crossing1 = true;
      } else {
        crossing1 = (o1s == COLLINEAR);
      }
    }

    Orientation o2s = orientation(seg1.source(), seg1.target(), seg2.source());
    Orientation o2e = orientation(seg1.source(), seg1.target(), seg2.target());
    if (o2s == RIGHT_TURN) {
      crossing2 = (o2e != RIGHT_TURN);
    } else {
      if (o2e != LEFT_TURN) {
        if (o2s == COLLINEAR && o2e == COLLINEAR)
          return internal::squared_distance_parallel(seg1, seg2, k);
        crossing2 = true;
      } else {
        crossing2 = (o2s == COLLINEAR);
      }
    }

    if (crossing1) {
      if (crossing2)
        return (FT)0;

      c2s = CGAL::abs(wcross(seg1.source(), seg1.target(), seg2.source(), k));
      c2e = CGAL::abs(wcross(seg1.source(), seg1.target(), seg2.target(), k));
      Comparison_result dm = compare(c2s,c2e);

      if (dm == SMALLER) {
        return internal::squared_distance(seg2.source(), seg1, k);
      } else {
        if (dm == LARGER) {
          return internal::squared_distance(seg2.target(), seg1, k);
        } else {
          // parallel, should not happen (no crossing)
          return internal::squared_distance_parallel(seg1, seg2, k);
        }
      }
    } else {
      c1s = CGAL::abs(wcross(seg2.source(), seg2.target(), seg1.source(), k));
      c1e = CGAL::abs(wcross(seg2.source(), seg2.target(), seg1.target(), k));
      Comparison_result dm = compare(c1s,c1e);
      if (crossing2) {
        if (dm == SMALLER) {
          return internal::squared_distance(seg1.source(), seg2, k);
        } else {
          if (dm == LARGER) {
            return internal::squared_distance(seg1.target(), seg2, k);
          } else {
            // parallel, should not happen (no crossing)
            return internal::squared_distance_parallel(seg1, seg2, k);
          }
        }
      } else {
        FT min1, min2;

        if (dm == EQUAL)
          return internal::squared_distance_parallel(seg1, seg2, k);
        min1 = (dm == SMALLER) ?
                 internal::squared_distance(seg1.source(), seg2, k):
                 internal::squared_distance(seg1.target(), seg2, k);

        c2s = CGAL::abs(wcross(seg1.source(), seg1.target(), seg2.source(), k));
        c2e = CGAL::abs(wcross(seg1.source(), seg1.target(), seg2.target(), k));
        dm = compare(c2s,c2e);

        if (dm == EQUAL)  // should not happen.
          return internal::squared_distance_parallel(seg1, seg2, k);
        min2 = (dm == SMALLER) ?
                 internal::squared_distance(seg2.source(), seg1, k):
                 internal::squared_distance(seg2.target(), seg1, k);
        return (min1 < min2) ? min1 : min2;
      }
    }
  }

  template <class K>
  typename K::FT
  squared_distance(const typename K::Segment_2 &seg1,
                   const typename K::Segment_2 &seg2,
                   const K& k,
                   const Homogeneous_tag&)
  {
    typedef typename K::RT RT;
    typedef typename K::FT FT;
    bool crossing1, crossing2;
    RT c1s, c1e, c2s, c2e;
    if (seg1.source() == seg1.target())
      return internal::squared_distance(seg1.source(), seg2, k);
    if (seg2.source() == seg2.target())
      return internal::squared_distance(seg2.source(), seg1, k);
    c1s = wcross(seg2.source(), seg2.target(), seg1.source(), k);
    c1e = wcross(seg2.source(), seg2.target(), seg1.target(), k);
    c2s = wcross(seg1.source(), seg1.target(), seg2.source(), k);
    c2e = wcross(seg1.source(), seg1.target(), seg2.target(), k);
    if (c1s < RT(0)) {
      crossing1 = (c1e >= RT(0));
    } else {
      if (c1e <= RT(0)) {
        if (c1s == RT(0) && c1e == RT(0))
          return internal::squared_distance_parallel(seg1, seg2, k);
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
          return internal::squared_distance_parallel(seg1, seg2, k);
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
        return internal::squared_distance(seg2.source(), seg1, k);
      } else {
        if (dm > RT(0)) {
          return internal::squared_distance(seg2.target(), seg1, k);
        } else {
          // parallel, should not happen (no crossing)
          return internal::squared_distance_parallel(seg1, seg2, k);
        }
      }
    } else {
      if (crossing2) {
        RT dm;
        dm =
          _distance_measure_sub<K>(c1s, c1e,seg1.source(),seg1.target());
        if (dm < RT(0)) {
          return internal::squared_distance(seg1.source(), seg2, k);
        } else {
          if (dm > RT(0)) {
            return internal::squared_distance(seg1.target(), seg2, k);
          } else {
            // parallel, should not happen (no crossing)
            return internal::squared_distance_parallel(seg1, seg2, k);
          }
        }
      } else {

        FT min1, min2;
        RT dm = _distance_measure_sub<K>(
                                      c1s, c1e, seg1.source(), seg1.target());
        if (dm == RT(0))
          return internal::squared_distance_parallel(seg1, seg2, k);
        min1 = (dm < RT(0)) ?
          internal::squared_distance(seg1.source(), seg2, k):
          internal::squared_distance(seg1.target(), seg2, k);
        dm = _distance_measure_sub<K>(
                                   c2s, c2e, seg2.source(), seg2.target());
        if (dm == RT(0))  // should not happen.
          return internal::squared_distance_parallel(seg1, seg2, k);
        min2 = (dm < RT(0)) ?
          internal::squared_distance(seg2.source(), seg1, k):
          internal::squared_distance(seg2.target(), seg1, k);
        return (min1 < min2) ? min1 : min2;
      }
    }
  }


  template <class K>
  inline typename K::RT
  _distance_measure_sub(const typename K::RT &startwcross,
                        const typename K::RT &endwcross,
                        const typename K::Vector_2 &start,
                        const typename K::Vector_2 &end)
  {
    return  CGAL_NTS abs(wmult((K*)0, startwcross, end.hw())) -
      CGAL_NTS abs(wmult((K*)0, endwcross, start.hw()));
  }

  template <class K>
  typename K::FT
  squared_distance_parallel(const typename K::Segment_2 &seg,
                            const typename K::Ray_2 &ray,
                            const K& k)
  {
    typedef typename K::Vector_2 Vector_2;
    const Vector_2 &dir1 = seg.direction().vector();
    const Vector_2 &dir2 = ray.direction().vector();

    if (same_direction(dir1, dir2, k)) {
      if (!is_acute_angle(seg.source(), seg.target(), ray.source(), k))
        return internal::squared_distance(seg.target(), ray.source(), k);
    } else {
      if (!is_acute_angle(seg.target(), seg.source(), ray.source(), k))
        return internal::squared_distance(seg.source(), ray.source(), k);
    }
    return internal::squared_distance(ray.source(), seg.supporting_line(), k);
  }

  template <class K>
  typename K::FT
  squared_distance(const typename K::Segment_2 &seg,
                   const typename K::Ray_2 &ray,
                   const K& k)
  {
    typename K::Construct_vector_2 construct_vector;
    typedef typename K::RT RT;
    typedef typename K::FT FT;
    typedef typename K::Vector_2 Vector_2;
    const Vector_2 &raydir = ray.direction().vector();
    Vector_2 startvec(construct_vector(ray.source(), seg.source()));
    Vector_2 endvec(construct_vector(ray.source(), seg.target()));
    typename K::Orientation_2 orientation;

    bool crossing1, crossing2;
    RT c1s, c1e;
    if (seg.source() == seg.target())
      return internal::squared_distance(seg.source(), ray, k);
    c1s = wcross(raydir, startvec, k);
    c1e = wcross(raydir, endvec, k);
    if (c1s < RT(0)) {
      crossing1 = (c1e >= RT(0));
    } else {
      if (c1e <= RT(0)) {
        if (c1s == RT(0) && c1e == RT(0))
          return internal::squared_distance_parallel(seg, ray, k);
        crossing1 = true;
      } else {
        crossing1 = (c1s == RT(0));
      }
    }
    switch (orientation(seg.source(), seg.target(), ray.source())) {
    case LEFT_TURN:
      crossing2 = right_turn(construct_vector(seg.source(), seg.target()), raydir, k);
      break;
    case RIGHT_TURN:
      crossing2 = left_turn(construct_vector(seg.source(), seg.target()), raydir, k);
      break;
    default:
      crossing2 = true;
      break;
    }

    if (crossing1) {
      if (crossing2)
        return FT(0);
      return internal::squared_distance(ray.source(), seg, k);
    } else {
      if (crossing2) {
        RT dm;
        dm = _distance_measure_sub<K>(c1s, c1e, startvec, endvec);
        if (dm < RT(0)) {
          return internal::squared_distance(seg.source(), ray, k);
        } else {
          if (dm > RT(0)) {
            return internal::squared_distance(seg.target(), ray, k);
          } else {
            // parallel, should not happen (no crossing)
            return internal::squared_distance_parallel(seg, ray, k);
          }
        }
      } else {
        FT min1, min2;
        RT dm = _distance_measure_sub<K>(c1s, c1e, startvec, endvec);
        if (dm == RT(0))
          return internal::squared_distance_parallel(seg, ray, k);
        min1 = (dm < RT(0))
          ? internal::squared_distance(seg.source(), ray, k)
          : internal::squared_distance(seg.target(), ray, k);
        min2 = internal::squared_distance(ray.source(), seg, k);
        return (min1 < min2) ? min1 : min2;
      }
    }
  }

  template <class K>
  inline typename K::FT
  squared_distance(const typename K::Ray_2 &ray,
                   const typename K::Segment_2 &seg,
                   const K& k)
  {
    return internal::squared_distance(seg, ray, k);
  }

  template <class K>
  typename K::FT
  _sqd_to_line(const typename K::Vector_2 &diff,
               const typename K::RT & wcross,
               const typename K::Vector_2 &dir )
  {
    typedef typename K::RT RT;
    typedef typename K::FT FT;
    RT numerator = CGAL_NTS square(wcross);
    RT denominator = wmult((K*)0, RT(wdot(dir,dir, K())),
                           diff.hw(), diff.hw());
    return Rational_traits<FT>().make_rational(numerator, denominator);
  }

  template <class K>
  typename K::FT
  squared_distance(const typename K::Segment_2 &seg,
                   const typename K::Line_2 &line,
                   const K& k)
  {
    typename K::Construct_vector_2 construct_vector;
    typedef typename K::RT RT;
    typedef typename K::FT FT;
    typedef typename K::Vector_2 Vector_2;
    typedef typename K::Point_2  Point_2;
    const Vector_2 &linedir = line.direction().vector();
    const Point_2 &linepoint = line.point();
    Vector_2 startvec(construct_vector(linepoint, seg.source()));
    Vector_2 endvec(construct_vector(linepoint, seg.target()));

    bool crossing1;
    RT c1s, c1e;
    if (seg.source() == seg.target())
      return internal::squared_distance(seg.source(), line, k);
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
  squared_distance(const typename K::Line_2 &line,
                   const typename K::Segment_2 &seg,
                   const K& k)
  {
    return internal::squared_distance(seg, line, k);
  }

  template <class K>
  typename K::FT
  ray_ray_squared_distance_parallel(
    const typename K::Vector_2 &ray1dir,
    const typename K::Vector_2 &ray2dir,
    const typename K::Vector_2 &from1to2,
    const K& k)
  {
    typedef typename K::RT RT;
    typedef typename K::FT FT;
    if (!is_acute_angle(ray1dir, from1to2, k)) {
      if (!same_direction(ray1dir, ray2dir, k))
        return (typename K::FT)k.compute_squared_length_2_object()(from1to2);
    }
    RT wcr, w;
    wcr = wcross(ray1dir, from1to2, k);
    w = from1to2.hw();
    return (typename K::FT)(FT(wcr*wcr)
                            / FT(wmult((K*)0, RT(wdot(ray1dir, ray1dir, k)), w, w)));
  }

  template <class K>
  typename K::FT
  squared_distance(const typename K::Ray_2 &ray1,
                   const typename K::Ray_2 &ray2,
                   const K& k)
  {
    typename K::Construct_vector_2 construct_vector;
    typedef typename K::Vector_2 Vector_2;
    typedef typename K::FT FT;
    const Vector_2 &ray1dir = ray1.direction().vector();
    const Vector_2 &ray2dir = ray2.direction().vector();
    Vector_2 diffvec(construct_vector(ray1.source(),ray2.source()));

    bool crossing1, crossing2;
    switch (orientation(ray1dir, ray2dir, k)) {
    case COUNTERCLOCKWISE:
      crossing1 = !clockwise(diffvec, ray2dir, k);
      crossing2 = !counterclockwise(ray1dir, diffvec, k);
      break;
    case CLOCKWISE:
      crossing1 = !counterclockwise(diffvec, ray2dir, k);
      crossing2 = !clockwise(ray1dir, diffvec, k);
      break;
    default:
      return ray_ray_squared_distance_parallel(ray1dir,ray2dir,diffvec,k);
    }

    if (crossing1) {
      if (crossing2)
        return (FT)0;
      return internal::squared_distance(ray2.source(), ray1, k);
    } else {
      if (crossing2) {
        return internal::squared_distance(ray1.source(), ray2, k);
      } else {

        FT min1, min2;
        min1 = internal::squared_distance(ray1.source(), ray2, k);
        min2 = internal::squared_distance(ray2.source(), ray1, k);
        return (min1 < min2) ? min1 : min2;
      }
    }
  }

  template <class K>
  typename K::FT
  squared_distance(const typename K::Line_2 &line,
                   const typename K::Ray_2 &ray,
                   const K& k)
  {
    typename K::Construct_vector_2 construct_vector;
    typedef typename K::FT FT;
    typedef typename K::Vector_2 Vector_2;
    Vector_2 normalvec(line.a(), line.b());
    Vector_2 diff = construct_vector(line.point(), ray.source());
    FT sign_dist = k.compute_scalar_product_2_object()(diff,normalvec);
    if (sign_dist < FT(0)) {
      if (is_acute_angle(normalvec, ray.direction().vector(), k) )
        return (FT)0;
    } else {
      if (is_obtuse_angle(normalvec, ray.direction().vector(), k) )
        return (FT)0;
    }
    return (typename K::FT)((sign_dist*sign_dist)/k.compute_squared_length_2_object()(normalvec));
  }

  template <class K>
  inline typename K::FT
  squared_distance(const typename K::Ray_2 &ray,
                   const typename K::Line_2 &line,
                   const K& k)
  {
    return internal::squared_distance(line, ray, k);
  }

  template <class K>
  inline typename K::FT
  squared_distance(const typename K::Line_2 &line1,
                   const typename K::Line_2 &line2,
                   const K& k)
  {
    typedef typename K::FT FT;
    if (internal::parallel(line1, line2, k))
      return internal::squared_distance(line1.point(), line2, k);
    else
      return (FT)0;
  }

  template <class K>
  void
  distance_index(int &ind,
                 const typename K::Point_2 &pt,
                 const typename K::Ray_2 &ray,
                 const K& k)
  {
    typename K::Construct_vector_2 construct_vector;
    if (!is_acute_angle(ray.direction().vector(), construct_vector(ray.source(), pt), k)) {
      ind = 0;
      return;
    }
    ind = -1;
  }

  template <class K>
  void
  distance_index(int &ind,
                 const typename K::Point_2 &pt,
                 const typename K::Segment_2 &seg,
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
  squared_distance_indexed(const typename K::Point_2 &pt,
                           const typename K::Ray_2 &ray,
                           int ind,
                           const K& k)
  {
    if (ind == 0)
      return internal::squared_distance(pt, ray.source(), k);
    return internal::squared_distance(pt, ray.supporting_line(), k);
  }

  template <class K>
  inline typename K::FT
  squared_distance_indexed(const typename K::Point_2 &pt,
                           const typename K::Segment_2 &seg,
                           int ind,
                           const K& k)
  {
    if (ind == 0)
      return internal::squared_distance(pt, seg.source(), k);
    if (ind == 1)
      return internal::squared_distance(pt, seg.target(), k);
    return internal::squared_distance(pt, seg.supporting_line(), k);
  }

} // namespace internal

template <class K>
inline typename K::FT
squared_distance(const Point_2<K> & pt1, const Point_2<K> & pt2)
{
  return internal::squared_distance(pt1, pt2, K());
}


template <class K>
inline typename K::FT
squared_distance(const Point_2<K> &pt, const Line_2<K> &line)
{
  return internal::squared_distance(pt, line, K());
}


template <class K>
inline typename K::FT
squared_distance(const Line_2<K> & line, const Point_2<K> & pt)
{
    return squared_distance(pt, line);
}


template <class K>
inline typename K::FT
squared_distance(const Point_2<K> &pt, const Ray_2<K> &ray)
{
  return internal::squared_distance(pt, ray, K());
}

template <class K>
inline typename K::FT
squared_distance(const Ray_2<K> & ray, const Point_2<K> & pt)
{
    return squared_distance(pt, ray);
}


template <class K>
inline typename K::FT
squared_distance(const Point_2<K> &pt, const Segment_2<K> &seg)
{
  return internal::squared_distance(pt, seg, K());
}


template <class K>
inline typename K::FT
squared_distance(const Segment_2<K> & seg, const Point_2<K> & pt)
{
  return internal::squared_distance(pt, seg, K());
}


template <class K>
inline typename K::FT
squared_distance(const Segment_2<K> &seg1, const Segment_2<K> &seg2)
{
  typedef typename K::Kernel_tag Tag;
  return internal::squared_distance(seg1, seg2, K(), Tag());
}

template <class K>
inline typename K::FT
squared_distance(const Segment_2<K> &seg, const Ray_2<K> &ray)
{
  return internal::squared_distance(seg, ray, K());
}

template <class K>
inline typename K::FT
squared_distance(const Ray_2<K> & ray, const Segment_2<K> & seg)
{
  return internal::squared_distance(seg, ray, K());
}

template <class K>
inline typename K::FT
squared_distance(const Segment_2<K> &seg, const Line_2<K> &line)
{
  return internal::squared_distance(seg, line, K());
}

template <class K>
inline typename K::FT
squared_distance(const Line_2<K> & line, const Segment_2<K> & seg)
{
  return internal::squared_distance(seg, line, K());
}

template <class K>
inline typename K::FT
squared_distance(const Ray_2<K> &ray1, const Ray_2<K> &ray2)
{
  return internal::squared_distance(ray1, ray2, K());
}

template <class K>
inline typename K::FT
squared_distance(const Line_2<K> &line, const Ray_2<K> &ray)
{
  return internal::squared_distance(line, ray, K());
}

template <class K>
inline typename K::FT
squared_distance(const Ray_2<K> & ray, const Line_2<K> & line)
{
  return internal::squared_distance(line, ray, K());
}

template <class K>
inline typename K::FT
squared_distance(const Line_2<K> &line1, const Line_2<K> &line2)
{
  return internal::squared_distance(line1, line2, K());
}

} //namespace CGAL

#endif
