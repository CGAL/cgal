// Copyright (c) 1998-2021
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
// Author(s)     : Geert-Jan Giezeman, Andreas Fabri


#ifndef CGAL_DISTANCE_3_3_H
#define CGAL_DISTANCE_3_3_H

#include <CGAL/squared_distance_3_2.h>

#include <CGAL/Point_3.h>
#include <CGAL/Tetrahedron_3.h>

namespace CGAL {

namespace internal {

template <class K>
inline
typename K::FT
squared_distance(const typename K::Tetrahedron_3 & t,
                 const typename K::Point_3 & pt,
                 const K& k)
{
  bool on_bounded_side = true;
  const typename K::Point_3 t0 = t[0];
  const typename K::Point_3 t1 = t[1];
  const typename K::Point_3 t2 = t[2];
  const typename K::Point_3 t3 = t[3];

  bool dmin_initialized = false;
  typename K::FT dmin;
  bool inside = false;
  if(orientation(t0,t1,t2, pt) == NEGATIVE){
    on_bounded_side = false;
    dmin = squared_distance_to_triangle(pt, t0, t1, t2, inside, k);
    dmin_initialized = true;
    if(inside){
      return dmin;
    }
  }

  if(orientation(t0,t3,t1, pt) == NEGATIVE){
    on_bounded_side = false;
    const typename K::FT d = squared_distance_to_triangle(pt, t0, t3, t1, inside, k);
    if(inside){
      return d;
    }
    if(! dmin_initialized){
      dmin = d;
      dmin_initialized = true;
    }else{
      dmin = (std::min)(d,dmin);
    }
  }

  if(orientation(t1,t3,t2, pt) == NEGATIVE){
    on_bounded_side = false;
    const typename K::FT d = squared_distance_to_triangle(pt, t1, t3, t2, inside, k);
    if(inside){
      return d;
    }
    if(! dmin_initialized){
      dmin = d;
      dmin_initialized = true;
    }else{
      dmin = (std::min)(d,dmin);
    }
  }

  if(orientation(t2,t3,t0, pt) == NEGATIVE){
    on_bounded_side = false;
    const typename K::FT d = squared_distance_to_triangle(pt, t2, t3, t0, inside, k);
    if(inside){
      return d;
    }
    if(! dmin_initialized){
      dmin = d;
      dmin_initialized = true;
    }else{
      dmin = (std::min)(d,dmin);
    }
  }

  if(on_bounded_side){
    return typename K::FT(0);
  }
  return dmin;
}

} // namespace internal


template <class K>
typename K::FT
squared_distance(const Tetrahedron_3<K> & t,
                 const Point_3<K> & pt)
{
  return internal::squared_distance(t,pt,K());
}


template <class K>
typename K::FT
squared_distance(const Point_3<K> & pt,
                 const Tetrahedron_3<K> & t)
{
  return internal::squared_distance(t,pt,K());
}

} //namespace CGAL


#endif
