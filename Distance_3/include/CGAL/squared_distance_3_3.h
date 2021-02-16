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

  if (! t.has_on_unbounded_side(pt)){
    return K::FT(0);
  }

  const typename K::Triangle_3 t0 = {t[0], t[1], t[2]};
  const typename K::Triangle_3 t1 = {t[0], t[1], t[3]};
  const typename K::Triangle_3 t2 = {t[1], t[2], t[3]};
  const typename K::Triangle_3 t3 = {t[0], t[2], t[3]};

  const typename K::FT d0 = squared_distance(pt, t0, k);
  const typename K::FT d1 = squared_distance(pt, t1, k);
  const typename K::FT d2 = squared_distance(pt, t2, k);
  const typename K::FT d3 = squared_distance(pt, t3, k);

  return (std::min)({d0, d1, d2, d3});
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
