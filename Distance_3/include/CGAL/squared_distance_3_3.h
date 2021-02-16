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
squared_distance(const typename K::Tetrahedron_3 & tet,
                 const typename K::Point_3 & pt,
                 const K& k)
{
  CGAL_assertion(orientation(tet[0],tet[1],tet[2],tet[3]) == POSITIVE);

  if(orientation(tet[0],tet[1],tet[2],pt) != POSITIVE){
    const Triangle_3<K> tr = {tet[0], tet[1], tet[2]};
    return CGAL::squared_distance(tr, pt);
  }

  if(orientation(tet[1],tet[3],tet[2],pt) != POSITIVE){
    const Triangle_3<K> tr = {tet[1], tet[3], tet[2]};
    return CGAL::squared_distance(tr, pt);
  }

  if(orientation(tet[2],tet[3],tet[0],pt) != POSITIVE){
    const Triangle_3<K> tr = {tet[2], tet[3], tet[0]};
    return CGAL::squared_distance(tr, pt);
  }

  if(orientation(tet[0],tet[3],tet[1],pt) != POSITIVE){
    const Triangle_3<K> tr = {tet[0], tet[3], tet[1]};
    return CGAL::squared_distance(tr, pt);
  }

    return K::FT(0);
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
