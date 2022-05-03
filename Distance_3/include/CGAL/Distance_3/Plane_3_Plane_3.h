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

#ifndef CGAL_DISTANCE_3_PLANE_3_PLANE_3_H
#define CGAL_DISTANCE_3_PLANE_3_PLANE_3_H

#include <CGAL/Distance_3/internal/squared_distance_utils_3.h>

#include <CGAL/Plane_3.h>

namespace CGAL {
namespace internal {

template <class K>
inline typename K::FT
squared_distance(const typename K::Plane_3& plane1,
                 const typename K::Plane_3& plane2,
                 const K& k)
{
  typename K::Construct_orthogonal_vector_3 ortho_vec = k.construct_orthogonal_vector_3_object();
  typename K::Compute_squared_distance_3 sq_dist = k.compute_squared_distance_3_object();

  if(!is_null(wcross(ortho_vec(plane1), ortho_vec(plane2), k), k))
    return typename K::FT(0);
  else
    return sq_dist(plane1.point(), plane2);
}

} // namespace internal

template <class K>
inline
typename K::FT
squared_distance(const Plane_3<K>& plane1,
                 const Plane_3<K>& plane2)
{
  return K().compute_squared_distance_3_object()(plane1, plane2);
}

} // namespace CGAL

#endif // CGAL_DISTANCE_3_PLANE_3_PLANE_3_H
