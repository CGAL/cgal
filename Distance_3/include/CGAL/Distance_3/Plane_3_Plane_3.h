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
#include <CGAL/Distance_3/Point_3_Plane_3.h>

#include <CGAL/Plane_3.h>

namespace CGAL {

template <class K>
inline
typename K::FT
squared_distance(const Plane_3<K>& plane1,
                 const Plane_3<K>& plane2)
{
  K k;
  typename K::Construct_orthogonal_vector_3 ortho_vec = k.construct_orthogonal_vector_3_object();

  if (!internal::is_null(internal::wcross(ortho_vec(plane1), ortho_vec(plane2), k), k))
    return typename K::FT(0);
  else
    return internal::squared_distance(plane1.point(), plane2, k);
}

} // namespace CGAL

#endif // CGAL_DISTANCE_3_PLANE_3_PLANE_3_H
