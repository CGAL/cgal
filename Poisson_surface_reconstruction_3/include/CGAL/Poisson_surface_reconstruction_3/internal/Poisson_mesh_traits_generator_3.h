// Copyright (c) 2006-2007  INRIA Sophia-Antipolis (France).
// Copyright (c) 2024       GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau, Jane Tournois

// This file is a copy-paste-adaptation of Surface_mesher/include/CGAL/Surface_mesh_traits_generator_3.h
// Surface_mesher that has been deprecated and will be removed in the future.

#ifndef CGAL_POISSON_MESH_TRAITS_GENERATOR_3_H
#define CGAL_POISSON_MESH_TRAITS_GENERATOR_3_H

#include <CGAL/license/Poisson_surface_reconstruction_3.h>


#include <CGAL/Poisson_surface_reconstruction_3/internal/Poisson_sphere_oracle_3.h>

namespace CGAL {

template <class K>
class Sphere_3;

/** Default traits class.
 *  Partial specialization will be in other headers
*/
template <typename Surface>
struct Poisson_mesh_traits_generator_3
{
  typedef typename Surface::Surface_mesher_traits_3 Type;
  typedef Type type; // for Boost compatibility (meta-programming)
};

  // specialization for Kernel::Sphere_3
template <typename Kernel>
struct Poisson_mesh_traits_generator_3<CGAL::Sphere_3<Kernel> >
{
  typedef Surface_mesher::Poisson_sphere_oracle_3<Kernel> Type;
  typedef Type type; // for Boost compatibility (meta-programming)
};

} // end namespace CGAL

#endif // CGAL_POISSON_MESH_TRAITS_GENERATOR_3_H
