// Copyright (c) 2006-2007  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_MULTI_SURFACE_3_H
#define CGAL_MULTI_SURFACE_3_H

#include <CGAL/license/Surface_mesher.h>

#include <CGAL/disable_warnings.h>

namespace CGAL {

  template<
    typename Surface_a,
    typename Surface_b
    >
  class Multi_surface_3
  {
    const Surface_a& surf_a;
    const Surface_b& surf_b;
  public:
    Multi_surface_3(const Surface_a& surface_a, const Surface_b& surface_b)
      : surf_a(surface_a), surf_b(surface_b)
    {
    }

    const Surface_a& surface_a() const
    {
      return surf_a;
    }


    const Surface_b& surface_b() const
    {
      return surf_b;
    }
  };
} // end namespace CGAL

#include <CGAL/Surface_mesher/Combining_oracle.h>
#include <CGAL/Surface_mesh_traits_generator_3.h>

namespace CGAL {
  template <typename Surface_a, typename Surface_b>
  struct Surface_mesh_traits_generator_3 <
    Multi_surface_3<Surface_a, Surface_b>
  >
  {
    typedef typename Surface_mesh_traits_generator_3<Surface_a>::type Oracle_a;
    typedef typename Surface_mesh_traits_generator_3<Surface_b>::type Oracle_b;

    typedef typename Surface_mesher::Combining_oracle<
      Oracle_a,
      Oracle_b
      > Type;

    typedef Type type; // Boost meta-programming compatibility
  };
} // end namespace CGAL, second occurrence.

#include <CGAL/enable_warnings.h>

#endif // CGAL_MULTI_SURFACE_3_H
