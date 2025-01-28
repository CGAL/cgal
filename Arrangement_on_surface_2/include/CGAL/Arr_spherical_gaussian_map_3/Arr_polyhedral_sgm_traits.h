// Copyright (c) 2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Efi Fogel          <efif@post.tau.ac.il>

#ifndef CGAL_ARR_POLYHEDRAL_SGM_TRAITS_H
#define CGAL_ARR_POLYHEDRAL_SGM_TRAITS_H

#include <CGAL/license/Arrangement_on_surface_2.h>


#include <CGAL/basic.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>

#if defined(CGAL_ARR_TRACING_TRAITS)
#include "CGAL/Arr_tracing_traits_2.h"
#elif defined(CGAL_COUNTING_TRAITS)
#include "CGAL/Arr_counting_traits_2.h"
#endif

namespace CGAL {

/*! \file
 * A traits class-template for constructing and maintaining Gaussian maps
 * (also known as normal diagrams) of polytopes in 3D. It is parameterized
 * by a (linear) 3D geometry a kernel.
 * It serves two main purposes as follows:
 * 1. It models the concept ArrangementTraits, and as such, it handles arcs,
 *    the type of which are present in Gaussian maps of polytopes in 3D, that
 *    is, arcs of great circles on 2S, also known as geodesic arcs, and,
 * 2. It provides a type of point in 3D used to represent vertices of
 *    polytopes, and an operation that adds two points of this type.
 * This traits is used to instantiate a type that represents an arrangement
 * induced by geodesic arcs embedded on a sphere, where each face of the
 * arrangement stores the corresponding vertex of the polytope. This type,
 * in turn, represents a Gaussian map, which is a unique dual representation
 * of a polytope.
 */
template <typename Kernel_, int atan_x = -1, int atan_y = 0>
class Arr_polyhedral_sgm_traits :
#if defined(CGAL_ARR_TRACING_TRAITS)
  public Arr_tracing_traits_2<
    Arr_geodesic_arc_on_sphere_traits_2<Kernel_, atan_x, atan_y> >
#elif defined(CGAL_ARR_COUNTING_TRAITS)
  public Arr_counting_traits_2<
    Arr_geodesic_arc_on_sphere_traits_2<Kernel_, atan_x, atan_y> >
#else
  public Arr_geodesic_arc_on_sphere_traits_2<Kernel_, atan_x, atan_y>
#endif
{
public:
  typedef Kernel_                                                     Kernel;
  typedef typename Kernel::Point_3                                    Point_3;
  typedef typename Kernel::Vector_3                                   Vector_3;

protected:
#if defined(CGAL_ARR_TRACING_TRAITS)
  typedef Arr_tracing_traits_2<
    Arr_geodesic_arc_on_sphere_traits_2<Kernel, atan_x, atan_y> >     Base;
#elif defined(CGAL_ARR_COUNTING_TRAITS)
  typedef Arr_counting_traits_2<
    Arr_geodesic_arc_on_sphere_traits_2<Kernel, atan_x, atan_y> >     Base;
#else
  typedef Arr_geodesic_arc_on_sphere_traits_2<Kernel, atan_x, atan_y> Base;
#endif

public:
};

} //namespace CGAL

#endif
