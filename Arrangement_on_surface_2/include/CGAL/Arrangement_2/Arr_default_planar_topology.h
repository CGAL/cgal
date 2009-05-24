// Copyright (c) 2006 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s): Ron Wein          <wein@post.tau.ac.il>

#ifndef CGAL_ARR_DEFAULT_PLANAR_TOPOLOGY_H
#define CGAL_ARR_DEFAULT_PLANAR_TOPOLOGY_H

/*! \file
 * Definition of a helper template that uses partial specialization to
 * select the default planar topology-traits class best fit for a
 * given geometry-traits class.
 */

#include <CGAL/tags.h>
#include <CGAL/Arr_bounded_planar_topology_traits_2.h>
#include <CGAL/Arr_unb_planar_topology_traits_2.h>

CGAL_BEGIN_NAMESPACE

template <class GeomTraits, class Dcel, class Tag>
struct Default_planar_topology
{};

template <class GeomTraits, class Dcel>
struct Default_planar_topology<GeomTraits, Dcel, Arr_unbounded_boundary_tag>
{
  // A topology-traits class that supports unbounded curves:
  typedef Arr_unb_planar_topology_traits_2<GeomTraits, Dcel>        Traits;
};

template <class GeomTraits, class Dcel>
struct Default_planar_topology<GeomTraits, Dcel, Arr_no_boundary_tag>
{
  // A topology-traits class that supports only bounded curves:
  typedef Arr_bounded_planar_topology_traits_2<GeomTraits, Dcel>    Traits;

};

CGAL_END_NAMESPACE

#endif

