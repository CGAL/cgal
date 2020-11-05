// Copyright (c) 2019 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Ron Wein         <wein@post.tau.ac.il>
//             Ophir Setter     <ophirset@post.tau.ac.il>
//             Sebastien Loriot <sebastien.loriot@cgal.org>
//             Efi Fogel        <efifogel@gmail.com>

#ifndef CGAL_GRAPH_TRAITS_DUAL_ARRANGEMENT_2_H
#define CGAL_GRAPH_TRAITS_DUAL_ARRANGEMENT_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

/*! \file
 * Definition of:
 * 1. the specialized Dual<Arrangement_2> class,
 * 2. the specialized boost::graph_traits<Dual<Arrangement_2> >class,
 * 3. The free functions required by the various graph concepts.
 */

// include this to avoid a VC15 warning
#include <CGAL/boost/graph/Named_function_parameters.h>

#include <CGAL/Arrangement_2.h>
#include <CGAL/Arrangement_2/graph_traits_dual.h>
#include <CGAL/disable_warnings.h>

namespace CGAL {

// The specialized Dual<Arrangement_2... class template.
template <typename GeomTraits_2, typename Dcel>
class Dual<Arrangement_2<GeomTraits_2, Dcel> > :
    public Dual_arrangement_on_surface<Arrangement_2<GeomTraits_2, Dcel> >
{
public:
  typedef Arrangement_2<GeomTraits_2, Dcel>               Arrangement;
  typedef typename Arrangement::Geometry_traits_2         Geometry_traits_2;
  typedef typename Arrangement::Topology_traits           Topology_traits;

private:
  typedef Dual_arrangement_on_surface<Arrangement>        Base;

public:
  /*! Default constructor. */
  Dual() : Base() {}

  /*! Constructor from an arrangement. */
  Dual(const Arrangement& arr) : Base(arr) {}
};

}

namespace boost {

// The specialized graph_traits<CGAL::Dual<CGAL::Arrangement_2... class template.
template <typename GeomTraits_2, typename Dcel>
class graph_traits<CGAL::Dual<CGAL::Arrangement_2<GeomTraits_2, Dcel> > > :
    public CGAL::Graph_traits_dual_arr_on_surface_impl<CGAL::Arrangement_2
                                                       <GeomTraits_2, Dcel> >
{};

}

namespace CGAL {

// Templates of free functions that handle
//   graph_traits<CGAL::Dual<CGAL::Arrangement_2... class template.
CGAL_DUAL_ARRANGEMENT_2_OUT_DEGREE(Arrangement_2)
CGAL_DUAL_ARRANGEMENT_2_OUT_EDGES(Arrangement_2)
CGAL_DUAL_ARRANGEMENT_2_SOURCE(Arrangement_2)
CGAL_DUAL_ARRANGEMENT_2_TARGET(Arrangement_2)
CGAL_DUAL_ARRANGEMENT_2_IN_DEGREE(Arrangement_2)
CGAL_DUAL_ARRANGEMENT_2_IN_EDGES(Arrangement_2)
CGAL_DUAL_ARRANGEMENT_2_DEGREE(Arrangement_2)
CGAL_DUAL_ARRANGEMENT_2_NUM_VERTICES(Arrangement_2)
CGAL_DUAL_ARRANGEMENT_2_VERTICES(Arrangement_2)
CGAL_DUAL_ARRANGEMENT_2_NUM_EDGES(Arrangement_2)
CGAL_DUAL_ARRANGEMENT_2_EDGES(Arrangement_2)

}

#endif
