// Copyright (c) 2025 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef ELEMENT_TOPO_H
#define ELEMENT_TOPO_H

#include "Prism_and_pyramid_creation.h"

enum cell_topo
  {
    SQUARE = 0,
    TRIANGLE = 1,
    HEXAHEDRON = 2,
    TETRAHEDRON = 3,
    PRISM = 4,
    PYRAMID = 5,
    GENERIC_2D = 6,
    GENERIC_3D = 7,
    EDGE = 8,
    NO_TYPE = -1
  };

/**
 * @brief To get the type of dimD cell of the LCC of dimlcc dimension.
 */
template<typename LCC, unsigned int dim, unsigned int dimlcc=LCC::dimension>
struct Get_cell_topo
{
  static cell_topo run(LCC&, typename LCC::Dart_handle dh,
                       typename LCC::Dart_handle& starting_dart)
  {
    starting_dart=dh;
    return NO_TYPE;
  }
};

/**
 * @brief To get the type associated of an edge. For now only one type.
 */
template<typename LCC, unsigned int dimlcc>
struct Get_cell_topo<LCC, 1, dimlcc>
{
  static cell_topo run(LCC&, typename LCC::Dart_handle it,
                       typename LCC::Dart_handle& starting_dart)
  {
    starting_dart=it;
    return EDGE;
  }
};

/**
 * @brief To get the type of 2D cell of the LCC of dimlcc dimension.
 */
template<typename LCC, unsigned int dimlcc>
struct Get_cell_topo<LCC, 2, dimlcc>
{
  static cell_topo run(LCC& lcc, typename LCC::Dart_handle it,
                       typename LCC::Dart_handle& starting_dart)
  {
    starting_dart=it;

    if (lcc.is_face_combinatorial_polygon(it, 3))
      return TRIANGLE;

    else if (lcc.is_face_combinatorial_polygon(it, 4))
      return SQUARE;

    return GENERIC_2D;
  }
};

/**
 * @brief To get the type of 3D cell of the LCC of dimension 3.
 */
template<typename LCC>
struct Get_cell_topo<LCC, 3, 3>
{
  static cell_topo run(LCC& lcc, typename LCC::Dart_handle it,
                       typename LCC::Dart_handle& starting_dart)
  {
    starting_dart=it;

    if (lcc.is_volume_combinatorial_tetrahedron(it))
      return TETRAHEDRON;

    else if (lcc.is_volume_combinatorial_hexahedron(it))
      return HEXAHEDRON;

    // For non symetric object, we need to test all darts
    for (auto itv=lcc.template darts_of_cell<3>(it).begin(),
         itvend=lcc.template darts_of_cell<3>(it).end(); itv!=itvend; ++itv)
    {
      starting_dart=itv;

      if (is_volume_combinatorial_prism(lcc, itv))
        return PRISM;

      else if (is_volume_combinatorial_pyramid(lcc, itv))
        return PYRAMID;
    }

    return GENERIC_3D;
  }
};

template<typename LCC>
cell_topo get_cell_topo(LCC& lcc, typename LCC::Dart_handle it,
                        typename LCC::Dart_handle& starting_dart)
{ return Get_cell_topo<LCC, LCC::dimension>::run(lcc, it, starting_dart); }

template<typename LCC>
cell_topo get_cell_topo(LCC& lcc, typename LCC::Dart_handle it)
{
  typename LCC::Dart_handle dummy;
  return get_cell_topo(lcc, it, dummy);
}

#endif // ELEMENT_TOPO_H
