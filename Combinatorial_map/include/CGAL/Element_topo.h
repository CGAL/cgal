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
////////////////////////////////////////////////////////////////////////////////
#ifndef CMAP_ELEMENT_TOPO_H
#define CMAP_ELEMENT_TOPO_H

#include <string>

namespace CGAL {
  namespace CMap {
    namespace Element_topo {

enum cell_topo
  {
    SQUARE=0,
    TRIANGLE=1,
    HEXAHEDRON=2,
    TETRAHEDRON=3,
    PRISM=4,
    PYRAMID=5,
    GENERIC_2D=6,
    GENERIC_3D=7,
    EDGE=8,
    TETRAHEDRON10=9,
    PENTAGONAL_PRISM=10,
    HEXAGONAL_PRISM=11,
    NO_TYPE=-1
  };

inline
std::string topo_name(cell_topo t)
{
  switch(t)
  {
    case SQUARE: return "SQUARE";
    case TRIANGLE: return "TRIANGLE";
    case HEXAHEDRON: return "HEXAHEDRON";
    case TETRAHEDRON: return "TETRAHEDRON";
    case PRISM: return "PRISM";
    case PYRAMID: return "PYRAMID";
    case GENERIC_2D: return "GENERIC_2D";
    case GENERIC_3D: return "GENERIC_3D";
    case EDGE: return "EDGE";
    case TETRAHEDRON10: return "TETRAHEDRON10";
    case PENTAGONAL_PRISM: return "PENTAGONAL_PRISM";
    case HEXAGONAL_PRISM: return "HEXAGONAL_PRISM";
    case NO_TYPE: return "NO_TYPE";
  }
  return "UNKNOWN";
}

inline
cell_topo topo_from_name(const std::string& t)
{
  if (t=="SQUARE") return SQUARE;
  if (t=="TRIANGLE") return TRIANGLE;
  if (t=="HEXAHEDRON") return HEXAHEDRON;
  if (t=="TETRAHEDRON") return TETRAHEDRON;
  if (t=="PRISM") return PRISM;
  if (t=="PYRAMID") return PYRAMID;
  if (t=="GENERIC_2D") return GENERIC_2D;
  if (t=="GENERIC_3D") return GENERIC_3D;
  if (t=="EDGE") return EDGE;
  if (t=="TETRAHEDRON10") return TETRAHEDRON10;
  if (t=="PENTAGONAL_PRISM") return PENTAGONAL_PRISM;
  if (t=="HEXAGONAL_PRISM") return HEXAGONAL_PRISM;
  if (t=="NO_TYPE") return NO_TYPE;
  return NO_TYPE;
}

/**
 * @brief To get the type of dimD cell of the CMap of cmapdim dimension.
 */
template<typename CMap, unsigned int dimcell,
         unsigned int cmapdim=CMap::dimension>
struct Get_cell_topo
{
  static cell_topo run(CMap&, typename CMap::Dart_descriptor dh,
                       typename CMap::Dart_descriptor& starting_dart)
  {
    starting_dart=dh;
    return NO_TYPE;
  }
};

/**
 * @brief To get the type associated of an edge. For now only one type.
 */
template<typename CMap, unsigned int cmapdim>
struct Get_cell_topo<CMap, 1, cmapdim>
{
  static cell_topo run(CMap&, typename CMap::Dart_descriptor it,
                       typename CMap::Dart_descriptor& starting_dart)
  {
    starting_dart=it;
    return EDGE;
  }
};

/**
 * @brief To get the type of 2D cell of the CMap of cmapdim dimension.
 */
template<typename CMap, unsigned int cmapdim>
struct Get_cell_topo<CMap, 2, cmapdim>
{
  static cell_topo run(CMap& cmap, typename CMap::Dart_descriptor it,
                       typename CMap::Dart_descriptor& starting_dart)
  {
    starting_dart=it;

    if (cmap.is_face_combinatorial_polygon(it, 3))
    { return TRIANGLE; }

    else if (cmap.is_face_combinatorial_polygon(it, 4))
    { return SQUARE; }

    return GENERIC_2D;
  }
};

/**
 * @brief To get the type of 3D cell of the CMap of dimension 3.
 */
template<typename CMap>
struct Get_cell_topo<CMap, 3, 3>
{
  static cell_topo run(CMap& cmap, typename CMap::Dart_descriptor it,
                       typename CMap::Dart_descriptor& starting_dart)
  {
    starting_dart=it;

    if (cmap.is_volume_combinatorial_tetrahedron(it))
    { return TETRAHEDRON; }

    else if (cmap.is_volume_combinatorial_hexahedron(it))
    { return HEXAHEDRON; }

    else if(cmap.is_volume_combinatorial_tetrahedron10(it))
    { return TETRAHEDRON10; }

    // For non symetric object, we need to test all darts
    for (auto itv=cmap.template darts_of_cell<3>(it).begin(),
         itvend=cmap.template darts_of_cell<3>(it).end(); itv!=itvend; ++itv)
    {
      starting_dart=itv;

      if (cmap.is_volume_combinatorial_prism(itv))
      { return PRISM; }

      else if (cmap.is_volume_combinatorial_pentagonal_prism(itv))
      { return PENTAGONAL_PRISM; }

      else if (cmap.is_volume_combinatorial_pyramid(itv))
      { return PYRAMID; }

      else if (cmap.is_volume_combinatorial_hexagonal_prism(itv))
      { return HEXAGONAL_PRISM; }

    }

    return GENERIC_3D;
  }
};

template<unsigned int dimcell, typename CMap>
cell_topo get_cell_topo(CMap& cmap, typename CMap::Dart_descriptor it,
                        typename CMap::Dart_descriptor& starting_dart)
{ return Get_cell_topo<CMap, dimcell>::run(cmap, it, starting_dart); }

template<unsigned int dimcell, typename CMap>
cell_topo get_cell_topo(CMap& cmap, typename CMap::Dart_descriptor it)
{
  typename CMap::Dart_descriptor dummy;
  return get_cell_topo<dimcell, CMap>(cmap, it, dummy);
}

template<unsigned int dimcell, typename CMap>
cell_topo get_cell_topo(const CMap& cmap, typename CMap::Dart_const_descriptor it,
    typename CMap::Dart_const_descriptor& starting_dart)
{
  typename CMap::Dart_descriptor it2=const_cast<CMap&>(cmap).dart_descriptor
                                       (cmap.darts().index(it));
  typename CMap::Dart_descriptor sd2;
  cell_topo res=Get_cell_topo<CMap, dimcell>::run(const_cast<CMap&>(cmap),
                                                    it2, sd2);
  starting_dart=sd2;
  return res;
}

template<unsigned int dimcell, typename CMap>
cell_topo get_cell_topo(const CMap& cmap, typename CMap::Dart_const_descriptor it)
{
  typename CMap::Dart_descriptor it2=it;
  return Get_cell_topo<CMap, dimcell>::run(const_cast<CMap&>(cmap), it2);
}

} } } // namespace CGAL::CMap::Element_topo

#endif // CMAP_ELEMENT_TOPO_H
