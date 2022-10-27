// Copyright (c) 2022 GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois

#ifndef TETRAHEDRAL_REMESHING_CELL_SELECTOR_H
#define TETRAHEDRAL_REMESHING_CELL_SELECTOR_H

#include <CGAL/license/Tetrahedral_remeshing.h>

#include <CGAL/Mesh_3/min_dihedral_angle.h> //todo : remove dependency on Mesh_3

#include <CGAL/property_map.h>

#include <unordered_set>

namespace CGAL
{
/*!
* \ingroup PkgTetrahedralRemeshingRef
*
*
* \sa `CGAL::tetrahedral_isotropic_remeshing()`
* 
* \todo use sliver_value() from Mesh_cell_base_3.h when available
*/
template<typename Triangulation>
class Sliver_selection_property_map
{
  using Tr = Triangulation;
  using Cell_handle = typename Tr::Cell_handle;
  using FT = typename Tr::Geom_traits::FT;
  using Subdomain_index
    = typename Tr::Triangulation_data_structure::Cell::Subdomain_index;

  const Tr& m_tr;
  std::unordered_set<Cell_handle> m_slivers;

public:
  using key_type    = Cell_handle;
  using value_type  = bool;
  using reference   = bool;
  using category    = boost::read_write_property_map_tag;

  Sliver_selection_property_map(const Tr& tr, const FT& sliver_bound)
    : m_tr(tr)
    , m_slivers()
  {
    collect_slivers(sliver_bound);
  }

  friend value_type get(const Sliver_selection_property_map& map,
                        const key_type& c)
  {
    return map.m_slivers.find(c) != map.m_slivers.end();
  }
  friend void put(Sliver_selection_property_map& map,
                  const key_type& c,
                  const value_type b)
  {
    if (b)               map.m_slivers.insert(c);
    else if (get(map, c)) map.m_slivers.erase(c);
  }

private:
  template<typename CellsSet>
  void insert_neighbors(Cell_handle c, CellsSet& cells)
  {
    for (int i = 0; i < 4; ++i)
    {
      Cell_handle ni = c->neighbor(i);
      if (Subdomain_index() != ni->subdomain_index())
        cells.insert(ni);
    }
  }

  void collect_slivers(const FT& sliver_bound)
  {
    using CGAL::Mesh_3::minimum_dihedral_angle;
    auto cp = m_tr.geom_traits().construct_point_3_object();

    for (Cell_handle c : m_tr.finite_cell_handles())
    {
      FT a = minimum_dihedral_angle(cp(c->vertex(0)->point()),
                                    cp(c->vertex(1)->point()),
                                    cp(c->vertex(2)->point()),
                                    cp(c->vertex(3)->point()),
                                    m_tr.geom_traits());
      if (a < sliver_bound)
        m_slivers.insert(c);
    }

    std::unordered_set<Cell_handle> neighbors;
    for (Cell_handle c : m_slivers)
      insert_neighbors(c, neighbors);
    //todo : we probably need to enlarge the selected region
  }


};//end class Sliver_selection_property_map

/**
*/
template<typename Tr>
Sliver_selection_property_map<Tr>
sliver_selection_property_map(const Tr& tr,
                              const typename Tr::Geom_traits::FT& angle)
{
  return Sliver_selection_property_map<Tr>(tr, angle);
}


}//end namespace CGAL

#endif //TETRAHEDRAL_REMESHING_CELL_SELECTOR_H
