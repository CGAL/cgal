// Copyright (c) 2020 GeometryFactory (France) and Telecom Paris (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois, Noura Faraj, Jean-Marc Thiery, Tamy Boubekeur

#ifndef CGAL_INTERNAL_PEEL_SLIVERS_H
#define CGAL_INTERNAL_PEEL_SLIVERS_H

#include <CGAL/license/Tetrahedral_remeshing.h>

#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>

namespace CGAL
{
namespace Tetrahedral_remeshing
{

template<typename C3T3, typename CellSelector>
std::size_t peel_slivers(C3T3& c3t3,
                         const typename C3T3::Triangulation::Geom_traits::FT& sliver_angle,
                         const CellSelector& cell_selector)
{
  using FT = typename C3T3::Triangulation::Geom_traits::FT;
  using Cell_handle = typename C3T3::Triangulation::Cell_handle;
  using Surface_patch_index = typename C3T3::Surface_patch_index;

  auto& tr = c3t3.triangulation();

  std::size_t nb_slivers_peel = 0;
  std::vector<std::pair<Cell_handle, std::array<bool, 4> > > peelable_cells;

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
  FT mindh = FT(180);
#endif
  for (Cell_handle cit : c3t3.cells_in_complex())
  {
    if(!get(cell_selector, cit))
      continue;

    std::array<bool, 4> facets_on_surface;

    const FT dh = min_dihedral_angle(tr, cit);
    if (dh < sliver_angle && is_peelable(c3t3, cit, facets_on_surface))
      peelable_cells.push_back(std::make_pair(cit, facets_on_surface));

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
      mindh = (std::min)(dh, mindh);
#endif
  }

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
  std::cout << "Min dihedral angle : " << mindh << std::endl;
  std::cout << "Peelable cells : " << peelable_cells.size() << std::endl;
#endif

  for (auto c_i : peelable_cells)
  {
    Cell_handle c = c_i.first;
    const std::array<bool, 4>& f_on_surface = c_i.second;

    boost::optional<Surface_patch_index> patch;
    for (int i = 0; i < 4; ++i)
    {
      if (f_on_surface[i])
      {
        Surface_patch_index spi = c3t3.surface_patch_index(c, i);
        if (patch != boost::none && patch != spi)
        {
          patch = boost::none;
          break;
        }
        else
        {
          patch = spi;
        }
      }
    }
    if (patch == boost::none)
      continue;

    for (int i = 0; i < 4; ++i)
    {
      if (f_on_surface[i])
        c3t3.remove_from_complex(c, i);
      else
        c3t3.add_to_complex(c, i, patch.get());
    }

    c3t3.remove_from_complex(c);
    ++nb_slivers_peel;
  }

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
  mindh = FT(180);
  for (Cell_handle cit : c3t3.cells_in_complex())
  {
    const FT dh = min_dihedral_angle(tr, cit);
    mindh = (std::min)(dh, mindh);
  }
  std::cout << "Peeling done (removed " << nb_slivers_peel << " slivers, "
    << "min dihedral angle = " << mindh << ")." << std::endl;
#endif

  return nb_slivers_peel;
}

} // end namespace Tetrahedral_remeshing
} // end namespace CGAL

#endif // CGAL_INTERNAL_PEEL_SLIVERS_H
