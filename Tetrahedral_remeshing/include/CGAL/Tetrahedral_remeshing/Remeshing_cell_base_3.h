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

#ifndef CGAL_TET_ADAPTIVE_REMESHING_CELL_BASE_3_H
#define CGAL_TET_ADAPTIVE_REMESHING_CELL_BASE_3_H

#include <CGAL/license/Tetrahedral_remeshing.h>

#include <CGAL/Compact_simplicial_mesh_cell_base_3.h>

#include <CGAL/assertions.h>

namespace CGAL {
namespace Tetrahedral_remeshing {

/*!
\ingroup PkgTetrahedralRemeshingClasses

The class `Remeshing_cell_base_3` is a model of the concept `RemeshingCellBase_3`.
It is designed to serve as cell base class for the 3D triangulation
used in the tetrahedral remeshing process.

\tparam Gt is the geometric traits class.
It has to be a model of the concept `RemeshingTriangulationTraits_3`.

\tparam Cb is a cell base class from which `Remeshing_cell_base_3` derives.
It must be a model of the `SimplicialMeshCellBase_3` concept.

\cgalModels{RemeshingCellBase_3}

*/
template<typename Gt,
         typename Cb = CGAL::Compact_simplicial_mesh_cell_base_3<
                         int /*Subdomain_index*/,
                         int /*Surface_patch_index*/> >
class Remeshing_cell_base_3
  : public Cb
{
public:
  using FT = typename Gt::FT;

  using Vertex_handle = typename Cb::Vertex_handle;
  using Cell_handle = typename Cb::Cell_handle;

  using Geom_traits = Gt;

public:
  template <typename TDS2>
  struct Rebind_TDS
  {
    using Cb2 = typename Cb::template Rebind_TDS<TDS2>::Other;
    using Other = Remeshing_cell_base_3<Gt, Cb2>;
  };

public:
  using Cb::Cb; // constructors

public:
  void set_sliver_value(const FT value)
  {
    sliver_cache_validity_ = true;
    sliver_value_ = value;
  }

  FT sliver_value() const
  {
    CGAL_assertion(is_cache_valid());
    return sliver_value_;
  }

  bool is_cache_valid() const { return sliver_cache_validity_; }
  void reset_cache_validity() const { sliver_cache_validity_ = false; }

private:
  FT sliver_value_ = 0.;
  mutable bool sliver_cache_validity_ = false;
};

} // namespace Tetrahedral_remeshing
} // namespace CGAL

#endif //CGAL_TET_ADAPTIVE_REMESHING_CELL_BASE_3_H
