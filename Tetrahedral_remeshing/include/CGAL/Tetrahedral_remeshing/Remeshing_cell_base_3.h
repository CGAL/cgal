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

#include <CGAL/Simplicial_mesh_cell_base_3.h>

namespace CGAL
{
namespace Tetrahedral_remeshing
{

template<typename Subdomain_index,
         typename Surface_patch_index,
         typename TDS>
class Remeshing_cell_3
  : public CGAL::Simplicial_mesh_cell_3<Subdomain_index,
                                        Surface_patch_index,
                                        TDS>
{
private:
  using Base = CGAL::Simplicial_mesh_cell_3<Subdomain_index,
                                            Surface_patch_index,
                                            TDS>;

public:
  using Base::Base;

  void set_sliver_value(double value)
  {
    sliver_cache_validity_ = true;
    sliver_value_ = value;
  }

  double sliver_value() const
  {
    CGAL_assertion(is_cache_valid());
    return sliver_value_;
  }

  bool is_cache_valid() const { return sliver_cache_validity_; }
  void reset_cache_validity() const { sliver_cache_validity_ = false; }

private:
  double sliver_value_ = 0.;
  mutable bool sliver_cache_validity_ = false;
};

/*!
\ingroup PkgTetrahedralRemeshingClasses

The class `Remeshing_cell_base_3` is a model of the concept `SimplicialMeshCellBase_3`.
It is designed to serve as cell base class for the 3D triangulation
used in the tetrahedral remeshing process.

\tparam Subdomain_index Type of indices for subdomains of the discretized geometric domain.
Must be a model of `CopyConstructible`, `Assignable`, `DefaultConstructible`
and `EqualityComparable`. The default constructed value must match the label
of the exterior of the domain (which contains at least the unbounded component).
 It must match the `Subdomain_index` of the model
  of the `MeshDomain_3` concept when used for mesh generation.

\tparam Surface_patch_index Type of indices for surface patches (boundaries and interfaces)
of the discretized geometric domain.
Must be a model of `CopyConstructible`, `Assignable`, `DefaultConstructible`
and `EqualityComparable`. The default constructed value must be the index value
assigned to a non surface facet.
 It must match the `Surface_patch_index` of the model
  of the `MeshDomain_3` concept when used for mesh generation.

\cgalModels{RemeshingCellBase_3,SimplicialMeshCellBase_3}

*/
template<typename Subdomain_index = int,
         typename Surface_patch_index = int>
class Remeshing_cell_base_3
{
public:
  template <typename TDS2>
  struct Rebind_TDS
  {
    typedef Remeshing_cell_3<Subdomain_index,
                             Surface_patch_index,
                             TDS2> Other;
  };
};

}//end namespace Tetrahedral_remeshing
}//end namespace CGAL

#endif //CGAL_TET_ADAPTIVE_REMESHING_CELL_BASE_3_H
