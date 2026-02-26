// Copyright (c) 2006-2009 INRIA Sophia-Antipolis (France).
// Copyright (c) 2011      GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau, Stephane Tayeb


#ifndef CGAL_MESH_TRIANGULATION_3_H
#define CGAL_MESH_TRIANGULATION_3_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Mesh_3/config.h>

#include <CGAL/Kernel_traits.h>

#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Robust_weighted_circumcenter_filtered_traits_3.h>

#include <CGAL/Mesh_vertex_base_3.h>
#include <CGAL/Compact_mesh_cell_base_3.h>
#include <CGAL/SMDS_3/io_signature.h>

namespace CGAL {

namespace details {

template<typename K>
struct Mesh_geom_traits_generator
{
private:
  typedef Robust_weighted_circumcenter_filtered_traits_3<K>   Geom_traits;

public:
  typedef Geom_traits                                         type;
  typedef type                                                Type;
};  // end struct Mesh_geom_traits_generator

} // namespace details

/*!
\ingroup PkgMesh3MeshClasses

The class `Mesh_triangulation_3` is a class template which provides the triangulation
type to be used for the 3D triangulation embedding the mesh.

\tparam MD must be a model of `MeshDomain_3`.

\tparam GT must be a model of `MeshTriangulationTraits_3` or `Default`
and defaults to `Kernel_traits<MD>::%Kernel`.

\tparam ConcurrencyTag enables sequential versus parallel meshing and optimization algorithms.
                       Possible values are `Sequential_tag` (the default), `Parallel_tag`,
                       and `Parallel_if_available_tag`.

\tparam VertexBase must be a model of `MeshVertexBase_3` or `Default`
and defaults to `Mesh_vertex_base_3<GT, MD>`.

\tparam CellBase must be a model of `MeshCellBase_3` or `Default`
and defaults to `Compact_mesh_cell_base_3<GT, MD>`.

\warning To improve the robustness of the meshing process, the input traits `GT`
         is wrapped with the traits class `Robust_weighted_circumcenter_filtered_traits_3`.
         The class `Robust_weighted_circumcenter_filtered_traits_3<GT>` upgrades the functors
         models of `Kernel::ConstructWeightedCircumcenter_3`, `Kernel::ComputeSquaredRadius_3`,
         and `Kernel::ComputeSquaredRadiusSmallestOrthogonalSphere_3` that are
         provided by `GT` to use exact computations when the geometric configuration
         is close to degenerate (e.g. almost coplanar points). <br><br>
         Users should therefore be aware that the traits class of the triangulation
         will have type `Robust_weighted_circumcenter_filtered_traits_3<GT>`.

\sa `make_mesh_3()`
\sa `Mesh_complex_3_in_triangulation_3<Tr,CornerIndex,CurveIndex>`

*/
template<class MD,
         class GT = Default,
         class ConcurrencyTag = Sequential_tag,
         class VertexBase = Default,
         class CellBase = Default>
struct Mesh_triangulation_3
{
private:
  using K = typename Default::Lazy_get<GT, Kernel_traits<MD> >::type;

  using Geom_traits = typename details::Mesh_geom_traits_generator<K>::type;

  using Indices_tuple = Mesh_3::internal::Indices_tuple_t<MD>;
  using Vertex_base = typename Default::Get<
    VertexBase,
    Mesh_vertex_generator_3<Geom_traits,
                            Indices_tuple,
                            typename MD::Index> >::type;
  using Cell_base = typename Default::Get<
    CellBase,
    Compact_mesh_cell_generator_3<Geom_traits,
                                  typename MD::Subdomain_index,
                                  typename MD::Surface_patch_index,
                                  typename MD::Index> >::type;
  using Concurrency_tag =
      typename Default::Get<ConcurrencyTag, Sequential_tag>::type;
  struct Tds : public Triangulation_data_structure_3<Vertex_base, Cell_base,
                                                     Concurrency_tag> {};
  using Triangulation = Regular_triangulation_3<Geom_traits, Tds>;

public:
#ifndef DOXYGEN_RUNNING
  using type = Triangulation;
  using Type = type;
#else
  /// \name Types
  /// @{

  /*!
  The triangulation type to be used for the 3D triangulation embedding the mesh.
  This type is a `Regular_triangulation_3` type whose vertex and cell base classes are respectively
  `VertexBase` and `CellBase`.
  */
  typedef unspecified_type type;

  /// @}
#endif // DOXYGEN_RUNNING
};

} // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_MESH_TRIANGULATION_3_H
