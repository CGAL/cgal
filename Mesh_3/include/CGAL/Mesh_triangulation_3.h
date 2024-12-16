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

} // end namespace details

template<class Gt_, class Tds_>
class Mesh_3_regular_triangulation_3_wrapper
  : public Regular_triangulation_3<Gt_, Tds_>
{
public:
  typedef Regular_triangulation_3<Gt_, Tds_>                  Base;

  typedef typename Base::Geom_traits                          Geom_traits;

  typedef typename Geom_traits::FT                            FT;
  typedef typename Base::Bare_point                           Bare_point;
  typedef typename Base::Weighted_point                       Weighted_point;
  typedef typename Base::Triangle                             Triangle;

  typedef typename Base::Vertex_handle                        Vertex_handle;
  typedef typename Base::Facet                                Facet;
  typedef typename Base::Cell_handle                          Cell_handle;

  typedef typename Geom_traits::Vector_3                      Vector;

  using Base::geom_traits;
  using Base::point;
  using Base::triangle;

  static std::string io_signature() { return Get_io_signature<Base>()(); }

  // The undocumented, straightforward functions below are required for Mesh_3
  // because of Periodic_3_mesh_3: they are functions for which both triangulations
  // have fundamentally different implementations (usually, Periodic_3_mesh_3
  // does not know the offset of a point and must brute-force check for all
  // possibilities). To enable Periodic_3_mesh_3 to use Mesh_3's files,
  // each mesh triangulation implements its own version.

  const Bare_point& get_closest_point(const Bare_point& /*p*/, const Bare_point& q) const
  {
    return q;
  }

  Triangle get_incident_triangle(const Facet& f, const Vertex_handle) const
  {
    return triangle(f);
  }

  void set_point(const Vertex_handle v,
                 const Vector& /*move*/,
                 const Weighted_point& new_position)
  {
    v->set_point(new_position);
  }

  FT compute_power_distance_to_power_sphere(const Cell_handle c,
                                            const Vertex_handle v) const
  {
    typedef typename Geom_traits::Compute_power_distance_to_power_sphere_3 Critical_radius;

    Critical_radius critical_radius =
        geom_traits().compute_power_distance_to_power_sphere_3_object();

    return critical_radius(point(c, 0), point(c, 1), point(c, 2), point(c, 3), point(v));
  }

  // \pre c->neighbor(i) is finite
  FT compute_power_distance_to_power_sphere(const Cell_handle c, const int i) const
  {
    Cell_handle nc = c->neighbor(i);
    CGAL_precondition(!this->is_infinite(nc));
    Vertex_handle v = nc->vertex(nc->index(c));

    return compute_power_distance_to_power_sphere(c, v);
  }

  typename Geom_traits::FT min_squared_distance(const Bare_point& p, const Bare_point& q) const
  {
    return geom_traits().compute_squared_distance_3_object()(p, q);
  }
};


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
  using Triangulation =
      Mesh_3_regular_triangulation_3_wrapper<Geom_traits, Tds>;

public:
#ifndef DOXYGEN_RUNNING
  using type = Triangulation;
  using Type = type;
#else
  /// \name Types
  /// @{

  /*!
  The triangulation type to be used for the 3D triangulation embedding the mesh.
  This type is a wrapper around the type `CGAL::Regular_triangulation_3`, whose vertex
  and cell base classes are respectively `VertexBase` and `CellBase`.
  */
  typedef unspecified_type type;

  /// @}
#endif // DOXYGEN_RUNNING
};

} // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_MESH_TRIANGULATION_3_H
