// Copyright (c) 2009-2014 INRIA Sophia-Antipolis (France).
// Copyright (c) 2010-2013 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     :  Mael Rouxel-Labbé, Maxime Gimeno, Jane Tournois
//
//******************************************************************************
// File Description :
//******************************************************************************

#ifndef CGAL_SMDS_3_TETRAHEDRON_SOUP_TO_C3T3_H
#define CGAL_SMDS_3_TETRAHEDRON_SOUP_TO_C3T3_H

#include <CGAL/license/SMDS_3.h>

#include <CGAL/SMDS_3/tet_soup_to_c3t3.h>
#include <CGAL/SMDS_3/internal/SMDS_3_helper.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <vector>
#include <array>
#include <boost/unordered_map.hpp>

namespace CGAL {

  /** \ingroup PkgSMDS3Functions
  * builds a 3D triangulation from a soup of tetrahedra.
  *
  * @tparam TetrahedronRange a model of `Range` whose value type is
  * a `Tetrahedron_3`. The point type of the tetrahedra must be convertible
  * to `Triangulation::Point`
  * @tparam Triangulation a 3D triangulation class that has
  * a vertex base model of `SimplicialMeshVertexBase_3`
  * and a cell base model of `SimplicialMeshCellBase_3`
  *
  * @param tets the set of finite tetrahedra of a valid CGAL triangulation.
  * Each element in the range is the geometric description of the
  * corresponding cell in the triangulation.
  *
  * @returns the 3D triangulation built from \p tets
  *
  * @post the output triangulation must be a triangulation of the convex hull of `tets`
  *
  * @sa @ref SMDS_3/tetrahedron_soup_to_c3t3_example.cpp
  *
  */
  template<typename Triangulation, typename TetrahedronRange>
  Triangulation tetrahedron_soup_to_triangulation_3(const TetrahedronRange& tets)
  {
    using Tr    = Triangulation;
    using Point = typename Tr::Point;

    Triangulation tr;
    std::vector<Point> points;
    std::vector<std::array<int, 4> > finite_cells;
    boost::unordered_map<std::array<int, 3>, typename Tr::Cell::Surface_patch_index> border_facets;
    boost::unordered_map<Point, int> p2i;

    CGAL_assertion_code(
      typename Triangulation::Geom_traits::Orientation_3 orientation =
      tr.geom_traits().orientation_3_object();
    );

    for (typename TetrahedronRange::value_type tet : tets)
    {
      CGAL_assertion(tet.orientation() == CGAL::POSITIVE);
      std::array<int, 4> cell;

      for (int i = 0; i < 4; ++i)
      {
        const Point& pi = tet[i];
        if (p2i.find(pi) == p2i.end())
        {
          points.push_back(pi);
          int index = static_cast<int>(points.size() - 1);
          p2i.insert(std::make_pair(pi, index));
          cell[i] = index;
        }
        else
          cell[i] = p2i.at(pi);
      }

      CGAL_assertion(orientation(points[cell[0]],
        points[cell[1]], points[cell[2]], points[cell[3]]) == CGAL::POSITIVE);

      finite_cells.push_back(cell);
    }

    typename Tr::Cell::Subdomain_index default_si(1);
    CGAL::SMDS_3::build_triangulation_one_subdomain(tr, points, finite_cells, default_si, border_facets,
      /*verbose = */false, /*replace_domain_0 = */false, /*allow_non_manifold =*/false);

    CGAL_assertion(CGAL::SMDS_3::internal::is_convex(tr));

    return tr;
  }

  /** \ingroup PkgSMDS3Functions
  * builds a 3D triangulation from a soup of tetrahedra.
  *
  * @tparam PointRange a model of the concept `RandomAccessContainer`
  * whose value type is the point type.
  * The point type must be convertible to `Triangulation::Point`.
  * @tparam TetrahedronRange a model of the concept `RandomAccessContainer` whose
  * value type is a model of the concept `RandomAccessContainer` whose value type is `std::size_t`
  * @tparam Triangulation a 3D triangulation class that has
  * a vertex base model of `MeshVertexBase_3`
  * and a cell base model of `MeshCellBase_3`
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param points points of the soup of tetrahedra
  * @param tets the set of finite tetrahedra of a valid CGAL triangulation.
  * Each element in the range describes a tetrahedron using the indices of the points
  * in `points`.
  * It must
  * describe a non self-intersecting set of tetrahedra, that cover the convex hull of the
  * corresponding point set. The tetrahedra must form a valid triangulation with each
  * pair of neighboring cells sharing exactly one triangle. Combinatorial validity and
  * validity of the geometric embedding are required.
  *
  * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{surface_facets}
  *     \cgalParamDescription{each element in the range describes a surface facet using the indices of points
  *       in `points` (indices 0 to 2), and the associated `Surface_patch_index` (index 3)}
  *     \cgalParamType{a class model of `AssociativeContainer`
  *                    whose key type is model of `RandomAccessContainer` containing `int`
  *                    and mapped type is `Tr::Cell::Surface_patch_index`}
  *     \cgalParamDefault{an empty `boost::unordered_map<std::array<int, 3>, typename Tr::Cell::Surface_patch_index>`}
  *     \cgalParamExtra{to avoid copies of large data sets, this parameter can be passed using `std::cref`}
  *   \cgalParamNEnd
  *   \cgalParamNBegin{subdomain_indices}
  *     \cgalParamDescription{each element in the range gives the
  *       `Triangulation::Cell::Subdomain_index` corresponding to the tetrahedron (cell)
  *        of same index in `tets`}
  *     \cgalParamType{a class model of `RandomAccessContainer` whose value type
  *                    is `Triangulation::Cell::Subdomain_index`}
  *     \cgalParamDefault{each finite cell of the output triangulation is
  *                    set to have `1` as `Subdomain_index`}
  *     \cgalParamExtra{to avoid copies of large data sets, this parameter can be passed using `std::cref`}
  *   \cgalParamNEnd
  *\cond SKIP_IN_MANUAL
  *   \cgalParamNBegin{allow_non_manifold}
  *     \cgalParamDescription{allows the construction of a triangulation with non-manifold edges
  *       and non manifold vertices. The triangulation is invalid if this situation is met,
  *       so it should be used only in advanced cases, and the triangulation will be hardly usable.}
  *     \cgalParamType{bool}
  *     \cgalParamDefault{false}
  *   \cgalParamNEnd
  *\endcond
  * \cgalNamedParamsEnd
  *
  * @returns the 3D triangulation built from parameters
  *
  * @pre `points` contains each point only once
  * @post the output triangulation must be a triangulation of the convex hull of `points`
  * @post `is_valid()` returns `true` for the returned triangulation
  *
  * @sa \link Polygon_mesh_processing::polygon_soup_to_polygon_mesh() `CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh()` \endlink
  * @sa @ref SMDS_3/tetrahedron_soup_to_c3t3_example.cpp
  *
  */
  template<typename Triangulation,
           typename PointRange,
           typename TetrahedronRange,
           typename NamedParameters = parameters::Default_named_parameters>
  Triangulation tetrahedron_soup_to_triangulation_3(const PointRange& points,
                                                    const TetrahedronRange& tets,
                                                    const NamedParameters& np = parameters::default_values())
  {
    using parameters::choose_parameter;
    using parameters::get_parameter;
    using parameters::get_parameter_reference;

    using Default_facet_map
      = boost::unordered_map<std::array<int, 3>,
                             typename Triangulation::Cell::Surface_patch_index>;
    using Facet_map_ref_type
      = typename internal_np::Lookup_named_param_def<
                             internal_np::surface_facets_t,
                             NamedParameters,
                             Default_facet_map>::reference;
    using Subdomains_ref_type
      = typename internal_np::Lookup_named_param_def<
                              internal_np::subdomain_indices_t,
                              NamedParameters,
                              int>::reference;

    Triangulation tr;
    Default_facet_map empty_map;
    const Facet_map_ref_type& facets = choose_parameter(
          get_parameter_reference(np, internal_np::surface_facets),
          empty_map);
    const Subdomains_ref_type& subdomains = choose_parameter(
          get_parameter_reference(np, internal_np::subdomain_indices),
          1);
    const bool non_manifold = choose_parameter(
          get_parameter(np, internal_np::allow_non_manifold),
          false);

    CGAL::SMDS_3::build_triangulation_with_subdomains_range(tr, points, tets, subdomains, facets,
      /*verbose = */false, /*replace_domain_0 = */false, non_manifold);

    CGAL_assertion(CGAL::SMDS_3::internal::is_convex(tr));

    return tr;
  }

} //namespace CGAL


#endif // CGAL_SMDS_3_TETRAHEDRON_SOUP_TO_C3T3_H
