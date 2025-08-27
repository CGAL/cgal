// Copyright (c) 2023-2025  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_CDT_3_INTERNAL_READ_POLYGON_MESH_FOR_CDT_3_H
#define CGAL_CDT_3_INTERNAL_READ_POLYGON_MESH_FOR_CDT_3_H

#include <CGAL/license/Constrained_triangulation_3.h>

#include <CGAL/Constrained_triangulation_3/internal/cdt_debug_io.h>
#include <CGAL/Constrained_triangulation_3/internal/ostream_redirect_guard.h>

#include <CGAL/Polygon_mesh_processing/polygon_soup_self_intersections.h>
#include <tl/expected.hpp>
#include <sstream>

namespace CGAL {

  template <typename PolygonMesh>
  struct CDT_3_read_polygon_mesh_output {
    tl::expected<PolygonMesh, std::string> polygon_mesh;

    std::size_t nb_of_duplicated_points = 0;
    std::size_t nb_of_simplified_polygons = 0;
    std::size_t nb_of_new_polygons = 0;
    std::size_t nb_of_removed_invalid_polygons = 0;
    std::size_t nb_of_removed_duplicated_polygons = 0;
    std::size_t nb_of_removed_isolated_points = 0;

    bool polygon_soup_self_intersects = false;
    bool polygon_mesh_is_manifold = true;
  };

  template <typename PolygonMesh,
            typename Points,
            typename Faces,
            typename NamedParameters = parameters::Default_named_parameters>
  CDT_3_read_polygon_mesh_output<PolygonMesh> convert_polygon_soup_to_polygon_mesh_for_cdt_3(
      Points& points, Faces& faces, const NamedParameters& np = parameters::default_values())
  {
    CDT_3_read_polygon_mesh_output<PolygonMesh> result;
    using Traits = typename GetPolygonGeomTraits<Points, Faces, NamedParameters>::type;

    using parameters::choose_parameter;
    using parameters::get_parameter;

    auto traits = choose_parameter<Traits>(get_parameter(np, internal_np::geom_traits));
    auto verbose = choose_parameter(get_parameter(np, internal_np::verbose), false);
    bool do_repair = choose_parameter(get_parameter(np, internal_np::repair_polygon_soup), true);
    const auto& return_error = parameters::get_parameter_reference(np, internal_np::callback);

    namespace PMP = CGAL::Polygon_mesh_processing;
    namespace PMP_internal = PMP::internal;

    if (do_repair)
    {
      result.nb_of_duplicated_points = PMP::merge_duplicate_points_in_polygon_soup(points, faces, np);
      result.nb_of_simplified_polygons = PMP_internal::simplify_polygons_in_polygon_soup(points, faces, traits);
      result.nb_of_new_polygons = PMP_internal::split_pinched_polygons_in_polygon_soup(points, faces, traits);
      result.nb_of_removed_invalid_polygons = PMP_internal::remove_invalid_polygons_in_polygon_soup(points, faces);
      result.nb_of_removed_duplicated_polygons = PMP::merge_duplicate_polygons_in_polygon_soup(points, faces, np);
      result.nb_of_removed_isolated_points = PMP::remove_isolated_points_in_polygon_soup(points, faces);
    }

    result.polygon_soup_self_intersects = PMP::does_polygon_soup_self_intersect(points, faces, np);

    if (!PMP::orient_polygon_soup(points, faces))
    {
      result.polygon_mesh_is_manifold = false;
      if (verbose)
        std::cerr << "Some duplication happened during polygon soup orientation" << std::endl;
    }

    if (!PMP::is_polygon_soup_a_polygon_mesh(faces))
    {
      if (verbose)
        std::cerr << "Warning: polygon soup does not describe a polygon mesh" << std::endl;
      if constexpr (std::is_invocable_v<decltype(return_error)>)
        return return_error();
      else
        return result;
    }
    PMP::polygon_soup_to_polygon_mesh(points, faces, *result.polygon_mesh, parameters::default_values(), np);
    return result;
  }

  template <typename PolygonMesh, typename NamedParameters = parameters::Default_named_parameters>
  CDT_3_read_polygon_mesh_output<PolygonMesh>
  read_polygon_mesh_for_cdt_3(const std::string &fname,
                              const NamedParameters &np = parameters::default_values())
  {
    CDT_3_read_polygon_mesh_output<PolygonMesh> result;

    using VPM = typename CGAL::GetVertexPointMap<PolygonMesh, NamedParameters>::type;
    using Point = typename boost::property_traits<VPM>::value_type;

    using parameters::choose_parameter;
    using parameters::get_parameter;

    auto verbose = choose_parameter(get_parameter(np, internal_np::verbose), false);

    std::ostringstream local_verbose_output;
    {
      auto cerr_redirect = CGAL::internal::ostream_redirect_guard::redirect(std::cerr).to(local_verbose_output);

      auto return_error = [&]() {
        result.polygon_mesh = tl::unexpected(std::move(local_verbose_output).str());
        return result;
      };

      using Points = std::vector<Point>;
      using Face = std::vector<std::size_t>;
      using Faces = std::vector<Face>;
      Points points;
      Faces faces;
      if(!CGAL::IO::read_polygon_soup(fname, points, faces, CGAL::parameters::verbose(true))) {
        if(verbose)
          std::cerr << "Warning: cannot read polygon soup" << std::endl;
        return return_error();
      }
      result = convert_polygon_soup_to_polygon_mesh_for_cdt_3<PolygonMesh>(points, faces, np.callback(return_error));
    }
    if(verbose) {
      std::cerr << std::move(local_verbose_output).str();
    }
    return result;
  }

} // end namespace CGAL

#endif // CGAL_CDT_3_INTERNAL_READ_POLYGON_MESH_FOR_CDT_3_H
