// Copyright (c) 2023-2024  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_CDT_3_DEBUG_IO_H
#define CGAL_CDT_3_DEBUG_IO_H

#include <CGAL/license/Constrained_triangulation_3.h>

#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>

#include <ostream>
#include <sstream>

#include <tl/expected.hpp>

namespace CGAL {

  template <typename Tr, typename Facets>
  auto export_facets_to_surface_mesh(const Tr& tr, Facets&& facets_range) {
    using Point_3 = typename Tr::Geom_traits::Point_3;
    const auto size = std::distance(facets_range.begin(), facets_range.end());
    std::vector<Point_3> points;
    points.reserve(size * 3);
    std::vector<std::array<std::size_t, 3>> facets;
    facets.reserve(size);

    std::size_t i = 0;
    for(const auto& [cell, facet_index] : facets_range) {
      const auto v0 = cell->vertex(Tr::vertex_triple_index(facet_index, 0));
      const auto v1 = cell->vertex(Tr::vertex_triple_index(facet_index, 1));
      const auto v2 = cell->vertex(Tr::vertex_triple_index(facet_index, 2));
      points.push_back(tr.point(v0));
      points.push_back(tr.point(v1));
      points.push_back(tr.point(v2));
      facets.push_back({i, i+1, i+2});
      i += 3;
    }
    CGAL::Polygon_mesh_processing::merge_duplicate_points_in_polygon_soup(points, facets);
    CGAL::Polygon_mesh_processing::orient_polygon_soup(points, facets);
    Surface_mesh<Point_3> mesh;
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, facets, mesh);
    return mesh;
  }

  template <typename Tr, typename Facets>
  void write_facets(std::ostream& out, const Tr& tr, Facets&& facets_range) {
    const auto mesh = export_facets_to_surface_mesh(tr, std::forward<Facets>(facets_range));
    CGAL::IO::write_OFF(out, mesh);
  }

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

  template <typename PolygonMesh, typename NamedParameters = parameters::Default_named_parameters>
  CDT_3_read_polygon_mesh_output<PolygonMesh>
  read_polygon_mesh_for_cdt_3(const std::string &fname,
                              const NamedParameters &np = parameters::default_values())
  {
    CDT_3_read_polygon_mesh_output<PolygonMesh> result;

    namespace PMP = CGAL::Polygon_mesh_processing;
    namespace PMP_internal = PMP::internal;

    using VPM = typename CGAL::GetVertexPointMap<PolygonMesh, NamedParameters>::type;
    using Point = typename boost::property_traits<VPM>::value_type;

    using parameters::choose_parameter;
    using parameters::get_parameter;

    auto verbose = choose_parameter(get_parameter(np, internal_np::verbose), false);

    std::ostringstream local_verbose_output;
    auto *cerr_buff = std::cerr.rdbuf();
    std::cerr.rdbuf(local_verbose_output.rdbuf());
    auto restore_cerr = make_scope_exit([&]
                                        { std::cerr.rdbuf(cerr_buff); });

    auto return_error = [&]() {
      result.polygon_mesh = tl::unexpected(std::move(local_verbose_output).str());
      return result;
    };

    using Points = std::vector<Point>;
    using Face = std::vector<std::size_t>;
    using Faces = std::vector<Face>;
    Points points;
    Faces faces;
    if (!CGAL::IO::read_polygon_soup(fname, points, faces, CGAL::parameters::verbose(true)))
    {
      if (verbose)
        std::cerr << "Warning: cannot read polygon soup" << std::endl;
      return return_error();
    }
    using Traits = typename PMP_internal::GetPolygonGeomTraits<Points, Faces, NamedParameters>::type;

    auto traits = choose_parameter<Traits>(get_parameter(np, internal_np::geom_traits));

    bool do_repair = choose_parameter(get_parameter(np, internal_np::repair_polygon_soup), true);
    if (do_repair)
    {
      result.nb_of_duplicated_points = PMP::merge_duplicate_points_in_polygon_soup(points, faces, np);
      result.nb_of_simplified_polygons = PMP_internal::simplify_polygons_in_polygon_soup(points, faces, traits);
      result.nb_of_new_polygons = PMP_internal::split_pinched_polygons_in_polygon_soup(points, faces, traits);
      result.nb_of_removed_invalid_polygons = PMP_internal::remove_invalid_polygons_in_polygon_soup(points, faces);
      result.nb_of_removed_duplicated_polygons = PMP::merge_duplicate_polygons_in_polygon_soup(points, faces, np);
      result.nb_of_removed_isolated_points = PMP::remove_isolated_points_in_polygon_soup(points, faces);
    }

    // check if the polygon soup is pure triangles, and create a triangulated copy otherwise
    bool is_pure_triangles = std::all_of(faces.begin(), faces.end(), [](const Face &f) { return f.size() == 3; });

    { // Now, call does_triangle_soup_self_intersect.
      // ... but that function requires a triangulated soup (triangle soup)
      // So, if needed, create a copy of the range of faces, and triangulate it on the fly.

      // create a non-deleting pointer to `faces` (with a null deleter)
      using Deleter_function = void(Faces*);
      using Deleter = Deleter_function*; // function pointer type
      using Ptr = std::unique_ptr<Faces, Deleter>;
      auto null_deleter = +[](Faces *) {};
      Ptr triangle_faces_ptr{&faces, null_deleter};
      if (!is_pure_triangles)
      {
        auto faces_copy_ptr = new Faces(faces); // copy `faces`
        auto delete_function_ptr = +[](Faces* vector){ delete vector; };
        triangle_faces_ptr = Ptr(faces_copy_ptr, delete_function_ptr);
        PMP::triangulate_polygons(points, *triangle_faces_ptr, np);
      }

      result.polygon_soup_self_intersects = PMP::does_triangle_soup_self_intersect(points, *triangle_faces_ptr, np);
    }

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
      return return_error();
    }
    PMP::polygon_soup_to_polygon_mesh(points, faces, *result.polygon_mesh, parameters::default_values(), np);

    return result;
  }

} // end namespace CGAL

#endif // CGAL_CDT_3_DEBUG_IO_H
