// Copyright (c) 2023-2026 GeometryFactory and Claudio Mancinelli.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Claudio Mancinelli and Sébastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_BSURF_UTILS_H
#define CGAL_POLYGON_MESH_PROCESSING_BSURF_UTILS_H

#include <CGAL/license/Vector_graphics_on_surfaces.h>

#include <CGAL/Polygon_mesh_processing/locate.h>

namespace CGAL {
namespace Vector_graphics_on_surfaces {


/*!
 * \ingroup VGSMiscellaneous
 * computes the length of a path on a triangle mesh.
 * \tparam FT floating point number type (float or double)
 * \tparam TriangleMesh a model of `FaceGraph`
 * \param path a path described as a range of face locations, with the property that
               for two consecutive face locations, there exists a face in `tmesh` containing the two corresponding points.
 * \param tmesh the triangle mesh supporting the path
 * \todo add named parameters
 * \todo generic range
 */
template <class FT, class TriangleMesh>
FT path_length(const std::vector<CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, FT>>& path,
               const TriangleMesh &tmesh)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  std::size_t lpath = path.size();
  if(lpath<2)
    return 0;

  using VPM = typename boost::property_map<TriangleMesh, CGAL::vertex_point_t>::const_type;
  VPM vpm = get(CGAL::vertex_point, tmesh);

  FT len(0);

  for (std::size_t i=0; i<lpath-1; ++i)
    len += sqrt(squared_distance(PMP::construct_point(path[i],tmesh),
                                 PMP::construct_point(path[i+1],tmesh)));

  return len;
}

/*!
 * \ingroup VGSMiscellaneous
 * computes the length of a path on a triangle mesh.
 * \tparam FT floating point number type (float or double)
 * \tparam TriangleMesh a model of `FaceGraph`
 * \param src source of the path
 * \param tgt target of the path
 * \param path a path described as a range of edge locations, with the property that
               for two consecutive edge locations, there exists a face in `tmesh` containing the two corresponding points.
 * \param tmesh the triangle mesh supporting the path
 * \todo add named parameters
 * \todo generic range
 */
template <class FT, class TriangleMesh>
FT path_length(const std::vector<CGAL::Polygon_mesh_processing::Edge_location<TriangleMesh,FT>>& path,
               const CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, FT>& src,
               const CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, FT>& tgt,
               const TriangleMesh &tmesh)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  std::size_t lpath = path.size();
  if(lpath==0)
    return sqrt(squared_distance(PMP::construct_point(src,tmesh),
                                 PMP::construct_point(tgt,tmesh)));

  using VPM = typename boost::property_map<TriangleMesh, CGAL::vertex_point_t>::const_type;
  VPM vpm = get(CGAL::vertex_point, tmesh);

  FT len=sqrt(squared_distance(PMP::construct_point(src,tmesh),
                               PMP::construct_point(path[0],tmesh)));

  for (std::size_t i=0; i<lpath-1; ++i)
    len += sqrt(squared_distance(PMP::construct_point(path[i],tmesh),
                                 PMP::construct_point(path[i+1],tmesh)));

  len+=sqrt(squared_distance(PMP::construct_point(path.back(),tmesh),
                             PMP::construct_point(tgt,tmesh)));

  return len;
}

/*!
 * \ingroup VGSMiscellaneous
 * converts a path on a triangle mesh to the corresponding polyline of points.
 * If `path` contains identical consecutive vertices, only one point will be put in `poly_out` for this vertex.
 * \tparam FT floating point number type (float or double)
 * \tparam TriangleMesh a model of `FaceGraph`
 * \tparam OutputIterator an output iterator accepting points from `tmesh`
 * \param path a path described as a range of face locations, with the property that
               for two consecutive face locations, there exists a face in `tmesh` containing the two corresponding points.
 * \param tmesh the triangle mesh supporting the path
 * \param poly_out output iterator where points of the polyline are put.
 * \todo add named parameters
 * \todo generic range
 */
template<class TriangleMesh, class FT, class OutputIterator>
OutputIterator
convert_path_to_polyline(const std::vector<CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, FT>>& path,
                         const TriangleMesh& tmesh,
                         OutputIterator poly_out)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;

  vertex_descriptor vd = boost::graph_traits<TriangleMesh>::null_vertex();
  for (const PMP::Face_location<TriangleMesh, FT>& loc : path)
  {
    std::optional<vertex_descriptor> ov = PMP::vertex_descriptor_from_location(loc, tmesh);
    if (ov.has_value())
    {
      if (vd==*ov) continue;
      vd=*ov;
    }
    *poly_out++=PMP::construct_point(loc, tmesh);
  }
  return poly_out;
}

/*!
 * \ingroup VGSMiscellaneous
 * converts a path on a triangle mesh to the corresponding polyline of points.
 * If `path` contains identical consecutive vertices, only one point will be put in `poly_out` for this vertex.
 * \tparam FT floating point number type (float or double)
 * \tparam TriangleMesh a model of `FaceGraph`
 * \tparam OutputIterator an output iterator accepting points from `tmesh`
 * \param src source of the path
 * \param tgt target of the path
 * \param path a path described as a range of edge locations, with the property that
               for two consecutive edge locations, there exists an edge in `tmesh` containing the two corresponding points.
 * \param tmesh the triangle mesh supporting the path
 * \param poly_out output iterator where points of the polyline are put.
 * \todo add named parameters
 * \todo generic range
 */
template<class TriangleMesh, class FT, class OutputIterator>
OutputIterator
convert_path_to_polyline(const CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, FT>& src,
                         const std::vector<CGAL::Polygon_mesh_processing::Edge_location<TriangleMesh, FT>>& path,
                         const CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, FT>& tgt,
                         const TriangleMesh& tmesh,
                         OutputIterator poly_out)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;

  *poly_out++=PMP::construct_point(src, tmesh);
  vertex_descriptor vd = boost::graph_traits<TriangleMesh>::null_vertex();
  for (const PMP::Edge_location<TriangleMesh, FT>& loc : path)
  {
    std::optional<vertex_descriptor> ov = PMP::vertex_descriptor_from_location(loc, tmesh);
    if (ov.has_value())
    {
      if (vd==*ov) continue;
      vd=*ov;
    }
    *poly_out++=PMP::construct_point(loc, tmesh);
  }
  *poly_out++=PMP::construct_point(tgt, tmesh);
  return poly_out;
}

/*!
 * \ingroup VGSMiscellaneous
 * converts the coordinates of a range of points into polar coordinates with respect to a given center
 * \tparam K a model of `PolygonTraits_2`
 * \tparam PointRange_2 a model of the concept `RandomAccessContainer` with `K::Point_2` as value type
 * \param points the input points to convert.
 * \param center the point of reference for the polar coordinates.
 *               If omitted, then centroid of `polygon` will be used.
 * \return the polar coordinates of the points as a pair (distance, angle)
 */
template <class K, class PointRange_2>
std::vector<std::pair<typename K::FT, typename K::FT>>
convert_to_polar_coordinates(const PointRange_2& points,
                             std::optional<typename K::Point_2> center = std::nullopt)
{
  std::vector<std::pair<typename K::FT, typename K::FT>> result;
  // compute naively the center
  if (center==std::nullopt)
  {
    Bbox_2 bb = bbox_2(points.begin(), points.end());
    center = typename K::Point_2((bb.xmax()+bb.xmin())/2., (bb.ymax()+bb.ymin())/2);
  }

  bool is_closed = points.front()==points.back();
  std::size_t nbp=points.size();
  if (is_closed) --nbp;
  for (std::size_t i=0; i<nbp; ++i)
  {
    const typename K::Point_2& pt = points[i];
    typename K::FT d = approximate_sqrt(squared_distance(*center, pt));
    typename K::FT polar = std::atan2((pt.y()-center->y())/* /d */, (pt.x()-center->x())/* /d */);
    result.emplace_back(d, polar);
    // std::cout << center->x()+d*std::cos(polar) << " " << center->x()+d*std::sin(polar) << " 0\n";
  }
  if (is_closed) result.push_back(result.front());

  return result;
}

} } // end of CGAL::Vector_graphics_on_surfaces

#endif // CGAL_POLYGON_MESH_PROCESSING_BSURF_UTILS_H