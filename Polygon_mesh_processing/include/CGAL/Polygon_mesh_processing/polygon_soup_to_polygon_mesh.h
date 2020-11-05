// Copyright (c) 2009-2013 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau and Ilker O. Yaz

#ifndef CGAL_POLYGON_MESH_PROCESSING_POLYGON_SOUP_TO_POLYGON_MESH_H
#define CGAL_POLYGON_MESH_PROCESSING_POLYGON_SOUP_TO_POLYGON_MESH_H

#include <CGAL/license/Polygon_mesh_processing/repair.h>

#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>

#include <CGAL/algorithm.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/property_map.h>

#include <boost/dynamic_bitset.hpp>
#include <boost/range/size.hpp>
#include <boost/range/value_type.hpp>
#include <boost/range/reference.hpp>

#include <array>
#include <set>
#include <type_traits>
#include <vector>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

template <typename PM_Point, typename PS_Point>
PM_Point convert_to_pm_point(const PS_Point& p)
{
  CGAL_static_assertion((std::is_convertible<PS_Point, PM_Point>::value));
  return PM_Point(p);
}

// just for backward compatibility reasons
template <typename PM_Point, typename PS_FT>
PM_Point convert_to_pm_point(const std::array<PS_FT, 3>& p)
{
  return PM_Point(p[0], p[1], p[2]);
}

template <typename PointRange,
          typename PolygonRange,
          typename PointMap = typename CGAL::GetPointMap<PointRange>::const_type>
class PS_to_PM_converter
{
  typedef typename boost::range_value<PolygonRange>::type                 Polygon;

public:
  /**
  * The constructor for modifier object.
  * @param points points of the soup of polygons.
  * @param polygons each element in the range describes a polygon using the index of the points in the range.
  */
  PS_to_PM_converter(const PointRange& points,
                     const PolygonRange& polygons,
                     const PointMap pm = PointMap())
    : m_points(points),
      m_polygons(polygons),
      m_pm(pm)
  { }

  template <typename PolygonMesh, typename VertexPointMap>
  void operator()(PolygonMesh& pmesh,
                  VertexPointMap vpm,
                  const bool insert_isolated_vertices = true)
  {
    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor    vertex_descriptor;

    typedef typename boost::property_traits<VertexPointMap>::value_type     PM_Point;

    reserve(pmesh, static_cast<typename boost::graph_traits<PolygonMesh>::vertices_size_type>(m_points.size()),
            static_cast<typename boost::graph_traits<PolygonMesh>::edges_size_type>(2*m_polygons.size()),
            static_cast<typename boost::graph_traits<PolygonMesh>::faces_size_type>(m_polygons.size()));

    boost::dynamic_bitset<> not_isolated;
    if(!insert_isolated_vertices)
    {
      not_isolated.resize(m_points.size());
      for(std::size_t i = 0, end = m_polygons.size(); i < end; ++i)
      {
        const Polygon& polygon = m_polygons[i];
        const std::size_t size = polygon.size();
        for(std::size_t j = 0; j < size; ++j)
          not_isolated.set(polygon[j], true);
      }
    }

    std::vector<vertex_descriptor> vertices(m_points.size());
    for(std::size_t i = 0, end = m_points.size(); i < end; ++i)
    {
      if(!insert_isolated_vertices && !not_isolated.test(i))
        continue;

      vertices[i] = add_vertex(pmesh);
      PM_Point pi = convert_to_pm_point<PM_Point>(get(m_pm, m_points[i]));
      put(vpm, vertices[i], pi);
    }

    for(std::size_t i = 0, end = m_polygons.size(); i < end; ++i)
    {
      const Polygon& polygon = m_polygons[i];
      const std::size_t size = polygon.size();

      std::vector<vertex_descriptor> vr(size); //vertex range
      vr.resize(size);
      for(std::size_t j = 0; j < size; ++j)
        vr[j] = vertices[polygon[j] ];

      CGAL_assertion_code(typename boost::graph_traits<PolygonMesh>::face_descriptor fd =)
          CGAL::Euler::add_face(vr, pmesh);
      CGAL_assertion(fd != boost::graph_traits<PolygonMesh>::null_face());
    }
  }

  template <typename PolygonMesh>
  void operator()(PolygonMesh& pmesh,
                  const bool insert_isolated_vertices = true)
  {
    return operator()(pmesh, get(CGAL::vertex_point, pmesh), insert_isolated_vertices);
  }

private:
  const PointRange& m_points;
  const PolygonRange& m_polygons;
  const PointMap m_pm;
};

} // namespace internal

/**
* \ingroup PMP_repairing_grp
*
* returns `true` if the soup of polygons defines a valid polygon
* mesh that can be handled by
* `CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh()`.
* It checks that each edge has at most two incident faces and such an edge
* is visited in opposite direction along the two face boundaries,
* no polygon has twice the same vertex,
* and the polygon soup describes a manifold surface.
* This function does not require a range of points as an argument
* since the check is purely topological. To each vertex of the mesh
* is associated an index that is used in the description of the
* boundaries of the polygons provided in `polygons`.
*
* @tparam PolygonRange a model of the concept `RandomAccessContainer`
* whose value_type is a model of the concept `RandomAccessContainer`
* whose value_type is `std::size_t`.
*
* @param polygons each element in the range describes a polygon
* using the indices of the vertices.
*
* @sa `orient_polygon_soup()`
*/
template<typename PolygonRange>
bool is_polygon_soup_a_polygon_mesh(const PolygonRange& polygons)
{
  typedef typename boost::range_value<PolygonRange>::type           Polygon;
  typedef typename boost::range_value<Polygon>::type                V_ID;

  if(boost::begin(polygons) == boost::end(polygons))
    return true;

  //check there is no duplicated ordered edge, and
  //check there is no polygon with twice the same vertex
  std::set<std::pair<V_ID, V_ID> > edge_set;
  V_ID max_id = 0;
  for(const Polygon& polygon : polygons)
  {
    std::size_t nb_edges = boost::size(polygon);
    if(nb_edges < 3)
      return false;

    std::set<V_ID> polygon_vertices;
    V_ID prev = *std::prev(boost::end(polygon));
    for(V_ID id : polygon)
    {
      if(max_id<id)
        max_id = id;

      if(! edge_set.insert(std::pair<V_ID, V_ID>(prev,id)).second)
        return false;
      else
        prev = id;

      if(!polygon_vertices.insert(id).second)
        return false; // vertex met twice in the same polygon
    }
  }

  //check manifoldness
  typedef std::vector<V_ID>                                           PointRange;
  typedef internal::Polygon_soup_orienter<PointRange, PolygonRange>   Orienter;

  typename Orienter::Edge_map edges(max_id+1);
  typename Orienter::Marked_edges marked_edges;
  Orienter::fill_edge_map(edges, marked_edges, polygons);

  //returns false if duplication is necessary
  if(!marked_edges.empty())
    return false;

  return Orienter::has_singular_vertices(static_cast<std::size_t>(max_id+1),polygons,edges,marked_edges);
}

/**
* \ingroup PMP_repairing_grp
* builds a polygon mesh from a soup of polygons.
*
* @pre the input polygon soup describes a consistently oriented
* polygon mesh. This can be checked using the function
* \link CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh()
* `CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(polygons)` \endlink.
*
* @tparam PolygonMesh a model of `MutableFaceGraph`
* @tparam PointRange a model of the concept `RandomAccessContainer`
* whose value type is the point type
* @tparam PolygonRange a model of the concept `RandomAccessContainer` whose
* value type is a model of the concept `RandomAccessContainer` whose value type is `std::size_t`
* @tparam NamedParameters_PS a sequence of \ref bgl_namedparameters "Named Parameters"
* @tparam NamedParameters_PM a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param points points of the soup of polygons
* @param polygons each element in the vector describes a polygon using the indices of the points in `points`
* @param out the polygon mesh to be built
* @param np_ps an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{point_map}
*     \cgalParamDescription{a property map associating points to the elements of the range `points`}
*     \cgalParamType{a model of `ReadablePropertyMap` whose value type is a point type convertible to the point type
*                    of the vertex point map associated to the polygon mesh}
*     \cgalParamDefault{`CGAL::Identity_property_map`}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* @param np_pm an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `out`}
*     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
*                    as key type and `%Point_3` as value type}
*     \cgalParamDefault{`boost::get(CGAL::vertex_point, out)`}
*     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
*                     must be available in `PolygonMesh`.}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* \sa `CGAL::Polygon_mesh_processing::orient_polygon_soup()`
* \sa `CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh()`
* \sa `CGAL::Polygon_mesh_processing::polygon_mesh_to_polygon_soup()`
*/
template<typename PolygonMesh,
         typename PointRange, typename PolygonRange,
         typename NamedParameters_PS, typename NamedParameters_PM>
void polygon_soup_to_polygon_mesh(const PointRange& points,
                                  const PolygonRange& polygons,
                                  PolygonMesh& out,
                                  const NamedParameters_PS& np_ps,
                                  const NamedParameters_PM& np_pm)
{
  CGAL_precondition_msg(is_polygon_soup_a_polygon_mesh(polygons),
                        "Input soup needs to define a valid polygon mesh! See is_polygon_soup_a_polygon_mesh() for further information.");

  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename CGAL::GetPointMap<PointRange, NamedParameters_PS>::const_type    Point_map;
  Point_map pm = choose_parameter<Point_map>(get_parameter(np_ps, internal_np::point_map));

  typedef typename CGAL::GetVertexPointMap<PolygonMesh, NamedParameters_PM>::type   Vertex_point_map;
  Vertex_point_map vpm = choose_parameter(get_parameter(np_pm, internal_np::vertex_point),
                                          get_property_map(CGAL::vertex_point, out));

  internal::PS_to_PM_converter<PointRange, PolygonRange, Point_map> converter(points, polygons, pm);
  converter(out, vpm);
}

template<typename PolygonMesh,
         typename PointRange, typename PolygonRange,
         typename NamedParameters_PS>
void polygon_soup_to_polygon_mesh(const PointRange& points,
                                  const PolygonRange& polygons,
                                  PolygonMesh& out,
                                  const NamedParameters_PS& np_ps)
{
  return polygon_soup_to_polygon_mesh(points, polygons, out, np_ps, parameters::all_default());
}

template<typename PolygonMesh, typename PointRange, typename PolygonRange>
void polygon_soup_to_polygon_mesh(const PointRange& points,
                                  const PolygonRange& polygons,
                                  PolygonMesh& out)
{
  return polygon_soup_to_polygon_mesh(points, polygons, out, parameters::all_default(), parameters::all_default());
}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_POLYGON_SOUP_TO_POLYGON_MESH_H
