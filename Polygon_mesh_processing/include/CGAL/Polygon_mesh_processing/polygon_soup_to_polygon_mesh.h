// Copyright (c) 2009-2013 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Laurent Rineau and Ilker O. Yaz

#ifndef CGAL_POLYGON_MESH_PROCESSING_POLYGON_SOUP_TO_POLYGON_MESH
#define CGAL_POLYGON_MESH_PROCESSING_POLYGON_SOUP_TO_POLYGON_MESH

#include <CGAL/license/Polygon_mesh_processing/repair.h>


#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/property_map.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/algorithm.h>
#include <set>
#include <boost/dynamic_bitset.hpp>

#include <boost/range/size.hpp>
#include <boost/range/value_type.hpp>
#include <boost/range/reference.hpp>

namespace CGAL
{
namespace Polygon_mesh_processing
{
namespace internal
{
template <typename PM
        , typename Point
        , typename Polygon>
class Polygon_soup_to_polygon_mesh
{
  const std::vector<Point>& _points;
  const std::vector<Polygon>& _polygons;

  typedef typename boost::property_map<PM, CGAL::vertex_point_t>::type Vpmap;
  typedef typename boost::property_traits<Vpmap>::value_type Point_3;

  typedef typename boost::graph_traits<PM>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PM>::face_descriptor face_descriptor;

public:
  /**
  * The constructor for modifier object.
  * @param points points of the soup of polygons.
  * @param polygons each element in the vector describes a polygon using the index of the points in the vector.
  */
  Polygon_soup_to_polygon_mesh(const std::vector<Point>& points,
                               const std::vector<Polygon>& polygons)
    : _points(points),
      _polygons(polygons)
  { }

  void operator()(PM& pmesh, const bool insert_isolated_vertices = true)
  {
    Vpmap vpmap = get(CGAL::vertex_point, pmesh);

    boost::dynamic_bitset<> not_isolated;
    if (!insert_isolated_vertices)
    {
      not_isolated.resize(_points.size());
      for (std::size_t i = 0, end = _polygons.size(); i < end; ++i)
      {
        const Polygon& polygon = _polygons[i];
        const std::size_t size = polygon.size();
        for (std::size_t j = 0; j < size; ++j)
          not_isolated.set(polygon[j], true);
      }
    }

    std::vector<vertex_descriptor> vertices(_points.size());
    for (std::size_t i = 0, end = _points.size(); i < end; ++i)
    {
      if (!insert_isolated_vertices && !not_isolated.test(i))
        continue;

      Point_3 pi(_points[i][0], _points[i][1], _points[i][2]);
      vertices[i] = add_vertex(pmesh);
      put(vpmap, vertices[i], pi);
    }

    for (std::size_t i = 0, end = _polygons.size(); i < end; ++i)
    {
      const Polygon& polygon = _polygons[i];
      const std::size_t size = polygon.size();

      std::vector<vertex_descriptor> vr(size); //vertex range
      vr.resize(size);
      for (std::size_t j = 0; j < size; ++j)
        vr[j] = vertices[polygon[j] ];

      CGAL_assertion_code(face_descriptor fd = )
      CGAL::Euler::add_face(vr, pmesh);
      CGAL_assertion(fd != boost::graph_traits<PM>::null_face());
    }
  }
};
}//end namespace internal


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
  template<class PolygonRange>
  bool is_polygon_soup_a_polygon_mesh(const PolygonRange& polygons)
  {
    typedef typename boost::range_value<
      typename boost::range_value<
        PolygonRange>::type >::type V_ID;
    typedef typename boost::range_value<
      PolygonRange>::type           Polygon;

    //check there is no duplicated ordered edge, and
    //check there is no polygon with twice the same vertex
    std::set< std::pair<V_ID, V_ID> > edge_set;
    V_ID max_id=0;
    BOOST_FOREACH(const Polygon& polygon, polygons)
    {
      std::size_t nb_edges = boost::size(polygon);
      if (nb_edges<3) return false;

      std::set<V_ID> polygon_vertices;
      V_ID prev= *cpp11::prev(boost::end(polygon));
      BOOST_FOREACH(V_ID id, polygon)
      {
        if (max_id<id) max_id=id;
        if (! edge_set.insert(std::pair<V_ID, V_ID>(prev,id)).second )
          return false;
        else
          prev=id;

        if (!polygon_vertices.insert(id).second)
          return false;//vertex met twice in the same polygon
      }
    }

    //check manifoldness
    typedef std::vector<V_ID> PointRange;
    typedef internal::Polygon_soup_orienter<PointRange, PolygonRange> Orienter;
    typename Orienter::Edge_map edges;
    typename Orienter::Marked_edges marked_edges;
    Orienter::fill_edge_map(edges, marked_edges, polygons);
    //returns false if duplication is necessary
    if (!marked_edges.empty())
      return false;
    return Orienter::has_singular_vertices(static_cast<std::size_t>(max_id+1),polygons,edges,marked_edges);
  }

  /**
  * \ingroup PMP_repairing_grp
  * builds a polygon mesh from a soup of polygons.
  * @pre the input polygon soup describes a consistently oriented
  * polygon mesh.
  *
  * @tparam PolygonMesh a model of `MutableFaceGraph` with an internal point property map
  * @tparam Point a point type that has an operator `[]` to access coordinates
  * @tparam Polygon a `std::vector<std::size_t>` containing the indices
  *         of the points of the face
  *
  * @param points points of the soup of polygons
  * @param polygons each element in the vector describes a polygon using the index of the points in `points`
  * @param out the polygon mesh to be built
  *
  * @pre `CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(polygons)`
  *
  * \sa `CGAL::Polygon_mesh_processing::orient_polygon_soup()`
  * \sa `CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh()`
  *
  */
  template<class PolygonMesh, class Point, class Polygon>
  void polygon_soup_to_polygon_mesh(
    const std::vector<Point>& points,
    const std::vector<Polygon>& polygons,
    PolygonMesh& out)
  {
    CGAL_precondition_msg(is_polygon_soup_a_polygon_mesh(polygons),
                          "Input soup needs to be a polygon mesh!");

    internal::Polygon_soup_to_polygon_mesh<PolygonMesh, Point, Polygon>
      converter(points, polygons);
    converter(out);
  }

}//end namespace Polygon_mesh_processing

}// end namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_POLYGON_SOUP_TO_POLYGON_MESH
