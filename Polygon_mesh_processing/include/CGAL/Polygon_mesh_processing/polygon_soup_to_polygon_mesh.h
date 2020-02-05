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

#ifndef CGAL_POLYGON_MESH_PROCESSING_POLYGON_SOUP_TO_POLYGON_MESH
#define CGAL_POLYGON_MESH_PROCESSING_POLYGON_SOUP_TO_POLYGON_MESH

#include <CGAL/license/Polygon_mesh_processing/repair.h>

#include <CGAL/disable_warnings.h>

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

#include <set>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

// \ingroup PMP_repairing_grp
//
// Adds the vertices and faces of a mesh into a (possibly non-empty) soup.
//
// @tparam PolygonMesh a model of `FaceListGraph` with an internal point property map
// @tparam PointRange a model of the concepts `RandomAccessContainer` and
//                    `BackInsertionSequence` whose value type is the point type
// @tparam PolygonRange a model of the concepts `RandomAccessContainer` and `BackInsertionSequence` whose
//                      `value_type` is itself a model of the concepts `RandomAccessContainer` and
//                      `BackInsertionSequence` whose `value_type` is `std::size_t`.
//
// @param mesh the mesh whose faces are being put in the soup
// @param points points of the soup of polygons
// @param polygons each element in the vector describes a polygon using the index of the points in `points`
//
// \sa `CGAL::Polygon_mesh_processing::orient_polygon_soup()`
// \sa `CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh()`
// \sa `CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh()`
//
template<typename PolygonMesh, typename PointRange, typename PolygonRange, typename NamedParameters>
void polygon_mesh_to_polygon_soup(const PolygonMesh& mesh,
                                  PointRange& points,
                                  PolygonRange& polygons,
                                  const NamedParameters& np)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor              vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor            halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor                face_descriptor;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type      VPM;
  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(vertex_point, mesh));

  typedef CGAL::dynamic_vertex_property_t<std::size_t>                              Vertex_index;
  typedef typename boost::property_map<PolygonMesh, Vertex_index>::type             VIM;
  VIM vim = get(Vertex_index(), mesh);

  typedef typename boost::range_value<PolygonRange>                                 Polygon;

  std::size_t index = points.size(); // so that multiple meshes can be put into the same soup

  points.reserve(points.size() + vertices(mesh).size());
  polygons.reserve(polygons.size() + faces(mesh).size());

  for(const vertex_descriptor v : vertices(mesh))
  {
    points.push_back(get(vpm, v));
    put(vim, v, index++);
  }

  for(const face_descriptor f : faces(mesh))
  {
    Polygon polygon;
    for(halfedge_descriptor h : CGAL::halfedges_around_face(halfedge(f, mesh), mesh))
      polygon.push_back(get(vim, target(h, mesh)));

    polygons.push_back(polygon);
  }
}

template<typename PolygonMesh, typename PointRange, typename PolygonRange>
void polygon_mesh_to_polygon_soup(const PolygonMesh& mesh,
                                  PointRange& points,
                                  PolygonRange& polygons)
{
  return polygon_mesh_to_polygon_soup(mesh, points, polygons, CGAL::parameters::all_default());
}

// -------------------------------------------------------------------------------------------------

template <typename PM
        , typename PointRange
        , typename PolygonRange>
class Polygon_soup_to_polygon_mesh
{
  const PointRange& _points;
  const PolygonRange& _polygons;

  typedef typename boost::property_map<PM, CGAL::vertex_point_t>::type Vpmap;
  typedef typename boost::property_traits<Vpmap>::value_type Point_3;

  typedef typename boost::graph_traits<PM>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PM>::face_descriptor face_descriptor;
  typedef typename PolygonRange::value_type Polygon;
  typedef typename PointRange::value_type Point;

public:
  /**
  * The constructor for modifier object.
  * @param points points of the soup of polygons.
  * @param polygons each element in the range describes a polygon using the index of the points in the range.
  */
  Polygon_soup_to_polygon_mesh(const PointRange& points,
                               const PolygonRange& polygons)
    : _points(points),
      _polygons(polygons)
  { }

  void operator()(PM& pmesh, const bool insert_isolated_vertices = true)
  {
    reserve(pmesh, static_cast<typename boost::graph_traits<PM>::vertices_size_type>(_points.size()),
                   static_cast<typename boost::graph_traits<PM>::edges_size_type>(2*_polygons.size()),
                   static_cast<typename boost::graph_traits<PM>::faces_size_type>(_polygons.size()) );

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

    if(boost::begin(polygons) == boost::end(polygons)){
      return true;
    }
    //check there is no duplicated ordered edge, and
    //check there is no polygon with twice the same vertex
    std::set< std::pair<V_ID, V_ID> > edge_set;
    V_ID max_id=0;
    for(const Polygon& polygon : polygons)
    {
      std::size_t nb_edges = boost::size(polygon);
      if (nb_edges<3) return false;

      std::set<V_ID> polygon_vertices;
      V_ID prev= *std::prev(boost::end(polygon));
      for(V_ID id : polygon)
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
    typename Orienter::Edge_map edges(max_id+1);
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
  * @tparam PolygonMesh a model of `MutableFaceGraph` with an internal point 
  * property map
  * @tparam PointRange a model of the concepts `RandomAccessContainer` and 
  * `BackInsertionSequence` whose value type is the point type
  * @tparam PolygonRange a model of the concept `RandomAccessContainer` whose 
  * `value_type` is a model of the concept `RandomAccessContainer` whose `value_type` is `std::size_t`.

  *
  * @param points points of the soup of polygons
  * @param polygons each element in the vector describes a polygon using the index of the points in `points`
  * @param out the polygon mesh to be built
  *
  * @pre \link CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh()
  *            CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(polygons) \endlink
  *
  * \sa `CGAL::Polygon_mesh_processing::orient_polygon_soup()`
  * \sa `CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh()`
  *
  */
  template<class PolygonMesh, class PointRange, class PolygonRange>
  void polygon_soup_to_polygon_mesh(
    const PointRange& points,
    const PolygonRange& polygons,
    PolygonMesh& out)
  {
    CGAL_precondition_msg(is_polygon_soup_a_polygon_mesh(polygons),
                          "Input soup needs to be a polygon mesh!");

    internal::Polygon_soup_to_polygon_mesh<PolygonMesh, PointRange, PolygonRange>
      converter(points, polygons);
    converter(out);
  }

}//end namespace Polygon_mesh_processing

}// end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_POLYGON_MESH_PROCESSING_POLYGON_SOUP_TO_POLYGON_MESH
