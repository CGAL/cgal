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

#include <CGAL/Modifier_base.h>
#include <CGAL/IO/generic_print_polyhedron.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/property_map.h>

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
  : public CGAL::Modifier_base< PM >
{
  const std::vector<Point>& _points;
  const std::vector<Polygon>& _polygons;

  typedef typename PM::Point Point_3;

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

  void operator()(PM& pmesh)
  {
    typename boost::property_map<PM, CGAL::vertex_point_t>::type
      vpmap = get(CGAL::vertex_point, pmesh);

    std::vector<vertex_descriptor> vertices(_points.size());
    for (std::size_t i = 0, end = _points.size(); i < end; ++i)
    {
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


  /// \cond SKIP_IN_MANUAL
  /**
  * \ingroup PkgPolygonMeshProcessing
  * returns `true` if the soup of polygons defines a valid polygon mesh
  * that can be handled by `CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh()`.
  *
  * @tparam Polygon a `std::vector<std::size_t>` containing the indices
  *         of the points of the polygon face
  *
  * @param polygons each element in the vector describes a polygon using the index of the vertices
  *
  */
  template<class Polygon>
  bool is_polygon_soup_a_polygon_mesh(const std::vector<Polygon>& polygons)
  {
    typedef typename std::iterator_traits<
              typename Polygon::iterator >::value_type                   V_ID;

    std::set< std::pair<V_ID, V_ID> > edge_set;
    BOOST_FOREACH(const Polygon& polygon, polygons)
    {
      std::size_t nb_edges = polygon.size();
      if (nb_edges<3) return false;
      V_ID prev=polygon.back();
      BOOST_FOREACH(V_ID id, polygon)
        if (! edge_set.insert(std::pair<V_ID, V_ID>(prev,id)).second )
          return false;
        else
          prev=id;
    }

    return true;
  }
  /// \endcond

  /**
  * \ingroup PkgPolygonMeshProcessing
  * builds a polygon mesh from a soup of polygons.
  * @pre the input polygon soup describes consistently oriented
  * polygon mesh.
  *
  * @tparam PolygonMesh a model of `MutableFaceGraph`
  * @tparam Point a point type that has an operator `[]` to access coordinates
  * @tparam Polygon a `std::vector<std::size_t>` containing the indices
  *         of the points of the face
  *
  * @param points points of the soup of polygons
  * @param polygons each element in the vector describes a polygon using the index of the points in `points`
  * @param out the polygon mesh to be built
  *
  * \sa `CGAL::Polygon_mesh_processing::orient_polygon_soup()`
  *
  */
  template<class PolygonMesh, class Point, class Polygon>
  void polygon_soup_to_polygon_mesh(
    const std::vector<Point>& points,
    const std::vector<Polygon>& polygons,
    PolygonMesh& out)
  {
    internal::Polygon_soup_to_polygon_mesh<PolygonMesh, Point, Polygon>
      converter(points, polygons);
    converter(out);
  }

}//end namespace Polygon_mesh_processing

}// end namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_POLYGON_SOUP_TO_POLYGON_MESH
