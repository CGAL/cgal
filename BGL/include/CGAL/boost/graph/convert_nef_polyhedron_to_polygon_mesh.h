// Copyright (c) 2016 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Sebastien Loriot

#ifndef CGAL_BOOST_GRAPH_NEF_POLYHEDRON_TO_POLYGON_MESH_H
#define CGAL_BOOST_GRAPH_NEF_POLYHEDRON_TO_POLYGON_MESH_H

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/algorithm.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/circulator.h>
#include <CGAL/Cartesian_converter.h>
#include <boost/unordered_map.hpp>

namespace CGAL{

namespace nef_to_pm{

// Visitor used to collect and index vertices of a shell
template<class Nef_polyhedron, class Point_3, class Converter>
struct Shell_vertex_index_visitor
{
  typedef boost::unordered_map<typename Nef_polyhedron::Vertex_const_handle, std::size_t> Vertex_index_map;
  std::vector<Point_3>& points;
  Vertex_index_map vertex_indices;
  const Converter& converter;

  Shell_vertex_index_visitor(std::vector<Point_3>& points, const Converter& converter)
    :points(points), converter(converter)
  {}

  void visit(typename Nef_polyhedron::Vertex_const_handle vh)
  {
    std::pair<typename Vertex_index_map::iterator, bool> insert_res =
      vertex_indices.insert( std::make_pair(vh, points.size()) );
    if (insert_res.second)
      points.push_back( converter(vh->point()));
  }
  void visit(typename Nef_polyhedron::Halfedge_const_handle )
  {}
  void visit(typename Nef_polyhedron::Halffacet_const_handle )
  {}
  void visit(typename Nef_polyhedron::SHalfedge_const_handle )
  {}
  void visit(typename Nef_polyhedron::SHalfloop_const_handle )
  {}
  void visit(typename Nef_polyhedron::SFace_const_handle )
  {}
};

//Visitor used to collect polygons
template <class Nef_polyhedron>
struct Shell_polygons_visitor
{
  typedef boost::unordered_map<typename Nef_polyhedron::Vertex_const_handle, std::size_t> Vertex_index_map;
  Vertex_index_map& vertex_indices;
  std::vector< std::vector<std::size_t> >& polygons;

  Shell_polygons_visitor(Vertex_index_map& vertex_indices,
                         std::vector< std::vector<std::size_t> >& polygons)
    : vertex_indices( vertex_indices )
    , polygons(polygons)
  {}

  void visit(typename Nef_polyhedron::Halffacet_const_handle opposite_facet)
  {
    bool is_marked=opposite_facet->incident_volume()->mark();

    CGAL_assertion(Nef_polyhedron::Infi_box::is_standard(opposite_facet->plane()));

    typename Nef_polyhedron::SHalfedge_const_handle se;
    typename Nef_polyhedron::Halffacet_cycle_const_iterator fc;

    typename Nef_polyhedron::Halffacet_const_handle f = opposite_facet->twin();

    typename Nef_polyhedron::SHalfedge_around_facet_const_circulator
      sfc1(f->facet_cycles_begin());

    // create a new polygon
    polygons.push_back( std::vector<std::size_t>() );
    fc = f->facet_cycles_begin();
    se = typename Nef_polyhedron::SHalfedge_const_handle(fc);
    CGAL_assertion(se!=0);
    typename Nef_polyhedron::SHalfedge_around_facet_const_circulator
      hc_start(se), hc_end(hc_start);
    // insert the vertex indices of the new polygon
    CGAL_For_all(hc_start,hc_end)
      polygons.back().push_back(vertex_indices[hc_start->source()->center_vertex()]);
    if (!is_marked)
      std::reverse(polygons.back().begin(), polygons.back().end());
  }

  void visit(typename Nef_polyhedron::SFace_const_handle)
  {}
  void visit(typename Nef_polyhedron::Halfedge_const_handle)
  {}
  void visit(typename Nef_polyhedron::Vertex_const_handle)
  {}
  void visit(typename Nef_polyhedron::SHalfedge_const_handle)
  {}
  void visit(typename Nef_polyhedron::SHalfloop_const_handle)
  {}
};

template <class Point_3, class Nef_polyhedron, class Converter>
void collect_polygon_mesh_info(
  std::vector<Point_3>& points,
  std::vector< std::vector<std::size_t> >& polygons,
  Nef_polyhedron& nef,
  typename Nef_polyhedron::Shell_entry_const_iterator shell,
  const Converter& converter)
{
  // collect points and set vertex indices
  Shell_vertex_index_visitor<Nef_polyhedron, Point_3, Converter> vertex_index_visitor(points, converter);
  nef.visit_shell_objects(typename Nef_polyhedron::SFace_const_handle(shell), vertex_index_visitor);

  // collect polygons
  Shell_polygons_visitor<Nef_polyhedron> polygon_visitor(
    vertex_index_visitor.vertex_indices,
    polygons);
  nef.visit_shell_objects(typename Nef_polyhedron::SFace_const_handle(shell), polygon_visitor);
}

} //end of namespace nef_to_pm

/// \ingroup PkgBGL
/// Converts an objet of type `Nef_polyhedron_3` into a polygon mesh model of `MutableFaceGraph`.
/// Note that contrary to `Nef_polyhedron_3::convert_to_polyhedron()`, the output is not triangulated.
/// The polygon mesh can be triangulated using the function `triangulate_faces()`.
/// \pre `Polygon_mesh` must have an internal point property map with value type being `Nef_polyhedron_3::Point_3`.
/// \pre `nef.simple()`
template <class Nef_polyhedron, class Polygon_mesh>
void convert_nef_polyhedron_to_polygon_mesh(const Nef_polyhedron& nef, Polygon_mesh& pm)
{
  typedef typename Nef_polyhedron::Point_3 Point_3;
  typedef typename boost::property_traits<typename boost::property_map<Polygon_mesh, vertex_point_t>::type>::value_type PM_Point;
  typedef typename Kernel_traits<PM_Point>::Kernel PM_Kernel;
  typedef typename Kernel_traits<Point_3>::Kernel Nef_Kernel;
  typedef Cartesian_converter<Nef_Kernel, PM_Kernel> Converter;

  typename Nef_polyhedron::Volume_const_iterator vol_it = nef.volumes_begin(),
                                                 vol_end = nef.volumes_end();
  if ( Nef_polyhedron::Infi_box::extended_kernel() ) ++vol_it; // skip Infi_box
  CGAL_assertion ( vol_it != vol_end );
  ++vol_it; // skip unbounded volume

  std::vector<PM_Point> points;
  std::vector< std::vector<std::size_t> > polygons;
  Converter to_inexact;
  for (;vol_it!=vol_end;++vol_it)
    nef_to_pm::collect_polygon_mesh_info<PM_Point>(points,
                                                  polygons,
                                                  nef,
                                                  vol_it->shells_begin(),
                                                  to_inexact);

  Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, pm);
}

} //end of namespace CGAL



#endif // CGAL_BOOST_GRAPH_NEF_POLYHEDRON_TO_POLYGON_MESH_H
