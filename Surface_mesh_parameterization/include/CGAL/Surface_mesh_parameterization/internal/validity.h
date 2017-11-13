// Copyright (c) 2016  GeometryFactory (France).
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
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_INTERNAL_VALIDITY_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_INTERNAL_VALIDITY_H

#include <CGAL/license/Surface_mesh_parameterization.h>

#include <CGAL/Surface_mesh_parameterization/internal/Containers_filler.h>
#include <CGAL/Surface_mesh_parameterization/internal/kernel_traits.h>

#include <CGAL/Bbox_2.h>
#include <CGAL/box_intersection_d.h>

#include <CGAL/boost/graph/properties.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/intersections.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>

#include <boost/foreach.hpp>
#include <boost/function_output_iterator.hpp>

#include <vector>

namespace CGAL {

namespace Surface_mesh_parameterization {

namespace internal {

template<typename TriangleMesh, typename VertexUVMap>
bool has_flips(const TriangleMesh& mesh,
               typename boost::graph_traits<TriangleMesh>::halfedge_descriptor bhd,
               const VertexUVMap uvmap)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor      face_descriptor;

  // Kernel subtypes:
  typedef typename internal::Kernel_traits<TriangleMesh>::Kernel      Kernel;
  typedef typename Kernel::Point_2                                    Point_2;
  typedef typename Kernel::Point_3                                    Point_3;
  typedef typename Kernel::Vector_3                                   Vector_3;

  // Fill containers
  boost::unordered_set<vertex_descriptor> vertices;
  std::vector<face_descriptor> faces;

  internal::Containers_filler<TriangleMesh> fc(mesh, vertices, &faces);
  Polygon_mesh_processing::connected_component(
                                    face(opposite(bhd, mesh), mesh),
                                    mesh,
                                    boost::make_function_output_iterator(fc));

  Vector_3 first_triangle_normal(0., 0., 0.);
  bool is_normal_set = false;

  BOOST_FOREACH(face_descriptor fd, faces) {
    // Get 3 vertices of the facet
    halfedge_descriptor hd = halfedge(fd, mesh);
    vertex_descriptor vd0 = target(hd, mesh);
    vertex_descriptor vd1 = target(next(hd, mesh), mesh);
    vertex_descriptor vd2 = source(hd, mesh);

    // Get the 3 vertices position in 2D
    const Point_2& p0 = get(uvmap, vd0);
    const Point_2& p1 = get(uvmap, vd1);
    const Point_2& p2 = get(uvmap, vd2);

    // Compute the facet normal
    Point_3 p0_3D(p0.x(), p0.y(), 0.);
    Point_3 p1_3D(p1.x(), p1.y(), 0.);
    Point_3 p2_3D(p2.x(), p2.y(), 0.);
    Vector_3 v01_3D = p1_3D - p0_3D;
    Vector_3 v02_3D = p2_3D - p0_3D;
    Vector_3 normal = CGAL::cross_product(v01_3D, v02_3D);

    // Check that all normals are oriented the same way
    if (!is_normal_set) {
      first_triangle_normal = normal;
      is_normal_set = true;
    } else {
      if (first_triangle_normal * normal < 0)
        return true;
    }
  }
  return false;
}

template <typename TriangleMesh, typename VertexUVMap>
class Intersect_facets
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor      face_descriptor;

  // Kernel subtypes:
  typedef typename internal::Kernel_traits<TriangleMesh>::Kernel    Kernel;
  typedef typename Kernel::Point_2                                  Point_2;
  typedef typename Kernel::Segment_2                                Segment_2;
  typedef typename Kernel::Triangle_2                               Triangle_2;
  typedef typename Kernel::FT                                       NT;

  typename Kernel::Construct_segment_2 segment_functor;
  typename Kernel::Construct_triangle_2 triangle_functor;
  typename Kernel::Do_intersect_2 do_intersect_2_functor;

  typedef CGAL::Box_intersection_d::Box_with_info_d<NT, 2, face_descriptor> Box;

  const TriangleMesh& mesh;
  const VertexUVMap uvmap;
  unsigned int& self_intersection_counter;

public:
  unsigned int number_of_self_intersections() const { return self_intersection_counter; }

  // callback functor that reports all truly intersecting triangles
  void operator()(const Box* a, const Box* b) const
  {
    // Boxes intersect, need to check if there is actually an intersection
    // and if it is not a subface of the faces
    halfedge_descriptor h = halfedge(a->info(), mesh);
    halfedge_descriptor g = halfedge(b->info(), mesh);

    // check for shared egde
    if(face(opposite(h, mesh), mesh) == b->info() ||
       face(opposite(prev(h, mesh), mesh), mesh) == b->info() ||
       face(opposite(next(h, mesh), mesh), mesh) == b->info()) {
      // shared edge

      // intersection if the orientations are not identical
      if(CGAL::orientation(get(uvmap, target(h, mesh)),
                           get(uvmap, target(next(h, mesh), mesh)),
                           get(uvmap, source(h, mesh))) !=
         CGAL::orientation(get(uvmap, target(g, mesh)),
                           get(uvmap, target(next(g, mesh), mesh)),
                           get(uvmap, source(g, mesh)))) {
        ++self_intersection_counter;
      }
      return;
    }

    // check for shared vertex --> possible intersection
    halfedge_descriptor hd;

    if(target(h, mesh) == target(g, mesh))
      hd = g;
    if(target(h, mesh) == target(next(g, mesh), mesh))
      hd = next(g, mesh);
    if(target(h, mesh) == target(next(next(g, mesh), mesh), mesh))
      hd = next(next(g, mesh), mesh);

    if(target(h, mesh) == target(g, mesh))
      hd = g;
    if(target(h, mesh) == target(next(g, mesh), mesh))
      hd = next(g, mesh);
    if(target(h, mesh) == target(next(next(g, mesh), mesh), mesh))
      hd = next(next(g, mesh), mesh);

    if(hd == halfedge_descriptor()) {
      h = next(h, mesh);
      if(target(h, mesh) == target(g, mesh))
        hd = g;
      if(target(h, mesh) == target(next(g, mesh), mesh))
        hd = next(g, mesh);
      if(target(h, mesh) == target(next(next(g, mesh), mesh), mesh))
        hd = next(next(g, mesh), mesh);
      if(hd == halfedge_descriptor()) {
        h = next(h, mesh);
        if(target(h, mesh) == target(g, mesh))
          hd = g;
        if(target(h, mesh) == target(next(g, mesh), mesh))
          hd = next(g, mesh);
        if(target(h, mesh) == target(next(next(g, mesh), mesh), mesh))
          hd = next(next(g, mesh), mesh);
      }
    }

    if(hd != halfedge_descriptor()) {
      // shared vertex
      CGAL_assertion(target(h, mesh) == target(hd, mesh));

      // geometric check if the opposite segments intersect the triangles
      Triangle_2 t1 = triangle_functor(get(uvmap, target(h, mesh)),
                                       get(uvmap, target(next(h, mesh), mesh)),
                                       get(uvmap, target(next(next(h, mesh), mesh), mesh)));
      Triangle_2 t2 = triangle_functor(get(uvmap, target(hd, mesh)),
                                       get(uvmap, target(next(hd, mesh), mesh)),
                                       get(uvmap, target(next(next(hd, mesh), mesh), mesh)));

      Segment_2 s1 = segment_functor(get(uvmap, target(next(h, mesh), mesh)),
                                     get(uvmap, target(next(next(h, mesh), mesh), mesh)));
      Segment_2 s2 = segment_functor(get(uvmap, target(next(hd, mesh), mesh)),
                                     get(uvmap, target(next(next(hd, mesh), mesh), mesh)));

      if(do_intersect_2_functor(t1, s2)) {
        ++self_intersection_counter;
      } else if(do_intersect_2_functor(t2, s1)) {
        ++self_intersection_counter;
      }
      return;
    }

    // check for geometric intersection
    Triangle_2 t1 = triangle_functor(get(uvmap, target(h, mesh)),
                                     get(uvmap, target(next(h, mesh), mesh)),
                                     get(uvmap, target(next(next(h, mesh), mesh), mesh)));
    Triangle_2 t2 = triangle_functor(get(uvmap, target(g, mesh)),
                                     get(uvmap, target(next(g, mesh), mesh)),
                                     get(uvmap, target(next(next(g, mesh), mesh), mesh)));
    if(do_intersect_2_functor(t1, t2)) {
      ++self_intersection_counter;
    }
  }

  Intersect_facets(const TriangleMesh& mesh_,
                   const VertexUVMap uvmap_,
                   unsigned int& counter)
    :
      mesh(mesh_),
      uvmap(uvmap_),
      self_intersection_counter(counter)
  { }
};

/// Check if the 3D -> 2D mapping is one-to-one.
/// This function is stronger than "has_flips()" because the parameterized
/// surface can loop over itself without creating any flips.
template <typename TriangleMesh,
          typename Faces_Container,
          typename VertexUVMap>
bool is_one_to_one_mapping(const TriangleMesh& mesh,
                           const Faces_Container& faces,
                           const VertexUVMap uvmap)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor      face_descriptor;

  // Kernel subtypes:
  typedef typename internal::Kernel_traits<TriangleMesh>::Kernel      Kernel;
  typedef typename Kernel::FT                                         NT;
  typedef typename Kernel::Point_2                                    Point_2;

  typedef CGAL::Box_intersection_d::Box_with_info_d<NT, 2, face_descriptor> Box;

  // Create the corresponding vector of bounding boxes
  std::vector<Box> boxes;

  BOOST_FOREACH(face_descriptor fd, faces) {
    halfedge_descriptor hd = halfedge(fd, mesh);
    vertex_descriptor vd0 = target(hd, mesh);
    vertex_descriptor vd1 = target(next(hd, mesh), mesh);
    vertex_descriptor vd2 = source(hd, mesh);

    // Get the 3 vertices position in 2D
    const Point_2& p0 = get(uvmap, vd0);
    const Point_2& p1 = get(uvmap, vd1);
    const Point_2& p2 = get(uvmap, vd2);

    Bbox_2 b = p0.bbox();
    b += p1.bbox();
    b += p2.bbox();

    boxes.push_back(Box(b, fd));
  }

  std::vector<const Box*> boxes_ptr;
  boxes_ptr.reserve(boxes.size());

  BOOST_FOREACH(Box& b, boxes)
    boxes_ptr.push_back(&b);

  // Run the self intersection algorithm with all defaults
  unsigned int counter = 0;
  Intersect_facets<TriangleMesh, VertexUVMap> intersect_facets(mesh, uvmap, counter);
  std::ptrdiff_t cutoff = 2000;
  CGAL::box_self_intersection_d(boxes_ptr.begin(), boxes_ptr.end(),
                                intersect_facets, cutoff);
  return (counter == 0);
}

template <typename TriangleMesh,
          typename VertexUVMap>
bool is_one_to_one_mapping(const TriangleMesh& mesh,
                           typename boost::graph_traits<TriangleMesh>::halfedge_descriptor bhd,
                           const VertexUVMap uvmap)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor  vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor    face_descriptor;

  boost::unordered_set<vertex_descriptor> vertices;
  std::vector<face_descriptor> faces;

  internal::Containers_filler<TriangleMesh> fc(mesh, vertices, &faces);
  Polygon_mesh_processing::connected_component(
                                    face(opposite(bhd, mesh), mesh),
                                    mesh,
                                    boost::make_function_output_iterator(fc));

  return is_one_to_one_mapping(mesh, faces, uvmap);
}

} // namespace internal

} // namespace Surface_mesh_parameterization

} // namespace CGAL

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_INTERNAL_VALIDITY_H
