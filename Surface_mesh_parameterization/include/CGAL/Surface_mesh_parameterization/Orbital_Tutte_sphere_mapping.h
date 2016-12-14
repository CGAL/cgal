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
// $URL$
// $Id$
//
//
// Author(s)     :

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_ORBITAL_TUTTE_SPHERE_MAPPING_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_ORBITAL_TUTTE_SPHERE_MAPPING_H

#include <CGAL/Surface_mesh_parameterization/internal/orbital_cone_helper.h>
#include <CGAL/Surface_mesh_parameterization/internal/kernel_traits.h>

#include <CGAL/Arr_default_overlay_traits.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/intersection_2.h>
#include <CGAL/Triangulation_2.h>

#include <Eigen/Dense>

#include <boost/iterator/transform_iterator.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/unordered_map.hpp>

#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <vector>

namespace CGAL {

namespace Surface_mesh_parameterization {

namespace internal {

/// Output a container of triangles as a .off file.
template<typename Triangle_container>
void output_triangle_set(Triangle_container& triangles,
                         const char* filename)
{
  typedef typename Triangle_container::value_type      Triangle_2;
  typedef typename Triangle_2::R                       Kernel;
  typedef typename Kernel::Point_2                     Point_2;

  const std::size_t triangles_n = triangles.size();

  std::map<Point_2, int> points_to_indices;
  std::vector<Point_2> ids_to_print;
  ids_to_print.reserve(3 * triangles_n);

  std::ofstream out(filename);
  std::ostringstream faces_ss;

  int counter = 0;
  typename Triangle_container::const_iterator trit = triangles.begin(),
                                              trend = triangles.end();
  for(; trit!=trend; ++trit) {
    const Triangle_2& tr = *trit;
    std::vector<int> ids(3);

    for(std::size_t i=0; i<3; ++i) {
      std::pair<typename std::map<Point_2, int>::iterator, bool> is_insert_successful =
          points_to_indices.insert(std::make_pair(tr[i], counter));
      if(!is_insert_successful.second) {
        ids[i] = is_insert_successful.first->second;
      } else {
        ids[i] = counter;
        ids_to_print.push_back(tr[i]);
        ++counter;
      }
    }

    faces_ss << "3 " << ids[0] << " "<< ids[1] << " " << ids[2] << '\n';
  }

  // OFF header
  out << "OFF" << '\n';
  out << points_to_indices.size() << " " << triangles_n << " 0" << '\n';

  // outputs points
  typename std::vector<Point_2>::iterator it = ids_to_print.begin(),
                                          end = ids_to_print.end();
  for(; it!=end; ++it) {
    out << *it << " 0" << '\n';
  }

  // outputs faces
  out << faces_ss.str() << std::endl;
}

/// Outputs an arrangement as a .off file.
template<typename Arrangement>
void output_arrangement_to_off(const Arrangement& arr,
                               std::ofstream& out)
{
  typedef typename Arrangement::Vertex_const_iterator       Vertex_const_iterator;
  typedef typename Arrangement::Edge_const_iterator         Edge_const_iterator;
  typedef typename Arrangement::Point_2                     Exact_point_2;

  std::cout << "Arrangement with " << arr.number_of_vertices() << " vertices and "
                                   << arr.number_of_edges() << " edges" << std::endl;

  std::map<Exact_point_2, int> indices;

  std::ostringstream out_vertices;
  std::size_t vertices_counter = 0;

  // arrangement vertices
  Vertex_const_iterator vit = arr.vertices_begin(), vend = arr.vertices_end();
  for(; vit!=vend; ++vit) {
    std::pair<typename std::map<Exact_point_2, int>::iterator, bool> is_insert_successful =
                            indices.insert(std::make_pair(vit->point(), vertices_counter));
    if(is_insert_successful.second) {
      out_vertices << vit->point() << " 0" << '\n';
      ++vertices_counter;
    }
  }
  CGAL_assertion(arr.number_of_vertices() == vertices_counter);

  out << "OFF" << "\n" << std::endl;
  out << vertices_counter << " " << arr.number_of_edges() << " 0" << '\n';
  out << out_vertices.str();

  // arrangement edges
  Edge_const_iterator eit = arr.edges_begin(), eend = arr.edges_end();
  for(; eit!=eend; ++eit) {
    out << "3 " << indices[eit->curve().source()] << " "
                << indices[eit->curve().target()] << " "
                << indices[eit->curve().source()] << '\n';
  }

  out << std::endl;
}

/// A type to regroup all the info and avoid having to pass it all in each function.
template<typename SeamMesh_,
         typename ConeMap_,
         typename VertexIndexMap_,
         typename VertexUVMap_>
class Embedded_mesh
{
public:
  typedef SeamMesh_                   SeamMesh;
  typedef ConeMap_                    ConeMap;
  typedef VertexIndexMap_             VertexIndexMap;
  typedef VertexUVMap_                VertexUVMap;

  const SeamMesh& mesh;
  const ConeMap& cmap;
  const VertexIndexMap vimap;
  const VertexUVMap uvmap;
  const Orbifold_type orb_type;

  Embedded_mesh(const SeamMesh& mesh,
                const ConeMap& cmap,
                const VertexIndexMap vimap,
                const VertexUVMap uvmap,
                const Orbifold_type orb_type)
    :
      mesh(mesh),
      cmap(cmap),
      vimap(vimap),
      uvmap(uvmap),
      orb_type(orb_type)
  { }
};

/// Affine transformation to express A = T*B + V
// @fixme transf needs a translation component too (for orb type IV)
template<typename Kernel>
class Affine_transformation
{
public:
  typedef typename Kernel::FT                 NT;
  typedef typename Kernel::Point_2            Point_2;

  // rotation of 'angle' around 'center'
  Eigen::Matrix3d rot_matrix;

public:
  /// Apply the affine transformation to the point `p`.
  Point_2 transform(const Point_2& p) const
  {
    Eigen::Vector3d v, res;
    v(0) = p.x(); v(1) = p.y(); v(2) = 1;
    res = rot_matrix * v;
    Point_2 newp(res(0), res(1));
//    std::cout << "------------------------------  transformed " << p << " in " << newp << std::endl;
    return newp;
  }

  /// Multiply the current rotation by `new_rot` to get a single transformation matrix.
  void combine_rotation(const Eigen::Matrix3d& new_rot)
  {
    rot_matrix = rot_matrix * new_rot;
  }

  Affine_transformation(const NT angle, const Point_2& center)
  {
    const NT c = std::cos(angle);
    const NT s = std::sin(angle);
    const NT center_x = center.x();
    const NT center_y = center.y();
    rot_matrix <<  c, -s, center_x - c * center_x + s * center_y,
                   s,  c, center_y - s * center_x - c * center_y,
                   0, 0,                                       1;

    std::cout << "built matrix for " << center << " and angle: " << angle << std::endl;
    std::cout << rot_matrix << std::endl;
  }

  Affine_transformation()
  {
    rot_matrix << 1, 0, 0,
                  0, 1, 0,
                  0, 0, 1;
  }
};

/// halfedge that carries a transformation that is used to compute the coordinates
/// of the incident faces in the tiled space.
template<typename SeamMesh>
class halfedge_with_transformation
{
  typedef halfedge_with_transformation<SeamMesh>                          Self;

public:
  typedef SeamMesh                                                        Seam_mesh;
  typedef typename internal::Kernel_traits<Seam_mesh>::Kernel             Kernel;
  typedef Affine_transformation<Kernel>                                   Transformation;

  typedef typename Seam_mesh::TriangleMesh                                TriangleMesh;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor TM_halfedge_descriptor;

  TM_halfedge_descriptor tmhd;
  Transformation transformation;

  halfedge_with_transformation(TM_halfedge_descriptor tmhd,
                               const Transformation& transformation)
    :
      tmhd(tmhd), transformation(transformation)
  { }
};

} // namespace internal

template<typename Arrangement,
         typename EmbeddedMesh>
class Orbifold_sphere_mapper
{
  // Types related to the embedded mesh
  typedef typename EmbeddedMesh::ConeMap                                  ConeMap;
  typedef typename EmbeddedMesh::SeamMesh                                 Seam_mesh;
  typedef typename EmbeddedMesh::VertexIndexMap                           VertexIndexMap;
  typedef typename EmbeddedMesh::VertexUVMap                              VertexUVMap;

  // Graph types for the seam mesh
  typedef typename boost::graph_traits<Seam_mesh>::vertex_descriptor      vertex_descriptor;
  typedef typename boost::graph_traits<Seam_mesh>::halfedge_descriptor    halfedge_descriptor;
  typedef typename boost::graph_traits<Seam_mesh>::edge_descriptor        edge_descriptor;
  typedef typename boost::graph_traits<Seam_mesh>::face_descriptor        face_descriptor;
  typedef typename boost::graph_traits<Seam_mesh>::vertex_iterator        vertex_iterator;

  // Graph types for the underlying mesh
  typedef typename Seam_mesh::TriangleMesh                                TriangleMesh;
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor   TM_vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor TM_halfedge_descriptor;

  // Kernel types
  typedef typename internal::Kernel_traits<Seam_mesh>::Kernel             Kernel;
  typedef typename Kernel::FT                                             NT;
  typedef typename Kernel::Point_2                                        Point_2;
  typedef typename Kernel::Segment_2                                      Segment_2;
  typedef typename Kernel::Vector_2                                       Vector_2;
  typedef typename Kernel::Triangle_2                                     Triangle_2;

  // Arrangement types
  typedef typename Arrangement::Geometry_traits_2                         Arr_traits;
  typedef CGAL::Arr_default_overlay_traits<Arrangement>                   Overlay_traits;

  // Exact kernel types
  typedef typename Arr_traits::Kernel                                     Exact_kernel;
  typedef typename Arr_traits::X_monotone_curve_2                         Arr_segment_2;

  typedef typename Exact_kernel::FT                                       Exact_NT;
  typedef typename Exact_kernel::Point_2                                  Exact_point_2;
  typedef typename Exact_kernel::Segment_2                                Exact_segment_2;
  typedef CGAL::Cartesian_converter<Kernel, Exact_kernel>                 IK_to_EK;

  // Transformation types
  typedef internal::Affine_transformation<Kernel>                         Aff_tr;
  typedef internal::halfedge_with_transformation<Seam_mesh>               halfedge_w_tr;

  // Triangulation type
  typedef CGAL::Triangulation_2<Kernel>                                   Triangulation;

  IK_to_EK to_exact;

  /// For all seam edges, mark both halfedges with the transformation to apply
  /// to obtain the next tile's coordinates from the previous tile.
  template<typename Seam_halfedges_segment,
           typename SeamTransformationMap,
           typename BorderHalfedges>
  void tag_seam_segment(const EmbeddedMesh& emesh,
                        const Seam_halfedges_segment& seam_segment,
                        SeamTransformationMap& seam_transformations,
                        BorderHalfedges& border_halfedges,
                        NT ang)
  {
    std::cout << "Segment of length " << seam_segment.size() << std::endl;

    const Seam_mesh& mesh = emesh.mesh;
    const VertexIndexMap vimap = emesh.vimap;
    const VertexUVMap uvmap = emesh.uvmap;

    // s and t are the two sides of the seam
    vertex_descriptor vf_s = target(seam_segment.front().first, mesh);
    vertex_descriptor vf_t = source(seam_segment.front().second, mesh);
    vertex_descriptor ve_s = source(seam_segment.back().first, mesh);
    vertex_descriptor ve_t = target(seam_segment.back().second, mesh);

    bool is_reversed = (get(vimap, ve_s) == get(vimap, ve_t));

    if(is_reversed) {
      ang *= -1;
    }

    Point_2 pb = is_reversed ? get(uvmap, ve_s) : get(uvmap, vf_s);
    Point_2 pe = is_reversed ? get(uvmap, ve_t) : get(uvmap, vf_t);

    Aff_tr trs(+ 2 * CGAL_PI / ang /*angle*/, pb /*center*/);
    Aff_tr trt(- 2 * CGAL_PI / ang /*angle*/, pe /*center*/);
    Aff_tr identity;

    typename Seam_halfedges_segment::const_iterator it = seam_segment.begin(),
                                                    end = seam_segment.end();
    for(; it!=end; ++it) {
      // Associate to the halfedges on the seam a transformation that can be
      // thought of as:
      // "if I cross this border, which affine transformation must I apply to
      // a triangle of the initial mesh to obtain the coordinates of this triangle
      // in the next tile"
      halfedge_descriptor hd1 = it->first;
      TM_halfedge_descriptor base_hd1 = hd1.tmhd;
      seam_transformations[base_hd1] = trs;

      halfedge_descriptor hd2 = it->second;
      TM_halfedge_descriptor base_hd2 = hd2.tmhd;
      seam_transformations[base_hd2] = trt;

      // The initial border edges are by definition in the initial mesh so they
      // do not have a transformation at the beginning
      halfedge_w_tr hdt1(base_hd1, identity);
      halfedge_w_tr hdt2(base_hd2, identity);
      border_halfedges.push_back(hdt1);
      border_halfedges.push_back(hdt2);
    }
  }

  /// Go through the seam edges to collect the initial border halfedges and mark the
  /// transformations on the seam.
  template<typename SeamTransformationMap,
           typename BorderHalfedges>
  void compute_initial_border(const EmbeddedMesh& emesh,
                              SeamTransformationMap& seam_transformations,
                              BorderHalfedges& border_halfedges)
  {
    const Seam_mesh& mesh = emesh.mesh;

    // the angles at the cones
    typedef std::vector<NT>                           Angle_container;
    const Angle_container& angs = get_angles_at_cones<Angle_container>(emesh.orb_type);

    // Find the cone tagged 'First_unique_cone'
    int start_cone_index = -1; // index of the beginning of the seam (useless here)
    vertex_descriptor start_cone;
    internal::find_start_cone(emesh.cmap, emesh.vimap, start_cone, start_cone_index);

    // By property of the seam mesh, the canonical halfedge that points to start_cone
    // is on the seam, and is not on the border
    halfedge_descriptor hd = halfedge(start_cone, mesh);
    CGAL_precondition(mesh.has_on_seam(hd));
    halfedge_descriptor bhd = opposite(hd, mesh);
    CGAL_precondition(is_border(bhd, mesh));

    // To count the segments of the seam (a segment is a set of edges between two cones)
    // there are 3 segment max in orb types I, II, III and IV
    std::size_t segment_index = 0;

    // The segments between cones
    typedef std::vector<std::pair<halfedge_descriptor, halfedge_descriptor> >  Seam_halfedges_segment;
    Seam_halfedges_segment seam_segment;

    // Walk the border till we find the first cone again, while tagging the seam edges
    // with the correct transformation (such that if we go over the seam, we unfold
    // properly in the plane
    while(true) { // breaking at the last cone
      // get the two halfedges that are not on the border
      halfedge_descriptor hd1 = opposite(bhd, mesh);

      // The non-border halfedge with same vertices (in the underlying mesh of the seam
      // mesh) as bhd is simply bhd with the 'seam' boolean set to false
      halfedge_descriptor hd2 = bhd;
      hd2.seam = false;

      CGAL_assertion(segment_index < angs.size());
      NT ang = angs[segment_index];

      seam_segment.push_back(std::make_pair(hd1, hd2));

      // Check if we have reached the end of the seam
      vertex_descriptor bhd_target = target(bhd, mesh);
      typename ConeMap::const_iterator is_in_map = emesh.cmap.find(bhd_target);
      if(is_in_map != emesh.cmap.end()) {
        tag_seam_segment<Seam_halfedges_segment,
                         SeamTransformationMap,
                         BorderHalfedges>(emesh,
                                          seam_segment,
                                          seam_transformations,
                                          border_halfedges,
                                          ang);

        // reached the end of the seam
        if(is_in_map->second == Second_unique_cone) {
          break;
        }

        std::cout << "-------------------------" << std::endl;
        seam_segment.clear();
        ++segment_index;
      }

      // move to the next halfedge couple (walking on the border of the seam)
      bhd = next(bhd, mesh);
      CGAL_postcondition(mesh.has_on_seam(bhd) && is_border(bhd, mesh));
    }

    std::cout << "built border & tags" << std::endl;
    std::cout << border_halfedges.size() << " border halfedges" << std::endl;
    CGAL_postcondition(border_halfedges.size() == 2 * mesh.number_of_seam_edges());
  }

  /// Insert all the base faces in the set of unique faces.
  template<typename Face_set,
           typename Face_container>
  void mark_base_mesh_faces(const EmbeddedMesh& emesh,
                            Face_set& inserted_faces,
                            Face_container& faces_to_draw)
  {
    const Seam_mesh& mesh = emesh.mesh;
    const VertexUVMap uvmap = emesh.uvmap;

    BOOST_FOREACH(face_descriptor fd, faces(mesh)) {
      halfedge_descriptor hd = halfedge(fd, mesh);
      vertex_descriptor vda = source(hd, mesh);
      vertex_descriptor vdb = target(hd, mesh);
      vertex_descriptor vdc = target(next(hd, mesh), mesh);

      const Point_2& tpa = get(uvmap, vda);
      const Point_2& tpb = get(uvmap, vdb);
      const Point_2& tpc = get(uvmap, vdc);

//      std::cout << "add " << tpa << " " << tpb << " " << tpc << std::endl;

      Triangle_2 tr(tpa, tpb, tpc);
      faces_to_draw.push_back(tr);

      Point_2 bar = CGAL::barycenter(tpa, 1., tpb, 1., tpc, 1.);

      // ugly hack to get rid of imprecisions @fixme
      // uv values are at the beginnig close to 1 so *1e7 is not too dangerous
      NT barx = 1e7 * bar.x();
      NT bary = 1e7 * bar.y();
      barx = static_cast<int>(barx >= 0 ? barx + 0.5 : barx - 0.5);
      bary = static_cast<int>(bary >= 0 ? bary + 0.5 : bary - 0.5);
      inserted_faces.insert(std::make_pair(barx, bary));
    }
  }

  template<typename BorderHalfedges,
           typename SeamTransformationMap,
           typename Face_container,
           typename Barycenter_set,
           typename Arr_segments>
  void expand_border(const EmbeddedMesh& source_embedded_mesh,
                     const Triangulation& tr,
                     const SeamTransformationMap& seam_transformations,
                     const BorderHalfedges& current_border_halfedges,
                     BorderHalfedges& next_border_halfedges,
                     Face_container& faces_to_draw,
                     Barycenter_set& inserted_faces,
                     Arr_segments& arr_segments)
  {
    const Seam_mesh& mesh = source_embedded_mesh.mesh;
    const VertexUVMap uvmap = source_embedded_mesh.uvmap;
    const TriangleMesh& tm = mesh.mesh();

    typename BorderHalfedges::const_iterator bhit = current_border_halfedges.begin(),
                                             bhend = current_border_halfedges.end();
    for(; bhit!=bhend; ++bhit) {
      // Get the next border edge
      halfedge_w_tr hdt = *bhit;
      TM_halfedge_descriptor hd = hdt.tmhd;

//        std::cout << "at halfedge " << source(hd, tm) << " " << target(hd, tm) << std::endl;
//        halfedge_descriptor shd(hd);
//        std::cout << "that is " << get(uvmap, (source(shd, mesh))) << " || "
//                                << get(uvmap, (target(shd, mesh))) << std::endl;

      // the transformation stocked in the halfedge
      Aff_tr hd_transformation = hdt.transformation;

//        std::cout << "base transformation: " << std::endl;
//        std::cout << hd_transformation.rot_matrix << std::endl;

      // check if [ab] with the current transformation intersects, otherwise stops
      const Point_2& tpa = hd_transformation.transform(get(uvmap, target(hd, mesh)));
      const Point_2& tpb = hd_transformation.transform(get(uvmap, source(hd, mesh)));
      const Segment_2 es_ab(tpa, tpb);
      if(!is_intersecting(tr, es_ab))
        continue;


      // There is an intersection and we must insert the face incident to ohd.
      // Get the opposite halfedge, which is incident to the triangle that we will add
      // (with transformed coordinates)
      TM_halfedge_descriptor ohd = opposite(hd, tm);
//        TM_vertex_descriptor tmvdc = target(next(ohd, tm), tm); // third point
//        std::cout << "third " << tmvdc;
//        std::cout<< " with pos: "  << get(uvmap, target(halfedge_descriptor(next(ohd, tm)), mesh)) << std::endl;

      // if the edge is on a seam, we combine the current transformation
      // with the transformation assigned to the seam halfedge to obtain the transformation
      if(mesh.has_on_seam(hd)) {
        typename SeamTransformationMap::const_iterator it = seam_transformations.find(hd);
        CGAL_assertion(it != seam_transformations.end());
        Aff_tr seam_transf = it->second;
        hd_transformation.combine_rotation(seam_transf.rot_matrix);
      }

      TM_halfedge_descriptor hd_bc = next(ohd, tm);
      vertex_descriptor vdc = target(halfedge_descriptor(hd_bc), mesh);
      const Point_2& tpc = hd_transformation.transform(get(uvmap, vdc));

      const Exact_point_2& etpa = to_exact(tpa);
      const Exact_point_2& etpb = to_exact(tpb);
      const Exact_point_2& etpc = to_exact(tpc);

      // Check if we have already inserted the face incident to ohd in the arrangement
      Point_2 bar = CGAL::barycenter(tpa, 1., tpb, 1., tpc, 1.);

      // ugly hack to get rid of imprecisions @fixme
      NT barx = 1e7 * bar.x();
      NT bary = 1e7 * bar.y();
      barx = static_cast<int>(barx >= 0 ? barx + 0.5 : barx - 0.5);
      bary = static_cast<int>(bary >= 0 ? bary + 0.5 : bary - 0.5);

      std::pair<typename Barycenter_set::iterator, bool> is_insert_successful =
                              inserted_faces.insert(std::make_pair(barx, bary));
      if(!is_insert_successful.second) // triangle has already been considered
        continue;

      std::cout.precision(20);
      std::cout << "new triangle with bar: " << barx << " || " << bary << std::endl;

      Triangle_2 tr(tpa, tpb, tpc);
      faces_to_draw.push_back(tr);

//        std::cout << "tps: " << tpa << " " << tpb << " " << tpc << std::endl;

      // Add the new segments
      Exact_segment_2 es_ac(etpa, etpc);
      Exact_segment_2 es_bc(etpb, etpc);
      // There is no need to insert ab because we have already inserted it
      // at previous steps. Leaving it here for clarity.
//      arr_segments.push_back(Exact_segment_2(etpa, etpb))
      arr_segments.push_back(es_ac);
      arr_segments.push_back(es_bc);

      // Insert the (other) edges of the face incident to ohd
      // there are two "new" border edges: ac and bc
      TM_halfedge_descriptor hd_ac = prev(ohd, tm);

      halfedge_w_tr hdt_ac(hd_ac, hd_transformation);
      next_border_halfedges.push_back(hdt_ac);

      halfedge_w_tr hdt_bc(hd_bc, hd_transformation);
      next_border_halfedges.push_back(hdt_bc);
    }
  }

  /// Check if a given segment intersects the target mesh (actually, a triangulation
  /// of the -- slightly grown -- convex hull of the target mesh)
  bool is_intersecting(const Triangulation& tr,
                       const Segment_2& segment)
  {
    typename Kernel::Construct_triangle_2 triangle_functor;
    typename Kernel::Do_intersect_2 do_intersect_2_functor;

    // Brute force, but the triangulation _should_ have very few faces
    typedef typename Triangulation::Finite_faces_iterator    Finite_faces_iterator;
    Finite_faces_iterator fit = tr.finite_faces_begin(), end = tr.finite_faces_end();
    for(; fit!=end; fit++) {
      Triangle_2 tr = triangle_functor(fit->vertex(0)->point(),
                                       fit->vertex(1)->point(),
                                       fit->vertex(2)->point());
      if(do_intersect_2_functor(segment, tr))
        return true;
    }

    return false;
  }

  /// Functor used to offer dereferencing vertex_iterator to UV values
  struct Vd_to_uv
  {
    typedef vertex_iterator                           argument_type;
    typedef Point_2                                   result_type;

    const VertexUVMap uvmap;

    result_type operator()(vertex_iterator vit) const
    {
      vertex_descriptor vd = *vit;
      return get(uvmap, vd);
    }

    result_type operator()(vertex_descriptor vd) const
    {
      return get(uvmap, vd);
    }

    Vd_to_uv(const VertexUVMap uvmap) : uvmap(uvmap) { }
  };

  /// Compute the convex hull of a mesh, and triangulate the convex hull.
  void triangulate_convex_hull(const EmbeddedMesh& emesh,
                               Triangulation& tr) const
  {
    const Seam_mesh& mesh = emesh.mesh;
    const VertexUVMap uvmap = emesh.uvmap;

    // compute the convex hull of the mesh
    std::vector<Point_2> convex_hull_pts;
    Vd_to_uv vd_to_uv(uvmap);
    CGAL::convex_hull_2(boost::make_transform_iterator(mesh.vertices_begin(), vd_to_uv),
                        boost::make_transform_iterator(mesh.vertices_end(), vd_to_uv),
                        std::back_inserter(convex_hull_pts));

    std::cout << convex_hull_pts.size() << " pts on the convex hull" << std::endl;

    /// Trick: slightly grow the convex hull so that the source mesh will completely
    /// englobe the target mesh

    // compute the center point (ish)
    NT weight = 1. / convex_hull_pts.size();
    NT bx = 0., by = 0.;
    for(std::size_t i=0; i<convex_hull_pts.size(); ++i) {
      bx += convex_hull_pts[i].x();
      by += convex_hull_pts[i].y();
    }

    bx *= weight;
    by *= weight;

    Point_2 bp(bx, by);

    // move points on the convex hull away from that center point
    for(std::size_t i=0; i<convex_hull_pts.size(); ++i) {
      Point_2& pci = convex_hull_pts[i];
      pci = pci + 0.01 * Vector_2(bp, pci);
    }

    // triangulate the (grown) convex hull
    tr.insert(convex_hull_pts.begin(), convex_hull_pts.end());
    CGAL_postcondition(tr.number_of_vertices() == convex_hull_pts.size());
    std::cout << tr.number_of_faces() << " faces in the triangulation" << std::endl;
  }

  /// Grow the source mesh until it covers the target mesh.
  Arrangement grow_source_embedding(const EmbeddedMesh& source_embedded_mesh,
                                    const EmbeddedMesh& target_embedded_mesh)
  {
    const Seam_mesh& mesh = source_embedded_mesh.mesh;
    const VertexUVMap uvmap = source_embedded_mesh.uvmap;

    // Seam tags (on the source mesh) indicate the transformation that must
    // be applied to the UV coordinates to obtain the coordinates of
    // the adjacent triangle in the tiled space when going over the seam edge
    typedef boost::unordered_map<TM_halfedge_descriptor, Aff_tr> SeamTransformationMap;
    SeamTransformationMap seam_transformations;

    // Border_halfedes is the border of the expending border of the source mesh.
    // We grow it until the border does not intersect the target mesh anymore.
    // note: ideally, it'd be a set, but I don't have a good hash function for it
    //       and in practice it's not slower to use a vector.
    typedef std::vector<halfedge_w_tr> BorderHalfedges;
    BorderHalfedges current_border_halfedges, next_border_halfedges;
    current_border_halfedges.reserve(2 * mesh.number_of_seam_edges());

    compute_initial_border(source_embedded_mesh, seam_transformations, current_border_halfedges);

    CGAL_postcondition(seam_transformations.size() == 2 * mesh.number_of_seam_edges());
    CGAL_postcondition(current_border_halfedges.size() == 2 * mesh.number_of_seam_edges());

    // Segments that will form the arrangement obtained by growing the source mesh
    std::vector<Arr_segment_2> arr_segments;
    arr_segments.reserve(num_edges(mesh));

    // Fill the arrangement with the edges of the source mesh
    BOOST_FOREACH(edge_descriptor ed, edges(mesh)) {
      vertex_descriptor vds = source(ed, mesh);
      vertex_descriptor vdt = target(ed, mesh);
      arr_segments.push_back(Exact_segment_2(to_exact(get(uvmap, vds)),
                                             to_exact(get(uvmap, vdt))));
    }

// -------------- INFO & OUTPUT
    std::cout << num_edges(mesh.mesh()) << " edges in the base mesh" << std::endl;
    std::cout << mesh.number_of_seam_edges() << " seams" << std::endl;
    std::cout << num_edges(mesh) << " edges in the source seam mesh" << std::endl;
    std::cout << arr_segments.size() << " segments in the base arrangement" << std::endl;
    CGAL_postcondition(arr_segments.size() == num_edges(mesh));

    Arrangement source_arrangement;
    CGAL::insert(source_arrangement, arr_segments.begin(), arr_segments.end());

    // Write the arrangement to a file.
    std::ofstream arr_out("arr_source.off");
    internal::output_arrangement_to_off(source_arrangement, arr_out);
// --------------------

    // Keep a set of faces. Can't sort by IDs, so using a set of barycenters instead...
    // Problem is that rotations introduce imprecisions so we can't just get a point
    // Some (super) ugly hack is to only consider part of the mantissa
    typedef std::set<std::pair<int, int> >                       Barycenter_set;
    Barycenter_set inserted_faces;
    std::vector<Triangle_2> faces_to_draw; // just for output
    mark_base_mesh_faces(source_embedded_mesh, inserted_faces, faces_to_draw);

    // The target mesh is a dense mesh, which makes it expensive to check for intersections.
    // Instead, build a triangulation of the convex hull, which is "close enough"
    // since the parameterisations give somewhat convex results
    Triangulation tr;
    triangulate_convex_hull(target_embedded_mesh, tr);

    // Expand the border until it doesn't intersect the target mesh
    std::cout << "initial border of length: " << current_border_halfedges.size() << std::endl;
    int counter = 0;
    while(true) {
      // For all halfedges in the current border (non-continuous), add the incident
      // triangle if the shared edge intersects the target mesh
      expand_border(source_embedded_mesh, tr, seam_transformations,
                    current_border_halfedges, next_border_halfedges,
                    faces_to_draw, inserted_faces, arr_segments);
      ++counter;

      std::cout << arr_segments.size() << " segments in the arrangement" << std::endl;
      std::cout << next_border_halfedges.size() << " in the next loop" << std::endl;

      if(next_border_halfedges.empty())
        break;

      current_border_halfedges = next_border_halfedges;
      next_border_halfedges.clear();
    }

    std::cout << "Finished growing!" << std::endl;

    // visualize the grown mesh
    internal::output_triangle_set<std::vector<Triangle_2> >(faces_to_draw, "arr_grown.off");

    Arrangement grown_arrangement;
    std::cout << "building grown arrangement with " << arr_segments.size() << " segments" << std::endl;
    CGAL::insert(grown_arrangement, arr_segments.begin(), arr_segments.end());

    // Write the arrangement to a file.
    std::ofstream arr_grown_out("arr_grown_old.off");
    internal::output_arrangement_to_off(grown_arrangement, arr_grown_out);

    return grown_arrangement;
  }

  /// Build an arrangement from a two-dimensional mesh.
  Arrangement get_arrangement_from_embedded_mesh(const EmbeddedMesh& emesh)
  {
    const Seam_mesh& mesh = emesh.mesh;
    const VertexUVMap uvmap = emesh.uvmap;

    std::vector<Exact_segment_2> arr_segments;
    arr_segments.reserve(num_edges(mesh));

    BOOST_FOREACH(edge_descriptor ed, edges(mesh)) {
      vertex_descriptor vds = source(ed, mesh);
      vertex_descriptor vdt = target(ed, mesh);
      arr_segments.push_back(Exact_segment_2(to_exact(get(uvmap, vds)),
                                             to_exact(get(uvmap, vdt))));
    }

    // Build the arrangement
    Arrangement arrangement;
    CGAL::insert(arrangement, arr_segments.begin(), arr_segments.end());

    // Write the arrangement to a file.
    std::ofstream arr_out("arr_target.off");
    internal::output_arrangement_to_off(arrangement, arr_out);

    return arrangement;
  }

  /// Overlay the source and target arrangements.
  void overlay_arrangements(const Arrangement& source_arrangement,
                            const Arrangement& target_arrangement)
  {
    std::cout << "The source arrangement size:" << std::endl
              << "   V = " << source_arrangement.number_of_vertices()
              << ",  E = " << source_arrangement.number_of_edges()
              << ",  F = " << source_arrangement.number_of_faces() << std::endl;

    std::cout << "The target arrangement size:" << std::endl
              << "   V = " << target_arrangement.number_of_vertices()
              << ",  E = " << target_arrangement.number_of_edges()
              << ",  F = " << target_arrangement.number_of_faces() << std::endl;

    Arrangement overlay_arrangement;
    Overlay_traits overlay_traits;
    CGAL::overlay(source_arrangement, target_arrangement,
                  overlay_arrangement, overlay_traits);

    // Write the arrangement to a file
    std::ofstream arr_out("arr_overlay.off");
    internal::output_arrangement_to_off(overlay_arrangement, arr_out);

    std::cout << "The overlaid arrangement size:" << std::endl
              << "   V = " << overlay_arrangement.number_of_vertices()
              << ",  E = " << overlay_arrangement.number_of_edges()
              << ",  F = " << overlay_arrangement.number_of_faces() << std::endl;
  }

public:
  /// Compute the map between the source and the target mesh.
  void compute_map_from_sphere_embeddings(const EmbeddedMesh& source_embedded_mesh,
                                          const EmbeddedMesh& target_embedded_mesh)
  {
    // Grow the source mesh until it covers the target mesh and convert it to arrangement
    const Arrangement& source_arrangement =
              grow_source_embedding(source_embedded_mesh, target_embedded_mesh);

    const Arrangement& target_arrangement =
                       get_arrangement_from_embedded_mesh(target_embedded_mesh);

    // Overlay the two arrangements
    overlay_arrangements(source_arrangement, target_arrangement);

    // Compute the one-to-one mapping
    // @todo
  }

  /// Constructor.
  Orbifold_sphere_mapper() : to_exact() { }
};

} // namespace Surface_mesh_parameterization

} // namespace CGAL

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_ORBITAL_TUTTE_SPHERE_MAPPING_H
