// Copyright (c) 2019 INRIA Sophia-Antipolis (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Kaimo Hu

#ifndef SRC_INTERNAL_MINANGLE_REMESHING_MESH_PROPERTIES_H_
#define SRC_INTERNAL_MINANGLE_REMESHING_MESH_PROPERTIES_H_

// C/C++
#include <limits>
#include <functional>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <utility>
#include <algorithm>
#include <fstream>
// CGAL
#include <CGAL/Timer.h>
#include <CGAL/Triangulation_hierarchy_2.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
// boost
#include <boost/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>
// local
#include "Bvd.h"
#include "Random.h"

// namespace definition
namespace PMP = CGAL::Polygon_mesh_processing;
namespace NP = CGAL::parameters;

// const data
const double DOUBLE_MAX = 1000000.0;
const double DOUBLE_MIN = -1000000.0;
const int MAX_VALUE = 10000;
const double MIN_VALUE = 0.0001;  // specified for numerical stability
const double SQUARED_MIN_VALUE = 0.00000001;

// numerical types
enum SampleNumberStrategy {
  // #samples per face is roughtly fixed (with respect to the sample strategy)
  k_fixed = 0,
  // #samples per face is variable with respect to size_of_faces()c
  k_variable
};

enum SampleStrategy {
  k_uniform = 0,  // #samples per face is proportional to its area
  k_adaptive      // #samples per face is roughly the same
};

enum OptimizeType {
  k_none = 0,
  k_input_to_remesh,
  k_remesh_to_input,
  k_both
};

enum OptimizeStrategy {
  k_approximation = 0,
  k_Interpolation
};

enum EdgeFlipStrategy {
  k_improve_valence = 0,
  k_improve_angle
};

enum RelocateStrategy {
  k_barycenter = 0,
  k_cvt_barycenter
};

enum VertexType {
  k_feature_vertex = 0,
  k_crease_vertex,
  k_smooth_vertex
};

struct NamedParameters {
  // general parameters
  double max_error_threshold;
  double min_angle_threshold;
  int max_mesh_complexity;
  double smooth_angle_delta;
  bool apply_edge_flip;
  EdgeFlipStrategy edge_flip_strategy;
  bool flip_after_split_and_collapse;
  bool relocate_after_local_operations;
  RelocateStrategy relocate_strategy;
  bool keep_vertex_in_one_ring;
  bool use_local_aabb_tree;
  int collapsed_list_size;
  bool decrease_max_errors;
  bool verbose_progress;
  bool apply_initial_mesh_simplification;
  bool apply_final_vertex_relocation;
  // sample parameters
  int samples_per_face_in;
  int samples_per_face_out;
  int max_samples_per_area;
  int min_samples_per_triangle;
  int bvd_iteration_count;
  SampleNumberStrategy sample_number_strategy;
  SampleStrategy sample_strategy;
  bool use_stratified_sampling;
  // feature function parameters
  double sum_theta;
  double sum_delta;
  double dihedral_theta;
  double dihedral_delta;
  double feature_difference_delta;
  double feature_control_delta;
  bool inherit_element_types;
  bool use_feature_intensity_weights;
  // vertex relocate_parameters
  int vertex_optimize_count;
  double vertex_optimize_ratio;
  int stencil_ring_size;
  OptimizeStrategy optimize_strategy;
  OptimizeType face_optimize_type;
  OptimizeType edge_optimize_type;
  OptimizeType vertex_optimize_type;
  bool optimize_after_local_operations;
};

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

template<typename Kernel>
class Mesh_properties {
 public:
  // Type definitions
  typedef CGAL::Bbox_3 Bbox;
  // 2D BVD
  typedef CGAL::Triangulation_vertex_base_2<Kernel> Vbb;
  typedef CGAL::Triangulation_hierarchy_vertex_base_2<Vbb> Vb;
  typedef CGAL::Triangulation_face_base_2<Kernel> Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
  // 3D BVD
  typedef CBvd<Kernel, Tds> Bvd;
  typedef typename Bvd::FT FT;
  typedef typename Bvd::Vector_3 Vector;
  typedef typename Bvd::Vector_3 Normal;
  typedef typename Bvd::Point_3 Point;
  typedef typename Bvd::Line_3 Line;
  typedef typename Bvd::Segment_3 Segment;
  typedef typename Bvd::Triangle_3 Triangle;
  typedef typename Bvd::Plane_3 Plane;
  // Point list
  typedef typename Bvd::Point_list Point_list;
  typedef typename Bvd::Point_iter Point_iter;
  typedef typename Bvd::Point_const_iter Point_const_iter;
  // Color list
  typedef typename Bvd::Color_list Color_list;
  typedef typename Bvd::Color_iter Color_iter;
  typedef typename Bvd::Color_const_iter Color_const_iter;
  // Local link related
  typedef std::pair<Point, Point> Point_pair;        // for out links
  typedef std::pair<FT, Point_pair> Link;
  typedef std::list<Link> Link_list;
  typedef typename Link_list::iterator Link_list_iter;
  typedef typename Link_list::const_iterator Link_list_const_iter;
  typedef std::list<Link_list_iter> Link_iter_list;  // for in links
  typedef typename Link_iter_list::iterator Link_iter_list_iter;
  typedef typename Link_iter_list::const_iterator Link_iter_list_const_iter;
  typedef std::list<Link*> Link_pointer_list;
  typedef typename std::list<Link*>::iterator Link_pointer_iter;
  typedef typename std::list<Link*>::const_iterator Link_pointer_const_iter;
  // Surface_mesh related
  typedef CGAL::Surface_mesh<Point> Mesh;
  typedef typename boost::graph_traits<Mesh>::halfedge_descriptor
                                              halfedge_descriptor;
  typedef typename boost::graph_traits<Mesh>::edge_descriptor edge_descriptor;
  typedef typename boost::graph_traits<Mesh>::vertex_descriptor
                                              vertex_descriptor;
  typedef typename boost::graph_traits<Mesh>::face_descriptor face_descriptor;
  typedef CGAL::Halfedge_around_target_circulator<Mesh>
                         Halfedge_around_target_circulator;
  typedef CGAL::Halfedge_around_face_circulator<Mesh>
                         Halfedge_around_face_circulator;
  typedef CGAL::Face_around_target_circulator<Mesh>
                         Face_around_target_circulator;
  // Element list
  typedef std::list<halfedge_descriptor> Halfedge_list;
  typedef typename std::list<halfedge_descriptor>::iterator Halfedge_iter;
  typedef typename std::list<halfedge_descriptor>::const_iterator
                                                   Halfedge_const_iter;
  typedef std::list<edge_descriptor> Edge_list;
  typedef typename std::list<edge_descriptor>::iterator Edge_iter;
  typedef typename std::list<edge_descriptor>::const_iterator Edge_const_iter;
  typedef std::list<vertex_descriptor> Vertex_list;
  typedef typename std::list<vertex_descriptor>::iterator Vertex_iter;
  typedef typename std::list<vertex_descriptor>::const_iterator
                                                 Vertex_const_iter;
  typedef std::list<face_descriptor> Face_list;
  typedef typename std::list<face_descriptor>::iterator Face_iter;
  typedef typename std::list<face_descriptor>::const_iterator Face_const_iter;
  // Property related
  typedef typename boost::property_map<Mesh,
      CGAL::dynamic_face_property_t<int>>::type Face_tags;         // faces
  typedef typename boost::property_map<Mesh,
      CGAL::dynamic_face_property_t<Normal>>::type Face_normals;
  typedef typename boost::property_map<Mesh,
      CGAL::dynamic_face_property_t<FT>>::type Face_max_errors;
  typedef typename boost::property_map<Mesh,
      CGAL::dynamic_face_property_t<Link_list>>::type Face_link_list;
  typedef typename boost::property_map<Mesh,
      CGAL::dynamic_face_property_t<Link_iter_list>>::type
      Face_link_iter_list;
  typedef typename boost::property_map<Mesh,
      CGAL::dynamic_face_property_t<Link_pointer_list>>::type
      Face_link_pointer_list;
  typedef typename boost::property_map<Mesh,
      CGAL::dynamic_halfedge_property_t<int>>::type Halfedge_tags;  // edges
  typedef typename boost::property_map<Mesh,
      CGAL::dynamic_halfedge_property_t<FT>>::type Halfedge_normal_dihedrals;
  typedef typename boost::property_map<Mesh,
      CGAL::dynamic_halfedge_property_t<bool>>::type Halfedge_are_creases;
  typedef typename boost::property_map<Mesh,
      CGAL::dynamic_halfedge_property_t<Link_list>>::type Halfedge_link_list;
  typedef typename boost::property_map<Mesh,
      CGAL::dynamic_vertex_property_t<int>>::type Vertex_tags;     // vertices
  typedef typename boost::property_map<Mesh,
      CGAL::dynamic_vertex_property_t<FT>>::type Vertex_max_dihedral;
  typedef typename boost::property_map<Mesh,
      CGAL::dynamic_vertex_property_t<FT>>::type Vertex_gaussian_curvature;
  typedef typename boost::property_map<Mesh,
      CGAL::dynamic_vertex_property_t<Link>>::type Vertex_link;
  // AABB tree
  typedef CGAL::AABB_face_graph_triangle_primitive<Mesh>
      Face_primitive;
  typedef CGAL::AABB_traits<Kernel, Face_primitive> Face_traits;
  typedef AABB_tree<Face_traits> Face_tree;
  typedef typename Face_tree::Point_and_primitive_id Point_and_primitive_id;
  // Dynamic priority queues
  typedef boost::bimap<boost::bimaps::set_of<halfedge_descriptor>,
      boost::bimaps::multiset_of<FT, std::greater<FT>>> DPQueue_halfedge_long;
  typedef typename DPQueue_halfedge_long::value_type    Halfedge_long;
  typedef typename boost::bimap<boost::bimaps::set_of<halfedge_descriptor>,
      boost::bimaps::multiset_of<FT, std::less<FT>>>    DPQueue_halfedge_short;
  typedef typename DPQueue_halfedge_short::value_type   Halfedge_short;
  typedef typename boost::bimap<boost::bimaps::set_of<vertex_descriptor>,
      boost::bimaps::multiset_of<FT, std::greater<FT>>> DPQueue_vertex_long;
  typedef typename DPQueue_vertex_long::value_type      Vertex_long;
  typedef typename boost::bimap<boost::bimaps::set_of<vertex_descriptor>,
      boost::bimaps::multiset_of<FT, std::less<FT>>>    DPQueue_vertex_short;
  typedef typename DPQueue_vertex_short::value_type     Vertex_short;
  typedef typename boost::bimap<boost::bimaps::set_of<face_descriptor>,
      boost::bimaps::multiset_of<FT, std::greater<FT>>> DPQueue_face_long;
  typedef typename DPQueue_face_long::value_type        Face_long;
  typedef typename boost::bimap<boost::bimaps::set_of<face_descriptor>,
      boost::bimaps::multiset_of<FT, std::less<FT>>>    DPQueue_face_short;
  typedef typename DPQueue_face_short::value_type       Face_short;

  // Struct definitions
  struct Point_Comp {
    bool operator () (const Point &a, const Point &b) const {
      Vector vec = a - b;
      if (vec * vec < SQUARED_MIN_VALUE) {
        return false;
      } else {
        return a < b;
      }
    }
  };

 public:
  // 1) life cycles
  explicit Mesh_properties(Mesh *mesh)
    : mesh_(*mesh) {
    // face related properties
    face_tags_ = get(CGAL::dynamic_face_property_t<int>(), mesh_);
    face_normals_ = get(CGAL::dynamic_face_property_t<Normal>(), mesh_);
    face_max_squared_errors_ = get(CGAL::dynamic_face_property_t<FT>(), mesh_);
    face_out_links_ = get(CGAL::dynamic_face_property_t<Link_list>(), mesh_);
    face_in_links_ =
        get(CGAL::dynamic_face_property_t<Link_iter_list>(), mesh_);
    edge_in_links_ =
        get(CGAL::dynamic_face_property_t<Link_iter_list>(), mesh_);
    vertex_in_links_ =
        get(CGAL::dynamic_face_property_t<Link_pointer_list>(), mesh_);
    for (typename Mesh::Face_range::const_iterator fi = mesh_.faces().begin();
        fi != mesh_.faces().end(); ++fi) {
      reset_face_properties(*fi, get_null_face());
    }
    // halfedge related properties
    halfedge_tags_ = get(CGAL::dynamic_halfedge_property_t<int>(), mesh_);
    halfedge_normal_dihedrals_ =
        get(CGAL::dynamic_halfedge_property_t<FT>(), mesh_);
    halfedge_are_creases_ =
        get(CGAL::dynamic_halfedge_property_t<bool>(), mesh_);
    halfedge_out_links_ =
        get(CGAL::dynamic_halfedge_property_t<Link_list>(), mesh_);
    for (typename Mesh::Halfedge_range::const_iterator hi = mesh_.halfedges().begin();
        hi != mesh_.halfedges().end(); ++hi) {
      reset_halfedge_properties(*hi, get_null_halfedge());
    }
    // vertex related properties
    vertex_tags_ = get(CGAL::dynamic_vertex_property_t<int>(), mesh_);
    vertex_max_dihedrals_ = get(CGAL::dynamic_vertex_property_t<FT>(), mesh_);
    vertex_gaussian_curvatures_ =
        get(CGAL::dynamic_vertex_property_t<FT>(), mesh_);
    vertex_out_link_ = get(CGAL::dynamic_vertex_property_t<Link>(), mesh_);
    for (typename Mesh::Vertex_range::const_iterator vi = mesh_.vertices().begin();
        vi != mesh_.vertices().end(); ++vi) {
      reset_vertex_properties(*vi, get_null_vertex());
    }
  }

  virtual ~Mesh_properties() {}

  // 2) elements access
  inline const Mesh& get_mesh() const { return mesh_; }
  inline face_descriptor get_face(halfedge_descriptor hd) const
      { return mesh_.face(hd); }
  inline halfedge_descriptor get_opposite(halfedge_descriptor hd) const
      { return mesh_.opposite(hd); }
  inline edge_descriptor get_edge(halfedge_descriptor hd) const
    { return mesh_.edge(hd); }
  inline Point& get_point(vertex_descriptor vd) { return mesh_.point(vd); }
  inline const Point& get_point(vertex_descriptor vd) const
      { return mesh_.point(vd); }
  inline face_descriptor get_null_face() const { return mesh_.null_face(); }
  inline halfedge_descriptor get_null_halfedge() const
      { return mesh_.null_halfedge(); }
  inline edge_descriptor get_null_edge() const { return mesh_.null_edge(); }
  inline vertex_descriptor get_null_vertex() const
      { return mesh_.null_vertex(); }
  vertex_descriptor get_opposite_vertex(halfedge_descriptor hd) {
    assert(!is_border(hd));
    if (is_border(hd)) {
      std::cout << "Error definition of opposite vertex for border halfedge."
        << std::endl;
      return get_null_vertex();
    } else {
      return mesh_.target(mesh_.next(hd));
    }
  }
  vertex_descriptor get_opposite_vertex(halfedge_descriptor hd) const {
    if (is_border(hd)) {
      std::cout << "Error definition of opposite vertex for border halfedge."
        << std::endl;
      return get_null_vertex();
    } else {
      return mesh_.target(mesh_.next(hd));
    }
  }

  // 3) properties access
  int& get_face_tag(face_descriptor fd) { return get(face_tags_, fd); }
  const int& get_face_tag(face_descriptor fd) const
      { return get(face_tags_, fd); }
  void set_face_tag(face_descriptor fd, int value)
      { put(face_tags_, fd, value); }

  Normal& get_face_normal(face_descriptor fd)
      { return get(face_normals_, fd); }
  const Normal& get_face_normal(face_descriptor fd) const
      { return get(face_normals_, fd); }
  void set_face_normal(face_descriptor fd, Vector value)
      { put(face_normals_, fd, value); }

  FT& get_face_max_squared_error(face_descriptor fd)
      { return get(face_max_squared_errors_, fd); }
  const FT& get_face_max_squared_error(face_descriptor fd) const
      { return get(face_max_squared_errors_, fd); }
  void set_face_max_squared_error(face_descriptor fd, FT value)
      { put(face_max_squared_errors_, fd, value); }

  Link_list& get_face_out_links(face_descriptor fd)
      { return get(face_out_links_, fd); }
  const Link_list& get_face_out_links(face_descriptor fd) const
      { return get(face_out_links_, fd); }
  void set_face_out_links(face_descriptor fd, const Link_list &value)
      { put(face_out_links_, fd, value); }

  Link_iter_list& get_face_in_links(face_descriptor fd)
      { return get(face_in_links_, fd); }
  const Link_iter_list& get_face_in_links(face_descriptor fd) const
      { return get(face_in_links_, fd); }
  void set_face_in_links(face_descriptor fd, const Link_iter_list &value)
      { put(face_in_links_, fd, value); }

  Link_iter_list& get_edge_in_links(face_descriptor fd)
      { return get(edge_in_links_, fd); }
  const Link_iter_list& get_edge_in_links(face_descriptor fd) const
      { return get(edge_in_links_, fd); }
  void set_edge_in_links(face_descriptor fd, const Link_iter_list &value)
      { put(edge_in_links_, fd, value); }

  Link_pointer_list& get_vertex_in_links(face_descriptor fd)
      { return get(vertex_in_links_, fd); }
  const Link_pointer_list& get_vertex_in_links(face_descriptor fd) const
      { return get(vertex_in_links_, fd); }
  void set_vertex_in_links(face_descriptor fd, const Link_pointer_list &value)
      { put(vertex_in_links_, fd, value); }

  int& get_halfedge_tag(halfedge_descriptor hd)
      { return get(halfedge_tags_, hd); }
  const int& get_halfedge_tag(halfedge_descriptor hd) const
      { return get(halfedge_tags_, hd); }
  void set_halfedge_tag(halfedge_descriptor hd, int value)
      { put(halfedge_tags_, hd, value); }

  FT& get_halfedge_normal_dihedral(halfedge_descriptor hd)
      { return get(halfedge_normal_dihedrals_, hd); }
  const FT& get_halfedge_normal_dihedral(halfedge_descriptor hd) const
      { return get(halfedge_normal_dihedrals_, hd); }
  void set_halfedge_normal_dihedral(halfedge_descriptor hd, FT value)
      { put(halfedge_normal_dihedrals_, hd, value); }

  bool get_halfedge_is_crease(halfedge_descriptor hd)
      { return get(halfedge_are_creases_, hd); }
  bool get_halfedge_is_crease(halfedge_descriptor hd) const
      { return get(halfedge_are_creases_, hd); }
  void set_halfedge_is_crease(halfedge_descriptor hd, bool value)
      { put(halfedge_are_creases_, hd, value); }

  Link_list& get_halfedge_out_links(halfedge_descriptor hd)
      { return get(halfedge_out_links_, hd); }
  const Link_list& get_halfedge_out_links(halfedge_descriptor hd) const
      { return get(halfedge_out_links_, hd); }
  void set_halfege_out_links(halfedge_descriptor hd, const Link_list &value)
      { put(halfedge_out_links_, hd, value); }

  int& get_vertex_tag(vertex_descriptor vd) { return get(vertex_tags_, vd); }
  const int& get_vertex_tag(vertex_descriptor vd) const
      { return get(vertex_tags_, vd); }
  void set_vertex_tag(vertex_descriptor vd, int value)
      { put(vertex_tags_, vd, value); }

  FT& get_vertex_max_dihedral(vertex_descriptor vd)
      { return get(vertex_max_dihedrals_, vd); }
  const FT& get_vertex_max_dihedral(vertex_descriptor vd) const
      { return get(vertex_max_dihedrals_, vd); }
  void set_vertex_max_dihedral(vertex_descriptor vd, FT value)
      { put(vertex_max_dihedrals_, vd, value); }

  FT& get_vertex_gaussian_curvature(vertex_descriptor vd)
      { return get(vertex_gaussian_curvatures_, vd); }
  const FT& get_vertex_gaussian_curvature(vertex_descriptor vd) const
      { return get(vertex_gaussian_curvatures_, vd); }
  void set_vertex_gaussian_curvature(vertex_descriptor vd, FT value)
      { put(vertex_gaussian_curvatures_, vd, value); }

  Link& get_vertex_out_link(vertex_descriptor vd)
      { return get(vertex_out_link_, vd); }
  const Link& get_vertex_out_link(vertex_descriptor vd) const
      { return get(vertex_out_link_, vd); }
  void set_vertex_out_link(vertex_descriptor vd, const Link &value)
      { put(vertex_out_link_, vd, value); }

  // 4) queues (for isotropic remeshing)
  void fill_collapse_candidate_edges(FT max_error_threshold_value,
    DPQueue_halfedge_long *large_error_queue,
    DPQueue_halfedge_short *collapse_candidate_queue,
    const NamedParameters &np) const {
    // step 1: fill the large error queue if necessary
    if (np.decrease_max_errors) {
      FT max_se_threshold = std::pow(max_error_threshold_value, 2);
      halfedge_descriptor max_error_halfedge;
      for (typename Mesh::Face_range::const_iterator fi = mesh_.faces().begin();
        fi != mesh_.faces().end(); ++fi) {
        FT max_se = get_face_max_squared_error(*fi);
        if (max_se >= max_se_threshold) {
          max_error_halfedge = get_longest_halfedge(*fi);
          large_error_queue->insert(Halfedge_long(max_error_halfedge, max_se));
        }
      }
    }
    // step 2: fill the collapse_candidate queue
    for (typename Mesh::Halfedge_range::const_iterator hi = mesh_.halfedges().begin();
      hi != mesh_.halfedges().end(); ++hi) {
      if (!is_border(*hi)) {
        FT len = length(*hi);
        FT radian = calculate_opposite_radian(*hi);
        collapse_candidate_queue->insert(Halfedge_short(*hi, len * radian));
      }
    }
  }

  void fill_small_radian_edges(FT max_error_threshold_value,
    DPQueue_halfedge_long *large_error_queue,
    DPQueue_halfedge_short *small_radian_queue,
    const NamedParameters &np) const {
    // step 1: fill the large error queue if necessary
    if (np.decrease_max_errors) {
      FT max_se_threshold = std::pow(max_error_threshold_value, 2);
      halfedge_descriptor max_error_halfedge;
      for (typename Mesh::Face_range::const_iterator fi = mesh_.faces().begin();
        fi != mesh_.faces().end(); ++fi) {
        FT max_se = get_face_max_squared_error(*fi);
        if (max_se >= max_se_threshold) {
          max_error_halfedge = get_longest_halfedge(*fi);
          large_error_queue->insert(Halfedge_long(max_error_halfedge, max_se));
        }
      }
    }
    // step 2: fill the small radian queue
    FT min_radian_threshold = to_radian(np.min_angle_threshold);
    for (typename Mesh::Halfedge_range::const_iterator hi = mesh_.halfedges().begin();
      hi != mesh_.halfedges().end(); ++hi) {
      if (!is_border(*hi)) {
        FT radian = calculate_opposite_radian(*hi);
        if (radian < min_radian_threshold) {
          small_radian_queue->insert(Halfedge_short(*hi, radian));
        }
      }
    }
  }

  void fill_relocate_candidate_vertices(
    DPQueue_vertex_short *relocate_candidate_queue) const {
    for (typename Mesh::Vertex_range::const_iterator vi = mesh_.vertices().begin();
      vi != mesh_.vertices().end(); ++vi) {
      FT min_radian = calculate_minimal_radian_incident_to_vertex(*vi);
      relocate_candidate_queue->insert(Vertex_short(*vi, min_radian));
    }
  }

  void update_relocate_candidate_vertices(vertex_descriptor vd,
    DPQueue_vertex_short *relocate_candidate_queue) const {
    std::set<vertex_descriptor> incident_vertices;
    collect_incident_vertices(vd, &incident_vertices);
    incident_vertices.insert(vd);
    for (auto it = incident_vertices.begin();
      it != incident_vertices.end(); ++it) {
      vertex_descriptor v = *it;
      FT min_radian = calculate_minimal_radian_incident_to_vertex(v);
      relocate_candidate_queue->insert(Vertex_short(v, min_radian));
    }
  }

  // 5) sample weights and capacities
  FT get_max_sample_weight() const {
    FT max_weight = 0.0;
    max_weight = CGAL::max(max_weight, get_max_face_sample_weight());
    max_weight = CGAL::max(max_weight, get_max_edge_sample_weight());
    max_weight = CGAL::max(max_weight, get_max_vertex_sample_weight());
    return max_weight;
  }

  FT get_max_face_sample_weight() const {
    FT max_weight = 0.0;
    for (typename Mesh::Face_range::const_iterator fi = mesh_.faces().begin();
      fi != mesh_.faces().end(); ++fi) {
      const Link_list &face_out_links = get_face_out_links(*fi);
      for (auto it = face_out_links.begin();
        it != face_out_links.end(); ++it) {
        const Link &link = *it;
        max_weight = CGAL::max(max_weight, link.first);
      }
    }
    return max_weight;
  }

  FT get_max_edge_sample_weight() const {
    FT max_weight = 0.0;
    for (typename Mesh::Edge_range::const_iterator ei = mesh_.edges().begin();
      ei != mesh_.edges().end(); ++ei) {
      halfedge_descriptor hd = mesh_.halfedge(*ei);
      if (get_halfedge_normal_dihedral(hd) == -1.0) {
        hd = get_opposite(hd);
      }
      const Link_list &edge_out_links = get_halfedge_out_links(hd);
      for (auto it = edge_out_links.begin();
        it != edge_out_links.end(); ++it) {
        const Link &link = *it;
        max_weight = CGAL::max(max_weight, link.first);
      }
    }
    return max_weight;
  }

  FT get_max_vertex_sample_weight() const {
    FT max_weight = 0.0;
    for (typename Mesh::Vertex_range::const_iterator vi = mesh_.vertices().begin();
      vi != mesh_.vertices().end(); ++vi) {
      const Link &link = get_vertex_out_link(*vi);
      max_weight = CGAL::max(max_weight, link.first);
    }
    return max_weight;
  }

  FT get_max_sample_capacity() const {
    FT max_capacity = 0.0;
    max_capacity = CGAL::max(max_capacity, get_max_face_sample_capacity());
    max_capacity = CGAL::max(max_capacity, get_max_edge_sample_capacity());
    max_capacity = CGAL::max(max_capacity, get_max_vertex_sample_capacity());
  }

  FT get_max_face_sample_capacity() const {
    FT max_capacity = 0.0;
    for (typename Mesh::Face_range::const_iterator fi = mesh_.faces().begin();
      fi != mesh_.faces().end(); ++fi) {
      const Link_list &face_out_links = get_face_out_links(*fi);
      FT capacity = area(*fi) / face_out_links.size();
      max_capacity = CGAL::max(max_capacity, capacity);
    }
    return max_capacity;
  }

  FT get_max_vertex_sample_capacity() const {
    FT max_capacity = 0.0, weight = 0.0, feature_intensity = 0.0;
    for (typename Mesh::Vertex_range::const_iterator vi = mesh_.vertices().begin();
      vi != mesh_.vertices().end(); ++vi) {
      const Link &link = get_vertex_out_link(*vi);
      weight = link.first;
      feature_intensity = calculate_feature_intensity(*vi);
      max_capacity = CGAL::max(max_capacity, weight / feature_intensity);
    }
    return max_capacity;
  }

  FT get_max_edge_sample_capacity() const {
    FT max_capacity = 0.0;
    for (typename Mesh::Edge_range::const_iterator ei = mesh_.edges().begin();
      ei != mesh_.edges().end(); ++ei) {
      halfedge_descriptor hd = mesh_.halfedge(*ei);
      if (get_halfedge_normal_dihedral(hd) == -1.0) {
        hd = get_opposite(hd);
      }
      FT face_area = area(get_face(hd));
      if (!is_border(get_opposite(hd))) {
        face_area += area(get_face(get_opposite(hd)));
      }
      const Link_list &edge_out_links = get_halfedge_out_links(hd);
      FT capacity = face_area / (3 * edge_out_links.size());
      max_capacity = CGAL::max(max_capacity, capacity);
    }
    return max_capacity;
  }

  // 6) aabb tree and links
  void build_face_tree(Face_tree *face_tree) const {
    face_tree->rebuild(mesh_.faces().begin(), mesh_.faces().end(), mesh_);
    face_tree->accelerate_distance_queries();
  }

  void generate_out_links(const Face_tree &face_tree,
    int samples_per_face_value, int bvd_iteration_count_value,
    Mesh_properties *mesh_properties, const NamedParameters &np) {
    std::string type = (mesh_properties == NULL ? "out" : "in");
    // step 1: collect the edges and faces
    Face_list faces;
    Edge_list edges;
    get_all_faces(&faces);
    get_all_edges(&edges);
    // step 2: calculate the number of samples per face
    calculate_nb_samples_per_face(samples_per_face_value, faces, np);
    // step 3: generate links
    CGAL::Timer timer;
    // edge links
    timer.start();
    std::cout << "Generating edge " << type << " links...";
    generate_edge_links(face_tree, mesh_properties, edges, np);
    std::cout << "Done, count: " << get_edge_out_link_count()
      << " (" << timer.time() << " s)" << std::endl;
    // vertex links
    timer.reset();
    std::cout << "Generating vertex " << type << " links...";
    generate_vertex_links(face_tree, mesh_properties,
      np.use_stratified_sampling);
    std::cout << "Done, count: " << get_vertex_out_link_count()
      << " (" << timer.time() << " s)" << std::endl;
    // face links
    timer.reset();
    std::cout << "Generating face " << type << " links...";
    generate_face_links(face_tree, mesh_properties, bvd_iteration_count_value,
      np.use_stratified_sampling, faces);
    std::cout << "Done, count: " << get_face_out_link_count()
      << " (" << timer.time() << " s)" << std::endl;
    // reset tags
    reset_face_tags(0, faces);
  }

  void generate_local_links(const Face_tree &face_tree,
    bool reset_normal_dihedral, const Link_iter_list &face_in_links,
    const Link_iter_list &edge_in_links,
    const Link_pointer_list &vertex_in_links, vertex_descriptor vd,
    const std::set<face_descriptor> &in_link_faces,
    const NamedParameters &np) {
    // step 1: update the local feature_intensity around vd
    update_local_feature_intensity(vd, reset_normal_dihedral, np);
    // step 2: clear local links
    clear_local_links(vd, in_link_faces);
    // step 3: generate local out links
    generate_local_out_links(face_tree, vd, np);
    // step 4: generate local in links
    generate_local_in_links(in_link_faces, face_in_links, edge_in_links,
      vertex_in_links, np.use_local_aabb_tree);
  }

  void restore_local_in_links(const Point_list &face_in_end_points,
    const Link_iter_list &face_in_links, const Point_list &edge_in_end_points,
    const Link_iter_list &edge_in_links, const Point_list &vertex_in_end_points,
    const Link_pointer_list &vertex_in_links) {
    // only for simulation purpose
    Link_iter_list_const_iter lit;
    Point_const_iter pit;
    Link_pointer_const_iter link_iter;
    for (lit = face_in_links.begin(), pit = face_in_end_points.begin();
      lit != face_in_links.end(); ++lit, ++pit) {
      Link_list_iter it = *lit;
      it->second.second = *pit;
    }
    for (lit = edge_in_links.begin(), pit = edge_in_end_points.begin();
      lit != edge_in_links.end(); ++lit, ++pit) {
      Link_list_iter it = *lit;
      it->second.second = *pit;
    }
    for (link_iter = vertex_in_links.begin(),
      pit = vertex_in_end_points.begin(); link_iter != vertex_in_links.end();
      ++link_iter, ++pit) {
      Link *link = *link_iter;
      link->second.second = *pit;
    }
  }

  void backup_local_in_links(const std::set<face_descriptor> &extended_faces,
    Link_iter_list *face_in_links, Point_list *face_in_end_points,
    Link_iter_list *edge_in_links, Point_list *edge_in_end_points,
    Link_pointer_list *vertex_in_links,
    Point_list *vertex_in_end_points) const {
    // step 1: backup the in link iterations
    backup_local_in_links(extended_faces, face_in_links, edge_in_links,
      vertex_in_links);
    // step 2: store the in link end points (for simulation purpose)
    for (Link_iter_list_const_iter it = face_in_links->begin();
      it != face_in_links->end(); ++it) {
      Link_list_const_iter llit = *it;
      face_in_end_points->push_back(llit->second.second);
    }
    for (Link_iter_list_const_iter it = edge_in_links->begin();
      it != edge_in_links->end(); ++it) {
      Link_list_const_iter llit = *it;
      edge_in_end_points->push_back(llit->second.second);
    }
    for (Link_pointer_const_iter it = vertex_in_links->begin();
      it != vertex_in_links->end(); ++it) {
      const Link* link = *it;
      vertex_in_end_points->push_back(link->second.second);
    }
  }

  void clear_out_links() {
    // step 1: clear out links in faces
    for (typename Mesh::Face_range::const_iterator fi = mesh_.faces().begin();
      fi != mesh_.faces().end(); ++fi) {
      Link_list &face_out_links = get_face_out_links(*fi);
      face_out_links.clear();
    }
    // step 2: clear out links in halfedges
    for (typename Mesh::Edge_range::const_iterator ei = mesh_.edges().begin();
      ei != mesh_.edges().end(); ++ei) {
      halfedge_descriptor hd = mesh_.halfedge(*ei);
      if (get_halfedge_normal_dihedral(hd) == -1.0) {
        hd = get_opposite(hd);
      }
      Link_list &edge_out_links = get_halfedge_out_links(hd);
      edge_out_links.clear();
    }
    // step 3: clear out links in vertices (do not need to do anything)
  }

  void clear_in_link_iterators() {
    for (typename Mesh::Face_range::const_iterator fi = mesh_.faces().begin();
      fi != mesh_.faces().end(); ++fi) {
      Link_iter_list &face_in_links = get_face_in_links(*fi);
      face_in_links.clear();
      Link_iter_list &edge_in_links = get_edge_in_links(*fi);
      edge_in_links.clear();
      Link_pointer_list &vertex_in_links = get_vertex_in_links(*fi);
      vertex_in_links.clear();
    }
  }

  // 7) face collections
  void collect_all_faces(std::set<face_descriptor> *faces) const {
    for (typename Mesh::Face_range::const_iterator fi = mesh_.faces().begin();
      fi != mesh_.faces().end(); ++fi) {
      faces->insert(*fi);
    }
  }

  void extend_faces(const std::set<face_descriptor> &one_ring_faces,
    int stencil_ring_size, std::set<face_descriptor> *extended_faces) const {
    extended_faces->clear();
    extended_faces->insert(one_ring_faces.begin(), one_ring_faces.end());
    for (int i = 0; i < stencil_ring_size; ++i) {
      extend_faces_by_one_ring(extended_faces);
    }
  }

  void collect_one_ring_faces_incident_to_edge(halfedge_descriptor hd,
    std::set<face_descriptor> *faces) const {
    vertex_descriptor vp = get_source_vertex(hd);
    vertex_descriptor vq = get_target_vertex(hd);
    collect_one_ring_faces_incident_to_vertex(vp, faces);
    collect_one_ring_faces_incident_to_vertex(vq, faces);
  }

  // 8) max errors
  void calculate_max_squared_errors() {
    // precondition: mesh_ has been sampled
    reset_face_max_squared_errors(-1.0);
    update_edge_in_max_squared_errors();
    update_edge_out_max_squared_errors();
    update_vertex_in_max_squared_errors();
    update_vertex_out_max_squared_errors();
  }

  void calculate_max_squared_errors(std::set<face_descriptor> *faces) {
    // precondition: mesh_ has been sampled
    for (auto it = faces->begin(); it != faces->end(); ++it) {
      FT max_se = 0.0, se = 0.0;
      face_descriptor fd = *it;
      // face in links
      const Link_iter_list &face_in_links = get_face_in_links(fd);
      for (Link_iter_list_const_iter it = face_in_links.begin();
        it != face_in_links.end(); ++it) {
        Link_list_const_iter llit = *it;
        const Link &link = *llit;
        se = CGAL::squared_distance(link.second.first, link.second.second);
        max_se = CGAL::max(max_se, se);
      }
      // edge in links
      const Link_iter_list &edge_in_links = get_edge_in_links(fd);
      for (Link_iter_list_const_iter it = edge_in_links.begin();
        it != edge_in_links.end(); ++it) {
        Link_list_const_iter llit = *it;
        const Link &link = *llit;
        se = CGAL::squared_distance(link.second.first, link.second.second);
        max_se = CGAL::max(max_se, se);
      }
      // vertex in links
      const Link_pointer_list &vertex_in_links = get_vertex_in_links(fd);
      for (Link_pointer_const_iter it = vertex_in_links.begin();
        it != vertex_in_links.end(); ++it) {
        const Link *link = *it;
        se = CGAL::squared_distance(link->second.first, link->second.second);
        max_se = CGAL::max(max_se, se);
      }
      // face out links
      const Link_list &face_out_links = get_face_out_links(fd);
      for (Link_list_const_iter it = face_out_links.begin();
        it != face_out_links.end(); ++it) {
        const Link &link = *it;
        se = CGAL::squared_distance(link.second.first, link.second.second);
        max_se = CGAL::max(max_se, se);
      }
      halfedge_descriptor hd = mesh_.halfedge(fd);
      for (int i = 0; i <= 2; ++i) {
        // edge out links
        halfedge_descriptor hi = hd;
        if (get_halfedge_normal_dihedral(hi) == -1.0) {
          hi = get_opposite(hi);
        }
        const Link_list &edge_out_links = get_halfedge_out_links(hi);
        for (Link_list_const_iter it = edge_out_links.begin();
          it != edge_out_links.end(); ++it) {
          const Link &link = *it;
          se = CGAL::squared_distance(link.second.first, link.second.second);
          max_se = CGAL::max(max_se, se);
        }
        // vervex out links
        vertex_descriptor v = mesh_.target(hd);
        const Link &link = get_vertex_out_link(v);
        se = CGAL::squared_distance(link.second.first, link.second.second);
        max_se = CGAL::max(max_se, se);
        hd = mesh_.next(hd);
      }
      set_face_max_squared_error(fd, max_se);
    }
  }

  halfedge_descriptor calculate_local_maximal_error(
    const std::set<face_descriptor> &faces, FT *max_error) const {
    assert(!faces.empty());
    typename std::set<face_descriptor>::const_iterator cit = faces.cbegin();
    face_descriptor max_error_fd = *cit;
    FT max_se = get_face_max_squared_error(max_error_fd), se = 0.0;
    ++cit;
    for (; cit != faces.cend(); ++cit) {
      se = get_face_max_squared_error(*cit);
      if (max_se < se) {
        max_se = se;
        max_error_fd = *cit;
      }
    }
    *max_error = CGAL::sqrt(max_se);
    return get_longest_halfedge(max_error_fd);
  }

  halfedge_descriptor calculate_maximal_error(FT *max_error) const {
    typename Mesh::Face_range::const_iterator fi = mesh_.faces().begin();
    face_descriptor max_error_face = *fi;
    FT max_se = get_face_max_squared_error(max_error_face), se = 0.0;
    ++fi;
    for (; fi != mesh_.faces().end(); ++fi) {
      se = get_face_max_squared_error(*fi);
      if (se > max_se) {
        max_se = se;
        max_error_face = *fi;
      }
    }
    *max_error = CGAL::sqrt(max_se);
    return get_longest_halfedge(max_error_face);
  }

  // 9) elements properties
  Triangle triangle(face_descriptor fd) const {
    halfedge_descriptor hd = mesh_.halfedge(fd);
    const Point &a = get_point(get_target_vertex(hd));
    const Point &b = get_point(get_opposite_vertex(hd));
    const Point &c = get_point(get_source_vertex(hd));
    return Triangle(a, b, c);
  }

  FT area(face_descriptor fd) const {
    FT smallest_radian = calculate_smallest_radian(fd);
    if (smallest_radian < MIN_VALUE) {
      return 0.0;
    } else {
      return CGAL::sqrt(triangle(fd).squared_area());
    }
  }

  Point centroid(face_descriptor fd) const {
    halfedge_descriptor hd = mesh_.halfedge(fd);
    const Point &a = get_point(get_target_vertex(hd));
    const Point &b = get_point(get_opposite_vertex(hd));
    const Point &c = get_point(get_source_vertex(hd));
    return CGAL::centroid(a, b, c);
  }

  halfedge_descriptor get_longest_halfedge(face_descriptor fd) const {
    halfedge_descriptor hd = mesh_.halfedge(fd);
    halfedge_descriptor longest_hd = hd;
    FT longest_sl = squared_length(longest_hd);
    FT sl = squared_length(mesh_.next(hd));   // check the next halfedge
    if (sl > longest_sl) {
      longest_sl = sl;
      longest_hd = mesh_.next(hd);
    }
    sl = squared_length(mesh_.prev(hd));      // check the previous halfedge
    if (sl > longest_sl) {
      longest_hd = mesh_.prev(hd);
    }
    return longest_hd;
  }

  FT squared_length(halfedge_descriptor hd) const {
    const Point &a = get_point(get_source_vertex(hd));
    const Point &b = get_point(get_target_vertex(hd));
    return CGAL::squared_distance(a, b);
  }

  inline Point midpoint(halfedge_descriptor hd) const {
    const Point &a = get_point(get_source_vertex(hd));
    const Point &b = get_point(get_target_vertex(hd));
    return CGAL::midpoint(a, b);
  }

  halfedge_descriptor longest_side_propagation(halfedge_descriptor hd) const {
    // precondition: hd is the longest edge in its incient face
    if (is_border(hd) || is_border(get_opposite(hd))) {
      return hd;
    }
    face_descriptor fd = get_face(get_opposite(hd));
    halfedge_descriptor longest_in_neighbor = get_longest_halfedge(fd);
    if (longest_in_neighbor == get_opposite(hd)) {
      return hd;
    } else {
      return longest_side_propagation(longest_in_neighbor);
    }
  }

  FT calculate_feature_intensity(vertex_descriptor vd) const {
    // range: range: [1, (PI + 1) ^ 2]
    FT gaussian_curvature = get_vertex_gaussian_curvature(vd);
    FT max_halfedge_dihedral = get_vertex_max_dihedral(vd);
    return (gaussian_curvature + 1) * (max_halfedge_dihedral + 1);
  }

  VertexType get_vertex_type(vertex_descriptor vd,
    Halfedge_list *effective_edges, const NamedParameters &np) const {
    if (np.inherit_element_types) {
      collect_crease_edges_in_one_ring(vd, effective_edges);
      if (effective_edges->size() <= 1) {
        return VertexType::k_smooth_vertex;
      } else if (effective_edges->size() == 2) {
        return VertexType::k_crease_vertex;
      } else {
        return VertexType::k_feature_vertex;
      }
    } else {
      collect_effective_edges_in_one_ring(np.feature_control_delta, vd,
        effective_edges);
      if (effective_edges->size() == 0) {
        return VertexType::k_feature_vertex;
      } else if (effective_edges->size() == 2) {
        return VertexType::k_crease_vertex;
      } else {
        return VertexType::k_smooth_vertex;
      }
    }
  }

  FT calculate_minimal_radian_around_vertex(vertex_descriptor vd) const {
    FT minimal_radian = CGAL_PI;
    Halfedge_around_target_circulator hb(mesh_.halfedge(vd), mesh_), he(hb);
    do {
      if (!is_border(*hb)) {
        face_descriptor fd = get_face(*hb);
        FT radian = calculate_smallest_radian(fd);
        minimal_radian = CGAL::min(minimal_radian, radian);
      }
      ++hb;
    } while (hb != he);
    return minimal_radian;
  }

  // 10) surface mesh properties
  inline size_t size_of_vertices() const { return mesh_.number_of_vertices(); }

  void calculate_normals() {
    Normal normal;
    for (typename Mesh::Face_range::const_iterator fi = mesh_.faces().begin();
      fi != mesh_.faces().end(); ++fi) {
      normal = PMP::compute_face_normal(*fi, mesh_, PMP::parameters::
          vertex_point_map(mesh_.points()).geom_traits(Kernel()));
      set_face_normal(*fi, normal);
    }
  }

  void calculate_feature_intensities(const NamedParameters &np) {
    calculate_edge_feature_intensities(np);
    calculate_vertex_feature_intensities(np);
    if (np.inherit_element_types) {
      calculate_edge_classifications(np.feature_control_delta);
    }
  }

  halfedge_descriptor calculate_minimal_radian(FT *minimal_radian) const {
    // calculate the minimal radian of the mesh
    *minimal_radian = CGAL_PI;
    halfedge_descriptor minimal_radian_hd;
    for (typename Mesh::Face_range::const_iterator fi = mesh_.faces().begin();
      fi != mesh_.faces().end(); ++fi) {
      FT radian = calculate_smallest_radian(*fi);
      if (radian < *minimal_radian) {
        *minimal_radian = radian;
        minimal_radian_hd = get_shortest_halfedge(*fi);
      }
    }
    return minimal_radian_hd;
  }

  Bbox calculate_bounding_box() const {
    Bbox bbox = Bbox(DOUBLE_MAX, DOUBLE_MAX, DOUBLE_MAX,
      DOUBLE_MIN, DOUBLE_MIN, DOUBLE_MIN);
    typename Mesh::Vertex_range::const_iterator vi = mesh_.vertices().begin();
    bbox = get_point(*vi).bbox();
    for (; vi != mesh_.vertices().end(); ++vi) {
      bbox = bbox + get_point(*vi).bbox();
    }
    return bbox;
  }

  FT calculate_average_length() const {
    if (mesh_.number_of_edges() == 0) {
      return 0;
    } else {
      FT sum_length = 0.0;
      for (typename Mesh::Edge_range::const_iterator ei = mesh_.edges().begin();
        ei != mesh_.edges().end(); ++ei) {
        sum_length += length(mesh_.halfedge(*ei));
      }
      return sum_length / mesh_.number_of_edges();
    }
  }

  FT calculate_diagonal_length() const {
    Bbox bbox = calculate_bounding_box();
    if (bbox.xmax() <= bbox.xmin()) {   // invalid case
      return -1.0;
    } else {
      FT x_length = bbox.xmax() - bbox.xmin();
      FT y_length = bbox.ymax() - bbox.ymin();
      FT z_length = bbox.zmax() - bbox.zmin();
      return std::sqrt(x_length * x_length +
        y_length * y_length +
        z_length * z_length);
    }
  }

  // 11) mesh properties
  void trace_properties() {
    std::cout << mesh_.number_of_vertices() << " vertices" << std::endl;
    std::cout << mesh_.number_of_faces() << " faces" << std::endl;
    std::cout << mesh_.number_of_edges() << " edges" << std::endl;
    std::cout << nb_boundaries() << " boundary(ies)" << std::endl;
    std::cout << nb_components() << " component(s)" << std::endl;
    trace_edge_length();
  }

  void trace_additional_properties(FT diagonal_length) {
    if (mesh_.number_of_faces() == 0) {
      return;
    }
    FT min_face_area = DOUBLE_MAX;
    FT max_squared_error = 0.0;
    FT min_radian = CGAL_PI;
    FT max_radian = 0.0;
    FT avg_min_radian = 0.0;
    size_t min_out_link_count = MAX_VALUE;
    for (typename Mesh::Face_range::const_iterator fi = mesh_.faces().begin();
      fi != mesh_.faces().end(); ++fi) {
      FT face_area = area(*fi);
      min_face_area = CGAL::min(min_face_area, face_area);
      const Link_list &face_out_links = get_face_out_links(*fi);
      min_out_link_count = CGAL::min(min_out_link_count,
        face_out_links.size());
      max_squared_error = CGAL::max(max_squared_error,
        get_face_max_squared_error(*fi));
      FT smallest_radian = calculate_smallest_radian(*fi);
      min_radian = CGAL::min(min_radian, smallest_radian);
      avg_min_radian += smallest_radian;
      FT largest_radian = calculate_largest_radian(*fi);
      max_radian = CGAL::max(max_radian, largest_radian);
    }
    avg_min_radian /= mesh_.number_of_faces();
    FT rms_error = calculate_rms_distance();
    FT max_error = CGAL::sqrt(max_squared_error);
    // FT diagonal_length = get_diagonal_length();
    std::cout << "Minimal quality: " << calculate_min_quality() << std::endl;
    std::cout << "Average quality: " << calculate_avg_quality() << std::endl;
    std::cout << "Minimal angle: " << to_angle(min_radian) << std::endl;
    std::cout << "Average minimal angle: "
      << to_angle(avg_min_radian) << std::endl;
    std::cout << "Maximal angle: " << to_angle(max_radian) << std::endl;
    std::cout << "Maximal error: " << max_error << std::endl;
    std::cout << "Maximal error ratio: "
      << max_error / diagonal_length * 100 << "%" << std::endl;
    std::cout << "RMS face error " << rms_error << std::endl;
    std::cout << "RMS face error ratio: "
      << rms_error / diagonal_length * 100 << "%" << std::endl;
    std::cout << "Angles smaller than 30 degree: "
      << calculate_smaller_angle_ratio(30.0) * 100 << "%" << std::endl;
    std::cout << "Regular vertices ratio: "
      << calculate_regular_vertex_ratio() * 100 << "%" << std::endl;
    std::cout << "Minimal face area: " << min_face_area << std::endl;
    std::cout << "Minimal face out link count: "
      << min_out_link_count << std::endl;
  }

  // 12) static utilities
  static inline FT to_radian(FT angle) { return angle * CGAL_PI / 180.0; }

  static inline FT to_angle(FT radian) { return radian * 180.0 / CGAL_PI; }

  // 13) input operations
  int eliminate_degenerated_faces() {
    FT radian_threshold = MIN_VALUE;
    // step 1: fill all the degenerated faces in the dynamic priority queue
    DPQueue_face_long degenerated_faces;
    fill_degenerated_faces_queue(radian_threshold, &degenerated_faces);
    // step 2: eliminate these degenerations one by one
    int nb_eliminations = static_cast<int>(degenerated_faces.size());
    FT longest_squared_length = 0.0;
    typename DPQueue_face_long::right_map::iterator eit;
    while (!degenerated_faces.empty()) {
      eit = degenerated_faces.right.begin();
      face_descriptor fd = eit->second;
      degenerated_faces.right.erase(eit);
      halfedge_descriptor shortest_hd = get_shortest_halfedge(fd);
      // case 1: the triangle is acute, so we only need to collpase
      if (squared_length(shortest_hd) < SQUARED_MIN_VALUE) {
        // remove faces before collapse
        remove_incident_faces(shortest_hd, &degenerated_faces);
        // collapse the edge
        Point new_point = midpoint(shortest_hd);
        vertex_descriptor vd = collapse_short_edge(new_point, shortest_hd);
        if (vd == get_null_vertex()) {
          return -1;
        }
        // add incident face to the queue if it is still generated face
        add_circulator_degenerate_faces(vd, radian_threshold,
          &degenerated_faces);
      } else {
        // case 2: the triangle is obtuse, so we split the longest edge first
        halfedge_descriptor longest_hd = get_longest_halfedge(fd);
        longest_squared_length = squared_length(longest_hd);
        Point new_point = get_point(get_opposite_vertex(longest_hd));
        // remove faces before split
        degenerated_faces.left.erase(fd);
        if (!is_border(get_opposite(longest_hd))) {
          face_descriptor fd = get_face(get_opposite(longest_hd));
          degenerated_faces.left.erase(fd);
        }
        // split the longest halfedge
        halfedge_descriptor hnew = split_long_edge(new_point, longest_hd);
        // add the new degenerated faces after split
        add_circulator_degenerate_faces(get_target_vertex(hnew),
          radian_threshold, &degenerated_faces);
      }
    }
    return nb_eliminations;
  }

  int split_long_edges() {
    // step 1: fill dynamic priority queue with long edges in front
    FT average_length = calculate_average_length();
    FT max_length = 4.0 / 3.0 * average_length;
    FT max_squared_length = max_length * max_length;
    DPQueue_halfedge_long long_edges;
    fill_queue_with_long_edges(max_squared_length, &long_edges);
    // step 2: split long edges one by one
    int nb_split = 0;
    while (!long_edges.empty()) {
      // get a copy of the candidate edge
      typename DPQueue_halfedge_long::right_map::iterator eit =
        long_edges.right.begin();
      halfedge_descriptor hd = eit->second;
      // remove the hd and its opposite from queue
      long_edges.left.erase(hd);
      long_edges.left.erase(get_opposite(hd));
      // split the specified long edge
      Point new_point = midpoint(hd);
      halfedge_descriptor hnew = split_long_edge(new_point, hd);
      // update the queue
      add_circular_long_edges(max_squared_length, get_target_vertex(hnew),
        &long_edges);
      ++nb_split;
    }
    return nb_split;
  }

  // 14) remesh operations

  // 14.1) split
  vertex_descriptor split_edge(const Face_tree &input_face_tree,
    FT max_error_threshold_value, FT max_error, FT min_radian,
    bool reduce_complexity, DPQueue_halfedge_long *large_error_queue,
    DPQueue_halfedge_short *small_value_queue, halfedge_descriptor hd,
    const NamedParameters &np) {
    // max_error > 0 means to reduce error; othewiese improve radian
    // step 1: backup the original in_links and the edge types
    std::set<face_descriptor> one_ring_faces, extended_faces;
    one_ring_faces.insert(get_face(hd));
    if (!is_border(get_opposite(hd))) {
      one_ring_faces.insert(get_face(get_opposite(hd)));
    }
    extend_faces(one_ring_faces, np.stencil_ring_size, &extended_faces);
    Link_iter_list face_in_links, edge_in_links;
    Link_pointer_list vertex_in_links;
    backup_local_in_links(extended_faces, &face_in_links,
      &edge_in_links, &vertex_in_links);
    std::map<vertex_descriptor, bool> crease_map;
    if (np.inherit_element_types) {
      bool is_crease = get_halfedge_normal_dihedral(hd) == -1.0 ?
        get_halfedge_is_crease(get_opposite(hd)) : get_halfedge_is_crease(hd);
      crease_map[get_source_vertex(hd)] = is_crease;
      crease_map[get_target_vertex(hd)] = is_crease;
    }
    // step 2: remove from queue if necessary
    if (small_value_queue != NULL) {
      remove_small_value_edges_before_split(hd, large_error_queue,
        small_value_queue, np);
    }
    // step 3: split the edge (also update edge types)
    Point new_point = midpoint(hd);
    halfedge_descriptor hnew = split_long_edge(new_point, hd);
    vertex_descriptor vd = get_target_vertex(hnew);
    if (np.inherit_element_types) {
      Halfedge_around_target_circulator hb(mesh_.halfedge(vd), mesh_), he(hb);
      typename std::map<vertex_descriptor, bool>::const_iterator it;
      do {
        it = crease_map.find(get_source_vertex(*hb));
        if (it != crease_map.end()) {   // splitted halfedge
          halfedge_descriptor hd = *hb;
          if (is_border(hd)) {
            hd = get_opposite(hd);
          }
          set_halfedge_is_crease(hd, it->second);
          set_halfedge_normal_dihedral(hd, -0.5);   // temperary mark
          set_halfedge_normal_dihedral(get_opposite(hd), -1.0);
        }
        ++hb;
      } while (hb != he);
    }
    // step 4: update local links (optimize the position if necessary)
    one_ring_faces.clear();
    extended_faces.clear();
    collect_one_ring_faces_incident_to_vertex(vd, &one_ring_faces);
    extend_faces(one_ring_faces, np.stencil_ring_size, &extended_faces);
    generate_local_links(input_face_tree, false, face_in_links, edge_in_links,
                         vertex_in_links, vd, extended_faces, np);
    if (max_error > 0 || np.optimize_after_local_operations) {
      // if reduce error, we definitely optimize; otherweise it depends
      optimize_vertex_position(input_face_tree, face_in_links,
        edge_in_links, vertex_in_links, vd, extended_faces, np);
    }
    calculate_max_squared_errors(&extended_faces);
    // step 5: add to queue if necessary
    if (small_value_queue != NULL) {
      add_small_value_edges_after_split(vd, max_error_threshold_value,
        reduce_complexity, large_error_queue, small_value_queue, np);
    }
    // step 6: update normals
    calculate_local_normals(&one_ring_faces);
    // step 7: flip or relocate if necessary
    if (max_error < 0) {
      if (np.flip_after_split_and_collapse) {
        Halfedge_list halfedges;
        Halfedge_around_target_circulator hb(mesh_.halfedge(vd), mesh_), he(hb);
        do {
          if (!is_border(*hb)) {
            halfedges.push_back(mesh_.prev(*hb));
          }
          ++hb;
        } while (hb != he);
        flip_edges(input_face_tree, max_error_threshold_value, max_error,
          min_radian, reduce_complexity, large_error_queue,
          small_value_queue, &halfedges, np);
      } else if (np.relocate_after_local_operations) {
        std::set<vertex_descriptor> one_ring_vertices;
        collect_incident_vertices(vd, &one_ring_vertices);
        relocate_vertices(input_face_tree, max_error_threshold_value,
          max_error, min_radian, reduce_complexity, large_error_queue,
          small_value_queue, one_ring_vertices, np);
      }
    }
    return vd;
  }

  // 14.2) collapse
  bool is_collapsable(halfedge_descriptor hd) const {
    if (is_border(hd)) {
      return false;
    }
    halfedge_descriptor ho = get_opposite(hd);
    if (is_border(ho)) {    // border case
      vertex_descriptor vs = get_target_vertex(mesh_.next(ho));
      vertex_descriptor vt = get_source_vertex(mesh_.prev(ho));
      if (vs == vt) {
        return false;   // do not close a degree-3 hole
      }
    } else {                        // inner case
      vertex_descriptor vr = get_opposite_vertex(hd);
      vertex_descriptor vs = get_opposite_vertex(ho);
      if (vr == vs) {
        return false;
      }
      // inner edge but two boundary vertices
      if (is_on_boundary(get_target_vertex(hd)) &&
        is_on_boundary(get_target_vertex(ho))) {
        return false;
      }
    }
    return check_link_condition(hd);   // link condition
  }

  bool collapse_would_cause_wrinkle(const Halfedge_list &halfedges,
    const Point &new_point, halfedge_descriptor hd) const {
    for (auto it = halfedges.begin(); it != halfedges.end(); ++it) {
      halfedge_descriptor h = *it;
      const Point &start_point = get_point(get_source_vertex(h));
      const Point &end_point = get_point(get_target_vertex(h));
      const Point &old_point = get_point(get_opposite_vertex(h));
      FT radian = calculate_radian(end_point, old_point, start_point);
      if (radian < MIN_VALUE || CGAL_PI - radian < MIN_VALUE) {
        return true;    // degenerate cases
      }
      Plane plane(end_point, old_point, start_point);
      Point projection = plane.projection(new_point);
      if (same_side(start_point, end_point, old_point, projection) < 0) {
        return true;
      }
    }
    return false;
  }

  bool predict_faces_after_collapse(halfedge_descriptor hd,
    Halfedge_list *halfedges) const {
    // return whether faces compose a ring, halfedges are inserted in order
    // we use the halfedges to represent the faces
    vertex_descriptor vp = get_source_vertex(hd);
    vertex_descriptor vq = get_target_vertex(hd);
    if (is_border(get_opposite(hd))) {
      if (!is_border(get_opposite(mesh_.next(hd)))) {
        halfedge_descriptor h_add =
          mesh_.next(get_opposite(mesh_.prev(get_opposite(hd))));
        while (h_add != mesh_.prev(hd)) {
          halfedges->push_back(h_add);
          h_add = mesh_.next(get_opposite(mesh_.next(h_add)));
        }
      }
      if (!is_border(get_opposite(mesh_.prev(hd)))) {
        halfedge_descriptor h_add = mesh_.next(get_opposite(mesh_.prev(hd)));
        while (h_add !=
          mesh_.prev(get_opposite(mesh_.next(get_opposite(hd))))) {
          halfedges->push_back(h_add);
          h_add = mesh_.next(get_opposite(mesh_.next(h_add)));
        }
        halfedges->push_back(h_add);
      }
      return false;
    } else if (is_on_boundary(vp)) {
      Halfedge_around_target_circulator hb(mesh_.halfedge(vp), mesh_), he(hb);
      do {
        if (is_border(*hb)) {
          break;
        }
        ++hb;
      } while (hb != he);
      halfedge_descriptor h_add = mesh_.next(get_opposite(*hb));
      while (h_add != mesh_.prev(get_opposite(hd))) {
        halfedges->push_back(h_add);
        h_add = mesh_.next(get_opposite(mesh_.next(h_add)));
      }
      h_add = mesh_.next(get_opposite(h_add));
      while (h_add != mesh_.prev(hd)) {
        halfedges->push_back(h_add);
        h_add = mesh_.next(get_opposite(mesh_.next(h_add)));
      }
      if (!is_border(get_opposite(h_add))) {
        h_add = mesh_.next(get_opposite(h_add));
        while (!is_border(get_opposite(mesh_.next(h_add)))) {
          halfedges->push_back(h_add);
          h_add = mesh_.next(get_opposite(mesh_.next(h_add)));
        }
        halfedges->push_back(h_add);
      }
      return false;
    } else if (is_on_boundary(vq)) {
      Halfedge_around_target_circulator hb(mesh_.halfedge(vq), mesh_), he(hb);
      do {
        if (is_border(*hb)) {
          break;
        }
        ++hb;
      } while (hb != he);
      halfedge_descriptor h_add = mesh_.next(get_opposite(*hb));
      while (h_add != mesh_.prev(hd)) {
        halfedges->push_back(h_add);
        h_add = mesh_.next(get_opposite(mesh_.next(h_add)));
      }
      h_add = mesh_.next(get_opposite(h_add));
      while (h_add != mesh_.prev(get_opposite(hd))) {
        halfedges->push_back(h_add);
        h_add = mesh_.next(get_opposite(mesh_.next(h_add)));
      }
      if (!is_border(get_opposite(h_add))) {
        h_add = mesh_.next(get_opposite(h_add));
        while (!is_border(get_opposite(mesh_.next(h_add)))) {
          halfedges->push_back(h_add);
          h_add = mesh_.next(get_opposite(mesh_.next(h_add)));
        }
        halfedges->push_back(h_add);
      }
      return false;
    } else {
      halfedge_descriptor h_add = mesh_.next(get_opposite(mesh_.prev(hd)));
      while (h_add != mesh_.prev(get_opposite(hd))) {
        halfedges->push_back(h_add);
        h_add = mesh_.next(get_opposite(mesh_.next(h_add)));
      }
      h_add = mesh_.next(get_opposite(h_add));
      while (h_add != mesh_.prev(hd)) {
        halfedges->push_back(h_add);
        h_add = mesh_.next(get_opposite(mesh_.next(h_add)));
      }
      return true;
    }
  }

  Point calculate_initial_point_for_collapse(
    halfedge_descriptor hd, const NamedParameters &np) const {
    // If inherit_element_type is true, we use the vertex_type;
    // Otherwise, we use the qem to get intial point for collapse
    if (np.inherit_element_types) {  // use the inherited vertex types
      vertex_descriptor vp = get_source_vertex(hd);
      vertex_descriptor vq = get_target_vertex(hd);
      Halfedge_list effective_edges;
      VertexType vtp = get_vertex_type(vp, &effective_edges, np);
      VertexType vtq = get_vertex_type(vq, &effective_edges, np);
      if (vtp == vtq) {
        // option: use the qem instend
        return midpoint(hd);
      } else {
        if (vtp == VertexType::k_feature_vertex) {
          return get_point(vp);
        } else if (vtp == VertexType::k_crease_vertex) {
          return vtq == VertexType::k_feature_vertex ? get_point(vq) :
            get_point(vp);
        } else {
          return get_point(vq);
        }
      }
    } else {                        // use the qem instead
      return get_least_qem_point(hd);
    }
    // backup option:
    /*FT fi_p = vp->feature_intensity();
    FT fi_q = vq->feature_intensity();
    if (CGAL::abs(fi_p - fi_q) <
    CGAL::max(fi_p, fi_q) * feature_difference_delta_) {
    return midpoint(hh);
    } else {
    return fi_p > fi_q ? vp->point() : vq->point();
    }*/
  }

  vertex_descriptor construct_local_mesh(
    const std::set<face_descriptor> &one_ring_faces,
    const std::set<face_descriptor> &extended_faces,
    const Halfedge_list &halfedges, const Point &new_point, bool is_ring,
    Mesh *mesh) const {
    // construct the 2-manifold surface_mesh
    // step 1: build the point to index map
    std::map<Point, size_t, Point_Comp> points_map;
    // std::set<Point, size_t> points_map;
    points_map.insert(std::make_pair(new_point, points_map.size()));
    for (auto it = halfedges.begin(); it != halfedges.end(); ++it) {
      const Point &p = get_point(get_target_vertex(*it));
      points_map.insert(std::make_pair(p, points_map.size()));
    }
    if (!is_ring) {
      halfedge_descriptor hd = *(halfedges.begin());
      hd = get_opposite(hd);
      const Point &p = get_point(get_target_vertex(hd));
      points_map.insert(std::make_pair(p, points_map.size()));
    }
    std::set<face_descriptor> differ_faces;
    std::set_difference(extended_faces.begin(), extended_faces.end(),
      one_ring_faces.begin(), one_ring_faces.end(),
      std::inserter(differ_faces, differ_faces.end()));
    for (auto it = differ_faces.begin(); it != differ_faces.end(); ++it) {
      halfedge_descriptor hd = mesh_.halfedge(*it);
      Point points[3];
      points[0] = get_point(get_target_vertex(hd));
      points[1] = get_point(get_opposite_vertex(hd));
      points[2] = get_point(get_source_vertex(hd));
      for (int i = 0; i <= 2; ++i) {
        if (points_map.find(points[i]) == points_map.end()) {
          points_map.insert(std::make_pair(points[i], points_map.size()));
        }
      }
    }
    // step 2: build the polygon soups
    std::vector<size_t> polygon(3, 0);
    std::vector<std::vector<size_t>> polygons;
    for (auto it = halfedges.begin(); it != halfedges.end(); ++it) {
      const Point &p1 = get_point(get_source_vertex(*it));
      const Point &p2 = get_point(get_target_vertex(*it));
      polygon[1] = points_map[p1];
      polygon[2] = points_map[p2];
      polygons.push_back(polygon);
    }
    for (auto it = differ_faces.begin(); it != differ_faces.end(); ++it) {
      halfedge_descriptor hd = mesh_.halfedge(*it);
      const Point &p0 = get_point(get_target_vertex(hd));
      const Point &p1 = get_point(get_opposite_vertex(hd));
      const Point &p2 = get_point(get_source_vertex(hd));
      polygon[0] = points_map[p0];
      polygon[1] = points_map[p1];
      polygon[2] = points_map[p2];
      polygons.push_back(polygon);
    }
    // step 3: build the 2-manifold mesh
    std::vector<Point> points(points_map.size());
    for (auto it = points_map.begin(); it != points_map.end(); ++it) {
      points[it->second] = it->first;
    }
    // PMP::repair_polygon_soup(points, polygons);
    // PMP::orient_polygon_soup(points, polygons);
    PMP::polygon_soup_to_polygon_mesh(points, polygons, *mesh);
    for (typename Mesh::Vertex_range::const_iterator vi = mesh->vertices().begin();
      vi != mesh->vertices().end(); ++vi) {
      if (mesh->point(*vi) == new_point) {
        return *vi;
      }
    }
    return mesh->null_vertex();
  }

  vertex_descriptor collapse_edge(const Face_tree &input_face_tree,
    FT max_error_threshold_value, FT min_radian, bool reduce_complexity,
    DPQueue_halfedge_long *large_error_queue,
    DPQueue_halfedge_short *small_value_queue,
    const Link_iter_list &face_in_links, const Link_iter_list &edge_in_links,
    const Link_pointer_list &vertex_in_links, halfedge_descriptor hd,
    Point new_point, const NamedParameters &np) {
    // step 1: backup the edge type if necessary
    std::map<vertex_descriptor, bool> crease_map;
    if (np.inherit_element_types) {
      vertex_descriptor vd = get_opposite_vertex(hd);
      crease_map[vd] = is_crease_edge(mesh_.prev(hd)) ||
        is_crease_edge(mesh_.next(hd));
      if (!is_border(get_opposite(hd))) {
        vd = get_opposite_vertex(get_opposite(hd));
        crease_map[vd] = is_crease_edge(mesh_.prev(get_opposite(hd))) ||
          is_crease_edge(mesh_.next(get_opposite(hd)));
      }
    }
    // step 2: remove from queue if necessary
    if (small_value_queue != NULL) {
      remove_small_value_edges_before_collapse(hd, large_error_queue,
        small_value_queue, np);
    }
    // step 3: collapse the edge (also update edge types here)
    vertex_descriptor v_joined = collapse_short_edge(new_point, hd);
    if (v_joined == get_null_vertex()) {
      return get_null_vertex();
    }
    if (np.inherit_element_types) {
      Halfedge_around_target_circulator
        hb(mesh_.halfedge(v_joined), mesh_), he(hb);
      typename std::map<vertex_descriptor, bool>::const_iterator it;
      do {
        it = crease_map.find(get_source_vertex(*hb));
        if (it != crease_map.end()) {  // the merged halfedge
          halfedge_descriptor hd = *hb;
          if (is_border(hd)) {
            hd = get_opposite(hd);
          }
          set_halfedge_is_crease(hd, it->second);
          set_halfedge_normal_dihedral(hd, -0.5);   // temperary mark
          set_halfedge_normal_dihedral(get_opposite(hd), -1.0);
        }
        ++hb;
      } while (hb != he);
    }
    // step 4: update local links (optimize already performed in simulation)
    std::set<face_descriptor> one_ring_faces, extended_faces;
    collect_one_ring_faces_incident_to_vertex(v_joined, &one_ring_faces);
    extend_faces(one_ring_faces, np.stencil_ring_size, &extended_faces);
    generate_local_links(input_face_tree, true, face_in_links, edge_in_links,
      vertex_in_links, v_joined, extended_faces, np);
    calculate_max_squared_errors(&extended_faces);
    // step 5: add to queue if necessary
    if (small_value_queue != NULL) {
      add_small_value_edges_after_collapse(v_joined, max_error_threshold_value,
        reduce_complexity, large_error_queue, small_value_queue, np);
    }
    // step 6: update normals
    calculate_local_normals(&one_ring_faces);
    // step 7: flip or relocate if necessary
    if (min_radian > 0) {
      if (np.flip_after_split_and_collapse) {
        Halfedge_list halfedges;
        Halfedge_around_target_circulator hb(mesh_.halfedge(v_joined), mesh_),
                                          he(hb);
        do {
          if (!is_border(*hb)) {
            halfedges.push_back(*hb);
          }
          ++hb;
        } while (hb != he);
        int flip_num = flip_edges(input_face_tree, max_error_threshold_value,
            -1.0, min_radian, reduce_complexity, large_error_queue,
            small_value_queue, &halfedges, np);
      } else if (np.relocate_after_local_operations) {
        std::set<vertex_descriptor> one_ring_vertices;
        collect_incident_vertices(v_joined, &one_ring_vertices);
        relocate_vertices(input_face_tree, max_error_threshold_value, -1.0,
          min_radian, reduce_complexity, large_error_queue,
          small_value_queue, one_ring_vertices, np);
      }
    }
    return v_joined;
  }

  // 14.3) flip
  int flip_applied(const Face_tree &input_face_tree,
    FT max_error_threshold_value, FT max_error, FT min_radian,
    bool reduce_complexity, DPQueue_halfedge_long *large_error_queue,
    DPQueue_halfedge_short *small_radian_queue, halfedge_descriptor hd,
    const NamedParameters &np) {
    // step 1: construct the incident_edges
    Halfedge_list incident_edges;   // two incident edges will be included
    FT sd_prev = squared_length(mesh_.prev(hd));
    FT sd_next = squared_length(mesh_.next(hd));
    if (sd_prev >= sd_next) {
      incident_edges.push_back(mesh_.prev(hd));
      incident_edges.push_back(mesh_.next(hd));
    } else {
      incident_edges.push_back(mesh_.next(hd));
      incident_edges.push_back(mesh_.prev(hd));
    }
    // step 2: try to flip
    return flip_edges(input_face_tree, max_error_threshold_value, max_error,
      min_radian, reduce_complexity, large_error_queue, small_radian_queue,
      &incident_edges, np);
  }

  halfedge_descriptor flip_edge(const Face_tree &input_face_tree,
    FT max_error_threshold_value, FT max_error, FT min_radian,
    bool reduce_complexity, DPQueue_halfedge_long *large_error_queue,
    DPQueue_halfedge_short *small_value_queue, halfedge_descriptor hd,
    const NamedParameters &np) {
    // max_error > 0 means we want to reduce error; otherwise improve radian
    // step 1: topology constraints and min_radian constraints check
    if (!is_flippable(hd, np)) {
      return get_null_halfedge();
    }
    // step 2: backup the original local links
    std::set<face_descriptor> one_ring_faces, extended_faces;
    one_ring_faces.insert(get_face(hd));
    one_ring_faces.insert(get_face(get_opposite(hd)));
    extend_faces(one_ring_faces, np.stencil_ring_size, &extended_faces);
    Link_iter_list face_in_links, edge_in_links;
    Link_pointer_list vertex_in_links;
    backup_local_in_links(extended_faces, &face_in_links,
      &edge_in_links, &vertex_in_links);
    // step 3: remove from queue if necessary
    if (small_value_queue != NULL) {
      remove_small_value_edges_before_flip(hd, large_error_queue,
        small_value_queue, np);
    }
    // step 4: flip (do not optimize the vertex position)
    flip_inner_edge(hd);
    // step 5: update local links
    one_ring_faces.clear();
    extended_faces.clear();
    one_ring_faces.insert(get_face(hd));
    one_ring_faces.insert(get_face(get_opposite(hd)));
    extend_faces(one_ring_faces, np.stencil_ring_size, &extended_faces);
    generate_local_links(input_face_tree, false, face_in_links, edge_in_links,
      vertex_in_links, hd, extended_faces, np);
    calculate_max_squared_errors(&extended_faces);
    FT error_after;
    calculate_local_maximal_error(extended_faces, &error_after);
    // step 6: rollback if the flip violated the max_error constraints;
    bool flipped;
    if (max_error > 0) {
      flipped = error_after < max_error;
    } else {
      flipped = error_after < max_error_threshold_value;
    }
    if (!flipped) {
      flip_inner_edge(hd);
      one_ring_faces.clear();
      extended_faces.clear();
      one_ring_faces.insert(get_face(hd));
      one_ring_faces.insert(get_face(get_opposite(hd)));
      extend_faces(one_ring_faces, np.stencil_ring_size, &extended_faces);
      generate_local_links(input_face_tree, false, face_in_links,
        edge_in_links, vertex_in_links, hd, extended_faces, np);
      calculate_max_squared_errors(&extended_faces);
    }
    // step 7: add to queue if necessary
    if (small_value_queue != NULL) {
      add_small_value_edges_after_flip(hd, max_error_threshold_value,
        reduce_complexity, large_error_queue, small_value_queue, np);
    }
    // step 8: update the normals
    calculate_local_normals(&one_ring_faces);
    // step 9: relocate if necessary (only when we want to improve angle)
    if (max_error < 0) {
      if (flipped && np.relocate_after_local_operations) {
        std::set<vertex_descriptor> vertices;
        vertices.insert(get_source_vertex(hd));
        vertices.insert(get_target_vertex(hd));
        vertices.insert(get_opposite_vertex(hd));
        vertices.insert(get_opposite_vertex(get_opposite(hd)));
        std::set<face_descriptor> faces;  // used to calculate local_min_radian
        for (auto it = vertices.begin(); it != vertices.end(); ++it) {
          collect_one_ring_faces_incident_to_vertex(*it, &faces);
        }
        FT local_min_radian = calculate_local_minimal_radian(faces);
        relocate_vertices(input_face_tree, max_error_threshold_value,
          max_error, local_min_radian, reduce_complexity, large_error_queue,
          small_value_queue, vertices, np);
      }
    }
    return flipped ? hd : get_null_halfedge();
  }

  // 14.4) relocate
  int relocate_applied(const Face_tree &input_face_tree,
    FT max_error_threshold_value, FT max_error, FT min_radian,
    bool reduce_complexity, DPQueue_halfedge_long *large_error_queue,
    DPQueue_halfedge_short *small_value_queue, halfedge_descriptor hd,
    const NamedParameters &np) {
    // If max_error > 0, we reduce error; otherwise, we improve min_radian
    // If reduce_complexity, we update collapse_candidate_queue;
    //    otherwise, we update small_radian_queue.
    // step 1: construct the vertices
    std::set<vertex_descriptor> incident_vertices;
    incident_vertices.insert(get_opposite_vertex(hd));
    incident_vertices.insert(get_target_vertex(hd));
    incident_vertices.insert(get_source_vertex(hd));
    // step 2: relocate these incident vertices
    return relocate_vertices(input_face_tree, max_error_threshold_value,
      max_error, min_radian, reduce_complexity, large_error_queue,
      small_value_queue, incident_vertices, np);
  }

  bool relocate_vertex(const Face_tree &input_face_tree,
    FT max_error_threshold_value, FT max_error, FT min_radian,
    bool reduce_complexity, DPQueue_halfedge_long *large_error_queue,
    DPQueue_halfedge_short *small_value_queue, const Point &new_point,
    vertex_descriptor vd, const NamedParameters &np) {
    // If max_error > 0, we want to reduce error; otherwiese increase radian
    // If reduce_complexity, we update collapse_candidate_queue;
    //    otherwise, we update small_radian_queue.
    // step 1: geometry constraints check
    if (CGAL::squared_distance(get_point(vd), new_point) <
      SQUARED_MIN_VALUE) {
      return false;
    }
    if (np.keep_vertex_in_one_ring &&
      relocate_would_cause_wrinkle(new_point, vd)) {
      return false;
    }
    // step 2: backup the original local links
    std::set<face_descriptor> one_ring_faces, extended_faces;
    collect_one_ring_faces_incident_to_vertex(vd, &one_ring_faces);
    extend_faces(one_ring_faces, np.stencil_ring_size, &extended_faces);
    Link_iter_list face_in_links, edge_in_links;
    Link_pointer_list vertex_in_links;
    std::map<face_descriptor, Link_list> face_out_map;
    std::map<halfedge_descriptor, Link_list> edge_out_map;
    Link vertex_out;
    backup_local_links(extended_faces, &face_out_map, &edge_out_map,
      &vertex_out, &face_in_links, &edge_in_links, &vertex_in_links, vd);
    // step 3: remove from queue if necessary
    if (small_value_queue != NULL) {
      remove_small_value_edges_before_relocate(vd, large_error_queue,
        small_value_queue, np);
    }
    // step 4: relocate the vertex position
    Point old_point = get_point(vd);  // backup position in case rolling back
    relocate_vertex_point(vd, new_point);
    // step 5: update local links (optimized included if necessary)
    generate_local_links(input_face_tree, false, face_in_links,
      edge_in_links, vertex_in_links, vd, extended_faces, np);
    if (max_error > 0 || np.optimize_after_local_operations) {
      // if reduce error, we definitely optimize; otherweise it depends
      optimize_vertex_position(input_face_tree, face_in_links, edge_in_links,
        vertex_in_links, vd, extended_faces, np);
    }
    calculate_max_squared_errors(&extended_faces);
    FT error, radian;
    calculate_local_maximal_error(extended_faces, &error);
    radian = calculate_minimal_radian_around_vertex(vd);
    // step 6: rollback if the max_error or min_radian is violated
    FT radian_delta = to_radian(np.smooth_angle_delta);
    bool relocated;
    if (max_error > 0) {      // only reduce the error
      relocated = (error < max_error);
    } else {                  // improve the min_radian
      relocated = (error < max_error_threshold_value &&
        radian >= min_radian + radian_delta);
    }
    if (!relocated) {
      relocate_vertex_point(vd, old_point);   // restore back the position
      restore_local_links(face_out_map, edge_out_map, vertex_out,
        face_in_links, edge_in_links, vertex_in_links,
        vd, extended_faces, np);
      calculate_max_squared_errors(&extended_faces);
    }
    // step 7: add the queue if necesary
    if (small_value_queue != NULL) {
      add_small_value_edges_after_relocate(vd, max_error_threshold_value,
        reduce_complexity, large_error_queue, small_value_queue, np);
    }
    // step 8: update the normals
    calculate_local_normals(&one_ring_faces);
    return relocated;
  }

  void optimize_vertex_position(const Face_tree &input_face_tree,
    const Link_iter_list &face_in_links, const Link_iter_list &edge_in_links,
    const Link_pointer_list &vertex_in_links, vertex_descriptor vd,
    const std::set<face_descriptor> &in_link_faces,
    const NamedParameters &np) {
    // precondition: the samples have been generated
    for (int i = 0; i < np.vertex_optimize_count; ++i) {
      // step 1: optimize the position of vh
      if (!optimize_vertex_by_one_ring(input_face_tree, vd, np)) {
        break;
      }
      // step 2: update the samples
      generate_local_links(input_face_tree, false, face_in_links,
        edge_in_links, vertex_in_links, vd, in_link_faces, np);
    }
  }

  Point calculate_initial_point_for_relocate(const Face_tree &face_tree,
    vertex_descriptor vd, const NamedParameters &np) const {
    FT sum_area = 0.0;
    Halfedge_around_target_circulator hb(mesh_.halfedge(vd), mesh_), he(hb);
    do {
      if (!is_border(*hb)) {
        sum_area += area(get_face(*hb));
      }
      ++hb;
    } while (hb != he);
    Point initial_point;
    if (sum_area <= 1.0 / np.max_samples_per_area * mesh_.degree(vd)) {
      // if the area is too small, return its average center
      initial_point = average_center(vd);
    } else {
      if (np.relocate_strategy == RelocateStrategy::k_cvt_barycenter) {
        initial_point = cvt_barycenter(vd, np);
      } else {
        initial_point = barycenter(vd, np);
      }
    }
    // if it is a on a crease, then it is already on the face_tree
    return face_tree.closest_point(initial_point);
  }

 private:
  // 1) static functions
  static FT area(const Point &a, const Point &b, const Point &c) {
    if (CGAL::squared_distance(a, b) < SQUARED_MIN_VALUE ||
      CGAL::squared_distance(a, c) < SQUARED_MIN_VALUE ||
      CGAL::squared_distance(b, c) < SQUARED_MIN_VALUE) {
      return 0.0;
    }
    Triangle t(a, b, c);
    return CGAL::sqrt(t.squared_area());
  }

  static FT calculate_smallest_radian(
      const Point &a, const Point &b, const Point &c) {
    FT ab = CGAL::squared_distance(a, b);
    FT ac = CGAL::squared_distance(a, c);
    FT bc = CGAL::squared_distance(b, c);
    if (ab < ac) {
      return ab < bc ? calculate_radian(b, c, a) : calculate_radian(c, a, b);
    } else {
      return ac < bc ? calculate_radian(a, b, c) : calculate_radian(c, a, b);
    }
  }

  static FT calculate_radian(const Point &a, const Point &b, const Point &c) {
    if (CGAL::squared_distance(a, b) < SQUARED_MIN_VALUE ||
      CGAL::squared_distance(a, c) < SQUARED_MIN_VALUE ||
      CGAL::squared_distance(b, c) < SQUARED_MIN_VALUE) {
      return 0.0;   // degenerated case
    }
    Vector v1 = a - b;        // v1
    FT v1_length = std::sqrt(v1 * v1);
    if (v1_length < MIN_VALUE) {
      return 0.0;
    }
    v1 = v1 / v1_length;
    Vector v2 = c - b;        // v2
    FT v2_length = std::sqrt(v2 * v2);
    if (v2_length < MIN_VALUE) {
      return 0.0;
    }
    v2 = v2 / v2_length;
    FT cos_value = v1 * v2;
    if (cos_value > 1.0) {
      cos_value = 1.0;
    }
    if (cos_value < -1.0) {
      cos_value = -1.0;
    }
    return std::acos(cos_value);
  }

  static FT calculate_angle(const Point &a, const Point &b, const Point &c) {
    FT radian = get_radian(a, b, c);
    return to_angle(radian);
  }

  // static FT calculate_radian(const Vector &v1, const Vector &v2) {
  //   if (v1 == CGAL::NULL_VECTOR || v2 == CGAL::NULL_VECTOR) {
  //     return 0.0;     // degenerated case
  //   }
  //   FT v1_length = std::sqrt(v1.squared_length());
  //   FT v2_length = std::sqrt(v2.squared_length());
  //   if (v1_len < MIN_VALUE || v2_len < MIN_VALUE) {
  //     return 0.0;
  //   }
  //   FT inner_product = (v1 * v2) / (v1_length * v2_length);
  //   if (inner_product > 1.0) {
  //     inner_product = 1.0;
  //   }
  //   if (inner_product < -1.0) {
  //     inner_product = -1.0;
  //   }
  //   return std::acos(inner_product);
  // }

  static FT calculate_angle(const Vector &v1, const Vector &v2) {
    FT radian = get_radian(v1, v2);
    return to_angle(radian);
  }

  static int same_side(const Point &a, const Point &b, const Point &c,
                       const Point &p) {
    // precondition: a, b, c and p are on the same plane
    // determine whether c and p lays on the same side of ab
    // 1: the same side;
    // 0: c or p is on the line ab or projected on the line;
    // -1: the opposite side
    Vector ab = b - a;
    Vector ac = c - a;
    Vector ap = p - a;
    Vector v1 = CGAL::cross_product(ab, ac);
    Vector v2 = CGAL::cross_product(ab, ap);
    FT cp = v1 * v2;
    if (cp > 0) {
      return 1;
    } else if (cp == 0.0) {
      return 0;
    } else {
      return -1;
    }
  }

  static Point calculate_nearest_point(const Segment &s, const Point &p) {
    Line line = s.supporting_line();
    if (CGAL::squared_distance(s, p) == CGAL::squared_distance(line, p)) {
      return line.projection(p);
    } else {
      return CGAL::squared_distance(s.source(), p) <
        CGAL::squared_distance(s.target(), p) ? s.source() : s.target();
    }
  }

  static bool point_in_triangle(const Point &a, const Point &b, const Point &c,
                                const Point &p) {
    return same_side(a, b, c, p) > 0 &&
      same_side(b, c, a, p) > 0 &&
      same_side(c, a, b, p) > 0;
  }

  // 2) elements access
  inline vertex_descriptor get_target_vertex(halfedge_descriptor hd)
      { return mesh_.target(hd); }

  inline vertex_descriptor get_target_vertex(halfedge_descriptor hd) const
      { return mesh_.target(hd); }

  inline vertex_descriptor get_source_vertex(halfedge_descriptor hd)
      { return mesh_.source(hd); }

  inline vertex_descriptor get_source_vertex(halfedge_descriptor hd) const
      { return mesh_.source(hd); }

  void get_all_faces(Face_list *faces) const {
    for (typename Mesh::Face_range::const_iterator fi = mesh_.faces().begin();
      fi != mesh_.faces().end(); ++fi) {
      faces->push_back(*fi);
    }
  }

  void get_all_halfedges(Halfedge_list *halfedges) const {
    for (typename Mesh::Halfedge_range::const_iterator hi = mesh_.halfedges().begin();
      hi != mesh_.halfedges().end(); ++hi) {
      halfedges->push_back(*hi);
    }
  }

  void get_all_edges(Edge_list *edges) const {
    for (typename Mesh::Edge_range::const_iterator ei = mesh_.edges().begin();
      ei != mesh_.edges().end(); ++ei) {
      edges->push_back(*ei);
    }
  }

  void get_all_vertices(Vertex_list *vertices) const {
    for (typename Mesh::Vertex_range::const_iterator vi = mesh_.vertices().begin();
      vi != mesh_.vertices().end(); ++vi) {
      vertices->push_back(*vi);
    }
  }

  // 3) element properties

  // 3.1) vertex properties
  bool is_on_boundary(vertex_descriptor vd) const {
    Halfedge_around_target_circulator hb(mesh_.halfedge(vd), mesh_), he(hb);
    do {
      if (is_border(*hb)) {  // only one halfedge will answer true
        return true;
      }
      ++hb;
    } while (hb != he);
    return false;
  }

  bool is_regular(vertex_descriptor vd) const {
    size_t degree = mesh_.degree(vd);
    if (is_on_boundary(vd)) {
      return degree >= 3 && degree <= 5;
      // return degree == 4;
    } else {
      return degree >= 5 && degree <= 7;
      // return degree == 6;
    }
  }

  bool are_neighbors(vertex_descriptor vd1, vertex_descriptor vd2) const {
    Halfedge_around_target_circulator hb(mesh_.halfedge(vd1), mesh_), he(hb);
    do {
      if (get_source_vertex(*hb) == vd2) {
        return true;
      }
      ++hb;
    } while (hb != he);
    return false;
  }

  Point barycenter(vertex_descriptor vd, const NamedParameters &np) const {
    Vector vec = CGAL::NULL_VECTOR;
    Point pivot = get_point(vd);
    Halfedge_list effective_edgs;
    VertexType vertex_type = get_vertex_type(vd, &effective_edgs, np);
    if (vertex_type == VertexType::k_feature_vertex) {
      return pivot;
    } else if (vertex_type == VertexType::k_crease_vertex) {
      assert(effective_edgs.size() == 2);
      // step 1: get the initial_point
      for (Halfedge_const_iter cit = effective_edgs.begin();
        cit != effective_edgs.end(); ++cit) {
        halfedge_descriptor hd = *cit;
        const Point &neighbor = get_point(get_source_vertex(hd));
        vec = vec + (neighbor - pivot);
      }
      Point initial_point = pivot + vec / 2.0;
      // step 2: project the initial_point
      Point min_projection(DOUBLE_MAX, DOUBLE_MAX, DOUBLE_MAX);
      FT min_sd = DOUBLE_MAX;
      for (Halfedge_const_iter cit = effective_edgs.begin();
        cit != effective_edgs.end(); ++cit) {
        Segment segment(get_point(get_source_vertex(*cit)),
          get_point(get_target_vertex(*cit)));
        Point nearest_point = calculate_nearest_point(segment, initial_point);
        FT sd = CGAL::squared_distance(initial_point, nearest_point);
        if (sd < min_sd) {
          min_sd = sd;
          min_projection = nearest_point;
        }
      }
      return min_projection;
    } else {
      Halfedge_around_target_circulator hb(mesh_.halfedge(vd), mesh_), he(hb);
      do {
        const Point &neighbor = get_point(get_source_vertex(*hb));
        vec = vec + (neighbor - pivot);
        ++hb;
      } while (hb != he);
      return pivot + vec / mesh_.degree(vd);
    }
  }

  Point average_center(vertex_descriptor vd) const {
    // step 1: check whether it is a boundary vertex
    Halfedge_list boundary_edges, inner_edges;
    Halfedge_around_target_circulator hb(mesh_.halfedge(vd), mesh_), he(hb);
    do {
      if (!is_border(*hb) || !is_border(get_opposite(*hb))) {
        boundary_edges.push_back(*hb);
      }
      inner_edges.push_back(*hb);   // in case vd is not on boundary
      ++hb;
    } while (hb != he);
    // step 2: get the barycenter according to different cases
    const Point &pivot = get_point(vd);
    Vector vec = CGAL::NULL_VECTOR;
    if (!boundary_edges.empty()) {  // boundary case
      for (Halfedge_const_iter cit = boundary_edges.begin();
        cit != boundary_edges.end(); ++cit) {
        const Point &neighbor = get_point(get_source_vertex(*cit));
        vec = vec + (neighbor - pivot);
      }
      return pivot + vec / boundary_edges.size();
    } else {                        // inner case
      for (Halfedge_const_iter cit = inner_edges.begin();
        cit != inner_edges.end(); ++cit) {
        const Point &neighbor = get_point(get_source_vertex(*cit));
        vec = vec + (neighbor - pivot);
      }
      return pivot + vec / inner_edges.size();
    }
  }

  Point cvt_barycenter(vertex_descriptor vd, const NamedParameters &np) const {
    Vector vec = CGAL::NULL_VECTOR;
    Point pivot = get_point(vd);
    FT denominator = 0.0;
    Halfedge_list effective_edges;
    VertexType vertex_type = get_vertex_type(vd, &effective_edges, np);
    if (vertex_type == VertexType::k_feature_vertex) {
      return pivot;
    } else if (vertex_type == VertexType::k_crease_vertex) {
      assert(effective_edges.size() == 2);
      // step 1: get the initial_point
      for (Halfedge_const_iter cit = effective_edges.begin();
        cit != effective_edges.end(); ++cit) {
        vec = vec + (midpoint(*cit) - pivot);
      }
      Point initial_point = pivot + vec / 2;
      // step 2: project the initial_point
      Point min_projection(DOUBLE_MAX, DOUBLE_MAX, DOUBLE_MAX);
      FT min_sd = DOUBLE_MAX;
      for (Halfedge_const_iter cit = effective_edges.begin();
        cit != effective_edges.end(); ++cit) {
        Segment segment(get_point(get_source_vertex(*cit)),
          get_point(get_target_vertex(*cit)));
        Point nearest_point = calculate_nearest_point(segment, initial_point);
        FT sd = CGAL::squared_distance(initial_point, nearest_point);
        if (sd < min_sd) {
          min_sd = sd;
          min_projection = nearest_point;
        }
      }
      return min_projection;
    } else {    // vertex_type == VertexType::k_smooth_vertex
      Halfedge_around_target_circulator hb(mesh_.halfedge(vd), mesh_), he(hb);
      Vector vector;
      do {
        halfedge_descriptor h1 = *hb;
        halfedge_descriptor h2 = mesh_.next(h1);
        halfedge_descriptor h3 = mesh_.prev(h1);
        const Point &ph1 = get_point(get_target_vertex(h1));
        const Point &ph2 = get_point(get_target_vertex(h2));
        const Point &ph3 = get_point(get_target_vertex(h3));
        const Point p = CGAL::centroid(ph1, ph2, ph3);
        const Point p1 = midpoint(h1);
        const Point p2 = midpoint(h2);
        const Point c1 = CGAL::centroid(pivot, p, p1);
        const Point c2 = CGAL::centroid(pivot, p2, p);
        const FT area1 = area(pivot, p, p1);
        const FT area2 = area(pivot, p2, p);
        vector = area1 * (c1 - CGAL::ORIGIN) + area2 * (c2 - CGAL::ORIGIN);
        vec = vec + vector;
        denominator += (area1 + area2);
        ++hb;
      } while (hb != he);
      if (denominator < SQUARED_MIN_VALUE) {
        return pivot;
      } else {
        return CGAL::ORIGIN + vec / denominator;
      }
    }
  }

  FT calculate_vertex_capacity(vertex_descriptor vd) const {
    FT vertex_capacity = 0.0;
    const Point &p = get_point(vd);
    Halfedge_around_target_circulator hb(mesh_.halfedge(vd), mesh_), he(hb);
    do {
      if (!is_border(*hb)) {
        face_descriptor fd = get_face(*hb);
        const Point p1 = midpoint(*hb);
        const Point c = centroid(fd);  // option: use the circumcenter instead?
        const Point p2 = midpoint(mesh_.next(*hb));
        vertex_capacity += area(p1, p, c);
        vertex_capacity += area(c, p, p2);
      }
      ++hb;
    } while (hb != he);
    return vertex_capacity;
  }

  FT calculate_minimal_radian_incident_to_vertex(vertex_descriptor vd) const {
    FT minimal_radian = CGAL_PI;
    Halfedge_around_target_circulator hb(mesh_.halfedge(vd), mesh_), he(hb);
    do {
      if (!is_border(*hb)) {
        FT radian = calculate_radian(get_point(get_source_vertex(*hb)),
          get_point(get_target_vertex(*hb)),
          get_point(get_opposite_vertex(*hb)));
        minimal_radian = CGAL::min(minimal_radian, radian);
      }
      ++hb;
    } while (hb != he);
    return minimal_radian;
  }

  FT calculate_min_squared_distance_in_one_ring_faces(
      vertex_descriptor vd) const {
    FT min_sd = std::numeric_limits<double>::max();
    Halfedge_around_target_circulator hb(mesh_.halfedge(vd), mesh_), he(hb);
    do {
      if (!is_border(*hb)) {
        Line line(get_point(get_source_vertex(*hb)),
          get_point(get_opposite_vertex(*hb)));
        FT sd = CGAL::squared_distance(get_point(vd), line);
        if (sd < min_sd) {
          min_sd = sd;
        }
      }
      ++hb;
    } while (hb != he);
    return min_sd;
  }

  // 3.2) face properties
  halfedge_descriptor get_shortest_halfedge(face_descriptor fd) {
    halfedge_descriptor hd = mesh_.halfedge(fd);
    halfedge_descriptor shortest_hd = hd;
    FT shortest_sl = squared_length(shortest_hd);
    FT sl = squared_length(mesh_.next(hd));   // check the next halfedge
    if (sl < shortest_sl) {
      shortest_sl = sl;
      shortest_hd = mesh_.next(hd);
    }
    sl = squared_length(mesh_.prev(hd));      // check the previous halfedge
    if (sl < shortest_sl) {
      shortest_hd = mesh_.prev(hd);
    }
    return shortest_hd;
  }

  halfedge_descriptor get_shortest_halfedge(face_descriptor fd) const {
    halfedge_descriptor hd = mesh_.halfedge(fd);
    halfedge_descriptor shortest_hd = hd;
    FT shortest_sl = squared_length(shortest_hd);
    FT sl = squared_length(mesh_.next(hd));   // check the next halfedge
    if (sl < shortest_sl) {
      shortest_sl = sl;
      shortest_hd = mesh_.next(hd);
    }
    sl = squared_length(mesh_.prev(hd));      // check the previous halfedge
    if (sl < shortest_sl) {
      shortest_hd = mesh_.prev(hd);
    }
    return shortest_hd;
  }

  FT calculate_smallest_radian(face_descriptor fd) const {
    halfedge_descriptor shortest_hd = get_shortest_halfedge(fd);
    return calculate_opposite_radian(shortest_hd);
  }

  FT calculate_smallest_angle(face_descriptor fd) const {
    FT radian = get_smallest_radian(fd);
    return to_angle(radian);
  }

  FT calculate_largest_radian(face_descriptor fd) const {
    halfedge_descriptor longest_hd = get_longest_halfedge(fd);
    return calculate_opposite_radian(longest_hd);
  }

  FT calculate_largest_angle(face_descriptor fd) const {
    FT radian = get_largest_radian(fd);
    return to_angle(radian);
  }

  FT calculate_quality(face_descriptor fd) const {
    halfedge_descriptor hd = mesh_.halfedge(fd);
    FT s_t = area(fd);                          // area of the triangle
    FT h_t = length(get_longest_halfedge(fd));  // longest edge
    FT a = length(hd);                          // length of the first edge
    FT b = length(mesh_.next(hd));              // length of the second edge
    FT c = length(mesh_.prev(hd));              // length of the third edge
    FT p_t = (a + b + c) / 2.0;                 // halfedge perimeter
    if (h_t <= 0.0) {   // invalid case
      return 0.0;
    } else {
      return 6.0 * s_t / (CGAL::sqrt(3.0) * p_t * h_t);
    }
  }

  FT calculate_local_minimal_radian(
    const std::set<face_descriptor> &faces) const {
    FT minimal_radian = CGAL_PI;
    for (auto it = faces.begin(); it != faces.end(); ++it) {
      face_descriptor fd = *it;
      FT radian = calculate_smallest_radian(fd);
      minimal_radian = CGAL::min(minimal_radian, radian);
    }
    return minimal_radian;
  }

  FT calculate_sum_qem_value(const std::set<face_descriptor> &faces,
    const Point &point) const {
    FT sum_qem_value = 0.0;
    for (auto it = faces.begin(); it != faces.end(); ++it) {
      if (area(*it) >= SQUARED_MIN_VALUE) {  // ignore too small faces
        halfedge_descriptor hd = mesh_.halfedge(*it);
        Plane plane(get_point(get_target_vertex(hd)),
          get_point(get_opposite_vertex(hd)),
          get_point(get_source_vertex(hd)));
        sum_qem_value += CGAL::squared_distance(point, plane);
      }
    }
    return sum_qem_value;
  }

  FT calculate_sum_area(const Face_list &faces) const {
    FT sum_area = 0.0;
    for (auto it = faces.begin(); it != faces.end(); ++it) {
      sum_area += area(*it);
    }
    return sum_area;
  }

  void calculate_local_normals(std::set<face_descriptor> *faces) {
    Normal normal;
    for (auto fit = faces->begin(); fit != faces->end(); ++fit) {
      normal = PMP::compute_face_normal(*fit, mesh_, PMP::parameters::
          vertex_point_map(mesh_.points()).geom_traits(Kernel()));
      set_face_normal(*fit, normal);
    }
  }

  // 3.3) halfedge properties
  inline FT calculate_opposite_radian(halfedge_descriptor hd) const {
    if (is_border(hd)) {
      std::cout << "Error definition of opposite angle for a border."
        << std::endl;
      return -1;
    } else {
      return calculate_radian(get_point(get_target_vertex(hd)),
        get_point(get_opposite_vertex(hd)),
        get_point(get_source_vertex(hd)));
    }
  }

  inline FT calculate_opposite_angle(halfedge_descriptor hd) const {
    FT radian = get_opposite_radian(hd);
    return to_angle(radian);
  }

  inline bool is_border(halfedge_descriptor hd) const
      { return mesh_.is_border(hd); }

  inline FT length(halfedge_descriptor hd) const {
    return CGAL::sqrt(squared_length(hd));
  }

  inline bool is_crease_edge(halfedge_descriptor hd) const {
    if (get_halfedge_normal_dihedral(hd) == -1.0) {
      hd = get_opposite(hd);
    }
    return get_halfedge_is_crease(hd);
  }

  Point get_least_qem_point(halfedge_descriptor hd) const {
    std::set<face_descriptor> faces;
    collect_one_ring_faces_incident_to_edge(hd, &faces);
    const Point &start_point = get_point(get_source_vertex(hd));
    const Point mid_point = midpoint(hd);
    const Point &end_point = get_point(get_target_vertex(hd));
    FT qem_start = calculate_sum_qem_value(faces, start_point);
    FT qem_mid = calculate_sum_qem_value(faces, mid_point);
    FT qem_end = calculate_sum_qem_value(faces, end_point);
    if (qem_mid <= qem_start) {
      return qem_mid <= qem_end ? mid_point : end_point;
    } else {
      return qem_start <= qem_end ? start_point : end_point;
    }
  }

  FT calculate_normal_dihedral(halfedge_descriptor hd) const {
    // get the dihedral of the normals between the two incident faces
    CGAL_precondition(!is_border(hd) && !is_border(get_opposite(hd)));
    face_descriptor fd1 = get_face(hd);
    face_descriptor fd2 = get_face(get_opposite(hd));
    const Normal &n1 = get_face_normal(fd1);
    const Normal &n2 = get_face_normal(fd2);
    FT cos_value = n1 * n2;
    if (cos_value > 1.0) {
      cos_value = 1.0;
    } else if (cos_value < -1.0) {
      cos_value = -1.0;
    }
    return std::acos(cos_value);      // expressed in radians [0, PI]
  }

  // 4) samples and links
  void calculate_nb_samples_per_face(int nb_samples_per_face_value,
    const Face_list &faces, const NamedParameters &np) {
    // the nb_sample is recorded in face_tags_
    size_t nb_faces = faces.size();
    FT total_area = calculate_sum_area(faces);
    if (total_area < SQUARED_MIN_VALUE) {   // degenerated case
      for (auto it = faces.begin(); it != faces.end(); ++it) {
        set_face_tag(*it, 0);
      }
      return;
    }
    if (np.sample_strategy == SampleStrategy::k_uniform) {
      // the number of samples per face is propotional to its area
      size_t nb_samples = nb_samples_per_face_value * nb_faces;
      for (auto it = faces.begin(); it != faces.end(); ++it) {
        face_descriptor fd = *it;
        FT face_area = area(fd);
        int nb_face_samples =
          static_cast<int>(nb_samples * face_area / total_area);
        int nb_max_samples = np.max_samples_per_area * face_area;
        nb_face_samples = std::min(nb_face_samples, nb_max_samples);
        nb_face_samples = std::max(nb_face_samples,
          np.min_samples_per_triangle);
        set_face_tag(fd, nb_face_samples);
      }
    } else {
      // SampleStrategy::k_adaptive, samples per face is the same
      for (auto it = faces.begin(); it != faces.end(); ++it) {
        face_descriptor fd = *it;
        std::set<face_descriptor> incident_faces;
        collect_incident_faces(fd, &incident_faces);
        FT sum_area = 0.0;
        for (auto it2 = incident_faces.begin();
          it2 != incident_faces.end(); ++it2) {
          sum_area += area(*it2);
        }
        FT face_area = area(fd);
        int nb_face_samples = static_cast<int>(incident_faces.size() *
          nb_samples_per_face_value * face_area / sum_area);
        int nb_max_samples = np.max_samples_per_area * face_area;
        nb_face_samples = std::min(nb_face_samples, nb_max_samples);
        nb_face_samples = std::max(nb_face_samples,
          np.min_samples_per_triangle);
        set_face_tag(fd, nb_face_samples);
      }
    }
  }

  void clear_local_links(halfedge_descriptor hd,
    const std::set<face_descriptor> &in_link_faces) {
    // step 1: clear all in link faces
    clear_local_in_links(in_link_faces);
    // step 2: clear the out links in faces that are incident to hd
    if (get_halfedge_normal_dihedral(hd) == -1.0) {
      hd = get_opposite(hd);
    }
    Link_list &edge_out_links1 = get_halfedge_out_links(hd);
    edge_out_links1.clear();
    Link_list &edge_out_links2 = get_halfedge_out_links(get_opposite(hd));
    edge_out_links2.clear();
    if (!is_border(hd)) {
      Link_list &face_out_links = get_face_out_links(get_face(hd));
      face_out_links.clear();
    }
    if (!is_border(get_opposite(hd))) {
      Link_list &face_out_links =
        get_face_out_links(get_face(get_opposite(hd)));
      face_out_links.clear();
    }
  }

  void clear_local_links(vertex_descriptor vd,
    const std::set<face_descriptor> &in_link_faces) {
    // step 1: clear all in links in in_link_faces
    clear_local_in_links(in_link_faces);
    // step 2: clear the out links in one-ring faces of vd
    Halfedge_around_target_circulator hb(mesh_.halfedge(vd), mesh_), he(hb);
    do {
      halfedge_descriptor hd = *hb;
      if (!is_border(hd)) {         // clear face out links
        Link_list &face_out_links = get_face_out_links(get_face(hd));
        face_out_links.clear();
      }
      if (get_halfedge_normal_dihedral(hd) == -1.0) {
        hd = get_opposite(hd);
      }
      Link_list &edge_out_links = get_halfedge_out_links(hd);
      edge_out_links.clear();       // clear edge out links
      ++hb;
    } while (hb != he);
  }

  void clear_local_in_links(const std::set<face_descriptor> &faces) {
    for (auto it = faces.begin(); it != faces.end(); ++it) {
      // step 1: clear the links
      face_descriptor fd = *it;
      Link_iter_list &face_in_links = get_face_in_links(fd);
      Link_iter_list &edge_in_links = get_edge_in_links(fd);
      Link_pointer_list &vertex_in_inks = get_vertex_in_links(fd);
      face_in_links.clear();
      edge_in_links.clear();
      vertex_in_inks.clear();
      // step 2: update the max squared_error
      set_face_max_squared_error(fd, 0.0);
    }
  }

  void clear_face_links(face_descriptor fd) {
    set_face_out_links(fd, Link_list());
    set_face_in_links(fd, Link_iter_list());
    set_edge_in_links(fd, Link_iter_list());
    set_vertex_in_links(fd, Link_pointer_list());
  }

  void clear_halfedge_links(halfedge_descriptor hd) {
    set_halfege_out_links(hd, Link_list());
  }

  void clear_vertex_links(vertex_descriptor vd) {
    set_vertex_out_link(vd, Link());
  }

  void generate_face_links(const Face_tree &face_tree,
    Mesh_properties *mesh_properties, int bvd_iteration_count_value,
    bool use_stratified_sampling, const Face_list &faces) {
    Point_list inner_samples;
    Point_const_iter pit;
    std::list<double> feature_weights;
    std::list<double>::const_iterator fit;
    Point_and_primitive_id pp;
    for (auto it = faces.begin(); it != faces.end(); ++it) {
      face_descriptor fd = *it;
      int nb_samples = get_face_tag(fd);
      if (nb_samples > 0) {
        generate_random_samples(use_stratified_sampling, nb_samples,
          bvd_iteration_count_value, fd, &inner_samples, &feature_weights);
        FT capacity = use_stratified_sampling ? area(fd) / nb_samples : 1.0;
        Link_list &face_out_links = get_face_out_links(fd);
        for (pit = inner_samples.begin(), fit = feature_weights.begin();
          pit != inner_samples.end(); ++pit, ++fit) {
          pp = face_tree.closest_point_and_primitive(*pit);
          // 1) insert the sample in the source
          Link_list_iter it = face_out_links.insert(face_out_links.end(),
            std::make_pair(capacity * (*fit),
            std::make_pair(*pit, pp.first)));
          // 2) insert the samples in the target if necessary
          if (mesh_properties != NULL) {
            face_descriptor closest_fd = pp.second;   // closest fd
            mesh_properties->get_face_in_links(closest_fd).push_back(it);
          }
        }
      }
    }
  }

  void generate_random_samples(bool use_stratified_sampling, int nb_samples,
    int bvd_iteration_count_value, face_descriptor fd,
    Point_list *inner_samples, std::list<double> *feature_weights) const {
    // step 1: generate the samples
    inner_samples->clear();
    feature_weights->clear();
    halfedge_descriptor hd = mesh_.halfedge(fd);
    const Point &a = get_point(get_target_vertex(hd));
    const Point &b = get_point(get_opposite_vertex(hd));
    const Point &c = get_point(get_source_vertex(hd));
    if (nb_samples <= 0) {
      return;
    } else if (nb_samples == 1) {
      inner_samples->push_back(CGAL::centroid(a, b, c));
    } else if (nb_samples == 2) {
      FT bc = CGAL::squared_distance(b, c);
      FT ca = CGAL::squared_distance(c, a);
      FT ab = CGAL::squared_distance(a, b);
      Point d;
      int longest_edge_index = 0;
      if (bc > ca) {
        longest_edge_index = bc > ab ? 0 : 2;
      } else {
        longest_edge_index = ab < ca ? 1 : 2;
      }
      switch (longest_edge_index) {
      case 0:   // split in edge bc
        d = midpoint(mesh_.prev(hd));
        inner_samples->push_back(CGAL::centroid(a, b, d));
        inner_samples->push_back(CGAL::centroid(a, d, c));
        break;
      case 1:   // split in edge ca
        d = midpoint(mesh_.halfedge(fd));
        inner_samples->push_back(CGAL::centroid(b, c, d));
        inner_samples->push_back(CGAL::centroid(b, d, a));
        break;
      case 2:   // split in edge ab
        d = midpoint(mesh_.next(hd));
        inner_samples->push_back(CGAL::centroid(c, a, d));
        inner_samples->push_back(CGAL::centroid(c, d, b));
        break;
      default:
        break;
      }
    } else {
      // step 1.1: generate the unique inner samples
      std::set<Point, Point_Comp> samples;
      // std::set<Point> samples;
      Random random;
      srand(time(NULL));
      Vector ab = b - a, ac = c - a;  // edge vectors
      while (samples.size() < nb_samples) {
        FT u = 0.0, v = 0.0;
        u = random.random_value<double>(0.0, 1.0);
        u = 0.9 * u + 0.05;
        while (v == 0.0 || u + v == 1.0) {
          v = random.random_value<double>(0.0, 1.0);
          v = 0.9 * v + 0.05;
        }
        if (u + v > 1.0) {    // flip over diag if needed
          u = 1.0 - u, v = 1.0 - v;
        }
        samples.insert(a + u * ab + v * ac);
      }
      // step 1.2: convert to Point_list data structure
      Point_list all_samples(samples.begin(), samples.end());
      // step 3: relocate them using BVD iteration
      for (int i = 0; i < bvd_iteration_count_value; ++i) {
        Bvd bvd(triangle(fd));
        bvd.run(all_samples);
      }
      // step 1.3: pick the inner samples
      inner_samples->clear();
      for (auto it = all_samples.begin(); it != all_samples.end(); ++it) {
        inner_samples->push_back(*it);
      }
    }
    // step 2: calcualte the feature weights
    FT fi_a = calculate_feature_intensity(get_target_vertex(hd));
    FT fi_b = calculate_feature_intensity(get_opposite_vertex(hd));
    FT fi_c = calculate_feature_intensity(get_source_vertex(hd));
    FT feature_weight = 0.0;
    // version 1: each sample has the same feature weight (efficient)
    /*feature_weight = (fi_a + fi_b + fi_c) / 3.0;
    for (Point_iter it = inner_samples->begin();
    it != inner_samples->end(); ++it) {
    feature_weights->push_back(feature_weight);
    }*/
    // version 2: each sample has different feature weights (accurate)
    FT face_area = area(fd);
    FT weight_a, weight_b, weight_c;
    for (Point_const_iter it = inner_samples->begin();
      it != inner_samples->end(); ++it) {
      const Point &p = *it;
      weight_a = area(p, b, c) / face_area;
      weight_b = area(p, c, a) / face_area;
      weight_c = area(p, a, b) / face_area;
      feature_weight = weight_a * fi_a + weight_b * fi_b + weight_c * fi_c;
      feature_weights->push_back(feature_weight);
    }
  }

  void generate_edge_links(const Face_tree &face_tree,
    Mesh_properties *mesh_properties, const Edge_list &edges,
    const NamedParameters &np) {
    // precondition: the number of samples per face has been calculated
    // If mesh_properties is NULL, it means is_in_link is false
    for (auto it = edges.begin(); it != edges.end(); ++it) {
      halfedge_descriptor hd = mesh_.halfedge(*it);
      if (get_halfedge_normal_dihedral(hd) == -1.0) {
        hd = get_opposite(hd);
      }
      face_descriptor fd = get_face(hd);
      FT face_area = area(fd);
      int nb_edge_out_links = get_face_tag(fd);
      if (!is_border(get_opposite(hd))) {
        fd = get_face(get_opposite(hd));
        face_area += area(fd);
        nb_edge_out_links += get_face_tag(fd);
      }
      if (nb_edge_out_links > 0) {
        FT area_per_sample = face_area / nb_edge_out_links;
        FT diameter = 2 * CGAL::sqrt(area_per_sample / CGAL_PI);
        int nb_samples = length(hd) / diameter;
        nb_samples = std::max(nb_samples, np.min_samples_per_triangle);
        FT capacity = np.use_stratified_sampling ?
          face_area / (3 * nb_samples) : 1.0;
        sample_edge_links(face_tree, mesh_properties, capacity, nb_samples, hd);
      }
    }
  }

  void sample_edge_links(const Face_tree &face_tree,
    Mesh_properties *mesh_properties, FT capacity, int nb_samples,
    halfedge_descriptor hd) {
    CGAL_precondition(get_halfedge_normal_dihedral(hd) != -1.0);
    vertex_descriptor vp = get_source_vertex(hd), vq = get_target_vertex(hd);
    const Point &p = get_point(vp), &q = get_point(vq);
    Vector vector = q - p;
    FT step = 1.0 / (nb_samples + 1);
    Point_and_primitive_id pp;
    Link_list &edge_out_links = get_halfedge_out_links(hd);
    for (int i = 1; i <= nb_samples; ++i) {
      const Point sample = p + vector * step * i;
      pp = face_tree.closest_point_and_primitive(sample);
      FT feature_weight = calculate_feature_intensity(vq) * i +
        calculate_feature_intensity(vp) * (nb_samples + 1 - i);
      feature_weight /= (nb_samples + 1);   // interpolated feature intensity
      // 1) insert the sample in the source
      Link_list_iter it = edge_out_links.insert(edge_out_links.end(),
        std::make_pair(feature_weight * capacity,
        std::make_pair(sample, pp.first)));
      // 2) insert the iterator in the target if necessary
      if (mesh_properties != NULL) {
        face_descriptor fd = pp.second;  // closest fd
        mesh_properties->get_edge_in_links(fd).push_back(it);
      }
    }
  }

  void generate_vertex_links(const Face_tree &face_tree,
    Mesh_properties *mesh_properties, bool use_stratified_sampling) {
    for (typename Mesh::Vertex_range::const_iterator vi = mesh_.vertices().begin();
      vi != mesh_.vertices().end(); ++vi) {
      generate_vertex_link(face_tree, mesh_properties,
        use_stratified_sampling, *vi);
    }
  }

  void generate_vertex_link(const Face_tree &face_tree,
    Mesh_properties *mesh_properties, bool use_stratified_sampling,
    vertex_descriptor vd) {
    // precondition: the weights of the vertices have been calculated
    Point_and_primitive_id pp;
    FT capacity = use_stratified_sampling ? calculate_vertex_capacity(vd) : 1.0;
    pp = face_tree.closest_point_and_primitive(get_point(vd));
    // 1) insert the sample in the source
    Link& vertex_out_link = get_vertex_out_link(vd);
    vertex_out_link.first = capacity * calculate_feature_intensity(vd);
    vertex_out_link.second.first = get_point(vd);
    vertex_out_link.second.second = pp.first;
    // 2) insert the sample in the target if necessary
    if (mesh_properties != NULL) {
      face_descriptor fd = pp.second;  // closest fd
      mesh_properties->get_vertex_in_links(fd).push_back(&vertex_out_link);
    }
  }

  void generate_local_links(const Face_tree &face_tree,
    bool reset_normal_dihedral, const Link_iter_list &face_in_links,
    const Link_iter_list &edge_in_links,
    const Link_pointer_list &vertex_in_links,
    halfedge_descriptor hd, const std::set<face_descriptor> &in_link_faces,
    const NamedParameters &np) {
    // step 1: update the local feature_intensity around endpoints of hd
    vertex_descriptor vp = get_source_vertex(hd);
    vertex_descriptor vq = get_target_vertex(hd);
    update_local_feature_intensity(vp, reset_normal_dihedral, np);
    update_local_feature_intensity(vq, reset_normal_dihedral, np);
    // step 2: clear local links
    clear_local_links(hd, in_link_faces);
    // step 3: generate local out links
    generate_local_out_links(face_tree, hd, np);
    // step 4: generate local in links
    generate_local_in_links(in_link_faces, face_in_links, edge_in_links,
      vertex_in_links, np.use_local_aabb_tree);
  }

  void generate_local_out_links(const Face_tree &face_tree,
    vertex_descriptor vd, const NamedParameters &np) {
    // step 1: collect the one_ring faces and edges
    Face_list faces;
    Edge_list edges;
    collect_faces_incident_to_vertex(vd, &faces);
    collect_edges_incident_to_vertex(vd, &edges);
    // step 2: calculate the number of samples per face
    calculate_nb_samples_per_face(np.samples_per_face_out, faces, np);
    // step 3: generate out links (is_in_link is set to false)
    generate_edge_links(face_tree, NULL, edges, np);
    generate_vertex_links(face_tree, NULL, np.use_stratified_sampling);
    generate_face_links(face_tree, NULL, np.bvd_iteration_count,
      np.use_stratified_sampling, faces);
    // step 4: reset the face tags
    reset_face_tags(0, faces);
  }

  void generate_local_out_links(const Face_tree &face_tree,
    halfedge_descriptor hd, const NamedParameters &np) {
    // step 1: collect the faces and edges
    Face_list faces;
    Edge_list edges;
    collect_faces_incident_to_edge(hd, &faces);
    edges.push_back(mesh_.edge(hd));
    // step 2: calculate the number of samples per face
    calculate_nb_samples_per_face(np.samples_per_face_out, faces, np);
    // step 3: generate out links (is_in_links is set to false)
    generate_edge_links(face_tree, NULL, edges, np);
    generate_face_links(face_tree, NULL, np.bvd_iteration_count,
      np.use_stratified_sampling, faces);
    // step 4: reset the face tags
    reset_face_tags(0, faces);
  }

  void generate_local_in_links(const std::set<face_descriptor> &in_link_faces,
    const Link_iter_list &face_in_links, const Link_iter_list &edge_in_links,
    const Link_pointer_list &vertex_in_links, bool use_local_aabb_tree) {
    Point_and_primitive_id pp;
    face_descriptor closest_fd;
    if (use_local_aabb_tree) {
      // step 1: build the local aabb tree
      Face_tree local_face_tree;
      local_face_tree.rebuild(in_link_faces.begin(),
        in_link_faces.end(), mesh_);
      local_face_tree.accelerate_distance_queries();
      // step 2: update the in links
      for (Link_iter_list_const_iter it = face_in_links.begin();
        it != face_in_links.end(); ++it) {
        Link_list_iter llit = *it;
        pp = local_face_tree.closest_point_and_primitive(llit->second.first);
        llit->second.second = pp.first;   // update the closest point
        closest_fd = pp.second;
        get_face_in_links(closest_fd).push_back(llit);
      }
      for (Link_iter_list_const_iter it = edge_in_links.begin();
        it != edge_in_links.end(); ++it) {
        Link_list_iter llit = *it;
        pp = local_face_tree.closest_point_and_primitive(llit->second.first);
        llit->second.second = pp.first;   // update the closest point
        closest_fd = pp.second;
        get_edge_in_links(closest_fd).push_back(llit);
      }
      for (Link_pointer_const_iter it = vertex_in_links.begin();
        it != vertex_in_links.end(); ++it) {
        Link *link = *it;
        pp = local_face_tree.closest_point_and_primitive(link->second.first);
        link->second.second = pp.first;   // update the closest point
        closest_fd = pp.second;
        get_vertex_in_links(closest_fd).push_back(link);
      }
    } else {
      // for each sample, update its closest point, push back to the primitive
      for (Link_iter_list_const_iter it = face_in_links.begin();
        it != face_in_links.end(); ++it) {
        Link_list_iter llit = *it;
        pp = get_closest_point_and_primitive(in_link_faces,
          llit->second.first);
        llit->second.second = pp.first;   // update the closest point
        closest_fd = pp.second;
        get_face_in_links(closest_fd).push_back(llit);
      }
      for (Link_iter_list_const_iter it = edge_in_links.begin();
        it != edge_in_links.end(); ++it) {
        Link_list_iter llit = *it;
        pp = get_closest_point_and_primitive(in_link_faces,
          llit->second.first);
        llit->second.second = pp.first;   // update the closest point
        closest_fd = pp.second;
        get_edge_in_links(closest_fd).push_back(llit);
      }
      for (Link_pointer_const_iter it = vertex_in_links.begin();
        it != vertex_in_links.end(); ++it) {
        Link *link = *it;
        pp = get_closest_point_and_primitive(in_link_faces,
          link->second.first);
        link->second.second = pp.first;   // update the closest point
        closest_fd = pp.second;
        get_vertex_in_links(closest_fd).push_back(link);
      }
    }
  }

  void backup_local_links(const std::set<face_descriptor> &extended_faces,
    std::map<face_descriptor, Link_list> *face_out_map,
    std::map<halfedge_descriptor, Link_list> *edge_out_map,
    Link *vertex_out, Link_iter_list *face_in_links,
    Link_iter_list *edge_in_links, Link_pointer_list *vertex_in_links,
    vertex_descriptor vd) const {
    // step 1: backup the out links
    Halfedge_around_target_circulator hb(mesh_.halfedge(vd), mesh_), he(hb);
    do {
      halfedge_descriptor hd = *hb;
      if (!is_border(hd)) {
        face_descriptor fd = get_face(hd);
        const Link_list &face_out_links = get_face_out_links(fd);
        (*face_out_map)[fd].insert((*face_out_map)[fd].end(),
          face_out_links.begin(), face_out_links.end());
      }
      if (get_halfedge_normal_dihedral(hd) == -1.0) {
        hd = get_opposite(hd);
      }
      const Link_list &edge_out_links = get_halfedge_out_links(hd);
      (*edge_out_map)[hd].insert((*edge_out_map)[hd].end(),
        edge_out_links.begin(), edge_out_links.end());
      ++hb;
    } while (hb != he);
    (*vertex_out) = get_vertex_out_link(vd);
    // step 2: backup the in links
    backup_local_in_links(extended_faces, face_in_links, edge_in_links,
      vertex_in_links);
  }

  void backup_local_in_links(const std::set<face_descriptor> &extended_faces,
    Link_iter_list *face_in_links, Link_iter_list *edge_in_links,
    Link_pointer_list *vertex_in_links) const {
    for (auto it = extended_faces.begin(); it != extended_faces.end(); ++it) {
      face_descriptor fd = *it;
      // backup face in links
      const Link_iter_list &fd_face_in_links = get_face_in_links(fd);
      face_in_links->insert(face_in_links->end(),
        fd_face_in_links.begin(), fd_face_in_links.end());
      // backup edge in links
      const Link_iter_list &fd_edge_in_links = get_edge_in_links(fd);
      edge_in_links->insert(edge_in_links->end(),
        fd_edge_in_links.begin(), fd_edge_in_links.end());
      // backup vertex in links
      const Link_pointer_list &fd_vertex_in_links = get_vertex_in_links(fd);
      vertex_in_links->insert(vertex_in_links->end(),
        fd_vertex_in_links.begin(), fd_vertex_in_links.end());
    }
  }

  void restore_local_links(
    const std::map<face_descriptor, Link_list> &face_out_map,
    const std::map<halfedge_descriptor, Link_list> &edge_out_map,
    const Link &vertex_out, const Link_iter_list &face_in_links,
    const Link_iter_list &edge_in_links,
    const Link_pointer_list &vertex_in_links,
    vertex_descriptor vd, const std::set<face_descriptor> &in_link_faces,
    const NamedParameters &np) {
    // only for relocate purpose
    // step 1: update the local feature_intensity around vh
    update_local_feature_intensity(vd, false, np);
    // step 2: clear local_links
    clear_local_links(vd, in_link_faces);
    // step 3: restore local out links
    restore_local_out_links(face_out_map, edge_out_map, vertex_out, vd);
    // step 4: generate local in links
    generate_local_in_links(in_link_faces, face_in_links, edge_in_links,
      vertex_in_links, np.use_local_aabb_tree);
  }

  void restore_local_out_links(
    const std::map<face_descriptor, Link_list> &face_out_map,
    const std::map<halfedge_descriptor, Link_list> &edge_out_map,
    const Link &vertex_out, vertex_descriptor vd) {
    // step 1: restore the face out links
    for (auto it = face_out_map.begin(); it != face_out_map.end(); ++it) {
      face_descriptor fd = it->first;
      Link_list &face_out_links = get_face_out_links(fd);
      face_out_links.insert(face_out_links.end(),
        it->second.begin(), it->second.end());
    }
    // step 2: restore the edge out links
    for (auto it = edge_out_map.begin(); it != edge_out_map.end(); ++it) {
      halfedge_descriptor hd = it->first;
      Link_list &edge_out_links = get_halfedge_out_links(hd);
      edge_out_links.insert(edge_out_links.end(),
        it->second.begin(), it->second.end());
    }
    // step 3: restore the vertex out link
    get_vertex_out_link(vd) = vertex_out;
  }

  // 5) queues

  // 5.1) for eliminating degenerated faces
  void fill_degenerated_faces_queue(FT radian_threshold,
    DPQueue_face_long *queue) const {
    for (typename Mesh::Face_range::const_iterator fi = mesh_.faces().begin();
      fi != mesh_.faces().end(); ++fi) {
      if (calculate_smallest_radian(*fi) < radian_threshold) {
        halfedge_descriptor hd = get_longest_halfedge(*fi);
        queue->insert(Face_long(*fi, squared_length(hd)));
      }
    }
  }

  // 5.2) for splitting long edges
  void fill_queue_with_long_edges(FT max_squared_length,
    DPQueue_halfedge_long *queue) const {
    for (typename Mesh::Edge_range::const_iterator ei = mesh_.edges().begin();
      ei != mesh_.edges().end(); ++ei) {
      halfedge_descriptor hd = mesh_.halfedge(*ei);
      if (is_border(hd)) {
        hd = get_opposite(hd);
      }
      const FT sl = squared_length(hd);
      if (sl >= max_squared_length) {
        queue->insert(Halfedge_long(hd, sl));
      }
    }
  }

  void add_circular_long_edges(FT max_squared_length, vertex_descriptor vd,
    DPQueue_halfedge_long *queue) const {
    Halfedge_around_target_circulator hb(mesh_.halfedge(vd), mesh_), he(hb);
    do {
      halfedge_descriptor hd = *hb;
      if (is_border(hd)) {
        hd = get_opposite(hd);
      }
      const FT sl = squared_length(hd);
      if (sl >= max_squared_length) {
        queue->insert(Halfedge_long(hd, sl));
      }
      ++hb;
    } while (hb != he);
  }

  // 5.3) private queues
  void remove_small_value_edges_before_collapse(halfedge_descriptor hd,
    DPQueue_halfedge_long *large_error_queue,
    DPQueue_halfedge_short *small_value_queue,
    const NamedParameters &np) const {
    std::set<face_descriptor> one_ring_faces;
    collect_one_ring_faces_incident_to_edge(hd, &one_ring_faces);
    remove_small_value_edges(one_ring_faces, large_error_queue,
      small_value_queue, np);
  }

  void add_small_value_edges_after_collapse(vertex_descriptor vd,
    FT max_error_threshold_value, bool reduce_complexity,
    DPQueue_halfedge_long *large_error_queue,
    DPQueue_halfedge_short *small_value_queue,
    const NamedParameters &np) const {
    std::set<face_descriptor> one_ring_faces;
    collect_one_ring_faces_incident_to_vertex(vd, &one_ring_faces);
    add_small_value_edges(one_ring_faces, max_error_threshold_value,
      reduce_complexity, large_error_queue, small_value_queue, np);
  }

  void remove_small_value_edges_before_split(halfedge_descriptor hd,
    DPQueue_halfedge_long *large_error_queue,
    DPQueue_halfedge_short *small_value_queue,
    const NamedParameters &np) const {
    std::set<face_descriptor> one_ring_faces;
    collect_faces_incident_to_edge(hd, &one_ring_faces);
    remove_small_value_edges(one_ring_faces, large_error_queue,
      small_value_queue, np);
  }

  void add_small_value_edges_after_split(vertex_descriptor vd,
    FT max_error_threshold_value, bool reduce_complexity,
    DPQueue_halfedge_long *large_error_queue,
    DPQueue_halfedge_short *small_value_queue,
    const NamedParameters &np) const {
    std::set<face_descriptor> one_ring_faces;
    collect_one_ring_faces_incident_to_vertex(vd, &one_ring_faces);
    add_small_value_edges(one_ring_faces, max_error_threshold_value,
      reduce_complexity, large_error_queue, small_value_queue, np);
  }

  void remove_small_value_edges_before_flip(halfedge_descriptor hd,
    DPQueue_halfedge_long *large_error_queue,
    DPQueue_halfedge_short *small_value_queue,
    const NamedParameters &np) const {
    std::set<face_descriptor> one_ring_faces;
    collect_faces_incident_to_edge(hd, &one_ring_faces);
    remove_small_value_edges(one_ring_faces, large_error_queue,
      small_value_queue, np);
  }

  void add_small_value_edges_after_flip(halfedge_descriptor hd,
    FT max_error_threshold_value, bool reduce_complexity,
    DPQueue_halfedge_long *large_error_queue,
    DPQueue_halfedge_short *small_value_queue,
    const NamedParameters &np) const {
    std::set<face_descriptor> one_ring_faces;
    collect_faces_incident_to_edge(hd, &one_ring_faces);
    add_small_value_edges(one_ring_faces, max_error_threshold_value,
      reduce_complexity, large_error_queue, small_value_queue, np);
  }

  void remove_small_value_edges_before_relocate(vertex_descriptor vd,
    DPQueue_halfedge_long *large_error_queue,
    DPQueue_halfedge_short *small_value_queue,
    const NamedParameters &np) const {
    std::set<face_descriptor> one_ring_faces;
    collect_one_ring_faces_incident_to_vertex(vd, &one_ring_faces);
    remove_small_value_edges(one_ring_faces, large_error_queue,
      small_value_queue, np);
  }

  void add_small_value_edges_after_relocate(vertex_descriptor vd,
    FT max_error_threshold_value, bool reduce_complexity,
    DPQueue_halfedge_long *large_error_queue,
    DPQueue_halfedge_short *small_value_queue,
    const NamedParameters &np) const {
    std::set<face_descriptor> one_ring_faces;
    collect_one_ring_faces_incident_to_vertex(vd, &one_ring_faces);
    add_small_value_edges(one_ring_faces, max_error_threshold_value,
      reduce_complexity, large_error_queue, small_value_queue, np);
  }

  void remove_small_value_edges(
    const std::set<face_descriptor> &one_ring_faces,
    DPQueue_halfedge_long *large_error_queue,
    DPQueue_halfedge_short *small_value_queue,
    const NamedParameters &np) const {
    // small_value_queue may be small_radian_queue or collapse_candidate_queue
    // step 1: remove from small_value_queue
    typename std::set<face_descriptor>::iterator it;
    for (it = one_ring_faces.begin(); it != one_ring_faces.end(); ++it) {
      face_descriptor fd = *it;
      halfedge_descriptor hd = mesh_.halfedge(fd);
      small_value_queue->left.erase(hd);
      small_value_queue->left.erase(mesh_.next(hd));
      small_value_queue->left.erase(mesh_.prev(hd));
    }
    // step 2: remove from large_error_queue if necesary
    if (np.decrease_max_errors) {
      std::set<face_descriptor> extended_faces;
      extend_faces(one_ring_faces, np.stencil_ring_size, &extended_faces);
      for (it = extended_faces.begin(); it != extended_faces.end(); ++it) {
        face_descriptor fd = *it;
        halfedge_descriptor hd = mesh_.halfedge(fd);
        large_error_queue->left.erase(hd);
        large_error_queue->left.erase(mesh_.next(hd));
        large_error_queue->left.erase(mesh_.prev(hd));
      }
    }
  }

  void add_small_value_edges(const std::set<face_descriptor> &one_ring_faces,
    FT max_error_threshold_value, bool reduce_complexity,
    DPQueue_halfedge_long *large_error_queue,
    DPQueue_halfedge_short *small_value_queue,
    const NamedParameters &np) const {
    // small_value_queue may be small_radian_queue or collapse_candidate_queue
    // If is former, the prority is radian; otherwise it is length * radian
    // step 1: add to small_value_queue
    typename std::set<face_descriptor>::iterator it;
    if (reduce_complexity) {    // for min rdian improvement
      for (it = one_ring_faces.begin(); it != one_ring_faces.end(); ++it) {
        face_descriptor fd = *it;
        halfedge_descriptor hd = mesh_.halfedge(fd);
        for (int i = 0; i <= 2; ++i) {
          FT radian = calculate_opposite_radian(hd);
          FT len = length(hd);
          small_value_queue->insert(Halfedge_short(hd, radian * len));
          hd = mesh_.next(hd);
        }
      }
    } else {                    // for mesh complexity reduction
      FT min_radian_threshold = to_radian(np.min_angle_threshold);
      for (it = one_ring_faces.begin(); it != one_ring_faces.end(); ++it) {
        face_descriptor fd = *it;
        halfedge_descriptor hd = mesh_.halfedge(fd);
        for (int i = 0; i <= 2; ++i) {
          FT radian = calculate_opposite_radian(hd);
          if (radian < min_radian_threshold) {
            small_value_queue->insert(Halfedge_short(hd, radian));
          }
          hd = mesh_.next(hd);
        }
      }
    }
    // step 2: add to large_error_queue if necessary
    if (np.decrease_max_errors) {
      FT max_se_threshold = std::pow(max_error_threshold_value, 2);
      std::set<face_descriptor> extended_faces;
      extend_faces(one_ring_faces, np.stencil_ring_size, &extended_faces);
      for (it = extended_faces.begin(); it != extended_faces.end(); ++it) {
        face_descriptor fd = *it;
        FT max_se = get_face_max_squared_error(fd);
        if (max_se >= max_se_threshold) {
          halfedge_descriptor longest_hd = get_longest_halfedge(fd);
          large_error_queue->insert(Halfedge_long(longest_hd, max_se));
        }
      }
    }
  }

  void remove_incident_faces(halfedge_descriptor hd,
    DPQueue_face_long *queue) const {
    vertex_descriptor vp = get_source_vertex(hd);
    vertex_descriptor vq = get_target_vertex(hd);
    Halfedge_around_target_circulator hb1(mesh_.halfedge(vp), mesh_), he1(hb1);
    do {
      if (!is_border(*hb1)) {
        face_descriptor fd = get_face(*hb1);
        queue->left.erase(fd);
      }
      ++hb1;
    } while (hb1 != he1);
    Halfedge_around_target_circulator hb2(mesh_.halfedge(vq), mesh_), he2(hb2);
    do {
      if (!is_border(*hb2)) {
        face_descriptor fd = get_face(*hb2);
        queue->left.erase(fd);
      }
      ++hb2;
    } while (hb2 != he2);
  }

  void add_circulator_degenerate_faces(vertex_descriptor vd,
    FT radian_threshold, DPQueue_face_long *queue) const {
    Halfedge_around_target_circulator hb(mesh_.halfedge(vd), mesh_), he(hb);
    FT longest_squared_length = 0.0;
    do {
      halfedge_descriptor hd = *hb;
      if (!is_border(hd)) {
        face_descriptor fd = get_face(hd);
        halfedge_descriptor shortest_hd = get_shortest_halfedge(fd);
        FT smallest_radian = calculate_opposite_radian(shortest_hd);
        if (smallest_radian < radian_threshold) {
          // degenerate edge (we want to make degenerate edge first)
          if (squared_length(shortest_hd) < SQUARED_MIN_VALUE) {
            queue->insert(Face_long(fd, MAX_VALUE));
          } else {    // degenerate face
            halfedge_descriptor longest_hd = get_longest_halfedge(fd);
            queue->insert(Face_long(fd, squared_length(longest_hd)));
          }
        }
      }
      ++hb;
    } while (hb != he);
  }

  // 6) local operations

  // 6.1) split
  halfedge_descriptor split_long_edge(const Point &new_point,
    halfedge_descriptor hd) {
    // 1) apply the local operator(s)
    // 2) update property maps of new generated faces, halfedges and vertices
    halfedge_descriptor hnew = CGAL::Euler::split_edge(hd, mesh_);
    vertex_descriptor vd = get_target_vertex(hnew);
    halfedge_descriptor h = get_null_halfedge();
    get_point(vd) = new_point;
    reset_halfedge_properties(hnew, hd);
    reset_halfedge_properties(get_opposite(hnew), get_opposite(hd));
    reset_vertex_properties(vd, get_null_vertex());
    if (!is_border(hnew)) {
      h = CGAL::Euler::split_face(hnew, mesh_.next(hd), mesh_);
      reset_halfedge_properties(h, get_null_halfedge());
      reset_halfedge_properties(get_opposite(h), get_null_halfedge());
      reset_face_properties(get_face(h), get_null_face());
      reset_face_properties(get_face(get_opposite(h)), get_null_face());
    }
    if (!is_border(get_opposite(hnew))) {
      h = CGAL::Euler::split_face(mesh_.next(get_opposite(hnew)),
        get_opposite(hd), mesh_);
      reset_halfedge_properties(h, get_null_halfedge());
      reset_halfedge_properties(get_opposite(h), get_null_halfedge());
      reset_face_properties(get_face(h), get_null_face());
      reset_face_properties(get_face(get_opposite(h)), get_null_face());
    }
    return hnew;
  }

  // 6.2) collapse
  vertex_descriptor collapse_short_edge(const Point &new_point,
    halfedge_descriptor hd) {
    // 1) apply the local operator(s)
    // 2) update property maps of new generated faces, halfedges and vertices
    edge_descriptor ed = mesh_.edge(hd);
    vertex_descriptor vd = CGAL::Euler::collapse_edge(ed, mesh_);
    get_point(vd) = new_point;
    // since no new elements added, we only reset the sample links
    Halfedge_around_target_circulator hb(mesh_.halfedge(vd), mesh_), he(hb);
    do {
      if (!is_border(*hb)) {
        face_descriptor fd = get_face(*hb);
        clear_face_links(fd);
      }
      clear_halfedge_links(*hb);
      clear_halfedge_links(get_opposite(*hb));
      ++hb;
    } while (hb != he);
    clear_vertex_links(vd);
    return vd;
  }

  bool check_link_condition(halfedge_descriptor hd) const {
    // h = (p -> q)
    vertex_descriptor vp = get_source_vertex(hd);
    vertex_descriptor vq = get_target_vertex(hd);
    if (is_border(get_opposite(hd))) {
      // (p, q, r)
      vertex_descriptor vr = get_opposite_vertex(hd);
      Halfedge_around_target_circulator hb(mesh_.halfedge(vp), mesh_), he(hb);
      do {
        vertex_descriptor v = get_source_vertex(*hb);
        if (v == vq || v == vr) {
          ++hb;
          continue;
        }
        if (are_neighbors(v, vq)) {
          return false;
        }
        ++hb;
      } while (hb != he);
      return true;
    }
    // (p, q, r)
    vertex_descriptor vr = get_opposite_vertex(hd);
    // (q, p, s)
    vertex_descriptor vs = get_opposite_vertex(get_opposite(hd));
    Halfedge_around_target_circulator hb(mesh_.halfedge(vp), mesh_), he(hb);
    do {
      vertex_descriptor v = get_source_vertex(*hb);
      if (v == vq || v == vr || v == vs) {
        ++hb;
        continue;
      }
      if (are_neighbors(v, vq)) {
        return false;
      }
      ++hb;
    } while (hb != he);
    return true;
  }

  // 6.3) flip
  bool flip_would_cause_wrinkle(halfedge_descriptor hd) const {
    const Point &p = get_point(get_source_vertex(hd));
    const Point &q = get_point(get_target_vertex(hd));
    const Point &s = get_point(get_opposite_vertex(hd));
    const Point &t = get_point(get_opposite_vertex(get_opposite(hd)));
    Plane plane(t, s, p);
    Point projection = plane.projection(q);
    if (same_side(t, s, p, projection) >= 0) {  // on the same side
      return true;
    } else {
      return false;
    }
  }

  bool is_flippable(halfedge_descriptor hd, const NamedParameters &np) const {
    // step 1: check whether it is border edge or the new_edge already exists
    if (is_border(hd) || is_border(get_opposite(hd))) {
      return false;
    }
    vertex_descriptor p = get_target_vertex(hd);
    vertex_descriptor q = get_source_vertex(hd);
    vertex_descriptor r = get_opposite_vertex(hd);
    vertex_descriptor s = get_opposite_vertex(get_opposite(hd));
    if (are_neighbors(r, s)) {
      return false;
    }
    // step 2: check whether we are flipping a crease edge
    if (np.inherit_element_types && is_crease_edge(hd)) {
      return false;
    }
    // step 3: check whether flip would cause wrinkle
    if (flip_would_cause_wrinkle(hd)) {
      return false;
    }
    // step 4: check whether flipping improves the valence or angle
    if (np.edge_flip_strategy == EdgeFlipStrategy::k_improve_valence) {
      return edge_flip_would_improve_valence(hd);
    } else {
      return edge_flip_would_improve_radian(hd);
    }
  }

  bool edge_flip_would_improve_valence(halfedge_descriptor hd) const {
    // incident vertices
    vertex_descriptor a = get_target_vertex(hd);
    vertex_descriptor b = get_opposite_vertex(hd);
    vertex_descriptor c = get_source_vertex(hd);
    vertex_descriptor d = get_opposite_vertex(get_opposite(hd));
    // target valences
    int tar_var_a = is_on_boundary(a) ? 4 : 6;
    int tar_var_b = is_on_boundary(b) ? 4 : 6;
    int tar_var_c = is_on_boundary(c) ? 4 : 6;
    int tar_var_d = is_on_boundary(d) ? 4 : 6;
    // real valence
    int var_a = static_cast<int>(mesh_.degree(a));
    int var_b = static_cast<int>(mesh_.degree(b));
    int var_c = static_cast<int>(mesh_.degree(c));
    int var_d = static_cast<int>(mesh_.degree(d));
    int deviation_pre = std::abs(var_a - tar_var_a) +
      std::abs(var_b - tar_var_b) +
      std::abs(var_c - tar_var_c) +
      std::abs(var_d - tar_var_d);
    // after flipping, the valences of a, b, c and d change
    --var_a, --var_c, ++var_b, ++var_d;
    int deviation_post = std::abs(var_a - tar_var_a) +
      std::abs(var_b - tar_var_b) +
      std::abs(var_c - tar_var_c) +
      std::abs(var_d - tar_var_d);
    return deviation_pre > deviation_post;
  }

  bool edge_flip_would_improve_radian(halfedge_descriptor hd) const {
    CGAL_precondition(!is_border(hd) && !is_border(get_opposite(hd)));
    const Point &p = get_point(get_source_vertex(hd));
    const Point &q = get_point(get_target_vertex(hd));
    const Point &s = get_point(get_opposite_vertex(hd));
    const Point &t = get_point(get_opposite_vertex(get_opposite(hd)));
    FT min_radian_before = calculate_smallest_radian(p, q, s);
    min_radian_before = CGAL::min(min_radian_before,
      calculate_smallest_radian(q, p, t));
    FT min_radian_after = calculate_smallest_radian(t, s, p);
    min_radian_after = CGAL::min(min_radian_after,
      calculate_smallest_radian(s, t, q));
    return min_radian_before < min_radian_after;
  }

  void flip_inner_edge(halfedge_descriptor hd) {
    // 1) apply the local operator(s)
    // 2) update property maps of new generated faces, halfedges and vertices
    CGAL::Euler::flip_edge(hd, mesh_);
    // since no element added, we only clear the links
    clear_halfedge_links(hd);
    clear_halfedge_links(get_opposite(hd));
    clear_face_links(get_face(hd));
    clear_face_links(get_face(get_opposite(hd)));
  }

  int flip_edges(const Face_tree &input_face_tree,
    FT max_error_threshold_value, FT max_error, FT min_radian,
    bool reduce_complexity, DPQueue_halfedge_long *large_error_queue,
    DPQueue_halfedge_short *small_radian_queue, Halfedge_list *halfedges,
    const NamedParameters &np) {
    int num = 0;
    for (Halfedge_iter it = halfedges->begin(); it != halfedges->end(); ++it) {
      halfedge_descriptor hd = *it;
      halfedge_descriptor hnew = flip_edge(input_face_tree,
        max_error_threshold_value, max_error, min_radian, reduce_complexity,
        large_error_queue, small_radian_queue, hd, np);
      if (hnew != get_null_halfedge()) {
        ++num;
      }
    }
    return num;
  }

  // 6.4) relocate
  bool relocate_would_cause_wrinkle(const Point &new_point,
    vertex_descriptor vd) const {
    Halfedge_around_target_circulator hb(mesh_.halfedge(vd), mesh_), he(hb);
    do {
      if (!is_border(*hb)) {
        const Point &start_point = get_point(get_opposite_vertex(*hb));
        const Point &end_point = get_point(get_source_vertex(*hb));
        const Point &old_point = get_point(get_target_vertex(*hb));
        FT radian = calculate_radian(end_point, old_point, start_point);
        if (radian < MIN_VALUE || CGAL_PI - radian < MIN_VALUE) {
          return true;    // degenerate case
        }
        Plane plane(end_point, old_point, start_point);
        Point projection = plane.projection(new_point);
        if (same_side(start_point, end_point, old_point, projection) < 0) {
          return true;
        }
      }
      ++hb;
    } while (hb != he);
    return false;
  }

  void relocate_vertex_point(vertex_descriptor vd, const Point &new_point) {
    // 1) apply the local operator(s)
    // 2) update property maps of new generated faces, halfedges and vertices
    get_point(vd) = new_point;
    // since no new element added, we only clear the links
    Halfedge_around_target_circulator hb(mesh_.halfedge(vd), mesh_), he(hb);
    do {
      if (!is_border(*hb)) {
        face_descriptor fd = get_face(*hb);
        clear_face_links(fd);
      }
      clear_halfedge_links(*hb);
      clear_halfedge_links(get_opposite(*hb));
      ++hb;
    } while (hb != he);
    clear_vertex_links(vd);
  }

  int relocate_vertices(const Face_tree &input_face_tree,
    FT max_error_threshold_value, FT max_error, FT min_radian,
    bool reduce_complexity, DPQueue_halfedge_long *large_error_queue,
    DPQueue_halfedge_short *small_radian_queue,
    const std::set<vertex_descriptor> &vertices, const NamedParameters &np) {
    // step 1: construct the vertex_map to make bigger distance relocate first
    std::map<FT, vertex_descriptor> vertex_map;
    for (auto it = vertices.begin(); it != vertices.end(); ++it) {
      vertex_descriptor vd = *it;
      Point p = calculate_initial_point_for_relocate(input_face_tree, vd, np);
      vertex_map[CGAL::squared_distance(get_point(vd), p)] = vd;
    }
    // step 2: relocate these vertices in order
    int num = 0;
    for (auto it = vertex_map.rbegin(); it != vertex_map.rend(); ++it) {
      vertex_descriptor vd = it->second;
      Point initial_point = calculate_initial_point_for_relocate(
        input_face_tree, vd, np);
      num += relocate_vertex(input_face_tree, max_error_threshold_value,
        max_error, min_radian, reduce_complexity, large_error_queue,
        small_radian_queue, initial_point, vd, np);
    }
    return num;
  }

  // 6.5) optimize
  Vector calculate_optimize_vector(bool use_face_in_links,
    bool use_face_out_links, bool use_edge_in_links, bool use_edge_out_links,
    bool use_vertex_in_links, bool use_vertex_out_links,
    bool use_feature_intensity_weights, vertex_descriptor vd,
    FT *denominator) const {
    Vector sum_vec = CGAL::NULL_VECTOR;
    FT weight = 0.0;
    Halfedge_around_target_circulator hb(mesh_.halfedge(vd), mesh_), he(hb);
    do {
      // vertex out links
      if (use_vertex_out_links) {
        const Point &out_start = get_vertex_out_link(vd).second.first;
        const Point &out_end = get_vertex_out_link(vd).second.second;
        weight = CGAL::sqrt(CGAL::squared_distance(out_start, out_end));
        if (use_feature_intensity_weights) {
          weight *= get_vertex_out_link(vd).first;
        }
        sum_vec = sum_vec + weight * (out_end - CGAL::ORIGIN);
        *denominator += weight;
      }
      // edge out links
      if (use_edge_out_links) {
        halfedge_descriptor hd = *hb;
        const Point &v0 = get_point(get_target_vertex(hd));
        const Point &v1 = get_point(get_source_vertex(hd));
        if (get_halfedge_normal_dihedral(hd) == -1.0) {
          hd = get_opposite(hd);
        }
        FT len = length(hd);
        if (len >= MIN_VALUE) {    // we ignore too short edges
          const Link_list &edge_out_links = get_halfedge_out_links(hd);
          for (auto it = edge_out_links.begin();
            it != edge_out_links.end(); ++it) {
            const Point &out_start = it->second.first;
            const Point &out_end = it->second.second;
            weight = CGAL::sqrt(CGAL::squared_distance(out_start, out_end));
            if (use_feature_intensity_weights) {
              weight *= it->first;
            }
            FT alpha =
              CGAL::sqrt(CGAL::squared_distance(out_start, v1)) / len;
            sum_vec = sum_vec + weight * alpha * (out_end - CGAL::ORIGIN) -
              weight * alpha * (1.0 - alpha) * (v1 - CGAL::ORIGIN);
            *denominator += weight * alpha * alpha;
          }
        }
      }
      halfedge_descriptor hd = *hb;
      if (is_border(hd)) {
        ++hb;
        continue;
      }
      // face out links
      face_descriptor fd = get_face(hd);
      const FT face_area = area(fd);
      if (face_area < SQUARED_MIN_VALUE) {  // we ignore too small faces
        ++hb;
        continue;
      }
      const Point &v0 = get_point(get_target_vertex(hd));
      const Point &v1 = get_point(get_opposite_vertex(hd));
      const Point &v2 = get_point(get_source_vertex(hd));
      if (use_face_out_links) {
        const Link_list &face_out_links = get_face_out_links(fd);
        for (auto it = face_out_links.begin();
          it != face_out_links.end(); ++it) {
          const Point &out_start = it->second.first;
          const Point &out_end = it->second.second;
          weight = CGAL::sqrt(CGAL::squared_distance(out_start, out_end));
          if (use_feature_intensity_weights) {
            weight *= it->first;
          }
          FT alpha_0 = area(v1, v2, out_start) / face_area;
          FT alpha_1 = area(v2, v0, out_start) / face_area;
          FT alpha_2 = area(v0, v1, out_start) / face_area;
          sum_vec = sum_vec + weight * alpha_0 * (out_end - CGAL::ORIGIN)
            - weight * alpha_0 * alpha_1 * (v1 - CGAL::ORIGIN)
            - weight * alpha_0 * alpha_2 * (v2 - CGAL::ORIGIN);
          *denominator += weight * alpha_0 * alpha_0;
        }
      }
      // vertex in links
      if (use_vertex_in_links) {
        const Link_pointer_list &vertex_in_links = get_vertex_in_links(fd);
        for (auto it = vertex_in_links.begin();
          it != vertex_in_links.end(); ++it) {
          Link *link = *it;
          const Point &in_start = link->second.first;
          const Point &in_end = link->second.second;
          weight = CGAL::sqrt(CGAL::squared_distance(in_start, in_end));
          if (use_feature_intensity_weights) {
            weight *= link->first;
          }
          FT alpha_0 = area(v1, v2, in_end) / face_area;
          FT alpha_1 = area(v2, v0, in_end) / face_area;
          FT alpha_2 = area(v0, v1, in_end) / face_area;
          sum_vec = sum_vec + weight * alpha_0 * (in_start - CGAL::ORIGIN)
            - weight * alpha_0 * alpha_1 * (v1 - CGAL::ORIGIN)
            - weight * alpha_0 * alpha_2 * (v2 - CGAL::ORIGIN);
          *denominator += weight * alpha_0 * alpha_0;
        }
      }
      // edge in links
      if (use_edge_in_links) {
        const Link_iter_list &edge_in_links = get_edge_in_links(fd);
        for (auto it = edge_in_links.begin();
          it != edge_in_links.end(); ++it) {
          Link_list_const_iter lit = *it;
          const Point &in_start = lit->second.first;
          const Point &in_end = lit->second.second;
          weight = CGAL::sqrt(CGAL::squared_distance(in_start, in_end));
          if (use_feature_intensity_weights) {
            weight *= lit->first;
          }
          FT alpha_0 = area(v1, v2, in_end) / face_area;
          FT alpha_1 = area(v2, v0, in_end) / face_area;
          FT alpha_2 = area(v0, v1, in_end) / face_area;
          sum_vec = sum_vec + weight * alpha_0 * (in_start - CGAL::ORIGIN)
            - weight * alpha_0 * alpha_1 * (v1 - CGAL::ORIGIN)
            - weight * alpha_0 * alpha_2 * (v2 - CGAL::ORIGIN);
          *denominator += weight * alpha_0 * alpha_0;
        }
      }
      // face in links
      if (use_face_in_links) {
        const Link_iter_list &face_in_links = get_face_in_links(fd);
        for (auto it = face_in_links.begin();
          it != face_in_links.end(); ++it) {
          Link_list_const_iter lit = *it;
          const Point &in_start = lit->second.first;
          const Point &in_end = lit->second.second;
          weight = CGAL::sqrt(CGAL::squared_distance(in_start, in_end));
          if (use_feature_intensity_weights) {
            weight *= lit->first;
          }
          FT alpha_0 = area(v1, v2, in_end) / face_area;
          FT alpha_1 = area(v2, v0, in_end) / face_area;
          FT alpha_2 = area(v0, v1, in_end) / face_area;
          sum_vec = sum_vec + weight * alpha_0 * (in_start - CGAL::ORIGIN)
            - weight * alpha_0 * alpha_1 * (v1 - CGAL::ORIGIN)
            - weight * alpha_0 * alpha_2 * (v2 - CGAL::ORIGIN);
          *denominator += weight * alpha_0 * alpha_0;
        }
      }
      ++hb;
    } while (hb != he);
    return sum_vec;
  }

  Point calculate_optimized_position(vertex_descriptor vd,
    const NamedParameters &np) const {
    FT denominator = 0.0;
    Vector vec_sum = CGAL::NULL_VECTOR;
    bool use_face_in_links = false, use_face_out_links = false;
    bool use_edge_in_links = false, use_edge_out_links = false;
    bool use_vertex_in_links = false, use_vertex_out_links = false;
    // step 1: collect the face samples
    use_face_in_links =
      np.face_optimize_type == OptimizeType::k_input_to_remesh ||
      np.face_optimize_type == OptimizeType::k_both;
    use_face_out_links =
      np.face_optimize_type == OptimizeType::k_remesh_to_input ||
      np.face_optimize_type == OptimizeType::k_both;
    // step 2: collect the edge samples
    use_edge_in_links =
      np.edge_optimize_type == OptimizeType::k_input_to_remesh ||
      np.edge_optimize_type == OptimizeType::k_both;
    use_edge_out_links =
      np.edge_optimize_type == OptimizeType::k_remesh_to_input ||
      np.edge_optimize_type == OptimizeType::k_both;
    // step 3: collect the vertex samples
    use_vertex_in_links =
      np.vertex_optimize_type == OptimizeType::k_input_to_remesh ||
      np.vertex_optimize_type == OptimizeType::k_both;
    use_vertex_out_links =
      np.vertex_optimize_type == OptimizeType::k_remesh_to_input ||
      np.vertex_optimize_type == OptimizeType::k_both;
    // step 4: get the optimized vector
    vec_sum = calculate_optimize_vector(use_face_in_links, use_face_out_links,
      use_edge_in_links, use_edge_out_links, use_vertex_in_links,
      use_vertex_out_links, np.use_feature_intensity_weights, vd,
      &denominator);
    // step 5: get the vertex position if valid
    const Point &p = get_point(vd);
    if (denominator > SQUARED_MIN_VALUE) {
      vec_sum = vec_sum / denominator;
      Point new_point = CGAL::ORIGIN + vec_sum;
      return p + np.vertex_optimize_ratio * (new_point - p);
    } else {
      return p;
    }
  }

  bool optimize_vertex_by_one_ring(const Face_tree &input_face_tree,
    vertex_descriptor vd, const NamedParameters &np) {
    // step 1: check whether it is too small to optimize
    size_t nb_out_links = 0, nb_in_links = 0;
    Halfedge_around_target_circulator hb(mesh_.halfedge(vd), mesh_), he(hb);
    do {
      if (!is_border(*hb)) {
        face_descriptor fd = get_face(*hb);
        nb_out_links += get_face_out_links(fd).size();
        nb_in_links += get_face_in_links(fd).size();
        nb_in_links += get_edge_in_links(fd).size();
        nb_in_links += get_vertex_in_links(fd).size();
      }
      ++hb;
    } while (hb != he);
    Point &p = get_point(vd);
    if (nb_out_links == 0 || nb_in_links == 0) {
      p = barycenter(vd, np);
      return false;   // the area is too small
    }
    // step 2: calculate the optimized position
    Point new_point = calculate_optimized_position(vd, np);
    if (CGAL::squared_distance(new_point, p) < SQUARED_MIN_VALUE) {
      return false;   // the optimization is too little
    }
    if (np.keep_vertex_in_one_ring) {
      FT min_sd = calculate_min_squared_distance_in_one_ring_faces(vd);
      FT sd = CGAL::squared_distance(new_point, p);
      if (sd > min_sd) {
        p = p + (new_point - p) * CGAL::sqrt(min_sd) / CGAL::sqrt(sd);
      }
    } else {
      p = new_point;
    }
    if (np.optimize_strategy == OptimizeStrategy::k_Interpolation) {
      p = input_face_tree.closest_point(p);
    }
    return true;
  }

  // 7) max errors
  void reset_face_max_squared_errors(FT value) {
    for (typename Mesh::Face_range::const_iterator fi = mesh_.faces().begin();
      fi != mesh_.faces().end(); ++fi) {
      set_face_max_squared_error(*fi, value);
    }
  }

  void update_face_in_max_squared_errors() {
    FT se = 0.0;    // squared error
    for (typename Mesh::Face_range::const_iterator fi = mesh_.faces().begin();
      fi != mesh_.faces().end(); ++fi) {
      Link_iter_list &face_in_links = get_face_in_links(*fi);
      for (Link_iter_list_iter it = face_in_links.begin();
        it != face_in_links.end(); ++it) {
        Link_list_iter llit = *it;
        const Link &link = *llit;
        se = CGAL::squared_distance(link.second.first, link.second.second);
        set_face_max_squared_error(*fi, CGAL::max(get_face_max_error(*fi), se));
      }
    }
  }

  void update_face_out_max_squared_errors() {
    FT se = 0.0;    // squared error
    for (typename Mesh::Face_range::const_iterator fi = mesh_.faces().begin();
      fi != mesh_.faces().end(); ++fi) {
      Link_list &face_out_links = get_face_out_links(*fi);
      for (Link_list_iter it = face_out_links.begin();
        it != face_out_links.end(); ++it) {
        const Link &link = *it;
        se = CGAL::squared_distance(link.second.first, link.second.second);
        set_face_max_squared_error(*fi,
                                   CGAL::max(get_face_max_error(*fi), se));
      }
    }
  }

  void update_edge_in_max_squared_errors() {
    FT max_se = 0, se = 0.0;
    for (typename Mesh::Face_range::const_iterator fi = mesh_.faces().begin();
      fi != mesh_.faces().end(); ++fi) {
      Link_iter_list &edge_in_links = get_edge_in_links(*fi);
      for (Link_iter_list_iter it = edge_in_links.begin();
        it != edge_in_links.end(); ++it) {
        Link_list_iter llit = *it;
        const Link &link = *llit;
        se = CGAL::squared_distance(link.second.first, link.second.second);
        max_se = CGAL::max(get_face_max_squared_error(*fi), se);
        set_face_max_squared_error(*fi, max_se);
      }
    }
  }

  void update_edge_out_max_squared_errors() {
    FT se = 0.0;
    for (typename Mesh::Edge_range::const_iterator ei = mesh_.edges().begin();
      ei != mesh_.edges().end(); ++ei) {
      halfedge_descriptor hd = mesh_.halfedge(*ei);
      if (get_halfedge_normal_dihedral(hd) == -1.0) {
        hd = get_opposite(hd);
      }
      FT max_se = 0.0;
      Link_list &edge_out_links = get_halfedge_out_links(hd);
      for (Link_list_iter it = edge_out_links.begin();
        it != edge_out_links.end(); ++it) {
        const Link &link = *it;
        se = CGAL::squared_distance(link.second.first, link.second.second);
        max_se = CGAL::max(max_se, se);
      }
      face_descriptor fd = get_face(hd);
      set_face_max_squared_error(fd,
          CGAL::max(get_face_max_squared_error(fd), max_se));
      if (!is_border(get_opposite(hd))) {
        fd = get_face(get_opposite(hd));
        set_face_max_squared_error(fd,
            CGAL::max(get_face_max_squared_error(fd), max_se));
      }
    }
  }

  void update_vertex_in_max_squared_errors() {
    FT se = 0.0;
    for (typename Mesh::Face_range::const_iterator fi = mesh_.faces().begin();
      fi != mesh_.faces().end(); ++fi) {
      Link_pointer_list &vertex_in_links = get_vertex_in_links(*fi);
      for (Link_pointer_iter it = vertex_in_links.begin();
        it != vertex_in_links.end(); ++it) {
        Link *link = *it;
        se = CGAL::squared_distance(link->second.first, link->second.second);
        set_face_max_squared_error(*fi,
            CGAL::max(get_face_max_squared_error(*fi), se));
      }
    }
  }

  void update_vertex_out_max_squared_errors() {
    FT se = 0.0;
    for (typename Mesh::Vertex_range::const_iterator vi = mesh_.vertices().begin();
      vi != mesh_.vertices().end(); ++vi) {
      const Link &link = get_vertex_out_link(*vi);
      se = CGAL::squared_distance(link.second.first, link.second.second);
      Face_around_target_circulator fb(mesh_.halfedge(*vi), mesh_), fe(fb);
      do {
        if (*fb != get_null_face()) {
          set_face_max_squared_error(*fb,
              CGAL::max(get_face_max_squared_error(*fb), se));
        }
        ++fb;
      } while (fb != fe);
    }
  }

  // 8) collections

  // vertex_descriptor collection
  void collect_vertices(const std::set<face_descriptor> &faces,
    std::set<vertex_descriptor> *vertices) const {
    for (auto it = faces.begin(); it != faces.end(); ++it) {
      halfedge_descriptor hd = mesh_.halfedge(*it);
      vertices->insert(get_target_vertex(hd));
      vertices->insert(get_opposite_vertex(hd));
      vertices->insert(get_source_vertex(hd));
    }
  }

  void collect_incident_vertices(vertex_descriptor vd,
    std::set<vertex_descriptor> *vertices) const {
    Halfedge_around_target_circulator hb(mesh_.halfedge(vd), mesh_), he(hb);
    do {
      vertices->insert(get_source_vertex(*hb));
      ++hb;
    } while (hb != he);
  }

  // halfedge_descriptor collection
  void collect_edges_incident_to_vertex(vertex_descriptor vd,
    Edge_list *edges) const {
    Halfedge_around_target_circulator hb(mesh_.halfedge(vd), mesh_), he(hb);
    do {
      halfedge_descriptor hd = *hb;
      edges->push_back(mesh_.edge(hd));
      ++hb;
    } while (hb != he);
  }

  void collect_crease_edges_in_one_ring(vertex_descriptor vd,
    Halfedge_list *crease_edges) const {
    Halfedge_around_target_circulator hb(mesh_.halfedge(vd), mesh_), he(hb);
    do {
      halfedge_descriptor hd = *hb;
      if (get_halfedge_normal_dihedral(hd) == -1.0) {
        hd = get_opposite(hd);
      }
      if (get_halfedge_is_crease(hd) == true) {
        crease_edges->push_back(*hb);
      }
      ++hb;
    } while (hb != he);
  }

  void collect_effective_edges_in_one_ring(FT feature_control_delta,
    vertex_descriptor vd, Halfedge_list *effective_edges) const {
    FT f1 = calculate_feature_intensity(vd) * feature_control_delta;
    FT f2 = (get_vertex_max_dihedral(vd) + 1) * feature_control_delta;
    Halfedge_around_target_circulator hb(mesh_.halfedge(vd), mesh_), he(hb);
    do {
      vertex_descriptor vs = get_source_vertex(*hb);
      bool b1 = calculate_feature_intensity(vs) >= f1;
      halfedge_descriptor hd = *hb;
      if (get_halfedge_normal_dihedral(hd) == -1.0) {
        hd = get_opposite(hd);
      }
      bool b2 = (get_halfedge_normal_dihedral(hd) + 1) >= f2;
      if (b1 && b2) {
        effective_edges->push_back(hd);
      }
      ++hb;
    } while (hb != he);
  }

  // face_descriptor collection
  void collect_one_ring_faces_incident_to_vertex(vertex_descriptor vd,
    std::set<face_descriptor> *faces) const {
    Halfedge_around_target_circulator hb(mesh_.halfedge(vd), mesh_), he(hb);
    do {
      if (!is_border(*hb)) {
        faces->insert(get_face(*hb));
      }
      ++hb;
    } while (hb != he);
  }

  void collect_faces_incident_to_edge(halfedge_descriptor hd,
    Face_list *faces) const {
    if (!is_border(hd)) {
      faces->push_back(get_face(hd));
    }
    if (!is_border(get_opposite(hd))) {
      faces->push_back(get_face(get_opposite(hd)));
    }
  }

  void collect_faces_incident_to_edge(halfedge_descriptor hd,
    std::set<face_descriptor> *faces) const {
    if (!is_border(hd)) {
      faces->insert(get_face(hd));
    }
    if (!is_border(get_opposite(hd))) {
      faces->insert(get_face(get_opposite(hd)));
    }
  }

  void collect_faces_incident_to_vertex(vertex_descriptor vd,
    Face_list *faces) const {
    Halfedge_around_target_circulator hb(mesh_.halfedge(vd), mesh_), he(hb);
    do {
      if (!is_border(*hb)) {
        faces->push_back(get_face(*hb));
      }
      ++hb;
    } while (hb != he);
  }

  void extend_faces_by_one_ring(std::set<face_descriptor> *faces) const {
    std::set<vertex_descriptor> vertices;
    collect_vertices(*faces, &vertices);
    for (auto it = vertices.begin(); it != vertices.end(); ++it) {
      collect_one_ring_faces_incident_to_vertex(*it, faces);
    }
  }

  void collect_incident_faces(face_descriptor fd,
    std::set<face_descriptor> *faces) const {
    halfedge_descriptor hd = mesh_.halfedge(fd);
    // the first vertex
    vertex_descriptor vd1 = get_target_vertex(hd);
    Halfedge_around_target_circulator hb1(mesh_.halfedge(vd1), mesh_), he1(hb1);
    do {
      if (!is_border(*hb1)) {
        faces->insert(get_face(*hb1));
      }
      ++hb1;
    } while (hb1 != he1);
    // the second vertex
    vertex_descriptor vd2 = get_opposite_vertex(hd);
    Halfedge_around_target_circulator hb2(mesh_.halfedge(vd2), mesh_), he2(hb2);
    do {
      if (!is_border(*hb2)) {
        faces->insert(get_face(*hb2));
      }
      ++hb2;
    } while (hb2 != he2);
    // the third vertex
    vertex_descriptor vd3 = get_source_vertex(hd);
    Halfedge_around_target_circulator hb3(mesh_.halfedge(vd3), mesh_), he3(hb3);
    do {
      if (!is_border(*hb3)) {
        faces->insert(get_face(*hb3));
      }
      ++hb3;
    } while (hb3 != he3);
  }

  // 9) feature intensities
  void update_local_feature_intensity(vertex_descriptor vd,
    bool reset_normal_dihedral, const NamedParameters &np) {
    Vertex_list incident_vertices;
    incident_vertices.push_back(vd);
    // step 1: calculate the edge feature intensities around it
    Halfedge_around_target_circulator hb(mesh_.halfedge(vd), mesh_), he(hb);
    do {
      edge_descriptor ed = mesh_.edge(*hb);
      calculate_edge_feature_intensity(ed, reset_normal_dihedral, np);
      incident_vertices.push_back(get_source_vertex(*hb));
      ++hb;
    } while (hb != he);
    // step 2: calculate the vertex feature intensity
    for (Vertex_const_iter cit = incident_vertices.begin();
      cit != incident_vertices.end(); ++cit) {
      calculate_vertex_feature_intensity(*cit, np);
    }
  }

  void calculate_vertex_feature_intensities(const NamedParameters &np) {
    // precondition: the halfedges' normal dihedrals have been updated
    for (typename Mesh::Vertex_range::const_iterator vi = mesh_.vertices().begin();
      vi != mesh_.vertices().end(); ++vi) {
      calculate_vertex_feature_intensity(*vi, np);
    }
  }

  void calculate_vertex_feature_intensity(vertex_descriptor vd,
    const NamedParameters &np) {
    bool on_boundary = false;
    FT sum_radian = 0.0;
    FT gaussian_curvature = 0.0;
    FT max_halfedge_dihedral = 0.0;
    halfedge_descriptor hd;
    Halfedge_around_target_circulator hb(mesh_.halfedge(vd), mesh_), he(hb);
    do {
      if (is_border(*hb)) {
        on_boundary = true;
        hd = get_opposite(*hb);
        max_halfedge_dihedral = CGAL::max(max_halfedge_dihedral,
          get_halfedge_normal_dihedral(hd));
      } else {
        hd = *hb;
        sum_radian += calculate_radian(get_point(get_source_vertex(hd)),
          get_point(get_target_vertex(hd)),
          get_point(get_opposite_vertex(hd)));
        if (get_halfedge_normal_dihedral(hd) == -1.0) {
          hd = get_opposite(hd);
        }
        max_halfedge_dihedral = CGAL::max(max_halfedge_dihedral,
          get_halfedge_normal_dihedral(hd));
      }
      ++hb;
    } while (hb != he);
    if (on_boundary) {
      gaussian_curvature = CGAL::abs(CGAL_PI - sum_radian);
    } else {
      gaussian_curvature = CGAL::abs(2 * CGAL_PI - sum_radian);
    }
    gaussian_curvature /= np.sum_delta;
    gaussian_curvature = CGAL::min(gaussian_curvature, np.sum_theta * CGAL_PI);
    set_vertex_gaussian_curvature(vd, gaussian_curvature);
    set_vertex_max_dihedral(vd, max_halfedge_dihedral);
  }

  void calculate_edge_feature_intensities(const NamedParameters &np) {
    for (typename Mesh::Edge_range::const_iterator ei = mesh_.edges().begin();
      ei != mesh_.edges().end(); ++ei) {
      calculate_edge_feature_intensity(*ei, false, np);
    }
  }

  void calculate_edge_feature_intensity(edge_descriptor ed,
    bool reset_normal_dihedral, const NamedParameters &np) {
    // 1) normalized diheral = min(diheral_theta, dihedral / dihedral_delta);
    // 2) default value: dihedral_theta = 1 (PI), dihedral_delta = 0.5
    // 3) the range of the normalized dihedral: [0, dihedral_theta * PI]
    halfedge_descriptor hd = mesh_.halfedge(ed);
    if (is_border(hd)) {
      hd = get_opposite(hd);
    }
    if (is_border(get_opposite(hd))) {  // ed is on the boundary
      set_halfedge_normal_dihedral(hd, np.dihedral_theta * CGAL_PI);
      set_halfedge_normal_dihedral(get_opposite(hd), -1.0);
      set_halfedge_is_crease(get_opposite(hd), false);
    } else {                            // ed is an inner edge
      FT normal_dihedral = calculate_normal_dihedral(hd);
      normal_dihedral /= np.dihedral_delta;  // normalize the normal_dihedral
      normal_dihedral = CGAL::min(normal_dihedral, np.dihedral_theta * CGAL_PI);
      if (reset_normal_dihedral) {
        if (get_halfedge_normal_dihedral(hd) == -1.0) {
          set_halfedge_is_crease(hd, get_halfedge_is_crease(get_opposite(hd)));
        }
        set_halfedge_normal_dihedral(hd, normal_dihedral);
      } else {
        if (get_halfedge_normal_dihedral(get_opposite(hd)) != -1.0) {
          hd = get_opposite(hd);
        }
        set_halfedge_normal_dihedral(hd, normal_dihedral);
      }
      set_halfedge_normal_dihedral(get_opposite(hd), -1.0);
      set_halfedge_is_crease(get_opposite(hd), false);
    }
  }

  void calculate_edge_classifications(FT feature_control_delta) {
    // precondition: edge and normal feature intensities has been computed
    for (typename Mesh::Vertex_range::const_iterator vi = mesh_.vertices().begin();
      vi != mesh_.vertices().end(); ++vi) {
      Halfedge_list effective_edges;
      collect_effective_edges_in_one_ring(feature_control_delta, *vi,
        &effective_edges);
      if (effective_edges.size() == 2) {
        for (auto it = effective_edges.begin();
          it != effective_edges.end(); ++it) {
          halfedge_descriptor hd = *it;
          if (get_halfedge_normal_dihedral(hd) == -1.0) {
            hd = get_opposite(hd);
          }
          set_halfedge_is_crease(hd, true);
        }
      }
    }
  }

  // 10) properties access
  inline void reset_face_tags(int value, const Face_list &faces) {
    for (auto it = faces.begin(); it != faces.end(); ++it) {
      face_descriptor fd = *it;
      set_face_tag(fd, value);
    }
  }

  void reset_face_properties(face_descriptor fd, face_descriptor fd_source) {
    // step 1: reset the basic properties
    if (fd_source == get_null_face()) {
      set_face_tag(fd, 0);
      set_face_normal(fd, CGAL::NULL_VECTOR);
      set_face_max_squared_error(fd, 0.0);
    } else {
      set_face_tag(fd, get_face_tag(fd_source));
      set_face_normal(fd, get_face_normal(fd_source));
      set_face_max_squared_error(fd, get_face_max_squared_error(fd_source));
    }
    // step 2: reset the sample links
    clear_face_links(fd);
  }

  void reset_halfedge_properties(halfedge_descriptor hd,
    halfedge_descriptor hd_source) {
    // step 1: reset the basic properties
    if (hd_source == get_null_halfedge()) {
      set_halfedge_tag(hd, 0);
      set_halfedge_normal_dihedral(hd, -1.0);
      set_halfedge_is_crease(hd, false);
    } else {
      set_halfedge_tag(hd, get_halfedge_tag(hd_source));
      set_halfedge_normal_dihedral(hd, get_halfedge_normal_dihedral(hd_source));
      set_halfedge_is_crease(hd, get_halfedge_is_crease(hd_source));
    }
    // step 2: reset the sample links
    clear_halfedge_links(hd);
  }

  void reset_vertex_properties(vertex_descriptor vd,
    vertex_descriptor vd_source) {
    // step 1: reset the basic properties
    if (vd_source == get_null_vertex()) {
      set_vertex_tag(vd, 0);
      set_vertex_max_dihedral(vd, -1.0);
      set_vertex_gaussian_curvature(vd, 0.0);
    } else {
      set_vertex_tag(vd, get_vertex_tag(vd_source));
      set_vertex_max_dihedral(vd, get_vertex_max_dihedral(vd_source));
      set_vertex_gaussian_curvature(vd,
          get_vertex_gaussian_curvature(vd_source));
    }
    // step 2: reset the sample links
    clear_vertex_links(vd);
  }

  void tag_component(face_descriptor seed_face,
    const int tag_free, const int tag_done) {
    set_face_tag(seed_face, tag_done);
    Face_list faces;
    faces.push_front(seed_face);
    while (!faces.empty()) {
      face_descriptor fd = faces.front();
      faces.pop_front();
      Halfedge_around_face_circulator fb(mesh_.halfedge(fd), mesh_), fe(fb);
      do {
        face_descriptor nfd = get_face(get_opposite(*fb));
        if (nfd != get_null_face() && get_face_tag(nfd) == tag_free) {
          set_face_tag(nfd, tag_done);
          faces.push_front(nfd);
        }
        ++fb;
      } while (fb != fe);
    }
  }

  void tag_faces(const int tag) {
    for (typename Mesh::Face_range::const_iterator fi = mesh_.faces().begin();
      fi != mesh_.faces().end(); ++fi) {
      set_face_tag(*fi, tag);
    }
  }

  void tag_halfedge(const int tag) {
    for (typename Mesh::Halfedge_range::const_iterator hi = mesh_.halfedges().begin();
      hi != mesh_.halfedges().end(); ++hi) {
      set_halfedge_tag(*hi, tag);
    }
  }

  // 11) surface mesh properties
  int get_face_out_link_count() const {
    size_t face_out_link_count = 0;
    for (typename Mesh::Face_range::const_iterator fi = mesh_.faces().begin();
      fi != mesh_.faces().end(); ++fi) {
      const Link_list &face_out_links = get_face_out_links(*fi);
      face_out_link_count += face_out_links.size();
    }
    return static_cast<int>(face_out_link_count);
  }

  int get_edge_out_link_count() const {
    size_t edge_out_link_count = 0;
    for (typename Mesh::Edge_range::const_iterator ei = mesh_.edges().begin();
      ei != mesh_.edges().end(); ++ei) {
      halfedge_descriptor hd = mesh_.halfedge(*ei);
      if (get_halfedge_normal_dihedral(hd) == -1.0) {
        hd = get_opposite(hd);
      }
      const Link_list &edge_out_links = get_halfedge_out_links(hd);
      edge_out_link_count += edge_out_links.size();
    }
    return static_cast<int>(edge_out_link_count);
  }

  int get_vertex_out_link_count() const {
    return static_cast<int>(mesh_.number_of_vertices());
  }

  FT calculate_avg_quality() const {
    FT sum_quality = 0.0;
    for (typename Mesh::Face_range::const_iterator fi = mesh_.faces().begin();
      fi != mesh_.faces().end(); ++fi) {
      sum_quality += calculate_quality(*fi);
    }
    if (mesh_.number_of_faces() == 0) {   // invalid case
      return -1.0;
    } else {
      return sum_quality / mesh_.number_of_faces();
    }
  }

  FT calculate_min_quality() const {
    FT min_quality = DOUBLE_MAX, quality = 0.0;
    for (typename Mesh::Face_range::const_iterator fi = mesh_.faces().begin();
      fi != mesh_.faces().end(); ++fi) {
      quality = calculate_quality(*fi);
      min_quality = CGAL::min(min_quality, quality);
    }
    return min_quality;
  }

  FT calculate_rms_distance() const {
    FT rms_distance = 0.0;
    size_t nb_samples = 0;
    for (typename Mesh::Face_range::const_iterator fi = mesh_.faces().begin();
      fi != mesh_.faces().end(); ++fi) {
      // face in links
      const Link_iter_list &face_in_links = get_face_in_links(*fi);
      for (Link_iter_list_const_iter cit = face_in_links.begin();
        cit != face_in_links.end(); ++cit) {
        Link_list_const_iter llcit = *cit;
        const Point &start = llcit->second.first;
        const Point &end = llcit->second.second;
        rms_distance += CGAL::squared_distance(start, end);
      }
      nb_samples += face_in_links.size();
      // edge in links
      const Link_iter_list &edge_in_links = get_edge_in_links(*fi);
      for (Link_iter_list_const_iter cit = edge_in_links.begin();
        cit != edge_in_links.end(); ++cit) {
        Link_list_const_iter llcit = *cit;
        const Point &start = llcit->second.first;
        const Point &end = llcit->second.second;
        rms_distance += CGAL::squared_distance(start, end);
      }
      nb_samples += edge_in_links.size();
      // vertex in links
      const Link_pointer_list &vertex_in_links = get_vertex_in_links(*fi);
      for (Link_pointer_const_iter cit = vertex_in_links.begin();
        cit != vertex_in_links.end(); ++cit) {
        const Link *link = *cit;
        const Point &start = link->second.first;
        const Point &end = link->second.second;
        rms_distance += CGAL::squared_distance(start, end);
      }
      nb_samples += vertex_in_links.size();
    }
    // face out links
    for (typename Mesh::Face_range::const_iterator fi = mesh_.faces().begin();
      fi != mesh_.faces().end(); ++fi) {
      const Link_list &face_out_links = get_face_out_links(*fi);
      for (Link_list_const_iter cit = face_out_links.begin();
        cit != face_out_links.end(); ++cit) {
        const Link &link = *cit;
        const Point &start = link.second.first;
        const Point &end = link.second.second;
        rms_distance += CGAL::squared_distance(start, end);
      }
      nb_samples += face_out_links.size();
    }
    // edge out links
    for (typename Mesh::Edge_range::const_iterator ei = mesh_.edges().begin();
      ei != mesh_.edges().end(); ++ei) {
      halfedge_descriptor hd = mesh_.halfedge(*ei);
      if (get_halfedge_normal_dihedral(hd) == -1.0) {
        hd = get_opposite(hd);
      }
      const Link_list &edge_out_links = get_halfedge_out_links(hd);
      for (Link_list_const_iter cit = edge_out_links.begin();
        cit != edge_out_links.end(); ++cit) {
        const Link &link = *cit;
        const Point &start = link.second.first;
        const Point &end = link.second.second;
        rms_distance += CGAL::squared_distance(start, end);
      }
      nb_samples += edge_out_links.size();
    }
    // vertex out links
    for (typename Mesh::Vertex_range::const_iterator vi = mesh_.vertices().begin();
      vi != mesh_.vertices().end(); ++vi) {
      const Link &link = get_vertex_out_link(*vi);
      const Point &start = link.second.first;
      const Point &end = link.second.second;
      rms_distance += CGAL::squared_distance(start, end);
    }
    nb_samples += mesh_.number_of_vertices();
    if (nb_samples == 0) {    // invalid case
      return -1.0;
    } else {
      return CGAL::sqrt(rms_distance / nb_samples);
    }
  }

  void trace_edge_length() const {
    if (mesh_.number_of_edges() == 0) {   // empty mesh
      std::cout << "Empty surface mesh" << std::endl;
      return;
    }
    FT sum_length = 0.0;
    FT min_length = MAX_VALUE;
    FT max_length = 0.0;
    int nb_edges = 0;
    for (typename Mesh::Edge_range::const_iterator ei = mesh_.edges().begin();
      ei != mesh_.edges().end(); ++ei) {
      const FT edge_length = length(mesh_.halfedge(*ei));
      sum_length += edge_length;
      min_length = CGAL::min(min_length, edge_length);
      max_length = CGAL::max(max_length, edge_length);
      ++nb_edges;
    }
    std::cout << "Min edge length: " << min_length << std::endl;
    std::cout << "Max edge length: " << max_length << std::endl;
    std::cout << "Average edge length: "
      << sum_length / nb_edges << std::endl;
  }

  FT calculate_regular_vertex_ratio() const {
    int nb_regular_vertices = 0;
    for (typename Mesh::Vertex_range::const_iterator vi = mesh_.vertices().begin();
      vi != mesh_.vertices().end(); ++vi) {
      nb_regular_vertices += is_regular(*vi);
    }
    size_t nb_vertices = mesh_.number_of_vertices();
    if (nb_vertices == 0) {  // invalid case
      return -1.0;
    } else {
      return static_cast<double>(nb_regular_vertices) / nb_vertices;
    }
  }

  FT calculate_smaller_angle_ratio(FT angle) const {
    FT radian = to_radian(angle);
    int nb_small_angle_faces = 0;
    for (typename Mesh::Face_range::const_iterator fi = mesh_.faces().begin();
      fi != mesh_.faces().end(); ++fi) {
      FT smallest_radian = calculate_smallest_radian(*fi);
      nb_small_angle_faces += smallest_radian < radian;
    }
    size_t nb_faces = mesh_.number_of_faces();
    if (nb_faces == 0) {    // invalid case
      return -1.0;
    } else {
      return static_cast<double>(nb_small_angle_faces) / nb_faces;
    }
  }

  unsigned int nb_boundaries() {     // count the number of boundaries
    unsigned int nb = 0;
    tag_halfedge(0);
    for (typename Mesh::Halfedge_range::const_iterator hi = mesh_.halfedges().begin();
      hi != mesh_.halfedges().end(); ++hi) {
      if (is_border(*hi) && get_halfedge_tag(*hi) == 0) {
        ++nb;
        halfedge_descriptor curr = *hi;
        do {
          set_halfedge_tag(curr, 1);
          curr = mesh_.next(curr);
        } while (curr != *hi);
      }
    }
    return nb;
  }

  unsigned int nb_components() {       // count the number of components
    unsigned int nb = 0;
    tag_faces(0);
    for (typename Mesh::Face_range::const_iterator fi = mesh_.faces().begin();
      fi != mesh_.faces().end(); ++fi) {
      if (get_face_tag(*fi) == 0) {
        ++nb;
        tag_component(*fi, 0, 1);
      }
    }
    return nb;
  }

  // 12) IO
  void save_as(const std::string &file_name) const {
    size_t pos = file_name.find_last_of('.');
    if (pos == std::string::npos) {
      std::cout << "Invalid file name." << std::endl;
      return;
    }
    std::string extension = file_name.substr(pos);
    std::transform(extension.begin(), extension.end(),
      extension.begin(), [](unsigned char c){ return std::tolower(c); } );
    if (extension == ".off") {
      save_as_off(file_name);
    } else {
      std::cout << "Invalid file name." << std::endl;
    }
  }

  void save_as_off(const std::string &file_name) const {
    std::ofstream ofs(file_name);
    if (!ofs) {
      return;
    }
    CGAL::set_ascii_mode(ofs);
    bool suc = CGAL::write_off(ofs, mesh_);
    ofs.close();
  }

  // 13) utilities
  Point_and_primitive_id get_closest_point_and_primitive(
      const std::set<face_descriptor> &in_link_faces,
      const Point &point) const {
    // find the closest point in the face set
    Point_and_primitive_id pp;
    FT min_sd = DOUBLE_MAX;
    Point nearest_point;
    for (auto it = in_link_faces.begin(); it != in_link_faces.end(); ++it) {
      face_descriptor fd = *it;
      FT sd = squared_distance(point, fd, &nearest_point);
      if (sd < min_sd) {
        pp.first = nearest_point;
        pp.second = fd;
        min_sd = sd;
      }
    }
    return pp;
  }

  FT squared_distance(const Point &p, face_descriptor fd,
    Point *nearest_point) const {
    halfedge_descriptor hd = mesh_.halfedge(fd);
    const Point &a = get_point(get_target_vertex(hd));
    const Point &b = get_point(get_opposite_vertex(hd));
    const Point &c = get_point(get_source_vertex(hd));
    Plane plane(a, b, c);
    Point projection = plane.projection(p);
    if (point_in_triangle(a, b, c, projection)) {
      *nearest_point = projection;
      return CGAL::squared_distance(p, projection);
    } else {
      Segment ab(a, b), bc(b, c), ca(c, a);
      FT sd_ab = CGAL::squared_distance(p, ab);
      FT sd_bc = CGAL::squared_distance(p, bc);
      FT sd_ca = CGAL::squared_distance(p, ca);
      Segment closest_segment;
      if (sd_ab < sd_bc) {
        closest_segment = sd_ab < sd_ca ? ab : ca;
      } else {
        closest_segment = sd_bc < sd_ca ? bc : ca;
      }
      *nearest_point = calculate_nearest_point(closest_segment, p);
      return CGAL::min(sd_ab, CGAL::min(sd_bc, sd_ca));
    }
  }

 private:
  Mesh &mesh_;

  Face_tags face_tags_;                         // face related properties
  Face_normals face_normals_;
  Face_max_errors face_max_squared_errors_;
  Face_link_list face_out_links_;
  Face_link_iter_list face_in_links_;
  Face_link_iter_list edge_in_links_;
  Face_link_pointer_list vertex_in_links_;

  Halfedge_tags halfedge_tags_;                 // halfedge related properties
  Halfedge_normal_dihedrals halfedge_normal_dihedrals_;
  Halfedge_are_creases halfedge_are_creases_;
  Halfedge_link_list halfedge_out_links_;

  Vertex_tags vertex_tags_;                     // vertex related properties
  Vertex_max_dihedral vertex_max_dihedrals_;
  Vertex_gaussian_curvature vertex_gaussian_curvatures_;
  Vertex_link vertex_out_link_;
};

}   // namespace internal
}   // namespace Polygon_mesh_processing
}   // namespace CGAL

#endif  // SRC_INTERNAL_MINANGLE_REMESHING_MESH_PROPERTIES_H_
