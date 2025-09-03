// Copyright (c) 2025 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Qijia Huang

#ifndef CGAL_VARIATIONAL_MEDIAL_AXIS_SAMPLING_H
#define CGAL_VARIATIONAL_MEDIAL_AXIS_SAMPLING_H

#include <CGAL/license/Variational_medial_axis.h>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_sphere_primitive_3.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Compact_container.h>
#include <CGAL/Compact_container_with_index.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <Eigen/Dense>
#include <algorithm>
#include <iterator>
#include <tuple>

#include "C3t3_type.h"

#ifdef CGAL_LINKED_WITH_TBB
#include <functional>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#endif

namespace CGAL {

#ifndef DOXYGEN_RUNNING

template <typename TriangleMesh, typename GT> class Medial_Sphere_Mesh;

template <typename TriangleMesh, typename GT> class Medial_Sphere
{
public:
  using FT = typename GT::FT;
  using Point_3 = typename GT::Point_3;
  using Sphere_3 = typename GT::Sphere_3;
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
  using MSMesh = Medial_Sphere_Mesh<TriangleMesh, GT>;
  using Sphere_ID = typename MSMesh::Sphere_ID;

  // Specalization for Compact_container_with_index
  std::size_t for_compact_container() const { return static_cast<std::size_t>(id_); }
  void for_compact_container(std::size_t idx) { id_ = Sphere_ID(idx); }

  Medial_Sphere(const Sphere_3& s)
      : sphere_(s)
      , split_vertex_(boost::graph_traits<TriangleMesh>::null_vertex())
      , error_(FT(0))
      , cluster_area_(FT(0))
      , do_not_split_(false)
      , id_(MSMesh::INVALID_SPHERE_ID) {}

  Medial_Sphere(Sphere_3&& s)
      : sphere_(std::move(s))
      , split_vertex_(boost::graph_traits<TriangleMesh>::null_vertex())
      , error_(FT(0))
      , cluster_area_(FT(0))
      , do_not_split_(false)
      , id_(MSMesh::INVALID_SPHERE_ID) {}

  void reset() {
    error_ = FT(0);
    split_vertex_ = boost::graph_traits<TriangleMesh>::null_vertex();
    cluster_area_ = FT(0);
    neighbors_.clear();
    cluster_vertices_.clear();
  }

  Sphere_3 get_sphere() const { return sphere_; }
  FT get_radius() const { return CGAL::approximate_sqrt(sphere_.squared_radius()); }
  Point_3 get_center() const { return sphere_.center(); }
  FT get_area() const { return cluster_area_; }
  Sphere_ID get_id() const { return id_; }
  const std::unordered_set<Sphere_ID>& get_neighbors() const { return neighbors_; }
  vertex_descriptor get_split_vertex() const { return split_vertex_; }
  FT get_error() const { return error_; }

  void set_center(const Point_3& p) { sphere_ = Sphere_3(p, get_radius()); }
  void set_radius(FT r) { sphere_ = Sphere_3(get_center(), r * r); }
  void set_cluster_area(FT area) { cluster_area_ = area; }
  void set_error(FT e) { error_ = e; }
  void set_split_vertex(vertex_descriptor v) { split_vertex_ = v; }
  void set_id(Sphere_ID id) { id_ = id; }
  void set_do_not_split(bool value) { do_not_split_ = value; }
  std::vector<vertex_descriptor>& get_cluster_vertices() { return cluster_vertices_; }
  void add_cluster_vertex(vertex_descriptor v) { cluster_vertices_.push_back(v); }
  void accumulate_cluster_area(FT area) { cluster_area_ += area; }
  void add_neighbor(Sphere_ID id) { neighbors_.insert(id); }
  bool can_split() const { return !do_not_split_ && split_vertex_ != boost::graph_traits<TriangleMesh>::null_vertex(); }

private:
  Sphere_3 sphere_;
  vertex_descriptor split_vertex_; // vertex chosen for split
  FT error_;
  FT cluster_area_;
  bool do_not_split_ = false; // flag to indicate if the sphere should not be split
  std::unordered_set<Sphere_ID> neighbors_;
  std::vector<vertex_descriptor> cluster_vertices_;
  Sphere_ID id_;
};

template <typename TriangleMesh, typename GT> class Medial_Sphere_Mesh
{
public:
  using FT = typename GT::FT;
  using Point_3 = typename GT::Point_3;
  using Sphere_3 = typename GT::Sphere_3;
  using MSphere = Medial_Sphere<TriangleMesh, GT>;
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;

  using Sphere_container = CGAL::Compact_container_with_index<
      MSphere,
      CGAL_ALLOCATOR(MSphere),
      Multiply_by_two_policy_for_cc_with_size<64>>;

  using Sphere_ID = typename Sphere_container::Index;

  inline static const Sphere_ID INVALID_SPHERE_ID = Sphere_container::null_descriptor;

public:
  Sphere_ID add_sphere(const Sphere_3& s) {
    Sphere_ID id = spheres_.emplace(s);
    spheres_[id].set_id(id);
    return id;
  }

  Sphere_ID add_sphere(Sphere_3&& s) {
    Sphere_ID id = spheres_.emplace(std::move(s));
    spheres_[id].set_id(id);
    return id;
  }
  MSphere& get_sphere(Sphere_ID idx) {
    CGAL_precondition(spheres_.owns(idx));
    return spheres_[idx];
  }

  const MSphere& get_sphere(Sphere_ID idx) const {
    CGAL_precondition(spheres_.owns(idx));
    return spheres_[idx];
  }

  void remove(Sphere_ID idx) {
    if(spheres_.owns(idx)) {
      spheres_.erase(idx);
    }
  }
  void reset() {
    for(auto& sphere : spheres_) {
      sphere.reset();
    }
  }
  std::size_t nb_spheres() const { return spheres_.size(); }
  bool empty() const { return spheres_.empty(); }
  void clear() { spheres_.clear(); }
  Sphere_container& spheres() { return spheres_; }
  const Sphere_container& spheres() const { return spheres_; }

private:
  Sphere_container spheres_;
};
#endif // DOXYGEN_RUNNING

/**  \ingroup PkgVMASRef
 * \brief Class representing the medial skeleton of a shape.
 *
 * This class provides methods to manage and export the medial skeleton of a shape.
 * It stores the medial spheres, edges, and faces of the skeleton.
 *
 * @tparam TriangleMesh The type of the triangle mesh representing the shape.
 * @tparam GT The geometric traits class used for geometric computations.
 *         <b>%Default:</b>
 * \code
 *     CGAL::Kernel_traits<
 *       boost::property_traits<
 *          boost::property_map<TriangleMesh, CGAL::vertex_point_t>::type
 *        >::value_type
 *      >::Kernel
 * \endcode
 */
template <typename TriangleMesh_, typename GeomTraits_ = Default> class Medial_Skeleton
{
  using GT = typename Default::Get<
      GeomTraits_,
      typename Kernel_traits<typename boost::property_traits<
          typename boost::property_map<TriangleMesh_, vertex_point_t>::type>::value_type>::Kernel>::type;
  using Sphere_3 = typename GT::Sphere_3;
  using Point_3 = typename GT::Point_3;
  using FT = typename GT::FT;
  using MSMesh = Medial_Sphere_Mesh<TriangleMesh_, GT>;
  using Sphere_ID = typename MSMesh::Sphere_ID;

public:
  /**
   * loads a medial skeleton from a PLY file.
   *
   * @param filepath
   *     Filepath to the PLY file containing the medial skeleton data.
   * @return
   *     True if the skeleton was successfully loaded, false otherwise.
   *
   * Note: The file format is :
   * ```
   * ply
   * format ascii 1.0
   * element vertex N
   * property float x
   * property float y
   * property float z
   * property float radius
   * element edge M
   * property int vertex1
   * property int vertex2
   * element face K
   * property list uchar int vertex_indices
   * end_header
   * x1 y1 z1 r1
   * ...          // N vertices
   * xn yn zn rn
   * vx vy
   * ...          // M edges
   * vz vw
   * 3 v1 v2 v3
   * ...          // K faces
   * 3 vx vy vz
   * ```
   */
  bool load_skeleton_from_ply(std::string& filepath) {
    clear();
    std::ifstream ifs(filepath);
    if(!ifs) {
      std::cerr << "Error opening file: " << filepath << std::endl;
      return false;
    }

    std::string line;
    std::size_t num_vertices = 0, num_edges = 0, num_faces = 0;
    bool in_header = true;

    while(std::getline(ifs, line) && in_header) {
      std::istringstream iss(line);
      std::string token;
      iss >> token;

      if(token == "ply") {
        continue;
      } else if(token == "format") {
        std::string format_type;
        iss >> format_type;
        if(format_type != "ascii") {
          std::cerr << "Error: Only ASCII PLY format is supported" << std::endl;
          return false;
        }
      } else if(token == "element") {
        std::string element_type;
        std::size_t count;
        iss >> element_type >> count;

        if(element_type == "vertex") {
          num_vertices = count;
        } else if(element_type == "edge") {
          num_edges = count;
        } else if(element_type == "face") {
          num_faces = count;
        }
      } else if(token == "end_header") {
        break;
      }
    }

    vertices_.reserve(num_vertices);
    for(std::size_t i = 0; i < num_vertices; ++i) {
      if(!std::getline(ifs, line)) {
        std::cerr << "Error: Unexpected end of file while reading vertices" << std::endl;
        return false;
      }

      std::istringstream iss(line);
      double x, y, z, radius;
      if(!(iss >> x >> y >> z >> radius)) {
        std::cerr << "Error: Invalid vertex data at line " << i + 1 << std::endl;
        return false;
      }

      Point_3 center(x, y, z);
      FT squared_radius = FT(radius * radius);
      vertices_.emplace_back(center, squared_radius);
    }

    edges_.reserve(num_edges);
    for(std::size_t i = 0; i < num_edges; ++i) {
      if(!std::getline(ifs, line)) {
        std::cerr << "Error: Unexpected end of file while reading edges" << std::endl;
        return false;
      }

      std::istringstream iss(line);
      std::size_t v1, v2;
      if(!(iss >> v1 >> v2)) {
        std::cerr << "Error: Invalid edge data at line " << i + 1 << std::endl;
        return false;
      }

      if(v1 >= num_vertices || v2 >= num_vertices) {
        std::cerr << "Error: Edge references invalid vertex indices" << std::endl;
        return false;
      }

      edges_.emplace_back(v1, v2);
    }

    faces_.reserve(num_faces);
    for(std::size_t i = 0; i < num_faces; ++i) {
      if(!std::getline(ifs, line)) {
        std::cerr << "Error: Unexpected end of file while reading faces" << std::endl;
        return false;
      }

      std::istringstream iss(line);
      std::size_t vertex_count;
      if(!(iss >> vertex_count)) {
        std::cerr << "Error: Invalid face data at line " << i + 1 << std::endl;
        return false;
      }

      if(vertex_count != 3) {
        std::cerr << "Error: Only triangular faces are supported" << std::endl;
        return false;
      }

      std::size_t v1, v2, v3;
      if(!(iss >> v1 >> v2 >> v3)) {
        std::cerr << "Error: Invalid face vertex indices" << std::endl;
        return false;
      }

      if(v1 >= num_vertices || v2 >= num_vertices || v3 >= num_vertices) {
        std::cerr << "Error: Face references invalid vertex indices" << std::endl;
        return false;
      }

      faces_.push_back({v1, v2, v3});
    }

    ifs.close();

    std::cout << "Successfully loaded skeleton from " << filepath << std::endl;
    std::cout << "Vertices: " << num_vertices << ", Edges: " << num_edges << ", Faces: " << num_faces << std::endl;

    return true;
  }
#ifndef DOXYGEN_RUNNING
  void build_skeleton_from_medial_sphere_mesh(const MSMesh& sphere_mesh) {
    clear();

    std::unordered_map<Sphere_ID, std::size_t> id_to_index;

    // Convert spheres to vertices (as Sphere_3 objects)
    std::size_t vertex_idx = 0;
    for(const auto& sphere : sphere_mesh.spheres()) {
      vertices_.push_back(sphere.get_sphere()); // Store the complete Sphere_3
      id_to_index[sphere.get_id()] = vertex_idx++;
    }

    // Convert sphere adjacencies to edges
    std::set<std::pair<std::size_t, std::size_t>> edge_set;
    for(const auto& sphere : sphere_mesh.spheres()) {
      std::size_t a = id_to_index[sphere.get_id()];
      for(Sphere_ID neighbor_id : sphere.get_neighbors()) {
        if(id_to_index.find(neighbor_id) != id_to_index.end()) {
          std::size_t b = id_to_index[neighbor_id];
          if(a < b) {
            edge_set.emplace(a, b);
          }
        }
      }
    }
    edges_.assign(edge_set.begin(), edge_set.end());

    // Convert triangular adjacencies to faces
    std::set<std::array<std::size_t, 3>> face_set;
    for(const auto& sphere : sphere_mesh.spheres()) {
      Sphere_ID a_id = sphere.get_id();
      const auto& neighbors_a = sphere.get_neighbors();

      for(Sphere_ID b_id : neighbors_a) {
        if(b_id <= a_id)
          continue;
        if(id_to_index.find(b_id) == id_to_index.end())
          continue;

        const auto& neighbors_b = sphere_mesh.get_sphere(b_id).get_neighbors();

        for(Sphere_ID c_id : neighbors_a) {
          if(c_id <= b_id)
            continue;
          if(id_to_index.find(c_id) == id_to_index.end())
            continue;

          if(neighbors_b.find(c_id) != neighbors_b.end()) {
            std::array<std::size_t, 3> tri = {id_to_index[a_id], id_to_index[b_id], id_to_index[c_id]};
            std::sort(tri.begin(), tri.end());
            face_set.insert(tri);
          }
        }
      }
    }
    faces_.assign(face_set.begin(), face_set.end());
  }

  MSMesh build_medial_sphere_mesh_from_skeleton() const {
    MSMesh sphere_mesh;
    std::vector<Sphere_ID> index_to_id_map;
    index_to_id_map.reserve(vertices_.size());

    for(const auto& v_sphere : vertices_) {
      index_to_id_map.push_back(sphere_mesh.add_sphere(v_sphere));
    }
    for(const auto& edge : edges_) {
      if(edge.first < index_to_id_map.size() && edge.second < index_to_id_map.size()) {
        Sphere_ID id1 = index_to_id_map[edge.first];
        Sphere_ID id2 = index_to_id_map[edge.second];
        sphere_mesh.get_sphere(id1).add_neighbor(id2);
        sphere_mesh.get_sphere(id2).add_neighbor(id1);
      }
    }

    return sphere_mesh;
  }
#endif // DoXYGEN_RUNNING

  /// \name Accessor methods
  /// @{
  /**
   * returns the container of vertices, where each vertex is a medial sphere (`Sphere_3`).
   */
  const std::vector<Sphere_3>& vertices() const { return vertices_; }
  /**
   * returns the container of edges, where each edge is represented as a pair of indices in the vertices vector.
   */
  const std::vector<std::pair<std::size_t, std::size_t>>& edges() const { return edges_; }
  /**
   * returns the container of faces, where each face is represented as an array of three indices in the vertices vector.
   */
  const std::vector<std::array<std::size_t, 3>>& faces() const { return faces_; }
  /**
   * returns the number of vertices in the medial skeleton.
   */
  std::size_t number_of_vertices() const { return vertices_.size(); }
  /**
   * returns the number of edges in the medial skeleton.
   */
  std::size_t number_of_edges() const { return edges_.size(); }
  /**
   * returns the number of faces in the medial skeleton.
   */
  std::size_t number_of_faces() const { return faces_.size(); }

  /**
   * clears the data for the medial skeleton.
   */
  void clear() {
    vertices_.clear();
    edges_.clear();
    faces_.clear();
  }
  /**
   * sets the data for the medial skeleton.
   * @param vertices A vector of `Sphere_3` representing the vertices (medial spheres).
   * @param edges A vector of pairs representing the edges, where each pair contains indices of vertices.
   * @param faces A vector of arrays representing the faces, where each array contains three indices of vertices.
   */
  void set_data(std::vector<Sphere_3>&& vertices,
                std::vector<std::pair<std::size_t, std::size_t>>&& edges,
                std::vector<std::array<std::size_t, 3>>&& faces) {
    vertices_ = std::move(vertices);
    edges_ = std::move(edges);
    faces_ = std::move(faces);
  }

  /// @}
private:
  // Data members
  std::vector<Sphere_3> vertices_; // Each vertex is a complete medial sphere
  std::vector<std::pair<std::size_t, std::size_t>> edges_;
  std::vector<std::array<std::size_t, 3>> faces_;
};

template <class TriangleMesh_, class GeomTraits_ = Default> class Medial_skeleton_offset_function
{
  using GT = typename Default::Get<
      GeomTraits_,
      typename Kernel_traits<typename boost::property_traits<
          typename boost::property_map<TriangleMesh_, vertex_point_t>::type>::value_type>::Kernel>::type;
  using FT = typename GT::FT;
  using Point_3 = typename GT::Point_3;
  using Vector_3 = typename GT::Vector_3;
  using Sphere_3 = typename GT::Sphere_3;
  using MSkeleton = Medial_Skeleton<TriangleMesh_, GeomTraits_>;
  // AABB tree types over spheres
  using Iterator = typename std::vector<Sphere_3>::const_iterator;
  using Primitive = CGAL::AABB_sphere_primitive_3<GT, Iterator>;
  using Traits = CGAL::AABB_traits_3<GT, Primitive>;
  using Tree = CGAL::AABB_tree<Traits>;


public:
  Medial_skeleton_offset_function(const MSkeleton& skeleton)
      : skeleton_(skeleton)
      , spheres_(&skeleton_.vertices())
      , faces_(&skeleton_.faces())
      , tree_(spheres_->begin(), spheres_->end()) {
    const std::size_t n = skeleton_.number_of_vertices();
    radii_.reserve(n);
    for(const auto& sph : *spheres_)
      radii_.push_back(CGAL::approximate_sqrt(sph.squared_radius()));

    // Build adjacency for cones (edges) and incident face lists (slabs)
    edge_adj_.resize(n);
    face_incident_.resize(n);

    const auto& edges = skeleton_.edges();
    for(const auto& e : edges) {
      const std::size_t i = e.first, j = e.second;
      if(i < n && j < n) {
        edge_adj_[i].push_back(j);
        edge_adj_[j].push_back(i);
      }
    }
    for(std::size_t fid = 0; fid < faces_->size(); ++fid) {
      const auto& f = (*faces_)[fid];
      const std::size_t a = f[0], b = f[1], c = f[2];
      if(a < n)
        face_incident_[a].push_back(fid);
      if(b < n)
        face_incident_[b].push_back(fid);
      if(c < n)
        face_incident_[c].push_back(fid);
    }

    tree_.accelerate_distance_queries();
  }

  FT operator()(const Point_3& p) const {
    FT dmin = (std::numeric_limits<FT>::max)();

    // Find closest sphere with AABB tree
    auto pp = tree_.closest_point_and_primitive(p);
    const Iterator it = pp.second; // iterator into spheres_
    const std::size_t i = static_cast<std::size_t>(std::distance(spheres_->begin(), it));

    eval_sphere(i, p, dmin);
    for(std::size_t j : edge_adj_[i]) {
      eval_sphere(j, p, dmin);
    }

    // Cones (edges) incident to the closest sphere
    for(std::size_t j : edge_adj_[i]) {
      eval_cone(i, j, p, dmin);
    }

    // Slabs (faces) incident to the closest sphere
    for(std::size_t fid : face_incident_[i]) {
      const auto& tri = (*faces_)[fid];
      eval_slab(tri[0], tri[1], tri[2], p, dmin);
    }

    return dmin;
  }

private:
  int solve_quadric(const FT a, const FT b, const FT c, FT& t1, FT& t2) const {
    if(a == FT(0)) {
      if(b == FT(0))
        return 0;
      t1 = -c / b;
      return 1;
    } else {
      FT delta = b * b - 4 * a * c;
      if(delta < FT(0))
        return 0;
      else if(delta == FT(0)) {
        t1 = -b / (FT(2) * a);
        return 1;
      } else {
        const FT sqrt_delta = CGAL::approximate_sqrt(delta);
        t1 = (-b - sqrt_delta) / (FT(2) * a);
        t2 = (-b + sqrt_delta) / (FT(2) * a);
        return 2;
      }
    }
  }

  FT sphere_distance(const Point_3& p, const Point_3& c, const FT r) const {
    return CGAL::approximate_sqrt((p - c).squared_length()) - r;
  }

  FT cone_distance(const Point_3& p, const Point_3& c1, const Point_3& c2, const FT r1, const FT r2, FT t) const {
    const Point_3 c = Point_3((FT(1) - t) * c1.x() + t * c2.x(), (FT(1) - t) * c1.y() + t * c2.y(), (FT(1) - t) * c1.z() + t * c2.z());
    const FT r = (FT(1) - t) * r1 + t * r2;
    return sphere_distance(p, c, r);
  }

  FT slab_distance(const Point_3& p,
                   const Point_3& c1,
                   const Point_3& c2,
                   const Point_3& c3,
                   const FT& r1,
                   const FT& r2,
                   const FT& r3,
                   FT t1,
                   FT t2) const {
    const Point_3 c = Point_3(t1 * c1.x() + t2 * c2.x() + (FT(1) - t1 - t2) * c3.x(),
                              t1 * c1.y() + t2 * c2.y() + (FT(1) - t1 - t2) * c3.y(),
                              t1 * c1.z() + t2 * c2.z() + (FT(1) - t1 - t2) * c3.z());
    const FT r = t1 * r1 + t2 * r2 + (FT(1) - t1 - t2) * r3;
    return sphere_distance(p, c, r);
  }

  void eval_sphere(std::size_t idx, const Point_3& p, FT& dmin) const {
    const Point_3& c = spheres_->at(idx).center();
    const FT r = radii_[idx];
    const FT d = CGAL::approximate_sqrt((p - c).squared_length()) - r;
    if(d < dmin)
      dmin = d;
  }

  void eval_cone(std::size_t i, std::size_t j, const Point_3& p, FT& dmin) const {
    const Point_3& c1 = spheres_->at(i).center();
    const Point_3& c2 = spheres_->at(j).center();
    const FT r1 = radii_[i];
    const FT r2 = radii_[j];

    const Vector_3 c21 = c2 - c1;
    const Vector_3 c1p = p - c1;
    const FT r21 = r2 - r1;
    const FT a = CGAL::scalar_product(c21, c21);
    const FT b = CGAL::scalar_product(c21, c1p);
    const FT c = CGAL::scalar_product(c1p, c1p);
    const FT A = a * (a - r21 * r21);
    const FT B = 2 * b * (r21 * r21 - a);
    const FT C = b * b - r21 * r21 * c;
    FT t1 = 0, t2 = 0, dist1, dist2;
    int root_nb = solve_quadric(A, B, C, t1, t2);
    if(root_nb != 0) {
      t1 = std::clamp(t1, FT(0), FT(1));
      dist1 = cone_distance(p, c1, c2, r1, r2, t1);
      dmin = std::min(dmin, dist1);
      if(root_nb != 1) {
        t2 = std::clamp(t2, FT(0), FT(1));
        dist2 = cone_distance(p, c1, c2, r1, r2, t2);
        dmin = std::min(dmin, dist2);
      }
    } else {
      dmin = std::min(dmin, cone_distance(p, c1, c2, r1, r2, FT(0)));
      dmin = std::min(dmin, cone_distance(p, c1, c2, r1, r2, FT(1)));
    }
  }

  void eval_slab(std::size_t ia, std::size_t ib, std::size_t ic, const Point_3& p, FT& dmin) const {
    const Point_3& c1 = spheres_->at(ia).center();
    const Point_3& c2 = spheres_->at(ib).center();
    const Point_3& c3 = spheres_->at(ic).center();
    const FT r1 = radii_[ia];
    const FT r2 = radii_[ib];
    const FT r3 = radii_[ic];

    // c(t1, t2) = t1 * c1 + t2 * c2+ (1-t1-t2) * c3 = (c1-c3) * t1 + (c2-c3) * t2 + c3
    // r(t1, t2) = t1 * r1 + t2 * r2 + (1-t1-t2) * r3 = (r1-r3) * t1 + (r2-r3) * t2 + r3
    // f(t1, t2) = ||c(t1, t2) - p || - r(t1, t2)
    // obj: argmin_(t1,t2) f(t1, t2)
    const Vector_3 c13 = c1 - c3;
    const Vector_3 c23 = c2 - c3;
    const Vector_3 c3p = c3 - p;
    const FT r13 = r1 - r3;
    const FT r23 = r2 - r3;
    // let x = c(t1, t2) - p = c13 * t1 + c23 * t2 + c3p
    //  f(t1, t2) = ||x|| - (r1-r3) * t1 + (r2-r3) * t2 + r3
    // df(t1,t2) / dt1 = (x*c13)/||x|| - r13 = 0 ===> x*c13 = ||x||* r13
    // df(t1,t2) / dt2 = (x*c23)/||x|| - r23 = 0 ===> x*c23 = ||x||* r23
    // x^2 = ||c13||^2 *t1^2 + ||c23||^2 *t2^2 + 2*(c13*c23)*t1*t2 + 2*(c13*c3p)*t1 + 2*(c23*c3p)*t2 + ||c3p||^2
    const FT a = CGAL::scalar_product(c13, c13);
    const FT b = CGAL::scalar_product(c23, c23);
    const FT c = CGAL::scalar_product(c13, c23);
    const FT d = CGAL::scalar_product(c13, c3p);
    const FT e = CGAL::scalar_product(c23, c3p);
    const FT f = CGAL::scalar_product(c3p, c3p);
    // x^2 = a*t1^2 + b*t2^2 + 2*c*t1*t2 + 2*d*t1 + 2*e*t2 + f
    FT t1 = 0, t2 = 0, dist1 = 0, dist2 = 0;

    // three spheres have the same radius
    if(r13 == FT(0) && r23 == FT(0)) {
      //     x*c13                                         = 0
      //===> ||c13||^2 * t1 + (c13 *c23)* t2 + (c13 * c3p) = 0
      //===> a*t1 + c*t2 + d                               = 0

      //     x*c23                                         = 0
      //===> ||c23||^2 * t2 + (c23 *c13)* t1 + (c23 * c3p) = 0
      //===> b*t2 + c*t1 + e                               = 0

      const FT denom = a * b - c * c;
      t1 = (c * e - b * d) / denom;
      t2 = (c * d - a * e) / denom;
    } else if(r13 == FT(0) && r23 != FT(0)) {
      // x * c13 = 0 ===> a*t1 + c*t2 + d = 0 ===> t1 = -(c/a)*t2 - d/a
      const FT h = -c / a;
      const FT k = -d / a;

      //         df(t1,t2) / dt2 = 0
      //===>               x*c23 = ||x||* r23
      //===>         (x * c23)^2 = ||x||^2 * r23^2
      //===> (b*t2 + c*t1 + e)^2 = ||x||^2 * r23^2
      const FT A = (b + c * h) * (b + c * h) - r23 * r23 * (a * h * h + b + FT(2) * c * h);
      const FT B =
          (FT(2) * (b + c * h) * (c * k + e) - r23 * r23 * (FT(2) * a * h * k + FT(2) * c * k + FT(2) * d * h + FT(2) * e));
      const FT C = (c * k + e) * (c * k + e) - r23 * r23 * (a * k * k + FT(2) * d * k + f);
      FT t1_1 = 0, t1_2 = 0, t2_1 = 0, t2_2 = 0;
      solve_quadric(A, B, C, t2_1, t2_2);
      t1_1 = h * t2_1 + k;
      t1_2 = h * t2_2 + k;
      dist1 = slab_distance(p, c1, c2, c3, r1, r2, r3, t1_1, t2_1);
      dist2 = slab_distance(p, c1, c2, c3, r1, r2, r3, t1_2, t2_2);
      if(dist2 < dist1) {
        t1 = t1_2;
        t2 = t2_2;
      } else {
        t1 = t1_1;
        t2 = t2_1;
      }
    } else if(r13 != FT(0) && r23 == FT(0)) {
      // x * c23 = 0 ===> b*t2 + c*t1 + e = 0 ===> t1 = -(b/c)*t2 - e/c
      const FT h = -b / c;
      const FT k = -e / c;

      //         df(t1,t2) / dt1 = 0
      //===>               x*c13 = ||x||* r13
      //===>         (x * c13)^2 = ||x||^2 * r13^2
      //===> (a*t1 + c*t2 + d)^2 = ||x||^2 * r13^2
      const FT A = (c + a * h) * (c + a * h) - r13 * r13 * (a * h * h + b + FT(2) * c * h);
      const FT B =
          (2 * (c + a * h) * (a * k + d) - r13 * r13 * (FT(2) * a * h * k + FT(2) * c * k + FT(2) * d * h + FT(2) * e));
      const FT C = (a * k + d) * (a * k + d) - r13 * r13 * (a * k * k + FT(2) * d * k + f);
      FT t1_1 = 0, t1_2 = 0, t2_1 = 0, t2_2 = 0;
      solve_quadric(A, B, C, t2_1, t2_2);
      t1_1 = h * t2_1 + k;
      t1_2 = h * t2_2 + k;
      dist1 = slab_distance(p, c1, c2, c3, r1, r2, r3, t1_1, t2_1);
      dist2 = slab_distance(p, c1, c2, c3, r1, r2, r3, t1_2, t2_2);
      if(dist2 < dist1) {
        t1 = t1_2;
        t2 = t2_2;
      } else {
        t1 = t1_1;
        t2 = t2_1;
      }
    } else {
      // both spheres have different radii
      // x * c13 = ||x|| * r13
      // x * c23 = ||x|| * r23
      // ===> (x* c13)/r13   = (x * c23)/r23
      // ===> r23 * (x* c13) = r13 * (x * c23)
      // ===> (r23 * a- r13* c)*t1 + (r23*c-r13*b) *t2- (r23 * d - r13 * e) = 0
      const FT u = r23 * a - r13 * c;
      const FT v = r23 * c - r13 * b;
      const FT w = r23 * d - r13 * e;
      // ===> u*t1 + v*t2 + w = 0
      if(u == 0 && v != 0) {
        t2 = -w / v;
        const FT A = a * a - r13 * r13 * a;
        const FT B = FT(2) * a * (c * t2 + d) - r13 * r13 * (FT(2) * c * t2 + FT(2) * d);
        const FT C = (c * t2 + d) * (c * t2 + d) - r13 * r13 * (b * t2 * t2 + FT(2) * e * t2 + f);
        FT t1_1 = 0, t1_2 = 0;
        solve_quadric(A, B, C, t1_1, t1_2);
        dist1 = slab_distance(p, c1, c2, c3, r1, r2, r3, t1_1, t2);
        dist2 = slab_distance(p, c1, c2, c3, r1, r2, r3, t1_2, t2);
        if(dist2 < dist1) {
          t1 = t1_2;
        } else {
          t1 = t1_1;
        }
      } else if(u != 0 && v == 0) {
        t1 = -w / u;
        const FT A = b * b - r23 * r23 * b;
        const FT B = FT(2) * b * (c * t1 + e) - r23 * r23 * (FT(2) * c * t1 + FT(2) * e);
        const FT C = (c * t1 + e) * (c * t1 + e) - r23 * r23 * (a * t1 * t1 + FT(2) * d * t1 + f);
        FT t2_1 = 0, t2_2 = 0;
        solve_quadric(A, B, C, t2_1, t2_2);
        dist1 = slab_distance(p, c1, c2, c3, r1, r2, r3, t1, t2_1);
        dist2 = slab_distance(p, c1, c2, c3, r1, r2, r3, t1, t2_2);
        if(dist2 < dist1) {
          t2 = t2_2;
        } else {
          t2 = t2_1;
        }
      } else {
        const FT h = -u / v;
        const FT k = -w / v;
        const FT A = (b + c * h) * (b + c * h) - r23 * r23 * (a * h * h + b + FT(2) * c * h);
        const FT B = (2 * (b + c * h) * (c * k + e) -
                      r23 * r23 * (FT(2) * a * h * k + FT(2) * c * k + FT(2) * d * h + FT(2) * e));
        const FT C = (c * k + e) * (c * k + e) - r23 * r23 * (a * k * k + FT(2) * d * k + f);
        FT t1_1 = 0, t1_2 = 0, t2_1 = 0, t2_2 = 0;
        solve_quadric(A, B, C, t2_1, t2_2);
        t1_1 = h * t2_1 + k;
        t1_2 = h * t2_2 + k;
        dist1 = slab_distance(p, c1, c2, c3, r1, r2, r3, t1_1, t2_1);
        dist2 = slab_distance(p, c1, c2, c3, r1, r2, r3, t1_2, t2_2);
        if(dist2 < dist1) {
          t1 = t1_2;
          t2 = t2_2;
        } else {
          t1 = t1_1;
          t2 = t2_1;
        }
      }
    }
    if((t1 + t2) < FT(1) && t1 > FT(0) && t2 > FT(0) && t1 < FT(1) && t2 < FT(1)) {
      dmin = std::min(dmin, slab_distance(p, c1, c2, c3, r1, r2, r3, t1, t2));
      return;
    }
    eval_cone(ia, ib, p, dmin);
    eval_cone(ia, ic, p, dmin);
    eval_cone(ib, ic, p, dmin);
  }

private:
  const MSkeleton& skeleton_;
  const std::vector<Sphere_3>* spheres_;
  const std::vector<std::array<std::size_t, 3>>* faces_;
  std::vector<FT> radii_;

  // adjacency
  std::vector<std::vector<std::size_t>> edge_adj_;
  std::vector<std::vector<std::size_t>> face_incident_;

  Tree tree_;
};
/// \ingroup PkgVMASRef
/// \brief Algorithm class for extracting a variational medial skeleton from a triangulated surface mesh.
///
/// This class takes as input a triangulated surface mesh and iteratively samples the medial axis
/// by optimizing the placement of medial spheres. The result is a non-manifold triangle mesh that
/// captures the medial structure of the shape, consisting of both curve segments (1D) and triangle face patches (2D),
/// derived from the adjacency of clusters associated with each medial sphere.
///
/// The process terminates when either the desired number of spheres is reached, or a maximum number
/// of iterations is exceeded.
///
/// \note This method is designed to generate coarse approximations of the medial axis. The number of
/// spheres that can be reliably generated depends on the density of the input surface sampling. For best
/// results, we recommend keeping the total number of medial spheres under nb_vertices/100.
/// In addition, the user can explicitly request cluster merging or splitting operations to locally simplify or refine
/// the skeleton structure.
///
/// @tparam TriangleMesh_
///         a model of `FaceListGraph`
///
/// @tparam ConcurrencyTag_
///         a tag indicating whether the algorithm should run sequentially or in parallel.
///         This determines the execution mode at compile time.
///         <b>%Default:</b> `CGAL::Sequential_tag`<br>
///         <b>%Valid values:</b> `CGAL::Sequential_tag`, `CGAL::Parallel_tag`, `CGAL::Parallel_if_available_tag`
///
///  @tparam AccelerationType_
///         a tag indicating whether the algorithm should use Kd-tree or BVH as acceleration structure.
///         <b>%Default:</b> `CGAL::KD_tree_tag`<br>
///         <b>%Valid values:</b> `CGAL::KD_tree_tag`, `CGAL::BVH_tag`,
///
/// @tparam GeomTraits_
///         a model of `Kernel`<br>
///         <b>%Default:</b>
/// \code
///     CGAL::Kernel_traits<
///       boost::property_traits<
///          boost::property_map<TriangleMesh_, CGAL::vertex_point_t>::type
///        >::value_type
///      >::Kernel
/// \endcode
///
/// @tparam VertexPointMap_
///         a model of `ReadWritePropertyMap`
///         with `boost::graph_traits<TriangleMesh_>::%vertex_descriptor` as key and
///         `GeomTraits_::Point_3` as value type.<br>
///         <b>%Default:</b>
/// \code
///   boost::property_map<TriangleMesh_, CGAL::vertex_point_t>::const_type.
/// \endcode
///
///

template <typename TriangleMesh_,
          typename ConcurrencyTag_ = Sequential_tag,
          typename AccelerationType_ = KD_tree_tag,
          typename GeomTraits_ = Default,
          typename VertexPointMap_ = Default>
class Variational_medial_axis
{
  using VPM = typename Default::Get<VertexPointMap_,
                                    typename boost::property_map<TriangleMesh_, vertex_point_t>::const_type>::type;
  using GT = typename Default::Get<
      GeomTraits_,
      typename Kernel_traits<typename boost::property_traits<
          typename boost::property_map<TriangleMesh_, vertex_point_t>::type>::value_type>::Kernel>::type;

  using FT = typename GT::FT;
  using Point_3 = typename GT::Point_3;
  using Vector_3 = typename GT::Vector_3;
  using Sphere_3 = typename GT::Sphere_3;
  using MSMesh = Medial_Sphere_Mesh<TriangleMesh_, GT>;
  using MSphere = typename MSMesh::MSphere;
  using Sphere_ID = typename MSMesh::Sphere_ID;

  using Tree = AABB_tree<AABB_traits_3<GT, AABB_face_graph_triangle_primitive<TriangleMesh_, VPM>>>;

  using face_descriptor = typename boost::graph_traits<TriangleMesh_>::face_descriptor;
  using edge_descriptor = typename boost::graph_traits<TriangleMesh_>::edge_descriptor;
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh_>::vertex_descriptor;

  // Property map types
  using Vertex_normal_tag = CGAL::dynamic_vertex_property_t<Vector_3>;
  using Vertex_normal_map = typename boost::property_map<TriangleMesh_, Vertex_normal_tag>::const_type;
  using Vertex_area_tag = CGAL::dynamic_vertex_property_t<FT>;
  using Vertex_area_map = typename boost::property_map<TriangleMesh_, Vertex_area_tag>::const_type;
  using Face_normal_tag = CGAL::dynamic_face_property_t<Vector_3>;
  using Face_normal_map = typename boost::property_map<TriangleMesh_, Face_normal_tag>::const_type;
  using Face_area_tag = CGAL::dynamic_face_property_t<FT>;
  using Face_area_map = typename boost::property_map<TriangleMesh_, Face_area_tag>::const_type;
  using Vertex_error_tag = CGAL::dynamic_vertex_property_t<FT>;
  using Vertex_error_map = typename boost::property_map<TriangleMesh_, Vertex_error_tag>::const_type;
  using Vertex_cluster_sphere_tag = CGAL::dynamic_vertex_property_t<Sphere_ID>;
  using Vertex_cluster_sphere_map = typename boost::property_map<TriangleMesh_, Vertex_cluster_sphere_tag>::const_type;
  using Vertex_medial_sphere_pos_tag = CGAL::dynamic_vertex_property_t<Point_3>;
  using Vertex_medial_sphere_pos_map =
      typename boost::property_map<TriangleMesh_, Vertex_medial_sphere_pos_tag>::const_type;
  using Vertex_medial_sphere_radius_tag = CGAL::dynamic_vertex_property_t<FT>;
  using Vertex_medial_sphere_radius_map =
      typename boost::property_map<TriangleMesh_, Vertex_medial_sphere_radius_tag>::const_type;

public:
  /// \name Constructor
  ///@{
  ///
  /// The constructor of a vmas object.
  ///
  /// @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  ///
  /// @param tmesh
  ///        The input triangle mesh without borders.
  /// @param gt
  ///        The geometric traits class used for geometric computations.
  /// @param np
  ///        An optional sequence of \ref bgl_namedparameters "Named Parameters", listed below:
  ///
  /// \cgalNamedParamsBegin
  ///   \cgalParamNBegin{number_of_spheres}
  ///     \cgalParamDescription{The desired number of medial spheres in the resulting skeleton.}
  ///     \cgalParamType{unsigned int}
  ///     \cgalParamDefault{100}
  /// \cgalParamNBegin{max_iteration_number}
  ///    \cgalParamDescription{The maximum number of iterations for the optimization process.}
  ///    \cgalParamType{int}
  ///    \cgalParamDefault{1000}
  ///  \cgalParamNEnd
  ///   \cgalParamNBegin{lambda}
  ///     \cgalParamDescription{A weight balancing the two energy terms (SQEM and Euclidean). Smaller values encourage
  ///     the skeleton to extend deeper into local geometric features of the shape.}
  ///     \cgalParamType{FT}
  ///     \cgalParamDefault{FT(0.2)}
  ///     \cgalParamExtra{The range of this parameter is (0,1].}
  ///   \cgalParamNEnd
  ///   \cgalParamNBegin{vertex_point_map}
  ///     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
  ///     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
  ///                    as key type and `%Point_3` as value type}
  ///     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
  ///   \cgalParamNEnd
  ///   \cgalParamNBegin{verbose}
  ///     \cgalParamDescription{If true, the algorithm will print detailed information about its progress.}
  ///     \cgalParamType{bool}
  ///     \cgalParamDefault{false}
  ///   \cgalParamNEnd
  /// \cgalNamedParamsEnd
  ///
  template <class NamedParameters = parameters::Default_named_parameters>
  Variational_medial_axis(const TriangleMesh_& tmesh,
                          const NamedParameters& np = parameters::default_values(),
                          const GT& gt = GT())
      : tmesh_(tmesh)
      , traits_(gt) {

    using parameters::choose_parameter;
    using parameters::get_parameter;

    desired_number_of_spheres_ = choose_parameter(get_parameter(np, internal_np::number_of_spheres), 100);
    lambda_ = choose_parameter(get_parameter(np, internal_np::lambda), FT(0.2));
    max_iteration_ = choose_parameter(get_parameter(np, internal_np::number_of_iterations), 1000);
    verbose_ = choose_parameter(get_parameter(np, internal_np::verbose), false);
    vpm_ = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, tmesh));
#ifndef CGAL_LINKED_WITH_TBB
    static_assert(!std::is_same_v<ConcurrencyTag_, Parallel_tag>, "Parallel_tag is enabled but TBB is unavailable.");
#endif
    init();
  }
  ///
  ///@}

  /**
   *
   * \brief computes a static skeleton based on the Variational Medial Axis Sampling method.
   *
   * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
   *
   * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters", listed below:
   *
   * \cgalNamedParamsBegin
   *   \cgalParamNBegin{number_of_spheres}
   *     \cgalParamDescription{The desired number of medial spheres in the resulting skeleton.}
   *     \cgalParamType{unsigned int}
   *     \cgalParamDefault{100}
   *   \cgalParamNEnd
   *  \cgalParamNBegin{max_iteration_number}
   *     \cgalParamDescription{The maximum number of iterations for the optimization process.}
   *     \cgalParamType{int}
   *     \cgalParamDefault{1000}
   *   \cgalParamNEnd
   *   \cgalParamNBegin{lambda}
   *     \cgalParamDescription{A weight balancing the two energy terms (SQEM and Euclidean). Smaller values encourage
   * the skeleton to extend deeper into local geometric features of the shape.}
   *     \cgalParamType{FT}
   *     \cgalParamDefault{FT(0.2)}
   *     \cgalParamExtra{The range of this parameter is (0,1].}
   *   \cgalParamNEnd
   *   \cgalParamNBegin{verbose}
   *     \cgalParamDescription{If true, the algorithm will print detailed information about its progress on the standard
   * output.}
   *     \cgalParamType{bool}
   *     \cgalParamDefault{false}
   *   \cgalParamNEnd
   * \cgalNamedParamsEnd
   *
   *
   */
  template <class NamedParameters = parameters::Default_named_parameters>
  bool compute_variational_medial_axis_sampling(const NamedParameters& np = parameters::default_values()) {

    using parameters::choose_parameter;
    using parameters::get_parameter;
    // Extract algorithm parameters from named parameters
    desired_number_of_spheres_ =
        choose_parameter(get_parameter(np, internal_np::number_of_spheres), desired_number_of_spheres_);
    lambda_ = choose_parameter(get_parameter(np, internal_np::lambda), lambda_);
    max_iteration_ = choose_parameter(get_parameter(np, internal_np::number_of_iterations), max_iteration_);
    verbose_ = choose_parameter(get_parameter(np, internal_np::verbose), verbose_);
    bool success = false;
    reset_algorithm_state();
    sphere_mesh_->spheres().reserve(desired_number_of_spheres_);
    // Initialize with one sphere
    Sphere_3 init_sphere(Point_3(0., 0., 0.), FT(1.0));
    sphere_mesh_->add_sphere(init_sphere);
    if constexpr(std::is_same_v<ConcurrencyTag_, Parallel_tag>) {
#if CGAL_LINKED_WITH_TBB
      // Compute the shrinking balls in parallel
      compute_shrinking_balls_parallel();
#endif
    } else {
      compute_shrinking_balls();
    }

    if(verbose_) {
      std::cout << "Starting variational medial axis computation..." << std::endl;
      std::cout << "Target number of spheres: " << desired_number_of_spheres_ << std::endl;
      std::cout << "Lambda: " << lambda_ << std::endl;
      std::cout << "Max iterations: " << max_iteration_ << std::endl;
    }
    // Main algorithm loop
    while(iteration_count_ < max_iteration_) {
      bool converged = update_single_step(true);
      if(converged) {
        success = true;
        break;
      }
      ++iteration_count_;
    }

    // Final neighbor update
    update_sphere_neighbors();
    if(verbose_) {
      if(success) {
        std::cout << "Algorithm completed after " << iteration_count_ << " iterations" << std::endl;
        std::cout << "Final number of spheres: " << sphere_mesh_->nb_spheres() << std::endl;
        std::cout << "Final total error: " << total_error_ << std::endl;
      } else {
        std::cout << "Final number of spheres: " << sphere_mesh_->nb_spheres() << std::endl;
        std::cout << "Final total error: " << total_error_ << std::endl;
        std::cout << "Consider decreasing the target number of spheres." << std::endl;
      }
    }
    return success;
  }
  /**
   * \brief updates the medial spheres by performing a single step of the algorithm.
   *
   * This function performs one iteration of the algorithm, updating sphere positions and computing errors.
   * It can optionally enable sphere splitting based on convergence criteria.
   *
   * \param enable_split If true, allows sphere splitting based on convergence criteria.
   * \return True if the algorithm has converged, false otherwise.
   */

  bool update_single_step(bool enable_split = false) {
    // Clean data
    sphere_mesh_->reset();

    if constexpr(std::is_same_v<ConcurrencyTag_, Parallel_tag>) {
#ifdef CGAL_LINKED_WITH_TBB
      // Compute the cluster sphere for each vertex
      assign_vertices_to_clusters_parallel();
      // Update the sphere by optimizing the combined metric
      optimize_sphere_positions_parallel(true);
      total_error_ = compute_sphere_errors_parallel();
#endif
    } else {
      // Compute the cluster sphere for each vertex
      assign_vertices_to_clusters();
      // Update the sphere by optimizing the combined metric
      optimize_sphere_positions(true);
      // Compute error of each sphere
      total_error_ = compute_sphere_errors();
    }
    total_error_diff_ = std::abs(total_error_ - last_total_error_);
    last_total_error_ = total_error_;
    if(verbose_) {
      std::cout << "Iteration " << iteration_count_ << ": spheres=" << sphere_mesh_->nb_spheres()
                << ", error=" << total_error_ << ", error_diff=" << total_error_diff_ << std::endl;
    }
    // Check convergence
    if((sphere_mesh_->nb_spheres() >= desired_number_of_spheres_ && total_error_diff_ < converged_threshold_)) {
      if(verbose_) {
        std::cout << "Converged: reached target number of spheres with low error change" << std::endl;
      }
      return true;
    }

    // Split spheres periodically or when converged
    if(enable_split && (total_error_diff_ < converged_threshold_ || iteration_count_ % 10 == 0)) {
      update_sphere_neighbors();
      split_spheres();
    }

    return false;
  }
  /**
   * \brief performs a specified number of algorithm iterations.
   *
   * This function allows manual control over the algorithm execution,
   * performing only position optimization without sphere splitting.
   *
   * \param nb_iteration Number of iterations to perform
   *
   */
  void update(std::size_t nb_iteration) {
    for(std::size_t i = 0; i < nb_iteration; i++) {
      update_single_step();
    }
    update_sphere_neighbors();
  }

  /**
   * \brief adds spheres by iteratively splitting existing spheres.
   *
   * This function attempts to add the specified number of spheres
   * by running the algorithm with sphere splitting enabled until
   * either the target number is reached or maximum iterations are exceeded.
   *
   * \param nb_sphere Number of spheres to add
   *
   * \pre nb_sphere must be positive
   * \pre At least one sphere must already exist
   * \return True if spheres were added successfully, false otherwise.
   */
  bool add_spheres(int nb_sphere) {
    if(nb_sphere == 0) {
      return false;
    }
    iteration_count_ = 0;
    int max_iteration = std::min(1000, std::max(100, 10 * nb_sphere));
    desired_number_of_spheres_ = sphere_mesh_->nb_spheres() + nb_sphere;
    bool converged = false;
    while(!converged && iteration_count_ < max_iteration) {
      converged = update_single_step(true);
      iteration_count_++;
    }
    if(verbose_) {
      if(converged) {
        std::cout << "Added " << nb_sphere << " spheres successfully." << std::endl;
      } else {
        std::cout << "Failed to add " << nb_sphere << " spheres within the maximum iterations." << std::endl;
      }
    }
    update_sphere_neighbors();
    return converged;
  }

  /**
   * adds a new sphere by splitting sphere with the id `sphere_id`.
   *
   * This function is aimed to be called during interactive sessions, where the use
   * can specify a sphere to split and add a new sphere based on the split vertex.
   * \param sphere_id
   *    The ID of the sphere to split.
   * \param nb_iteration
   *    Number of optimization iterations to perform after adding the sphere (default: 10).
   */
  void add_sphere_by_id(Sphere_ID sphere_id, int nb_iteration = 10) {
    auto sphere = sphere_mesh_->get_sphere(sphere_id);
    vertex_descriptor split_vertex = sphere.get_split_vertex();
    Point_3 center = get(vertex_medial_sphere_pos_map_, split_vertex);
    FT radius = get(vertex_medial_sphere_radius_map_, split_vertex);
    sphere_mesh_->add_sphere(Sphere_3(center, radius * radius));
    // update the spheres for nb_iteration times
    update(nb_iteration);
  }
  /**
   * removes a sphere by its sphere_id.
   * This function is aimed to be called during interactive sessions, where the user
   * can specify a sphere to remove.
   * \param sphere_id
   *    The ID of the sphere to remove.
   * \param nb_iteration
   *    Number of optimization iterations to perform after removing the sphere (default: 10).
   */
  void remove_sphere_by_id(Sphere_ID sphere_id, int nb_iteration = 10) {
    sphere_mesh_->remove(sphere_id);
    // update the spheres
    update(nb_iteration);
  }

  /**
   * \brief exports the medial skeleton as a `Medial_Skeleton` object.
   *
   * This function builds a `Medial_Skeleton` from the current state of the medial sphere mesh.
   * It extracts the vertices, edges, and faces from the medial sphere mesh and constructs
   * the medial skeleton accordingly.
   *
   * \return
   *     A `Medial_Skeleton` object containing the medial skeleton data.
   */
  Medial_Skeleton<TriangleMesh_> export_skeleton() const {
    Medial_Skeleton<TriangleMesh_> skeleton;
    skeleton.build_skeleton_from_medial_sphere_mesh(*sphere_mesh_);
    return skeleton;
  }

  /// \name Parameters
  /// @{
  /**
   * \brief gets the Lambda parameter for the algorithm.
   *
   * This parameter controls the balance between the SQEM and Euclidean energy terms.
   * Smaller values encourage the skeleton to extend deeper into local geometric features of the shape.
   */
  FT lambda_param() const { return lambda_; }

  /**
   * \brief sets function for `lambda_param()`.
   * Note: The lambda must be strictly positive; if set to zero, it will default to 0.2.
   */
  void set_lambda(FT lambda) {
    if(lambda <= FT(0)) {
      std::cerr << "Warning: Lambda must be strictly positive. Setting to default value of 0.2." << std::endl;
      lambda_ = FT(0.2);
    } else {
      lambda_ = lambda;
    }
  }
  /**
   * \brief Get the desired number of spheres for the algorithm.
   */
  int number_of_spheres() const { return desired_number_of_spheres_; }
  /**
   * \brief sets the expected number of spheres.
   */
  void set_number_of_spheres(int num) { desired_number_of_spheres_ = num; }

  /**
   * \brief gets the maximum number of iterations for the algorithm.
   *
   * In case the algorithm does not converge, it will stop after this number of iterations.
   */
  int max_iteration() const { return max_iteration_; }

  /**
   * \brief sets the maximum number of iterations.
   */
  void set_max_iteration(int max_iter) { max_iteration_ = max_iter; }
  ///@}
private:
  /// Initialization that compute some global variable for the algorithm.
  void init() {
    namespace PMP = CGAL::Polygon_mesh_processing;

    // Build AABB-tree
    tree_ = std::make_unique<Tree>(faces(tmesh_).begin(), faces(tmesh_).end(), tmesh_, vpm_);
    if constexpr(std::is_same_v<AccelerationType_, KD_tree_tag>) {
      tree_->accelerate_distance_queries(vertices(tmesh_).begin(), vertices(tmesh_).end(), vpm_);
    } else {
      tree_->accelerate_distance_queries();
    }
    // get bounding box of the mesh
    auto bbox = tree_->bbox();
    scale_ = std::max(bbox.xmax() - bbox.xmin(), std::max(bbox.ymax() - bbox.ymin(), bbox.zmax() - bbox.zmin()));
    converged_threshold_ = scale_ * scale_ * 1e-5;
    // Create property maps
    vertex_normal_map_ = get(Vertex_normal_tag(), tmesh_, Vector_3(0., 0., 0.));
    vertex_area_map_ = get(Vertex_area_tag(), tmesh_, 0.);
    vertex_error_map_ = get(Vertex_error_tag(), tmesh_, 0.);
    vertex_cluster_sphere_map_ = get(Vertex_cluster_sphere_tag(), tmesh_, MSMesh::INVALID_SPHERE_ID);
    vertex_medial_sphere_pos_map_ = get(Vertex_medial_sphere_pos_tag(), tmesh_, Point_3(0., 0., 0.));
    vertex_medial_sphere_radius_map_ = get(Vertex_medial_sphere_radius_tag(), tmesh_, FT(0.));
    face_normal_map_ = get(Face_normal_tag(), tmesh_, Vector_3(0., 0., 0.));
    face_area_map_ = get(Face_area_tag(), tmesh_, 0.);

    // Compute normals
    PMP::compute_vertex_normals(tmesh_, vertex_normal_map_);
    PMP::compute_face_normals(tmesh_, face_normal_map_);

    // Compute vertex areas
    for(face_descriptor f : faces(tmesh_)) {
      double area = PMP::face_area(f, tmesh_);
      put(face_area_map_, f, area);
      for(vertex_descriptor v : vertices_around_face(halfedge(f, tmesh_), tmesh_)) {
        put(vertex_area_map_, v, get(vertex_area_map_, v) + area / 3.0);
      }
    }
    sphere_mesh_ = std::make_unique<MSMesh>();

    // Algorithm variables
    iteration_count_ = 0;
    total_error_ = FT(0.0);
    total_error_diff_ = (std::numeric_limits<FT>::max)();
    last_total_error_ = total_error_;
  }

  inline FT cosine_angle(const Vector_3& v1, const Vector_3& v2) {
    FT norm_v1v2 = CGAL::approximate_sqrt(v1.squared_length() * v2.squared_length());
    FT res = norm_v1v2 > 1e-20 ? (v1 * v2) / norm_v1v2 : FT(1);
    return std::clamp(res, FT(-1), FT(1)); // Ensure the result is within [-1, 1]
  }

  inline FT compute_radius(const Point_3& p, const Vector_3& n, const Point_3& q) {
    Vector_3 qp = p - q;
    FT d = CGAL::approximate_sqrt(qp.squared_length());
    FT cos_angle = cosine_angle(qp, n);
    return d / (2 * cos_angle);
  }
  std::pair<Point_3, FT>
  shrinking_ball_algorithm_bvh(std::vector<face_descriptor> incident_faces,
                               const Point_3& p,                  // point on the surface
                               const Vector_3& n,                 // inverse of search direction
                               FT delta_convergence = FT(1e-5)) { // model has to be normalized in [0, 1]^3
    using face_descriptor = typename Tree::Primitive_id;
    auto approximate_angle = GT().compute_approximate_angle_3_object();
    delta_convergence *= scale_;
    const FT denoise_preserve = FT(20.0); // in degree
    const int iteration_limit = 30;
    int j = 0;
    FT r = FT(1.0) * scale_; // initial radius
    Point_3 c = p - (r * n);
    face_descriptor last_face;
    while(true) {
      // get the hint directly from the kd-tree only to illustrate the use of the internal kd-tree
      // auto hint = tree_->kd_tree().closest_point(c);
      // auto [q_next, closest_face] = tree_->closest_point_and_primitive(c, hint);
      auto [q_next, closest_face] = tree_->closest_point_and_primitive(c);
      FT squared_dist = (q_next - c).squared_length();
      if(squared_dist >= (r - delta_convergence) * (r - delta_convergence)) {
        break; // convergence
      }
      bool should_break = false;
      for(auto f : incident_faces) {
        if(f == closest_face) {
          should_break = true;
          break; // found the face
        }
      }
      if(should_break) {
        break; // closest face is incident to the vertex
      }
      if(j > 0 && closest_face == last_face) {
        break; // no change in closest face
      }
      FT r_next = compute_radius(p, n, q_next);
      if(!CGAL::is_finite(r_next) || r_next <= FT(0)) {
        std::cerr << "Invalid radius at iteration " << j << ": " << r_next << std::endl;
        break;
      }
      Point_3 c_next = p - (r_next * n);
      FT separation_angle = approximate_angle(p - c_next, q_next - c_next);
      if(j > 0 && separation_angle < denoise_preserve) {
        break; // denoise preserve angle achieved
      }
      c = c_next;
      r = r_next;
      last_face = closest_face;
      j++;
      if(j > iteration_limit)
        break;
    }

    return {c, r};
  }
  std::pair<Point_3, FT> shrinking_ball_algorithm_kdt(const Point_3& p,  // point on the surface
                                                      const Vector_3& n, // inverse of search direction
                                                      FT delta_convergence = FT(1e-5)) {

    auto approximate_angle = GT().compute_approximate_angle_3_object();
    delta_convergence *= scale_;
    const FT denoise_preserve = FT(20.0); // in degree
    const int iteration_limit = 30;
    int j = 0;
    FT r = FT(0.25) * scale_; // initial radius
    Point_3 c = p - (r * n);

    while(true) {

      auto hint = tree_->kd_tree().closest_point(c);
      Point_3 q_next = hint.first;

      FT squared_dist = (q_next - c).squared_length();
      if(squared_dist >= (r - delta_convergence) * (r - delta_convergence) || p == q_next) {
        // std::cout << "Convergence achieved at iteration " << j << ": " << CGAL::to_double(squared_dist) << std::endl;
        break; // convergence
      }
      FT r_next = compute_radius(p, n, q_next);
      if(!CGAL::is_finite(r_next) || r_next <= FT(0)) {
        std::cerr << "Invalid radius at iteration " << j << ": " << r_next << std::endl;
        break;
      }
      Point_3 c_next = p - (r_next * n);
      FT seperation_angle = approximate_angle(p - c_next, q_next - c_next);
      if(j > 0 && seperation_angle < denoise_preserve) {
        /* std::cout << "Denoise preserve angle achieved at iteration: " << j
                  << ", Angle: " << CGAL::to_double(seperation_angle) << " degrees, Center: " << c_next
                  << ", Radius: " << r_next << std::endl;*/
        break;
      }
      c = c_next;
      r = r_next;
      j++;
      if(j > iteration_limit)
        break;
    }

    return {c, r};
  }

  std::pair<Point_3, FT> compute_shrinking_ball_impl(const Point_3& p, const Vector_3& normal, CGAL::KD_tree_tag) {

    return shrinking_ball_algorithm_kdt(p, normal);
  }

  std::pair<Point_3, FT> compute_shrinking_ball_impl(const Point_3& p, const Vector_3& normal, CGAL::BVH_tag) {
    auto [q, closest_face] = tree_->closest_point_and_primitive(p);
    auto face_range = CGAL::faces_around_face(halfedge(closest_face, tmesh_), tmesh_);
    std::vector<face_descriptor> incident_faces(face_range.begin(), face_range.end());
    return shrinking_ball_algorithm_bvh(incident_faces, p, normal);
  }
  void compute_one_vertex_shrinking_ball(vertex_descriptor v) {

    Vector_3 normal = get(vertex_normal_map_, v);
    Point_3 p = get(vpm_, v);
    auto [center, radius] = compute_shrinking_ball_impl(p, normal, AccelerationType_{});
    put(vertex_medial_sphere_pos_map_, v, center);
    put(vertex_medial_sphere_radius_map_, v, radius);
  }

  void compute_shrinking_balls() {
    for(auto v : vertices(tmesh_)) {
      compute_one_vertex_shrinking_ball(v);
    }
  }

  void assign_vertices_to_clusters() {
    for(vertex_descriptor v : vertices(tmesh_)) {

      Point_3 p = get(vpm_, v);
      FT min_distance = (std::numeric_limits<FT>::max)();
      Vector_3 normal = get(vertex_normal_map_, v);
      Sphere_ID closest_sphere_id = 0;

      // Find the sphere with smallest distance to the vertex
      for(auto& sphere : sphere_mesh_->spheres()) {
        Point_3 center = sphere.get_center();
        FT radius = sphere.get_radius();

        // compute euclidean distance
        FT dist_eucl = CGAL::approximate_sqrt((p - center).squared_length()) - radius;
        dist_eucl *= dist_eucl;

        // compute sqem distance
        FT dist_sqem = CGAL::scalar_product(p - center, normal) - radius;
        dist_sqem *= dist_sqem;

        FT distance = dist_sqem + lambda_ * dist_eucl;
        if(distance < min_distance) {
          min_distance = distance;
          closest_sphere_id = sphere.get_id();
        }
      }
      FT area = get(vertex_area_map_, v);
      // Update the closest sphere
      put(vertex_cluster_sphere_map_, v, closest_sphere_id);
      sphere_mesh_->get_sphere(closest_sphere_id).accumulate_cluster_area(area);
      sphere_mesh_->get_sphere(closest_sphere_id).add_cluster_vertex(v);
    }
    std::vector<Sphere_ID> sphere_ids_to_remove;
    for(auto& sphere : sphere_mesh_->spheres()) {
      auto& cluster_vertices = sphere.get_cluster_vertices();
      if(cluster_vertices.size() <= 4) {
        for(vertex_descriptor v : cluster_vertices) {
          put(vertex_cluster_sphere_map_, v, MSMesh::INVALID_SPHERE_ID);
        }
        sphere_ids_to_remove.push_back(sphere.get_id());
      }
    }
    for(Sphere_ID id : sphere_ids_to_remove) {
      if(verbose_) {
        std::cout << "Removing sphere with ID: " << id << " due to small cluster size." << std::endl;
      }
      sphere_mesh_->remove(id); // remove spheres with small clusters
    }
  }

  void correct_sphere(MSphere& sphere, const Eigen::Matrix<FT, 4, 1>& optimized_sphere_params) {

    Side_of_triangle_mesh<TriangleMesh_, GT, VPM, Tree> side_of(*tree_, traits_);
    Point_3 optimal_center(optimized_sphere_params(0), optimized_sphere_params(1), optimized_sphere_params(2));
    auto [cp, closest_face] = tree_->closest_point_and_primitive(optimal_center);
    FT len = (optimal_center - cp).squared_length();

    Vector_3 normal = (cp - optimal_center) / CGAL::approximate_sqrt(len);

    if((side_of(optimal_center) == CGAL::ON_UNBOUNDED_SIDE))
      normal = -normal; // if the center is outside, flip the normal

    auto [c, r] = compute_shrinking_ball_impl(cp, normal, AccelerationType_{});

    sphere.set_center(c);
    sphere.set_radius(r);
  }

  void optimize_single_sphere(MSphere& sphere, bool use_shrinking_ball_correction = false) {
    using EMat = Eigen::Matrix<FT, Eigen::Dynamic, Eigen::Dynamic>;
    using EVec = Eigen::Matrix<FT, Eigen::Dynamic, 1>;
    using EVec3 = Eigen::Matrix<FT, 3, 1>;
    using EVec4 = Eigen::Matrix<FT, 4, 1>;
    using LDLTSolver = Eigen::LDLT<EMat>;

    auto& cluster_vertices = sphere.get_cluster_vertices();
    Point_3 center = sphere.get_center();
    FT radius = sphere.get_radius();

    auto nrows = 2 * cluster_vertices.size();
    EMat J(nrows, 4);
    J.setZero();
    EVec b(nrows);
    b.setZero();
    int idx = 0;
    EVec4 s;
    s << center.x(), center.y(), center.z(), radius;

    for(int i = 0; i < 10; i++) {
      idx = 0;

      for(vertex_descriptor v : cluster_vertices) {
        Point_3 p = get(vpm_, v);
        EVec3 pos(p.x(), p.y(), p.z());

        // compute sqem energy
        EVec4 lhs = EVec4::Zero();
        FT rhs = 0.0;
        for(face_descriptor f : faces_around_target(halfedge(v, tmesh_), tmesh_)) {
          Vector_3 normal = get(face_normal_map_, f);
          EVec3 normal_eigen(normal.x(), normal.y(), normal.z());
          EVec4 n4(normal_eigen(0), normal_eigen(1), normal_eigen(2), 1.0);
          FT area = CGAL::approximate_sqrt(get(face_area_map_, f) / 3.0);
          lhs += -n4 * area;
          rhs += -1.0 * ((pos - EVec3(s(0), s(1), s(2))).dot(normal_eigen) - s(3)) * area;
        }
        J.row(idx) = lhs;
        b(idx) = rhs;
        ++idx;

        // compute euclidean energy
        EVec3 d = pos - EVec3(s(0), s(1), s(2));
        FT l = d.norm();
        FT area = CGAL::approximate_sqrt(get(vertex_area_map_, v));
        J.row(idx) = EVec4(-(d[0] / l), -(d[1] / l), -(d[2] / l), -1.0) * area * lambda_;
        b(idx) = -(l - s(3)) * area * lambda_;
        ++idx;
      }

      LDLTSolver solver(J.transpose() * J);
      EVec4 delta_s = solver.solve(J.transpose() * b);
      s += delta_s;
      if(delta_s.norm() < 1e-8) {
        // std::cout << "Convergence achieved after " << i + 1 << " iterations." << std::endl;
        break; // convergence
      }
      sphere.set_do_not_split(false); // reset the split flag
    }
    if(use_shrinking_ball_correction) {
      correct_sphere(sphere, s);
    } else {
      // Update sphere with optimized parameters
      sphere.set_center(Point_3(s(0), s(1), s(2)));
      sphere.set_radius(s(3));
    }
  }

  FT compute_single_sphere_error(MSphere& sphere) {
    auto& cluster_vertices = sphere.get_cluster_vertices();
    Point_3 center = sphere.get_center();
    FT radius = sphere.get_radius();
    FT max_dist = (std::numeric_limits<FT>::min)();
    FT error = 0.0;
    FT area = sphere.get_area();
    if(area <= FT(0.0)) {
      std::cerr << "Warning: sphere with zero area encountered, skipping error computation." << std::endl;
      return FT(0.0);
    }

    for(vertex_descriptor v : cluster_vertices) {
      Point_3 p = get(vpm_, v);
      FT area = get(vertex_area_map_, v);
      FT dist_sqem = CGAL::scalar_product(p - center, get(vertex_normal_map_, v)) - radius;
      dist_sqem *= dist_sqem; // square the distance

      FT dist_eucl = CGAL::approximate_sqrt((p - center).squared_length()) - radius;
      dist_eucl *= dist_eucl; // square the distance

      FT distance = area * (dist_sqem + lambda_ * dist_eucl);
      error += distance;

      put(vertex_error_map_, v, distance);
      if(distance > max_dist) {
        max_dist = distance;
        sphere.set_split_vertex(v); // update the split vertex
      }
    }
    error = error / sphere.get_area(); // average error
    sphere.set_error(error);

    return error;
  }

#ifdef CGAL_LINKED_WITH_TBB

  void assign_vertices_to_clusters_parallel() {
    std::vector<vertex_descriptor> vertices_vector;
    vertices_vector.reserve(num_vertices(tmesh_));
    for(vertex_descriptor v : vertices(tmesh_)) {
      vertices_vector.push_back(v);
    }

    std::vector<std::pair<Sphere_ID, FT>> vertex_assignments(vertices_vector.size());

    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, vertices_vector.size()),
                      [&](const tbb::blocked_range<std::size_t>& range) {
                        for(std::size_t i = range.begin(); i != range.end(); ++i) {

                          vertex_descriptor v = vertices_vector[i];
                          FT area = get(vertex_area_map_, v);
                          Point_3 p = get(vpm_, v);
                          FT min_distance = (std::numeric_limits<FT>::max)();
                          Vector_3 normal = get(vertex_normal_map_, v);
                          Sphere_ID closest_sphere_id;

                          // Find the sphere with smallest distance to the vertex
                          for(const auto& sphere : sphere_mesh_->spheres()) {
                            Point_3 center = sphere.get_center();
                            FT radius = sphere.get_radius();

                            // compute euclidean distance
                            FT dist_eucl = CGAL::approximate_sqrt((p - center).squared_length()) - radius;
                            dist_eucl *= dist_eucl;

                            // compute sqem distance
                            FT dist_sqem = CGAL::scalar_product(p - center, normal) - radius;
                            dist_sqem *= dist_sqem;

                            FT distance = area * (dist_sqem + lambda_ * dist_eucl);
                            if(distance < min_distance) {
                              min_distance = distance;
                              closest_sphere_id = sphere.get_id();
                            }
                          }
                          vertex_assignments[i] = {closest_sphere_id, area};
                        }
                      });
    for(std::size_t i = 0; i < vertices_vector.size(); ++i) {
      vertex_descriptor v = vertices_vector[i];
      auto& [closest_sphere_id, area] = vertex_assignments[i];

      if(closest_sphere_id != MSMesh::INVALID_SPHERE_ID) {
        put(vertex_cluster_sphere_map_, v, closest_sphere_id);
        sphere_mesh_->get_sphere(closest_sphere_id).accumulate_cluster_area(area);
        sphere_mesh_->get_sphere(closest_sphere_id).add_cluster_vertex(v);
      }
    }
    std::vector<Sphere_ID> sphere_ids_to_remove;
    for(auto& sphere : sphere_mesh_->spheres()) {
      auto& cluster_vertices = sphere.get_cluster_vertices();
      if(cluster_vertices.size() <= 4) {
        for(vertex_descriptor v : cluster_vertices) {
          put(vertex_cluster_sphere_map_, v, MSMesh::INVALID_SPHERE_ID);
        }
        sphere_ids_to_remove.push_back(sphere.get_id());
      }
    }
    for(Sphere_ID id : sphere_ids_to_remove) {
      if(verbose_) {
        std::cout << "Removing sphere with ID: " << id << " due to small cluster size." << std::endl;
      }
      sphere_mesh_->remove(id); // remove spheres with small clusters
    }
  }

  void compute_shrinking_balls_parallel() {
    std::vector<vertex_descriptor> vertices_vector;
    vertices_vector.reserve(num_vertices(tmesh_));

    for(vertex_descriptor v : vertices(tmesh_)) {
      vertices_vector.push_back(v);
    }
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, vertices_vector.size()),
                      [&](const tbb::blocked_range<std::size_t>& range) {
                        for(std::size_t i = range.begin(); i != range.end(); i++) {
                          compute_one_vertex_shrinking_ball(vertices_vector[i]);
                        }
                      });
  }
  void optimize_sphere_positions_parallel(bool use_shrinking_ball_correction = false) {
    auto& spheres = sphere_mesh_->spheres();
    std::vector<std::reference_wrapper<MSphere>> sphere_refs(spheres.begin(), spheres.end());

    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, sphere_refs.size()),
                      [&](const tbb::blocked_range<std::size_t>& range) {
                        for(std::size_t i = range.begin(); i != range.end(); i++) {
                          optimize_single_sphere(sphere_refs[i], use_shrinking_ball_correction);
                        }
                      });
  }

  FT compute_sphere_errors_parallel() {
    auto& spheres = sphere_mesh_->spheres();
    std::vector<std::reference_wrapper<MSphere>> sphere_refs(spheres.begin(), spheres.end());
    FT total_error = tbb::parallel_reduce(
        tbb::blocked_range<std::size_t>(0, sphere_refs.size()), FT(0.0),
        [&](const tbb::blocked_range<std::size_t>& range, FT local_sum) -> FT {
          for(std::size_t i = range.begin(); i != range.end(); ++i) {
            FT sphere_error = compute_single_sphere_error(sphere_refs[i]);
            local_sum += sphere_error;
          }
          return local_sum;
        },
        std::plus<FT>());
    return total_error;
  }
#endif

  void optimize_sphere_positions(bool use_shrinking_ball_correction = false) {
    for(auto& sphere : sphere_mesh_->spheres()) {
      optimize_single_sphere(sphere, use_shrinking_ball_correction);
    }
  }

  FT compute_sphere_errors() {
    FT total_error = FT(0.0);
    for(auto& sphere : sphere_mesh_->spheres()) {
      FT error = compute_single_sphere_error(sphere);
      total_error += error;
    }
    return total_error;
  }

  void update_sphere_neighbors() {
    for(edge_descriptor e : edges(tmesh_)) {
      vertex_descriptor v1 = source(e, tmesh_);
      vertex_descriptor v2 = target(e, tmesh_);
      Sphere_ID s1 = get(vertex_cluster_sphere_map_, v1);
      Sphere_ID s2 = get(vertex_cluster_sphere_map_, v2);
      if(s1 != s2 && s1 != MSMesh::INVALID_SPHERE_ID && s2 != MSMesh::INVALID_SPHERE_ID) {
        auto& sphere1 = sphere_mesh_->get_sphere(s1);
        auto& sphere2 = sphere_mesh_->get_sphere(s2);

        if(sphere1.get_neighbors().find(s2) == sphere1.get_neighbors().end()) {
          sphere1.add_neighbor(s2);
          sphere2.add_neighbor(s1);
        }
      }
    }
  }

  void split_spheres() {
    if(sphere_mesh_->nb_spheres() >= desired_number_of_spheres_) {
      return;
    }
    if(verbose_) {
      std::cout << "Start splitting spheres" << std::endl;
    }
    std::vector<Sphere_ID> sorted_sphere_ids;
    sorted_sphere_ids.reserve(sphere_mesh_->nb_spheres());
    for(const auto& sphere : sphere_mesh_->spheres()) {
      sorted_sphere_ids.push_back(sphere.get_id());
    }
    std::sort(sorted_sphere_ids.begin(), sorted_sphere_ids.end(), [&](Sphere_ID a, Sphere_ID b) {
      return sphere_mesh_->get_sphere(a).get_error() > sphere_mesh_->get_sphere(b).get_error(); // sort by error
    });

    int to_split_max = std::min(int(std::ceil(sphere_mesh_->nb_spheres() * 0.2)), 10);
    for(auto& sphere_id : sorted_sphere_ids) {
      if(sphere_mesh_->nb_spheres() >= desired_number_of_spheres_)
        break;
      auto& sphere = sphere_mesh_->get_sphere(sphere_id);
      if(to_split_max > 0 && sphere.can_split()) {
        for(Sphere_ID neighbor_id : sphere.get_neighbors()) {
          sphere_mesh_->get_sphere(neighbor_id).set_do_not_split(true);
        }
        vertex_descriptor split_vertex = sphere.get_split_vertex();
        Point_3 center = get(vertex_medial_sphere_pos_map_, split_vertex);
        FT radius = get(vertex_medial_sphere_radius_map_, split_vertex);
        sphere_mesh_->add_sphere(Sphere_3(center, radius * radius));
        to_split_max--;
      }
    }
  }
  void reset_algorithm_state() {
    iteration_count_ = 0;
    total_error_ = FT(0.0);
    total_error_diff_ = (std::numeric_limits<FT>::max)();
    last_total_error_ = total_error_;
    sphere_mesh_->reset();
  }

private:
  const TriangleMesh_& tmesh_;

  GT traits_;
  VPM vpm_;
  FT lambda_;
  std::size_t desired_number_of_spheres_;
  int max_iteration_;
  int iteration_count_;
  FT total_error_diff_;
  FT last_total_error_;
  FT total_error_;
  FT converged_threshold_;
  std::unique_ptr<Tree> tree_;
  std::unique_ptr<MSMesh> sphere_mesh_;
  FT scale_;
  bool verbose_;
  // Property maps
  Vertex_normal_map vertex_normal_map_;
  Vertex_area_map vertex_area_map_;
  Vertex_error_map vertex_error_map_;
  Vertex_cluster_sphere_map vertex_cluster_sphere_map_;
  Vertex_medial_sphere_pos_map vertex_medial_sphere_pos_map_;
  Vertex_medial_sphere_radius_map vertex_medial_sphere_radius_map_;
  Face_normal_map face_normal_map_;
  Face_area_map face_area_map_;
};

namespace IO {
/**
 * @ingroup PkgVMASRef
 * @brief writes the medial skeleton to a PLY file.
 *
 * @tparam TriangleMesh
 *         a model of `FaceListGraph`
 * @tparam GeomTraits
 *         a model of `Kernel`
 * @tparam NamedParameters
 *         a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param skeleton The medial skeleton to write.
 * @param filepath The name of the file to write to.
 * @param np optional \ref bgl_namedparameters "Named Parameters" described below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{stream_precision}
 *   \cgalParamDescription{a parameter used to set the precision (i.e. how many digits are generated) of the output
 * stream}
 *   \cgalParamType{int}
 *   \cgalParamDefault{the precision of the stream `os`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * Note:  This function writes the medial skeleton to a PLY file in ASCII format. The format is :
 * ```
 * ply
 * format ascii 1.0
 * element vertex N
 * property float x
 * property float y
 * property float z
 * property float radius
 * element edge M
 * property int vertex1
 * property int vertex2
 * element face K
 * property list uchar int vertex_indices
 * end_header
 * x1 y1 z1 r1
 * ...        //N vertices
 * xN yN zN rN
 * v1 v2
 * ...        //M edges
 * vM vN
 * 3 v1 v2 v3
 * ...        //K faces
 * 3 vX vY vZ
 * ```
 *
 * @returns `true` if writing was successful, `false` otherwise.
 */
template <typename TriangleMesh, typename GeomTraits, class NamedParameters = parameters::Default_named_parameters>
bool write_PLY(const Medial_Skeleton<TriangleMesh, GeomTraits>& skeleton,
               const std::string& filepath,
               const NamedParameters& np = parameters::default_values()) {
  std::ofstream ofs(filepath);
  set_stream_precision_from_NP(ofs, np);
  if(!ofs) {
    std::cerr << "Error opening file: " << filepath << std::endl;
    return false;
  }

  // Write header
  ofs << "ply\nformat ascii 1.0\n";
  ofs << "element vertex " << skeleton.number_of_vertices() << "\n";
  ofs << "property float x\nproperty float y\nproperty float z\n";
  ofs << "property float radius\n";
  ofs << "element edge " << skeleton.number_of_edges() << "\n";
  ofs << "property int vertex1\nproperty int vertex2\n";
  ofs << "element face " << skeleton.number_of_faces() << "\n";
  ofs << "property list uchar int vertex_indices\n";
  ofs << "end_header\n";

  // Write vertices (sphere centers and radii)
  for(const auto& sphere : skeleton.vertices()) {
    const auto& center = sphere.center();
    auto radius = CGAL::sqrt(sphere.squared_radius());
    ofs << CGAL::to_double(center.x()) << " " << CGAL::to_double(center.y()) << " " << CGAL::to_double(center.z())
        << " " << CGAL::to_double(radius) << "\n";
  }

  // Write edges
  for(const auto& e : skeleton.edges()) {
    ofs << e.first << " " << e.second << "\n";
  }

  // Write faces
  for(const auto& f : skeleton.faces()) {
    ofs << "3 " << f[0] << " " << f[1] << " " << f[2] << "\n";
  }
  return true;
}
} // namespace IO

} // namespace CGAL
#endif
