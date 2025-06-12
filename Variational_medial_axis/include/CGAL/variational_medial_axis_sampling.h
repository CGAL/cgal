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

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <Eigen/Dense>

namespace CGAL {

template <typename TriangleMesh, typename GT> 
class MedialSphereMesh
{
public:
  using FT = typename GT::FT;
  using Point_3 = typename GT::Point_3;
  using Sphere_3 = typename GT::Sphere_3;
  using Sphere_ID = std::size_t;
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
  static const Sphere_ID INVALID_SPHERE_ID = std::numeric_limits<Sphere_ID>::max();
  struct MedialSphere
  {
    Sphere_ID id;
    Sphere_3 sphere;
    FT error;
    FT cluster_area;
    bool do_not_split = false;      // flag to indicate if the sphere should not be split
    vertex_descriptor split_vertex; // the vertex for split
    std::unordered_set<Sphere_ID> neighbors;
    std::vector<vertex_descriptor> cluster_vertices;

    MedialSphere(const Sphere_3& s, Sphere_ID i)
        : id(i)
        , sphere(s)
        , split_vertex(boost::graph_traits<TriangleMesh>::null_vertex())
        , error(FT(0))
        , cluster_area(FT(0)) {}

    void reset() {
      error = FT(0);
      split_vertex = boost::graph_traits<TriangleMesh>::null_vertex();
      cluster_area = FT(0);
      neighbors.clear();
      cluster_vertices.clear();
    }

    FT get_radius() const { return CGAL::approximate_sqrt(sphere.squared_radius()); }
    Point_3 get_center() const { return sphere.center(); }
    FT get_area() const { return cluster_area; }
    Sphere_ID get_id() const { return id; }
    void set_center(const Point_3& p) { sphere = Sphere_3(p, get_radius()); }
    void set_radius(FT r) { sphere = Sphere_3(get_center(), r * r); }
    void set_cluter_area(FT area) { cluster_area = area; }
    void accumulate_cluster_area(FT area) { cluster_area += area; }
    void add_neighbor(Sphere_ID id) { neighbors.insert(id); }
    bool can_split() const { return !do_not_split && split_vertex != boost::graph_traits<TriangleMesh>::null_vertex(); }
    void set_do_not_split(bool value) { do_not_split = value; }
    std::vector<vertex_descriptor>& get_cluster_vertices() { return cluster_vertices; }
    void add_cluster_vertex(vertex_descriptor v) { cluster_vertices.push_back(v); }
  };

private:
  std::vector<std::shared_ptr<MedialSphere>> spheres_; // TO DO: use Compact container?
  std::unordered_map<Sphere_ID, std::size_t> id_to_index_;
  Sphere_ID next_id_ = 0;

public:
  Sphere_ID add_sphere(const Sphere_3 s) {
    Sphere_ID id = next_id_++;
    spheres_.push_back(std::make_shared<MedialSphere>(std::move(s), id));
    id_to_index_[id] = spheres_.size() - 1;
    return id;
  
  }


  void remove(Sphere_ID id) {
    auto it = id_to_index_.find(id);
    if(it == id_to_index_.end())
      return;
    std::size_t index_to_remove = it->second;
    std::size_t last_index = spheres_.size() - 1;

    if(index_to_remove != last_index) {
      std::swap(spheres_[index_to_remove], spheres_[last_index]);
      Sphere_ID swapped_id = spheres_[index_to_remove]->get_id();
      id_to_index_[swapped_id] = index_to_remove;
    }
    spheres_.pop_back();
    id_to_index_.erase(it);
  }
  void reset() {
    for(auto& sphere : spheres_) {
      sphere->reset();
    }
  }
  std::size_t nb_spheres() { return id_to_index_.size(); }
  std::shared_ptr<MedialSphere> get_sphere(Sphere_ID id) { return spheres_[id_to_index_.at(id)]; }
  std::shared_ptr<const MedialSphere> get_sphere(Sphere_ID id) const { return spheres_[id_to_index_.at(id)]; }
  std::vector<std::shared_ptr<MedialSphere>>& spheres() { return spheres_; }
  
  void export_to_ply_with_edges_faces(std::ostream& os) const {
    using Point_3 = typename GT::Point_3;

    std::unordered_map<Sphere_ID, std::size_t> id_to_vertex_index;
    std::vector<Point_3> vertices;
    std::vector<FT> radii;
    std::vector<std::pair<std::size_t, std::size_t>> edges;
    std::vector<std::array<std::size_t, 3>> faces;

    // Assign indices and store positions/radii
    std::size_t vid = 0;
    for(const auto& sphere : spheres_) {
      id_to_vertex_index[sphere->get_id()] = vid++;
      vertices.push_back(sphere->get_center());
      radii.push_back(sphere->get_radius());

    }

    // Collect edges (undirected)
    std::set<std::pair<std::size_t, std::size_t>> edge_set;
    for(const auto& sphere_ptr : spheres_) {
      if(!sphere_ptr)
        continue;
      std::size_t a = id_to_vertex_index.at(sphere_ptr->get_id());
      for(Sphere_ID n_id : sphere_ptr->neighbors) {
        std::size_t b = id_to_vertex_index.at(n_id);
        if(a < b)
          edge_set.emplace(a, b);
      }
    }
    edges.assign(edge_set.begin(), edge_set.end());

    // Collect triangle faces (3 mutually adjacent spheres)
    std::set<std::array<std::size_t, 3>> face_set;
    for(const auto& sphere_ptr : spheres_) {
      if(!sphere_ptr)
        continue;
      Sphere_ID a_id = sphere_ptr->get_id();
      const auto& neighbors_a = sphere_ptr->neighbors;

      for(Sphere_ID b_id : neighbors_a) {
        if(b_id <= a_id)
          continue;
        const auto& neighbors_b = get_sphere(b_id)->neighbors;

        for(Sphere_ID c_id : neighbors_a) {
          if(c_id <= b_id)
            continue;
          if(neighbors_b.find(c_id) != neighbors_b.end()) {
            std::array<std::size_t, 3> tri = {id_to_vertex_index[a_id], id_to_vertex_index[b_id],
                                              id_to_vertex_index[c_id]};
            std::sort(tri.begin(), tri.end());
            face_set.insert(tri);
          }
        }
      }
    }
    faces.assign(face_set.begin(), face_set.end());

    // Write header
    os << "ply\nformat ascii 1.0\n";
    os << "element vertex " << vertices.size() << "\n";
    os << "property float x\nproperty float y\nproperty float z\n";
    os << "property float radius\n";
    os << "element edge " << edges.size() << "\n";
    os << "property int vertex1\nproperty int vertex2\n";
    os << "element face " << faces.size() << "\n";
    os << "property list uchar int vertex_indices\n";
    os << "end_header\n";

    // Write vertices
    for(std::size_t i = 0; i < vertices.size(); ++i) {
      const auto& p = vertices[i];
      os << p.x() << " " << p.y() << " " << p.z() << " " << CGAL::to_double(radii[i]) << "\n";
    }

   // Write edges
    for(const auto& e : edges)
      os << e.first << " " << e.second << "\n";

    // Write faces
   for(const auto& f : faces)
      os << "3 " << f[0] << " " << f[1] << " " << f[2] << "\n";
      

  }

  void write_to_ply_file(const std::string& filename) const {
    std::ofstream ofs(filename);
    if(!ofs) {
      std::cerr << "Error opening file: " << filename << std::endl;
      return;
    }
    export_to_ply_with_edges_faces(ofs);
    ofs.close();
  }
};



template <typename GT>
inline typename GT::FT cosine_angle(const typename GT::Vector_3& v1, const typename GT::Vector_3& v2) {
  GT::FT norm_v1v2 = CGAL::approximate_sqrt(v1.squared_length() * v2.squared_length());
  GT::FT res = norm_v1v2 > 1e-20 ? (v1 * v2) / norm_v1v2 : GT::FT(1);
  return std::clamp(res, GT::FT(-1), GT::FT(1)); // Ensure the result is within [-1, 1]
}


template <typename GT>
inline typename GT::FT
compute_radius(const typename GT::Point_3& p, const typename GT::Vector_3& n, const typename GT::Point_3& q) {
  using Vector_3 = typename GT::Vector_3;
  using FT = typename GT::FT;
  Vector_3 qp = p - q;
  FT d = CGAL::approximate_sqrt(qp.squared_length());
  FT cos_angle = cosine_angle<GT>(qp, n);
  return d / (2 * cos_angle);
}

template <typename TriangleMesh, typename GT, typename Tree, typename FaceNormalMap>
std::pair<typename GT::Point_3, typename GT::FT>
shrinking_ball_algorithm(const TriangleMesh& tmesh,
                         const FaceNormalMap& face_normal_map, // face normal map
                         const typename GT::Point_3& p,  // point on the surface
                         const typename GT::Vector_3& n, // inverse of reaserch direction
                         const Tree& tree) {
  using Vector_3 = typename GT::Vector_3;
  using face_descriptor = typename boost::graph_traits<TriangleMesh>::face_descriptor;

  using Point_3 = typename GT::Point_3;
  using FT = typename GT::FT;

  const FT denoise_preserve = FT(40.0);//in degree
  const FT delta_convergence = FT(1e-5);
  const zero_n3 = 1e-3;
  const zero_n4 = 1e-4; 
  const int iteration_limit = 30;
  int j = 0;
  FT r = FT(1.0); // initial radius
  Point_3 c = p - (r * n);
  Point_3 q = p - (2 * r * n);

  while(true) {
    //std::cout << "Iteration: " << j << ", Current center: " << c << ", Radius: " << r << std::endl;
    auto res = tree.closest_point_and_primitive(c);
    Point_3 q_next = res.first;
    face_descriptor closest_face = res.second;
    
    
    /* FT dist_cq = (c - q_next).squared_length();
    FT dist_pq = CGAL::approximate_sqrt((p - q_next).squared_length());
    FT eps = (r-denoise_preserve)*(r-denoise_preserve);
    if(eps < zero_n3) {
      //std::cout << "Epsilon is too small at iteration: " << j << ", Epsilon: " << eps << std::endl;
      break;
    }
    if(dist_cq >= eps || dist_pq < zero_n3) {
      //std::cout << "Distance is too small at iteration: " << j << ", Distance: " << dist_cq << std::endl;
      break;
    }*/
    FT dist = CGAL::approximate_sqrt((q - q_next).squared_length());
     if(std::abs(dist - r) <= delta_convergence ||
       (CGAL::approximate_sqrt((q - q_next).squared_length()) < delta_convergence))
    {
      // convergence
      //std::cout << "Delta Convergence achieved at iteration: " << j << ", Center: " << c << ", Radius: " << r << std::endl;

      break;
    }
    Vector_3 closest_face_normal = get(face_normal_map, closest_face);
    FT dist_to_face = std::abs(CGAL::scalar_product(q_next - c, closest_face_normal));
    FT ratio = std::abs(r - dist_to_face) / r;
    if(ratio <= FT(0.01)) {
      //std::cout << "Distance to face is too small at iteration: " << j << ", Ratio: " << ratio << std::endl;
      break;
    }
    FT r_next = compute_radius<GT>(p, n, q_next);
    if(!CGAL::is_finite(r_next) || r_next <= FT(0)) {
      //std::cerr << "Invalid radius at iteration " << j << ": " << r_next << std::endl;
      break;
    }
    
    Point_3 c_next = p - (r_next * n);
    FT seperation_angle = CGAL::approximate_angle<GT>(p - c_next, q_next - c_next);
    if(j > 0 && seperation_angle < denoise_preserve) {
      /* std::cout << "Denoise preserve angle achieved at iteration: " << j
                << ", Angle: " << CGAL::to_double(seperation_angle) << " degrees, Center: " << c_next
                << ", Radius: " << r_next << std::endl;*/
      break;
    }
    c = c_next;
    r = r_next;
    q = q_next;
    j++;
    if(j > iteration_limit)
      break;
  }
  
  return {c, r};
}



template <typename TriangleMesh, typename GT, typename VertexPointMap, typename FaceAreaMap, typename VertexAreaMap>
void compute_vertex_areas(const TriangleMesh& tmesh,
                          const VertexPointMap& vpm,
                          FaceAreaMap& face_area_map,
                          VertexAreaMap& vertex_area_map) {
  namespace PMP = CGAL::Polygon_mesh_processing;
  using face_descriptor = typename boost::graph_traits<TriangleMesh>::face_descriptor;
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;

  for(face_descriptor f : faces(tmesh)) {
    double area = PMP::face_area(f, tmesh, parameters::vertex_point_map(vpm));
    put(face_area_map, f, area);
    for(vertex_descriptor v : vertices_around_face(halfedge(f, tmesh), tmesh)) {
      put(vertex_area_map, v, get(vertex_area_map, v) + area / 3.0);
    }
  }
}

template <typename GT,
          typename TriangleMesh,
          typename Tree,
          typename FaceNormalMap,
          typename VertexPointMap,
          typename VertexNormalMap,
          typename VertexMedialSpherePosMap,
          typename VertexMedialSphereRadiusMap>
void compute_shrinking_balls(const TriangleMesh& tmesh,
                             const VertexPointMap& vpm,
                             const Tree& tree,
                             const FaceNormalMap& face_normal_map, // face normal map    
                             const VertexNormalMap& vertex_normal_map,
                             VertexMedialSpherePosMap& vertex_medial_sphere_pos_map,
                             VertexMedialSphereRadiusMap& vertex_medial_sphere_radius_map) {
  using Vector_3 = typename GT::Vector_3;
  using Point_3 = typename GT::Point_3;
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;

  std::cout << "Compute shrinking ball for each vertex" << std::endl;
  for(vertex_descriptor v : vertices(tmesh)) {
    Vector_3 normal = get(vertex_normal_map, v);
    Point_3 p = get(vpm, v);
    auto [center, radius] = shrinking_ball_algorithm<TriangleMesh, GT, Tree>(tmesh, face_normal_map, p, normal, tree);
    put(vertex_medial_sphere_pos_map, v, center);
    put(vertex_medial_sphere_radius_map, v, radius);
  }
}

template <typename TriangleMesh,
          typename GT, 
          typename Tree, 
          typename FaceNormalMap,
          typename VertexPointMap, 
          typename MedialSphereMesh,
          typename VertexNormalMap,
          typename VertexMedialSpherePosMap,
          typename VertexMedialSphereRadiusMap>
void compute_shrinking_balls_and_save_result(const TriangleMesh& tmesh,
                                             const VertexPointMap& vpm,
                                             const Tree& tree,
                                             const FaceNormalMap& face_normal_map,    
                                             const VertexNormalMap& vertex_normal_map,
                                             MedialSphereMesh& sphere_mesh,
                                             VertexMedialSpherePosMap& vertex_medial_sphere_pos_map,
                                             VertexMedialSphereRadiusMap& vertex_medial_sphere_radius_map) {
  using Vector_3 = typename GT::Vector_3;
  using Point_3 = typename GT::Point_3;
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;

  std::cout << "Compute shrinking ball for each vertex" << std::endl;
  for(vertex_descriptor v : vertices(tmesh)) {
    Vector_3 normal = get(vertex_normal_map, v);
    Point_3 p = get(vpm, v);
    auto [center, radius] = shrinking_ball_algorithm<TriangleMesh, GT, Tree>(tmesh,face_normal_map, p, normal, tree);
    put(vertex_medial_sphere_pos_map, v, center);
    put(vertex_medial_sphere_radius_map, v, radius);
    sphere_mesh.add_sphere(typename MedialSphereMesh::Sphere_3(center, radius*radius));
  }
  std::string filename = "medial_sphere_mesh.ply";
  sphere_mesh.write_to_ply_file(filename);
  std::cout << "Medial sphere mesh written to " << filename << "\n";
}

template <typename TriangleMesh,
          typename GT,
          typename VertexPointMap,
          typename VertexAreaMap,
          typename VertexNormalMap,
          typename VertexClusterSphereMap>
void assign_vertices_to_clusters(const TriangleMesh& tmesh,
                                 const VertexPointMap& vpm,
                                 MedialSphereMesh<TriangleMesh, GT>& sphere_mesh,
                                 typename GT::FT lambda,
                                 const VertexAreaMap& vertex_area_map,
                                 const VertexNormalMap& vertex_normal_map,
                                 VertexClusterSphereMap& vertex_cluster_sphere_map) {
  using FT = typename GT::FT;
  using Point_3 = typename GT::Point_3;
  using Vector_3 = typename GT::Vector_3;
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
  using MSMesh = MedialSphereMesh<TriangleMesh, GT>;
  using Sphere_ID = typename MSMesh::Sphere_ID;
  for(vertex_descriptor v : vertices(tmesh)) {
    FT area = get(vertex_area_map, v);
    Point_3 p = get(vpm, v);
    FT min_distance = std::numeric_limits<FT>::max();
    Vector_3 normal = get(vertex_normal_map, v);
    Sphere_ID closest_sphere_id;

    // Find the sphere with smallest distance to the vertex
    for(auto& sphere : sphere_mesh.spheres()) {
      Point_3 center = sphere->get_center();
      FT radius = sphere->get_radius();

      // compute euclidean distance
      FT dist_eucl = CGAL::approximate_sqrt((p - center).squared_length()) - radius;
      dist_eucl *= dist_eucl;
     
      // compute sqem distance
      FT dist_sqem = CGAL::scalar_product(p - center, normal) - radius;
      dist_sqem *= dist_sqem;

      FT distance = area * (dist_sqem + lambda * dist_eucl);
      if(distance < min_distance) {
        min_distance = distance;
        closest_sphere_id = sphere->get_id();
      }
    }
    // Update the closest sphere
    put(vertex_cluster_sphere_map, v, closest_sphere_id);
    sphere_mesh.get_sphere(closest_sphere_id)->accumulate_cluster_area(area);
    sphere_mesh.get_sphere(closest_sphere_id)->add_cluster_vertex(v);
  }
  std::vector<Sphere_ID> sphere_ids_to_remove;
    for(auto& sphere : sphere_mesh.spheres()) { 
      auto& cluster_vertices = sphere->get_cluster_vertices();
      if(cluster_vertices.size() <= 4) {
        for(vertex_descriptor v : cluster_vertices) {
          put(vertex_cluster_sphere_map, v, MSMesh::INVALID_SPHERE_ID);
        }
        sphere_ids_to_remove.push_back(sphere->get_id());
      }

  }
    for(Sphere_ID id : sphere_ids_to_remove) {
    std::cout << "Removing sphere with ID: " << id << " due to small cluster size." << std::endl;
      sphere_mesh.remove(id); // remove spheres with small clusters
    }
 
}

template <typename TriangleMesh, typename GT, typename Tree, typename FaceNormalMap, typename VertexPointMap>
void correct_sphere(
    const TriangleMesh& tmesh,
    const Tree& tree,
    const FaceNormalMap& face_normal_map,
    const VertexPointMap& vpm,
    const GT& traits,
    std::shared_ptr<typename MedialSphereMesh<TriangleMesh, GT>::MedialSphere> sphere,
    const Eigen::Matrix<typename GT::FT, 4, 1>& optimized_sphere_params) {
  using FT = typename GT::FT;
  using Point_3 = typename GT::Point_3;
  using Vector_3 = typename GT::Vector_3;

  Side_of_triangle_mesh<TriangleMesh, GT, VertexPointMap, Tree> side_of(tree, traits);

  Point_3 optimal_center(optimized_sphere_params(0), optimized_sphere_params(1), optimized_sphere_params(2));
  Point_3 cp = tree.closest_point(optimal_center);
  FT len = (optimal_center - cp).squared_length();

  if(len < 1e-12) {
    std::cerr << "Warning: optimal center is too close to the closest point, skipping shrinking ball optimization.\n";
    return;
  }

  Vector_3 normal = (cp - optimal_center) / CGAL::approximate_sqrt(len);

  if((side_of(optimal_center) == CGAL::ON_UNBOUNDED_SIDE))
    normal = -normal; // if the center is outside, flip the normal

  auto [c, r] = shrinking_ball_algorithm<TriangleMesh, GT, Tree>(tmesh, face_normal_map, cp, normal, tree);
  sphere->set_center(c);
  sphere->set_radius(r);

  /* std::cout << "Optimized center: " << optimized_sphere_params(0) << ", " << optimized_sphere_params(1) << ", "
            << optimized_sphere_params(2) << ", radius: " << optimized_sphere_params(3)
            << "\nCorrected center: " << c.x() << ", " << c.y() << ", " << c.z() << ", radius: " << r << std::endl;
*/
            }

template <typename TriangleMesh,
          typename GT,
          typename VertexPointMap,
          typename VertexAreaMap,
          typename FaceNormalMap,
          typename FaceAreaMap,
          typename Tree>
void optimize_sphere_positions(const TriangleMesh& tmesh,
                               const VertexPointMap& vpm,
                               MedialSphereMesh<TriangleMesh, GT>& sphere_mesh,
                               typename GT::FT lambda,
                               const VertexAreaMap& vertex_area_map,
                               const FaceNormalMap& face_normal_map,
                               const FaceAreaMap& face_area_map,
                               const Tree& tree,
                               const GT& traits,
                               bool use_shrinking_ball_correction = false) {
  using FT = typename GT::FT;
  using Point_3 = typename GT::Point_3;
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
  using face_descriptor = typename boost::graph_traits<TriangleMesh>::face_descriptor;
  using Sphere_ID = typename MedialSphereMesh<TriangleMesh, GT>::Sphere_ID;
  using EMat = Eigen::Matrix<FT, Eigen::Dynamic, Eigen::Dynamic>;
  using EVec = Eigen::Matrix<FT, Eigen::Dynamic, 1>;
  using EVec3 = Eigen::Matrix<FT, 3, 1>;
  using EVec4 = Eigen::Matrix<FT, 4, 1>;
  using LDLTSolver = Eigen::LDLT<EMat>;

  for(auto& sphere : sphere_mesh.spheres()) {
    Sphere_ID id = sphere->get_id();
    auto& cluster_vertices = sphere->get_cluster_vertices();
    Point_3 center = sphere->get_center();
    FT radius = sphere->get_radius();

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
        Point_3 p = get(vpm, v);
        EVec3 pos(p.x(), p.y(), p.z());

        // compute sqem energy
        EVec4 lhs = EVec4::Zero();
        FT rhs = 0.0;
        for(face_descriptor f : faces_around_target(halfedge(v, tmesh), tmesh)) {
          typename GT::Vector_3 normal = get(face_normal_map, f);
          EVec3 normal_eigen(normal.x(), normal.y(), normal.z());
          EVec4 n4(normal_eigen(0), normal_eigen(1), normal_eigen(2), 1.0);
          FT area = CGAL::approximate_sqrt(get(face_area_map, f) / 3.0);
          lhs += -n4 * area;
          rhs += -1.0 * ((pos - EVec3(s(0), s(1), s(2))).dot(normal_eigen) - s(3)) * area;
        }
        J.row(idx) = lhs;
        b(idx) = rhs;
        ++idx;

        // compute euclidean energy
        EVec3 d = pos - EVec3(s(0), s(1), s(2));
        FT l = d.norm();
        FT area = CGAL::approximate_sqrt(get(vertex_area_map, v));
        J.row(idx) = EVec4(-(d[0] / l), -(d[1] / l), -(d[2] / l), -1.0) * area * lambda;
        b(idx) = -(l - s(3)) * area * lambda;
        ++idx;
      }

      LDLTSolver solver(J.transpose() * J);
      EVec4 delta_s = solver.solve(J.transpose() * b);
      s += delta_s;
      if(delta_s.norm() < 1e-8) {
        //std::cout << "Convergence achieved after " << i + 1 << " iterations." << std::endl;
        break; // convergence
      }
      sphere->set_do_not_split(false); // reset the split flag
    }
    if(use_shrinking_ball_correction) {
      correct_sphere<TriangleMesh, GT, Tree>(tmesh, tree,face_normal_map, vpm, traits, sphere, s);
    } else {
      // Update sphere with optimized parameters
      sphere->set_center(Point_3(s(0), s(1), s(2)));
      sphere->set_radius(s(3));
    }
  }
  
}


template <typename TriangleMesh,
          typename GT,
          typename VertexPointMap,
          typename VertexAreaMap,
          typename VertexNormalMap,
          typename VertexErrorMap>
typename GT::FT compute_sphere_errors(const TriangleMesh& tmesh,
                                      const VertexPointMap& vpm,
                                      MedialSphereMesh<TriangleMesh, GT>& sphere_mesh,
                                      typename GT::FT lambda,
                                      const VertexAreaMap& vertex_area_map,
                                      const VertexNormalMap& vertex_normal_map,
                                      VertexErrorMap& vertex_error_map) {
  using FT = typename GT::FT;
  using Point_3 = typename GT::Point_3;
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;

  FT total_error = FT(0.0);

  for(auto& sphere : sphere_mesh.spheres()) {
    auto& cluster_vertices = sphere->get_cluster_vertices();
    Point_3 center = sphere->get_center();
    FT radius = sphere->get_radius();
    FT max_dist = std::numeric_limits<FT>::min();
    FT error = 0.0;
    FT area = sphere->get_area();
    if(area <= FT(0.0)) {
      std::cerr << "Warning: sphere with zero area encountered, skipping error computation." << std::endl;
      continue; // skip spheres with zero area
    }

    for(vertex_descriptor v : cluster_vertices) {
      Point_3 p = get(vpm, v);
      FT area = get(vertex_area_map, v);
      FT dist_sqem = CGAL::scalar_product(p - center, get(vertex_normal_map, v)) - radius;
      dist_sqem *= dist_sqem; // square the distance
      

      FT dist_eucl = CGAL::approximate_sqrt((p - center).squared_length()) - radius;
      dist_eucl *= dist_eucl; // square the distance
      
      FT distance =  area * (dist_sqem + lambda * dist_eucl);
      error += distance;

      put(vertex_error_map, v, distance);
      if(distance > max_dist) {
        max_dist = distance;
        sphere->split_vertex = v; // update the split vertex
      }
    }
    sphere->error = error / sphere->get_area(); // average error
    total_error += sphere->error;
  }

  return total_error;
}

template <typename TriangleMesh, typename GT, typename VertexClusterSphereMap>
void update_sphere_neighbors(const TriangleMesh& tmesh,
                             MedialSphereMesh<TriangleMesh, GT>& sphere_mesh,
                             const VertexClusterSphereMap& vertex_cluster_sphere_map) {
  using edge_descriptor = typename boost::graph_traits<TriangleMesh>::edge_descriptor;
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
  using Sphere_ID = typename MedialSphereMesh<TriangleMesh, GT>::Sphere_ID;
  using MedialSphereMesh = typename MedialSphereMesh<TriangleMesh, GT>;

  for(edge_descriptor e : edges(tmesh)) {
    vertex_descriptor v1 = source(e, tmesh);
    vertex_descriptor v2 = target(e, tmesh);
    Sphere_ID s1 = get(vertex_cluster_sphere_map, v1);
    Sphere_ID s2 = get(vertex_cluster_sphere_map, v2);
    if(s1 != s2 && s1 != MedialSphereMesh::INVALID_SPHERE_ID && s2 != MedialSphereMesh::INVALID_SPHERE_ID) {
      auto sphere1 = sphere_mesh.get_sphere(s1);
      auto sphere2 = sphere_mesh.get_sphere(s2);

      if(sphere1->neighbors.find(s2) == sphere1->neighbors.end()) {
        sphere1->add_neighbor(s2);
        sphere2->add_neighbor(s1);
        
      }
    }
  }
}

template <typename TriangleMesh, typename GT, typename VertexMedialSpherePosMap, typename VertexMedialSphereRadiusMap>
void split_spheres(MedialSphereMesh<TriangleMesh, GT>& sphere_mesh,
                   int desired_number_of_spheres,
                   const VertexMedialSpherePosMap& vertex_medial_sphere_pos_map,
                   const VertexMedialSphereRadiusMap& vertex_medial_sphere_radius_map) {
  using MedialSphere = typename MedialSphereMesh<TriangleMesh, GT>::MedialSphere;
  using Sphere_ID = typename MedialSphereMesh<TriangleMesh, GT>::Sphere_ID;
  using Sphere_3 = typename GT::Sphere_3;
  using Point_3 = typename GT::Point_3;
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
  using FT = typename GT::FT;
  if(sphere_mesh.nb_spheres() >= desired_number_of_spheres) {
    return;
  }

  std::cout << "Start Split spheres" << std::endl;
  std::vector<std::shared_ptr<MedialSphere>> sorted_sphere = sphere_mesh.spheres();
  std::sort(sorted_sphere.begin(), sorted_sphere.end(),
            [](const std::shared_ptr<MedialSphere>& a, const std::shared_ptr<MedialSphere>& b) {
              return a->get_radius() > b->get_radius(); // sort by radius
            });

  int to_split_max = std::min(int(std::ceil(sphere_mesh.nb_spheres() * 0.2)), 10);
  for(auto& sphere : sorted_sphere) {
    if(sphere_mesh.nb_spheres() >= desired_number_of_spheres)
      break;
    if(to_split_max > 0 && sphere->can_split()) {
      for(Sphere_ID neighbor_id : sphere->neighbors) {
          sphere_mesh.get_sphere(neighbor_id)->set_do_not_split(true);
      }
      vertex_descriptor split_vertex = sphere->split_vertex;
      Point_3 center = get(vertex_medial_sphere_pos_map, split_vertex);
      FT radius = get(vertex_medial_sphere_radius_map, split_vertex);
      sphere_mesh.add_sphere(Sphere_3(center, radius*radius));
      to_split_max--;
    }
  }
}


/**
 * \ingroup PkgVMASRef
 * computes a static skeleton based on variational medial axis sampling method.
 *
 * @tparam TriangleMesh a model of `HalfedgeListGraph`, `FaceListGraph`, and `MutableFaceGraph`
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param tmesh input triangulated surface mesh
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{concurrency_tag}
 *     \cgalParamDescription{a tag indicating if the task should be done using one or several threads.}
 *     \cgalParamType{Either `CGAL::Sequential_tag`, or `CGAL::Parallel_tag`, or `CGAL::Parallel_if_available_tag`}
 *     \cgalParamDefault{`CGAL::Sequential_tag`}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `PMPSelfIntersectionTraits`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tm`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with
 * `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     must be available in `TriangleMesh`.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 */
template <class TriangleMesh, class NamedParameters = parameters::Default_named_parameters>
void variational_medial_axis_sampling(const TriangleMesh& tmesh,
                                      const NamedParameters& np = parameters::default_values()) {
  namespace PMP = CGAL::Polygon_mesh_processing;
  using parameters::choose_parameter;
  using parameters::get_parameter;

  using GT = typename GetGeomTraits<TriangleMesh, NamedParameters>::type;
  using FT = typename GT::FT;
  using Point_3 = typename GT::Point_3;
  using Vector_3 = typename GT::Vector_3;
  using Sphere_3 = typename GT::Sphere_3;
  using VPM = typename CGAL::GetVertexPointMap<TriangleMesh, NamedParameters>::type;
  using Tree = AABB_tree<AABB_traits_3<GT, AABB_face_graph_triangle_primitive<TriangleMesh, VPM>>>;

  using face_descriptor = typename boost::graph_traits<TriangleMesh>::face_descriptor;
  using edge_descriptor = typename boost::graph_traits<TriangleMesh>::edge_descriptor;
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;

  using Sphere_ID = typename MedialSphereMesh<TriangleMesh, GT>::Sphere_ID;
  using MedialSphere = typename MedialSphereMesh<TriangleMesh, GT>::MedialSphere;
  using MedialSphereMesh = typename MedialSphereMesh<TriangleMesh, GT>;

  // Property map types
  using Vertex_normal_tag = CGAL::dynamic_vertex_property_t<Vector_3>;
  using Vertex_normal_map = typename boost::property_map<TriangleMesh, Vertex_normal_tag>::const_type;
  using Vertex_area_tag = CGAL::dynamic_vertex_property_t<FT>;
  using Vertex_area_map = typename boost::property_map<TriangleMesh, Vertex_area_tag>::const_type;
  using Vertex_error_tag = CGAL::dynamic_vertex_property_t<FT>;
  using Vertex_error_map = typename boost::property_map<TriangleMesh, Vertex_error_tag>::const_type;
  using Vertex_cluster_sphere_tag = CGAL::dynamic_vertex_property_t<Sphere_ID>;
  using Vertex_cluster_sphere_map = typename boost::property_map<TriangleMesh, Vertex_cluster_sphere_tag>::const_type;
  using Vertex_medial_sphere_pos_tag = CGAL::dynamic_vertex_property_t<Point_3>;
  using Vertex_medial_sphere_pos_map =
      typename boost::property_map<TriangleMesh, Vertex_medial_sphere_pos_tag>::const_type;
  using Vertex_medial_sphere_radius_tag = CGAL::dynamic_vertex_property_t<FT>;
  using Vertex_medial_sphere_radius_map =
      typename boost::property_map<TriangleMesh, Vertex_medial_sphere_radius_tag>::const_type;
  using Face_normal_tag = CGAL::dynamic_face_property_t<Vector_3>;
  using Face_normal_map = typename boost::property_map<TriangleMesh, Face_normal_tag>::const_type;
  using Face_area_tag = CGAL::dynamic_face_property_t<FT>;
  using Face_area_map = typename boost::property_map<TriangleMesh, Face_area_tag>::const_type;

  typedef typename internal_np::Lookup_named_param_def<internal_np::concurrency_tag_t, NamedParameters,
                                                       Sequential_tag>::type Concurrency_tag;

  constexpr bool parallel_execution = std::is_same_v<Parallel_tag, Concurrency_tag>;

  GT traits = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));
  VPM vpm =
      choose_parameter(get_parameter(np, internal_np::vertex_point), get_const_property_map(CGAL::vertex_point, tmesh));

  // Build AABB-tree
  Tree tree(faces(tmesh).begin(), faces(tmesh).end(), tmesh, vpm);
  tree.accelerate_distance_queries();
  CGAL::Bbox_3 bbox = tree.bbox();
  std::cout << "AABB-tree bounding box: " << bbox << std::endl;

  // Create property maps
  Vertex_normal_map vertex_normal_map = get(Vertex_normal_tag(), tmesh, Vector_3(0., 0., 0.));
  Vertex_area_map vertex_area_map = get(Vertex_area_tag(), tmesh, 0.);
  Vertex_error_map vertex_error_map = get(Vertex_error_tag(), tmesh, 0.);
  Vertex_cluster_sphere_map vertex_cluster_sphere_map =
      get(Vertex_cluster_sphere_tag(), tmesh, MedialSphereMesh::INVALID_SPHERE_ID);
  Vertex_medial_sphere_pos_map vertex_medial_sphere_pos_map =
      get(Vertex_medial_sphere_pos_tag(), tmesh, Point_3(0., 0., 0.));
  Vertex_medial_sphere_radius_map vertex_medial_sphere_radius_map =
      get(Vertex_medial_sphere_radius_tag(), tmesh, FT(0.));
  Face_normal_map face_normal_map = get(Face_normal_tag(), tmesh, Vector_3(0., 0., 0.));
  Face_area_map face_area_map = get(Face_area_tag(), tmesh, 0.);

  // Compute normals
  PMP::compute_vertex_normals(tmesh, vertex_normal_map, parameters::geom_traits(traits).vertex_point_map(vpm));
  PMP::compute_face_normals(tmesh, face_normal_map, parameters::geom_traits(traits).vertex_point_map(vpm));

  // Compute vertex areas
  compute_vertex_areas<TriangleMesh, GT>(tmesh, vpm, face_area_map, vertex_area_map);

  MedialSphereMesh sphere_mesh;
  int desired_number_of_spheres = 300; // TO DO: pass desired number of sphere by parameter
  FT lambda = FT(0.2);                 // TO DO: pass lambda by parameter
  int iteration_count = 0;             
  int max_iteration =1000;             // TO DO: pass iteration count by parameter  
  FT total_error = FT(0.0);
  FT total_error_diff = (std::numeric_limits<FT>::max)();
  FT last_total_error = total_error;
  Sphere_3 init_sphere(Point_3(0., 0., 0.), FT(1.0)); // initial sphere

  sphere_mesh.add_sphere(init_sphere);

  // Compute shrinking ball of each vertex
  compute_shrinking_balls<GT, TriangleMesh, Tree>(tmesh, vpm, tree, face_normal_map,
                                                   vertex_normal_map, vertex_medial_sphere_pos_map,
                                                  vertex_medial_sphere_radius_map);

  while(iteration_count< max_iteration) {
    // Clean data
    sphere_mesh.reset();

    // Compute the cluster sphere for each vertex
    assign_vertices_to_clusters<TriangleMesh, GT>(tmesh, vpm, sphere_mesh, lambda, vertex_area_map, vertex_normal_map,
                                                  vertex_cluster_sphere_map);

#ifdef CGAL_EIGEN3_ENABLED
    // Update the sphere by optimizing the combined metric
    optimize_sphere_positions<TriangleMesh, GT>(tmesh, vpm, sphere_mesh, lambda, vertex_area_map, face_normal_map,
                                                face_area_map, tree, traits, true);
#endif

    // Compute error of each sphere
    total_error = compute_sphere_errors<TriangleMesh, GT>(tmesh, vpm, sphere_mesh, lambda, vertex_area_map,
                                                          vertex_normal_map, vertex_error_map);

    total_error_diff = std::abs(total_error - last_total_error);
    last_total_error = total_error;
    if((sphere_mesh.nb_spheres() >= desired_number_of_spheres && total_error_diff < 1e-5)) {
      break;
    }
    std::cout << "Total error: " << total_error << ", Total error diff: " << total_error_diff
              << ", Iteration count: " << iteration_count << ", Sphere count: " << sphere_mesh.nb_spheres()
              << std::endl;
    if(total_error_diff < 1e-5 || iteration_count % 10 == 0) {
      
      update_sphere_neighbors<TriangleMesh, GT>(tmesh, sphere_mesh, vertex_cluster_sphere_map);
      
      split_spheres<TriangleMesh, GT>(sphere_mesh, desired_number_of_spheres, vertex_medial_sphere_pos_map,
                                      vertex_medial_sphere_radius_map);
    }
   
    iteration_count++;
  }

  // Write the medial sphere mesh to a file
  update_sphere_neighbors<TriangleMesh, GT>(tmesh, sphere_mesh, vertex_cluster_sphere_map);
  std::string filename = "medial_sphere_mesh.ply";
  sphere_mesh.write_to_ply_file(filename);
  std::cout << "Medial sphere mesh written to " << filename << "\n";
  /* compute_shrinking_balls_and_save_result<TriangleMesh, GT, Tree>(
      tmesh, vpm, tree, face_normal_map,vertex_normal_map, sphere_mesh, vertex_medial_sphere_pos_map, vertex_medial_sphere_radius_map);*/
}

} // namespace CGAL
#endif

  
  
// end of CGAL namespace
// get the closest point to the origin
// Point_3 query(0.,0.,0.);
// Point_3 cp = tree.closest_point(query);

// std::cout << "closest point is " << cp << "\n";

// face_descriptor closest_face = tree.closest_point_and_primitive(query).second;
// halfedge_descriptor h = halfedge(closest_face, tmesh);

// std::array<vertex_descriptor, 3> fvertices;
// fvertices[0] = target(h, tmesh);
// fvertices[1] = target(next(h, tmesh), tmesh);
// fvertices[2] = target(opposite(h, tmesh), tmesh); // or source(h, tmesh)

// std::cout << "closest triangle is " << get(vpm, fvertices[0]) << " "
//                                     << get(vpm, fvertices[1]) << " "
//                                     << get(vpm, fvertices[2]) << "\n";

// Side_of_triangle_mesh<TriangleMesh, GT, VPM, Tree> side_of(tree, traits);
// std::cout << "query inside bounded volume? " << (side_of(query)==CGAL::ON_BOUNDED_SIDE) << "\n";

// Sphere_3 sphere(query, 10.);
// std::cout << "big sphere intersects triangles? " << tree.do_intersect(sphere) << "\n";

// sphere = Sphere_3(cp, 0.01);
// std::cout << "small sphere intersects triangles? " << tree.do_intersect(sphere) << "\n";
 // CGAL_VARIATIONAL_MEDIAL_AXIS_SAMPLING_H
