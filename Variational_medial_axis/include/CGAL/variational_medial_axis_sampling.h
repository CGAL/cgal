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
#include <tuple>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>
#include <functional>
#endif

namespace CGAL {
namespace Internal {
template <typename TriangleMesh, typename GT> class MedialSphere
{
public:
  using FT = typename GT::FT;
  using Point_3 = typename GT::Point_3;
  using Sphere_3 = typename GT::Sphere_3;
  using Sphere_ID = std::size_t;
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;

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
  Sphere_3 get_sphere() const { return sphere; }
  FT get_radius() const { return CGAL::approximate_sqrt(sphere.squared_radius()); }
  Point_3 get_center() const { return sphere.center(); }
  FT get_area() const { return cluster_area; }
  Sphere_ID get_id() const { return id; }
  const std::unordered_set<Sphere_ID>& get_neighbors() const { return neighbors; }
  vertex_descriptor get_split_vertex() const { return split_vertex; }
  FT get_error() const { return error; }
  void set_center(const Point_3& p) { sphere = Sphere_3(p, get_radius()); }
  void set_radius(FT r) { sphere = Sphere_3(get_center(), r * r); }
  void set_cluter_area(FT area) { cluster_area = area; }
  void set_error(FT e) { error = e; }
  void set_split_vertex(vertex_descriptor v) { split_vertex = v; }
  void accumulate_cluster_area(FT area) { cluster_area += area; }
  void add_neighbor(Sphere_ID id) { neighbors.insert(id); }
  bool can_split() const { return !do_not_split && split_vertex != boost::graph_traits<TriangleMesh>::null_vertex(); }
  void set_do_not_split(bool value) { do_not_split = value; }
  std::vector<vertex_descriptor>& get_cluster_vertices() { return cluster_vertices; }
  void add_cluster_vertex(vertex_descriptor v) { cluster_vertices.push_back(v); }

private:
  Sphere_ID id;
  Sphere_3 sphere;
  vertex_descriptor split_vertex; // the vertex for split
  FT error;
  FT cluster_area;
  bool do_not_split = false;      // flag to indicate if the sphere should not be split
  std::unordered_set<Sphere_ID> neighbors;
  std::vector<vertex_descriptor> cluster_vertices;
};

template <typename TriangleMesh, typename GT>
class MedialSphereMesh
{
public:
  using FT = typename GT::FT;
  using Point_3 = typename GT::Point_3;
  using Sphere_3 = typename GT::Sphere_3;
  using MSphere = MedialSphere<TriangleMesh, GT>;
  using Sphere_ID = typename MSphere::Sphere_ID;
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
  static constexpr Sphere_ID INVALID_SPHERE_ID = (std::numeric_limits<Sphere_ID>::max)();

private:
  std::vector<std::shared_ptr<MSphere>> spheres_; // TODO: use Compact container?
  std::unordered_map<Sphere_ID, std::size_t> id_to_index_;
  Sphere_ID next_id_ = 0;

public:
  Sphere_ID add_sphere(const Sphere_3 s) {
    Sphere_ID id = next_id_++;
    spheres_.push_back(std::make_shared<MSphere>(std::move(s), id));
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
  std::shared_ptr<MSphere> get_sphere(Sphere_ID id) { return spheres_[id_to_index_.at(id)]; }
  std::shared_ptr<const MSphere> get_sphere(Sphere_ID id) const { return spheres_[id_to_index_.at(id)]; }
  std::vector<std::shared_ptr<MSphere>>& spheres() { return spheres_; }
  const std::vector<std::shared_ptr<MSphere>>& spheres() const { return spheres_; }

  // TODO: This function will be deleted. It is only used for debugging purpose.
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
      for(Sphere_ID n_id : sphere_ptr->get_neighbors()) {
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
      const auto& neighbors_a = sphere_ptr->get_neighbors();

      for(Sphere_ID b_id : neighbors_a) {
        if(b_id <= a_id)
          continue;
        const auto& neighbors_b = get_sphere(b_id)->get_neighbors();

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
  typename GT::FT norm_v1v2 = CGAL::approximate_sqrt(v1.squared_length() * v2.squared_length());
  typename GT::FT res = norm_v1v2 > 1e-20 ? (v1 * v2) / norm_v1v2 : typename GT::FT(1);
  return std::clamp(res, typename GT::FT(-1), typename GT::FT(1)); // Ensure the result is within [-1, 1]
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

template <typename GT, typename Tree>
std::pair<typename GT::Point_3, typename GT::FT>
shrinking_ball_algorithm(const typename GT::Point_3& p,        // point on the surface
                         const typename GT::Vector_3& n,       // inverse of reaserch direction
                         const Tree& tree) {
  using face_descriptor = typename Tree::Primitive_id;

  using Point_3 = typename GT::Point_3;
  using FT = typename GT::FT;
//  const FT denoise_radius = FT(1e-4);   // model has to be normalized in [0, 1]^3
  const FT denoise_preserve = FT(30.0); // in degree
  const FT delta_convergence = FT(1e-6);

  const int iteration_limit = 30;
  int j = 0;
  FT r = FT(1.0); // initial radius
  Point_3 c = p - (r * n);
  Point_3 q = p - (2 * r * n);
  face_descriptor last_face;
  while(true) {
    auto res = tree.closest_point_and_primitive(c);
    Point_3 q_next = res.first;
    face_descriptor closest_face = res.second;

    FT squared_dist = (q_next - c).squared_length();
    if(squared_dist >= (r - delta_convergence) * (r - delta_convergence)) {
      // std::cout << "Convergence achieved at iteration " << j << ": " << CGAL::to_double(squared_dist) << std::endl;
      break; // convergence
    }
    if((q_next - p).squared_length() <= delta_convergence) {
      // std::cout << "Convergence achieved at iteration " << j << ": " << CGAL::to_double((q_next -
      // p).squared_length())
      //           << std::endl;
      break; // convergence
    }
    if(j>0 && closest_face == last_face) {
      // std::cout << "No change in closest face at iteration " << j << std::endl;
      break; // no change in closest face
    }
    FT r_next = compute_radius<GT>(p, n, q_next);
    if(!CGAL::is_finite(r_next) || r_next <= FT(0)) {
      // std::cerr << "Invalid radius at iteration " << j << ": " << r_next << std::endl;
      break;
    }
    if((r_next - delta_convergence)
      *(r_next - delta_convergence)<1e-4) {
        // std::cout << "Denoise radius achieved at iteration: " << j << ", Radius: " << CGAL::to_double(r_next)
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
    last_face = closest_face;
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
          typename Vertex_Descriptor,
          typename FaceNormalMap,
          typename VertexPointMap,
          typename VertexNormalMap,
          typename VertexMedialSpherePosMap,
          typename VertexMedialSphereRadiusMap>
void compute_one_vertex_shrinking_ball(const TriangleMesh& tmesh,
                             const Tree& tree,
                             const VertexPointMap& vpm,
                             const Vertex_Descriptor v,
                             const FaceNormalMap& face_normal_map, // face normal map
                             const VertexNormalMap& vertex_normal_map,
                             VertexMedialSpherePosMap& vertex_medial_sphere_pos_map,
                             VertexMedialSphereRadiusMap& vertex_medial_sphere_radius_map) {
  using Vector_3 = typename GT::Vector_3;
  using Point_3 = typename GT::Point_3;

  Vector_3 normal = get(vertex_normal_map, v);
  Point_3 p = get(vpm, v);
  auto [center, radius] = shrinking_ball_algorithm<GT>(p, normal, tree);
  put(vertex_medial_sphere_pos_map, v, center);
  put(vertex_medial_sphere_radius_map, v, radius);

}

template <typename GT,
          typename TriangleMesh,
          typename Tree,
          typename FaceNormalMap,
          typename VertexPointMap,
          typename VertexNormalMap,
          typename VertexMedialSpherePosMap,
          typename VertexMedialSphereRadiusMap>
void compute_shrinking_ball(const TriangleMesh& tmesh,
                            const Tree& tree,
                            const VertexPointMap& vpm,
                            const FaceNormalMap& face_normal_map, // face normal map
                            const VertexNormalMap& vertex_normal_map,
                            VertexMedialSpherePosMap& vertex_medial_sphere_pos_map,
                            VertexMedialSphereRadiusMap& vertex_medial_sphere_radius_map) {
  for(auto v : vertices(tmesh)) {
    compute_one_vertex_shrinking_ball(tmesh, tree, vpm, v, face_normal_map, vertex_normal_map,
                                      vertex_medial_sphere_pos_map, vertex_medial_sphere_radius_map);
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
                             const Tree& tree,
                             const VertexPointMap& vpm,
                             const FaceNormalMap& /* face_normal_map */, // face normal map
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
    auto [center, radius] = shrinking_ball_algorithm<GT>(p, normal, tree);
    put(vertex_medial_sphere_pos_map, v, center);
    put(vertex_medial_sphere_radius_map, v, radius);
  }
}
// TODO: This function will be deleted. It is only used for debugging purposes.
template <typename GT,
          typename TriangleMesh,
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
  using Sphere_3 = typename GT::Sphere_3;
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;

  std::cout << "Compute shrinking ball for each vertex" << std::endl;
  for(vertex_descriptor v : vertices(tmesh)) {
    Vector_3 normal = get(vertex_normal_map, v);
    Point_3 p = get(vpm, v);
    auto [center, radius] = shrinking_ball_algorithm<GT>(p, normal, tree);
    put(vertex_medial_sphere_pos_map, v, center);
    put(vertex_medial_sphere_radius_map, v, radius);
    sphere_mesh.add_sphere(Sphere_3(center, radius * radius));
  }
  std::string filename = "medial_sphere_mesh.ply";
  sphere_mesh.write_to_ply_file(filename);
  std::cout << "Medial sphere mesh written to " << filename << "\n";
}

template <typename GT,
          typename TriangleMesh,
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
    Sphere_ID closest_sphere_id = 0.;

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

template <typename GT, typename TriangleMesh, typename Tree, typename FaceNormalMap, typename VertexPointMap>
void correct_sphere(const TriangleMesh& /* tmesh */,
                    const Tree& tree,
                    const FaceNormalMap& /* face_normal_map */,
                    const VertexPointMap& /* vpm */,
                    const GT& traits,
                    std::shared_ptr<MedialSphere<TriangleMesh, GT>> sphere,
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

  auto [c, r] = shrinking_ball_algorithm<GT>(cp, normal, tree);
  sphere->set_center(c);
  sphere->set_radius(r);

  /* std::cout << "Optimized center: " << optimized_sphere_params(0) << ", " << optimized_sphere_params(1) << ", "
            << optimized_sphere_params(2) << ", radius: " << optimized_sphere_params(3)
            << "\nCorrected center: " << c.x() << ", " << c.y() << ", " << c.z() << ", radius: " << r << std::endl;
*/
}


//template <typename GT,
//          typename TriangleMesh,
//          typename VertexPointMap,
//          typename VertexAreaMap,
//          typename FaceNormalMap,
//          typename FaceAreaMap,
//          typename Tree>
//void optimize_sphere_positions(const TriangleMesh& tmesh,
//                               const VertexPointMap& vpm,
//                               MedialSphereMesh<TriangleMesh, GT>& sphere_mesh,
//                               typename GT::FT lambda,
//                               const VertexAreaMap& vertex_area_map,
//                               const FaceNormalMap& face_normal_map,
//                               const FaceAreaMap& face_area_map,
//                               const Tree& tree,
//                               const GT& traits,
//                               bool use_shrinking_ball_correction = false) {
//  using FT = typename GT::FT;
//  using Point_3 = typename GT::Point_3;
//  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
//  using face_descriptor = typename boost::graph_traits<TriangleMesh>::face_descriptor;
//  using Sphere_ID = typename MedialSphereMesh<TriangleMesh, GT>::Sphere_ID;
//  using EMat = Eigen::Matrix<FT, Eigen::Dynamic, Eigen::Dynamic>;
//  using EVec = Eigen::Matrix<FT, Eigen::Dynamic, 1>;
//  using EVec3 = Eigen::Matrix<FT, 3, 1>;
//  using EVec4 = Eigen::Matrix<FT, 4, 1>;
//  using LDLTSolver = Eigen::LDLT<EMat>;
//
//  for(auto& sphere : sphere_mesh.spheres()) {
//    Sphere_ID id = sphere->get_id();
//    auto& cluster_vertices = sphere->get_cluster_vertices();
//    Point_3 center = sphere->get_center();
//    FT radius = sphere->get_radius();
//
//    auto nrows = 2 * cluster_vertices.size();
//    EMat J(nrows, 4);
//    J.setZero();
//    EVec b(nrows);
//    b.setZero();
//    int idx = 0;
//    EVec4 s;
//    s << center.x(), center.y(), center.z(), radius;
//
//    for(int i = 0; i < 10; i++) {
//      idx = 0;
//
//      for(vertex_descriptor v : cluster_vertices) {
//        Point_3 p = get(vpm, v);
//        EVec3 pos(p.x(), p.y(), p.z());
//
//         compute sqem energy
//        EVec4 lhs = EVec4::Zero();
//        FT rhs = 0.0;
//        for(face_descriptor f : faces_around_target(halfedge(v, tmesh), tmesh)) {
//          typename GT::Vector_3 normal = get(face_normal_map, f);
//          EVec3 normal_eigen(normal.x(), normal.y(), normal.z());
//          EVec4 n4(normal_eigen(0), normal_eigen(1), normal_eigen(2), 1.0);
//          FT area = CGAL::approximate_sqrt(get(face_area_map, f) / 3.0);
//          lhs += -n4 * area;
//          rhs += -1.0 * ((pos - EVec3(s(0), s(1), s(2))).dot(normal_eigen) - s(3)) * area;
//        }
//        J.row(idx) = lhs;
//        b(idx) = rhs;
//        ++idx;
//
//         compute euclidean energy
//        EVec3 d = pos - EVec3(s(0), s(1), s(2));
//        FT l = d.norm();
//        FT area = CGAL::approximate_sqrt(get(vertex_area_map, v));
//        J.row(idx) = EVec4(-(d[0] / l), -(d[1] / l), -(d[2] / l), -1.0) * area * lambda;
//        b(idx) = -(l - s(3)) * area * lambda;
//        ++idx;
//      }
//
//      LDLTSolver solver(J.transpose() * J);
//      EVec4 delta_s = solver.solve(J.transpose() * b);
//      s += delta_s;
//      if(delta_s.norm() < 1e-8) {
//         std::cout << "Convergence achieved after " << i + 1 << " iterations." << std::endl;
//        break; // convergence
//      }
//      sphere->set_do_not_split(false); // reset the split flag
//    }
//    if(use_shrinking_ball_correction) {
//      correct_sphere<GT>(tmesh, tree, face_normal_map, vpm, traits, sphere, s);
//    } else {
//       Update sphere with optimized parameters
//      sphere->set_center(Point_3(s(0), s(1), s(2)));
//      sphere->set_radius(s(3));
//    }
//  }
//}
template <typename GT,
          typename TriangleMesh,
          typename MedialSphere,
          typename VertexPointMap,
          typename VertexAreaMap,
          typename FaceNormalMap,
          typename FaceAreaMap,
          typename Tree>
void optimize_single_sphere(const TriangleMesh& tmesh,
                            const VertexPointMap& vpm,
                            std::shared_ptr<MedialSphere>& sphere,
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
  // using Sphere_ID = typename MedialSphere::Sphere_ID;
  using EMat = Eigen::Matrix<FT, Eigen::Dynamic, Eigen::Dynamic>;
  using EVec = Eigen::Matrix<FT, Eigen::Dynamic, 1>;
  using EVec3 = Eigen::Matrix<FT, 3, 1>;
  using EVec4 = Eigen::Matrix<FT, 4, 1>;
  using LDLTSolver = Eigen::LDLT<EMat>;

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
      // std::cout << "Convergence achieved after " << i + 1 << " iterations." << std::endl;
      break; // convergence
    }
    sphere->set_do_not_split(false); // reset the split flag
  }
  if(use_shrinking_ball_correction) {
    correct_sphere<GT>(tmesh, tree, face_normal_map, vpm, traits, sphere, s);
  } else {
    // Update sphere with optimized parameters
    sphere->set_center(Point_3(s(0), s(1), s(2)));
    sphere->set_radius(s(3));
  }
}
template <typename GT,
          typename TriangleMesh,
          typename VertexPointMap,
          typename MedialSphere,
          typename VertexAreaMap,
          typename VertexNormalMap,
          typename VertexErrorMap>
typename GT::FT compute_single_sphere_error(const TriangleMesh& /* tmesh */,
                                            const VertexPointMap& vpm,
                                            std::shared_ptr<MedialSphere>& sphere,
                                            typename GT::FT lambda,
                                            const VertexAreaMap& vertex_area_map,
                                            const VertexNormalMap& vertex_normal_map,
                                            VertexErrorMap& vertex_error_map) {
  using FT = typename GT::FT;
  using Point_3 = typename GT::Point_3;
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;

  auto& cluster_vertices = sphere->get_cluster_vertices();
  Point_3 center = sphere->get_center();
  FT radius = sphere->get_radius();
  FT max_dist = (std::numeric_limits<FT>::min)();
  FT error = 0.0;
  FT area = sphere->get_area();
  if(area <= FT(0.0)) {
    std::cerr << "Warning: sphere with zero area encountered, skipping error computation." << std::endl;
    return FT(0.0);
  }

  for(vertex_descriptor v : cluster_vertices) {
    Point_3 p = get(vpm, v);
    FT area = get(vertex_area_map, v);
    FT dist_sqem = CGAL::scalar_product(p - center, get(vertex_normal_map, v)) - radius;
    dist_sqem *= dist_sqem; // square the distance

    FT dist_eucl = CGAL::approximate_sqrt((p - center).squared_length()) - radius;
    dist_eucl *= dist_eucl; // square the distance

    FT distance = area * (dist_sqem + lambda * dist_eucl);
    error += distance;

    put(vertex_error_map, v, distance);
    if(distance > max_dist) {
      max_dist = distance;
      sphere->set_split_vertex(v); // update the split vertex
    }
  }
  error = error / sphere->get_area(); // average error
  sphere->set_error(error);

  return error;
}

#ifdef CGAL_LINKED_WITH_TBB
template <typename GT,
          typename TriangleMesh,
          typename VertexPointMap,
          typename VertexAreaMap,
          typename VertexNormalMap,
          typename VertexClusterSphereMap>
void assign_vertices_to_clusters_parallel(const TriangleMesh& tmesh,
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

  std::vector<vertex_descriptor> vertices_vector;
  vertices_vector.reserve(num_vertices(tmesh));
  for(vertex_descriptor v : vertices(tmesh)) {
    vertices_vector.push_back(v);
  }
  auto& spheres = sphere_mesh.spheres();
  std::vector<std::pair<Sphere_ID, FT>> vertex_assignments(vertices_vector.size());


  tbb::parallel_for(tbb::blocked_range<std::size_t>(0, vertices_vector.size()),
                    [&](const tbb::blocked_range<std::size_t>& range) {
                      for(std::size_t i = range.begin(); i != range.end(); ++i) {

                        vertex_descriptor v = vertices_vector[i];
                        FT area = get(vertex_area_map, v);
                        Point_3 p = get(vpm, v);
                        FT min_distance = std::numeric_limits<FT>::max();
                        Vector_3 normal = get(vertex_normal_map, v);
                        Sphere_ID closest_sphere_id;

                        // Find the sphere with smallest distance to the vertex
                        for(const auto& sphere : sphere_mesh.spheres()) {
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
                        vertex_assignments[i] = {closest_sphere_id, area};
                      }
                    });
  for(std::size_t i = 0; i < vertices_vector.size(); ++i) {
    vertex_descriptor v = vertices_vector[i];
    auto& [closest_sphere_id, area] = vertex_assignments[i];

    if(closest_sphere_id != MSMesh::INVALID_SPHERE_ID) {
      put(vertex_cluster_sphere_map, v, closest_sphere_id);
      sphere_mesh.get_sphere(closest_sphere_id)->accumulate_cluster_area(area);
      sphere_mesh.get_sphere(closest_sphere_id)->add_cluster_vertex(v);
    }
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


template <typename GT,
          typename TriangleMesh,
          typename Tree,
          typename FaceNormalMap,
          typename VertexPointMap,
          typename VertexNormalMap,
          typename VertexMedialSpherePosMap,
          typename VertexMedialSphereRadiusMap>
void compute_shrinking_balls_parallel(const TriangleMesh& tmesh,
                                      const Tree& tree,
                                      const VertexPointMap& vpm,
                                      const FaceNormalMap& face_normal_map,
                                      const VertexNormalMap& vertex_normal_map,
                                      VertexMedialSpherePosMap& vertex_medial_sphere_pos_map,
                                      VertexMedialSphereRadiusMap& vertex_medial_sphere_radius_map) {
  using Vector_3 = typename GT::Vector_3;
  using Point_3 = typename GT::Point_3;
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
  std::vector<vertex_descriptor> vertices_vector;
  vertices_vector.reserve(num_vertices(tmesh));

  for(vertex_descriptor v : vertices(tmesh)) {
    vertices_vector.push_back(v);
  }
  tbb::parallel_for(tbb::blocked_range<std::size_t>(0, vertices_vector.size()),[&](
      const tbb::blocked_range<std::size_t>& range) {
    for(std::size_t i = range.begin(); i != range.end(); i++) {
      compute_one_vertex_shrinking_ball<GT>(tmesh, tree, vpm, vertices_vector[i], face_normal_map, vertex_normal_map,
                                            vertex_medial_sphere_pos_map, vertex_medial_sphere_radius_map);
    }
  });
}

template <typename GT,
          typename TriangleMesh,
          typename VertexPointMap,
          typename VertexAreaMap,
          typename FaceNormalMap,
          typename FaceAreaMap,
          typename Tree>
void optimize_sphere_positions_parallel(const TriangleMesh& tmesh,
                                        const VertexPointMap& vpm,
                                        MedialSphereMesh<TriangleMesh, GT>& sphere_mesh,
                                        typename GT::FT lambda,
                                        const VertexAreaMap& vertex_area_map,
                                        const FaceNormalMap& face_normal_map,
                                        const FaceAreaMap& face_area_map,
                                        const Tree& tree,
                                        const GT& traits,
                                        bool use_shrinking_ball_correction = false) {
  auto& spheres = sphere_mesh.spheres();

  tbb::parallel_for(tbb::blocked_range<std::size_t>(0,spheres.size()), [&](const tbb::blocked_range<std::size_t>& range) {
        for(std::size_t i = range.begin(); i != range.end(); i++) {
          optimize_single_sphere(tmesh, vpm, spheres[i], lambda, vertex_area_map, face_normal_map, face_area_map, tree,
                                 traits, use_shrinking_ball_correction);
        }
  });
}

template <typename GT,
          typename TriangleMesh,
          typename VertexPointMap,
          typename VertexAreaMap,
          typename VertexNormalMap,
          typename VertexErrorMap>
typename GT::FT compute_sphere_errors_parallel(const TriangleMesh& tmesh,
                                      const VertexPointMap& vpm,
                                      MedialSphereMesh<TriangleMesh, GT>& sphere_mesh,
                                      typename GT::FT lambda,
                                      const VertexAreaMap& vertex_area_map,
                                      const VertexNormalMap& vertex_normal_map,
                                      VertexErrorMap& vertex_error_map) {
  using FT = typename GT::FT;
  using Point_3 = typename GT::Point_3;
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;

  auto& spheres = sphere_mesh.spheres();
  FT total_error = tbb::parallel_reduce(
      tbb::blocked_range<std::size_t>(0, spheres.size()),
      FT(0.0),
      [&](const tbb::blocked_range<std::size_t>& range, FT local_sum) -> FT {
        for(std::size_t i = range.begin(); i != range.end(); ++i) {
          FT sphere_error = compute_single_sphere_error<GT>(tmesh, vpm, spheres[i], lambda, vertex_area_map,
                                                            vertex_normal_map, vertex_error_map);
          local_sum += sphere_error;
        }
        return local_sum;
      },
      std::plus<FT>()
  );
  return total_error;
}
#endif
template <typename GT,
          typename TriangleMesh,
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
  for(auto& sphere : sphere_mesh.spheres()) {
    optimize_single_sphere<GT>(tmesh, vpm, sphere, lambda, vertex_area_map, face_normal_map, face_area_map, tree,
                               traits, use_shrinking_ball_correction);
  }
}

template <typename GT,
          typename TriangleMesh,
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

  FT total_error = FT(0.0);

  for(auto& sphere : sphere_mesh.spheres()) {
    FT error = compute_single_sphere_error<GT>(tmesh, vpm, sphere, lambda, vertex_area_map, vertex_normal_map,
                                               vertex_error_map);
    total_error += error;
  }
  return total_error;
}

template <typename GT, typename TriangleMesh, typename VertexClusterSphereMap>
void update_sphere_neighbors(const TriangleMesh& tmesh,
                             MedialSphereMesh<TriangleMesh, GT>& sphere_mesh,
                             const VertexClusterSphereMap& vertex_cluster_sphere_map) {
  using edge_descriptor = typename boost::graph_traits<TriangleMesh>::edge_descriptor;
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
  using Sphere_ID = typename MedialSphereMesh<TriangleMesh, GT>::Sphere_ID;
  using MedialSphereMesh = MedialSphereMesh<TriangleMesh, GT>;

  for(edge_descriptor e : edges(tmesh)) {
    vertex_descriptor v1 = source(e, tmesh);
    vertex_descriptor v2 = target(e, tmesh);
    Sphere_ID s1 = get(vertex_cluster_sphere_map, v1);
    Sphere_ID s2 = get(vertex_cluster_sphere_map, v2);
    if(s1 != s2 && s1 != MedialSphereMesh::INVALID_SPHERE_ID && s2 != MedialSphereMesh::INVALID_SPHERE_ID) {
      auto sphere1 = sphere_mesh.get_sphere(s1);
      auto sphere2 = sphere_mesh.get_sphere(s2);

      if(sphere1->get_neighbors().find(s2) == sphere1->get_neighbors().end()) {
        sphere1->add_neighbor(s2);
        sphere2->add_neighbor(s1);
      }
    }
  }
}

template <typename GT, typename TriangleMesh, typename VertexMedialSpherePosMap, typename VertexMedialSphereRadiusMap>
void split_spheres(MedialSphereMesh<TriangleMesh, GT>& sphere_mesh,
                   std::size_t desired_number_of_spheres,
                   const VertexMedialSpherePosMap& vertex_medial_sphere_pos_map,
                   const VertexMedialSphereRadiusMap& vertex_medial_sphere_radius_map) {
  using MSphere = MedialSphere<TriangleMesh, GT>;
  using Sphere_ID = typename MedialSphereMesh<TriangleMesh, GT>::Sphere_ID;
  using Sphere_3 = typename GT::Sphere_3;
  using Point_3 = typename GT::Point_3;
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
  using FT = typename GT::FT;
  if(sphere_mesh.nb_spheres() >= desired_number_of_spheres) {
    return;
  }

  std::cout << "Start Split spheres" << std::endl;
  std::vector<std::shared_ptr<MSphere>> sorted_sphere = sphere_mesh.spheres();
  std::sort(sorted_sphere.begin(), sorted_sphere.end(),
            [](const std::shared_ptr<MSphere>& a, const std::shared_ptr<MSphere>& b) {
              return a->get_error() > b->get_error(); // sort by error
            });

  int to_split_max = std::min(int(std::ceil(sphere_mesh.nb_spheres() * 0.2)), 10);
  for(auto& sphere : sorted_sphere) {
    if(sphere_mesh.nb_spheres() >= desired_number_of_spheres)
      break;
    if(to_split_max > 0 && sphere->can_split()) {
      for(Sphere_ID neighbor_id : sphere->get_neighbors()) {
        sphere_mesh.get_sphere(neighbor_id)->set_do_not_split(true);
      }
      vertex_descriptor split_vertex = sphere->get_split_vertex();
      Point_3 center = get(vertex_medial_sphere_pos_map, split_vertex);
      FT radius = get(vertex_medial_sphere_radius_map, split_vertex);
      sphere_mesh.add_sphere(Sphere_3(center, radius * radius));
      to_split_max--;
    }
  }
}
template <typename TriangleMesh, typename GT, typename VertexMedialSpherePosMap, typename VertexMedialSphereRadiusMap>
void insert_sphere(MedialSphereMesh<TriangleMesh, GT>& sphere_mesh,
                   const typename MedialSphereMesh<TriangleMesh, GT>::Sphere_ID id,
                   const VertexMedialSpherePosMap& vertex_medial_sphere_pos_map,
                   const VertexMedialSphereRadiusMap& vertex_medial_sphere_radius_map) {
//  using MSphere = MedialSphere<TriangleMesh, GT>;
//  using Sphere_ID = typename MedialSphereMesh<TriangleMesh, GT>::Sphere_ID;
  using Sphere_3 = typename GT::Sphere_3;
  using Point_3 = typename GT::Point_3;
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
  using FT = typename GT::FT;

  auto& sphere = sphere_mesh.get_sphere(id);
  vertex_descriptor split_vertex = sphere->get_split_vertex();
  Point_3 center = get(vertex_medial_sphere_pos_map, split_vertex);
  FT radius = get(vertex_medial_sphere_radius_map, split_vertex);
  sphere_mesh.add_sphere(Sphere_3(center, radius * radius));
}

template <typename TriangleMesh, typename GT>
void remove_sphere(MedialSphereMesh<TriangleMesh, GT>& sphere_mesh,
                   const typename MedialSphereMesh<TriangleMesh, GT>::Sphere_ID sphere_id) {
    sphere_mesh.remove(sphere_id);
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
template <class TriangleMesh,
          class GT,
          class VPM,
          class Tree,
          class VertexNormalMap,
          class VertexAreaMap,
          class FaceNormalMap,
          class FaceAreaMap,
          class VertexErrorMap,
          class VertexClusterSphereMap,
          class VertexMedialSpherePosMap,
          class VertexMedialSphereRadiusMap,
          class NamedParameters = parameters::Default_named_parameters>
MedialSphereMesh<TriangleMesh, GT>
variational_medial_axis_sampling(const TriangleMesh& tmesh,
                                 const GT& traits,
                                 const VPM& vpm,
                                 const Tree& tree,
                                 VertexNormalMap vertex_normal_map,
                                 VertexAreaMap vertex_area_map,
                                 FaceNormalMap face_normal_map,
                                 FaceAreaMap face_area_map,
                                 VertexErrorMap vertex_error_map,
                                 VertexClusterSphereMap vertex_cluster_sphere_map,
                                 VertexMedialSpherePosMap vertex_medial_sphere_pos_map,
                                 VertexMedialSphereRadiusMap vertex_medial_sphere_radius_map,
                                 const NamedParameters& np = parameters::default_values()) {
  using parameters::choose_parameter;
  using parameters::get_parameter;

  using FT = typename GT::FT;
  using Point_3 = typename GT::Point_3;
  using Sphere_3 = typename GT::Sphere_3;
  using MSMesh = typename Internal::MedialSphereMesh<TriangleMesh, GT>;

  // Extract algorithm parameters from named parameters
  std::size_t desired_number_of_spheres = choose_parameter(get_parameter(np, internal_np::number_of_spheres), 300);
  FT lambda = choose_parameter(get_parameter(np, internal_np::lambda), FT(0.2));
  int max_iteration = choose_parameter(get_parameter(np, internal_np::number_of_iterations), 1000);

  typedef typename internal_np::Lookup_named_param_def<internal_np::concurrency_tag_t, NamedParameters,
                                                       Sequential_tag>::type Concurrency_tag;

  constexpr bool parallel_execution = std::is_same_v<Parallel_tag, Concurrency_tag>;
#ifndef CGAL_LINKED_WITH_TBB
  static_assert(!parallel_execution, "Parallel_tag is enabled but TBB is unavailable.");
#endif
  MSMesh sphere_mesh;

  // Algorithm variables
  int iteration_count = 0;
  FT total_error = FT(0.0);
  FT total_error_diff = (std::numeric_limits<FT>::max)();
  FT last_total_error = total_error;

  // Initialize with one sphere
  Sphere_3 init_sphere(Point_3(0., 0., 0.), FT(1.0));
  sphere_mesh.add_sphere(init_sphere);
#if CGAL_LINKED_WITH_TBB
  if(parallel_execution)
  // Compute shrinking ball of each vertex
     compute_shrinking_balls_parallel<GT>(tmesh, tree, vpm, face_normal_map, vertex_normal_map,
                                       vertex_medial_sphere_pos_map,
                                     vertex_medial_sphere_radius_map);
  else
#endif
  // Compute shrinking ball of each vertex
  compute_shrinking_balls<GT>(tmesh, tree, vpm, face_normal_map, vertex_normal_map,
                                                       vertex_medial_sphere_pos_map, vertex_medial_sphere_radius_map);

  std::cout << "Starting variational medial axis computation..." << std::endl;
  std::cout << "Target number of spheres: " << desired_number_of_spheres << std::endl;
  std::cout << "Lambda: " << lambda << std::endl;
  std::cout << "Max iterations: " << max_iteration << std::endl;

  // Main algorithm loop
  while(iteration_count < max_iteration) {
    // Clean data
    sphere_mesh.reset();
#if CGAL_LINKED_WITH_TBB
    if(parallel_execution) {
    // Compute the cluster sphere for each vertex
    assign_vertices_to_clusters_parallel<GT>(tmesh, vpm, sphere_mesh, lambda, vertex_area_map, vertex_normal_map,
                                    vertex_cluster_sphere_map);

    // Update the sphere by optimizing the combined metric

      optimize_sphere_positions_parallel<GT>(tmesh, vpm, sphere_mesh, lambda, vertex_area_map, face_normal_map,
                                             face_area_map, tree, traits, true);
      total_error = compute_sphere_errors_parallel<GT>(tmesh, vpm, sphere_mesh, lambda, vertex_area_map,
                                                       vertex_normal_map, vertex_error_map);
    } else {
#endif
      // Compute the cluster sphere for each vertex
      assign_vertices_to_clusters<GT>(tmesh, vpm, sphere_mesh, lambda, vertex_area_map, vertex_normal_map,
                                      vertex_cluster_sphere_map);
      // Update the sphere by optimizing the combined metric
      optimize_sphere_positions<GT>(tmesh, vpm, sphere_mesh, lambda, vertex_area_map, face_normal_map, face_area_map,
                                    tree, traits, true);

      // Compute error of each sphere
      total_error = compute_sphere_errors<GT>(tmesh, vpm, sphere_mesh, lambda, vertex_area_map, vertex_normal_map,
                                              vertex_error_map);
#if CGAL_LINKED_WITH_TBB
  }
#endif
    total_error_diff = std::abs(total_error - last_total_error);
    last_total_error = total_error;

    std::cout << "Iteration " << iteration_count << ": spheres=" << sphere_mesh.nb_spheres()
              << ", error=" << total_error << ", error_diff=" << total_error_diff << std::endl;

    // Check convergence
    if((sphere_mesh.nb_spheres() >= desired_number_of_spheres && total_error_diff < 1e-5)) {
      std::cout << "Converged: reached target number of spheres with low error change" << std::endl;
      break;
    }

    // Split spheres periodically or when converged
    if(total_error_diff < 1e-5 || iteration_count % 10 == 0) {
      update_sphere_neighbors<GT>(tmesh, sphere_mesh, vertex_cluster_sphere_map);
      split_spheres<GT>(
          sphere_mesh, desired_number_of_spheres, vertex_medial_sphere_pos_map, vertex_medial_sphere_radius_map);
    }

    iteration_count++;
  }

  // Final neighbor update
  update_sphere_neighbors<GT>(tmesh, sphere_mesh, vertex_cluster_sphere_map);

  std::cout << "Algorithm completed after " << iteration_count << " iterations" << std::endl;
  std::cout << "Final number of spheres: " << sphere_mesh.nb_spheres() << std::endl;
  std::cout << "Final total error: " << total_error << std::endl;


  return sphere_mesh;
}
} // namespace Internal

template <typename TriangleMesh, typename GT>
class Skeleton
{
public:
  using Sphere_3 = typename GT::Sphere_3;
  using Point_3 = typename GT::Point_3;
  using FT = typename GT::FT;
  using Sphere_ID = std::size_t;
  using MSMesh = typename Internal::MedialSphereMesh<TriangleMesh, GT>;

  void build_skeleton_from_medial_sphere_mesh(const MSMesh& sphere_mesh) {
    vertices_.clear();
    edges_.clear();
    faces_.clear();

    std::unordered_map<Sphere_ID, std::size_t> id_to_index;

    // Convert spheres to vertices (as Sphere_3 objects)
    std::size_t vertex_idx = 0;
    for(const auto& sphere : sphere_mesh.spheres()) {
      vertices_.push_back(sphere->get_sphere()); // Store the complete Sphere_3
      id_to_index[sphere->get_id()] = vertex_idx++;
    }

    // Convert sphere adjacencies to edges
    std::set<std::pair<std::size_t, std::size_t>> edge_set;
    for(const auto& sphere : sphere_mesh.spheres()) {
      std::size_t a = id_to_index[sphere->get_id()];
      for(Sphere_ID neighbor_id : sphere->get_neighbors()) {
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
      Sphere_ID a_id = sphere->get_id();
      const auto& neighbors_a = sphere->get_neighbors();

      for(Sphere_ID b_id : neighbors_a) {
        if(b_id <= a_id)
          continue;
        if(id_to_index.find(b_id) == id_to_index.end())
          continue;

        const auto& neighbors_b = sphere_mesh.get_sphere(b_id)->get_neighbors();

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

    // Add spheres to the mesh
    for(std::size_t i = 0; i < vertices_.size(); ++i) {
      sphere_mesh.add_sphere(vertices_[i]);
    }

    // Rebuild adjacencies based on edges
    const auto& spheres = sphere_mesh.spheres();
    for(const auto& edge : edges_) {
      if(edge.first < spheres.size() && edge.second < spheres.size()) {
        auto sphere1 = spheres[edge.first];
        auto sphere2 = spheres[edge.second];
        sphere1->add_neighbor(sphere2->get_id());
        sphere2->add_neighbor(sphere1->get_id());
      }
    }

    return sphere_mesh;
  }

  void write_to_ply_file(const std::string& filename) const {
    std::ofstream ofs(filename);
    if(!ofs) {
      std::cerr << "Error opening file: " << filename << std::endl;
      return;
    }

    // Write header
    ofs << "ply\nformat ascii 1.0\n";
    ofs << "element vertex " << vertices_.size() << "\n";
    ofs << "property float x\nproperty float y\nproperty float z\n";
    ofs << "property float radius\n";
    ofs << "element edge " << edges_.size() << "\n";
    ofs << "property int vertex1\nproperty int vertex2\n";
    ofs << "element face " << faces_.size() << "\n";
    ofs << "property list uchar int vertex_indices\n";
    ofs << "end_header\n";

    // Write vertices (sphere centers and radii)
    for(const auto& sphere : vertices_) {
      const Point_3& center = sphere.center();
      FT radius = CGAL::sqrt(sphere.squared_radius());
      ofs << CGAL::to_double(center.x()) << " " << CGAL::to_double(center.y()) << " " << CGAL::to_double(center.z())
          << " " << CGAL::to_double(radius) << "\n";
    }

    // Write edges
    for(const auto& e : edges_) {
      ofs << e.first << " " << e.second << "\n";
    }

    // Write faces
    for(const auto& f : faces_) {
      ofs << "3 " << f[0] << " " << f[1] << " " << f[2] << "\n";
    }

    ofs.close();
  }

  // Accessor methods
  const std::vector<Sphere_3>& vertices() const { return vertices_; }
  const std::vector<std::pair<std::size_t, std::size_t>>& edges() const { return edges_; }
  const std::vector<std::array<std::size_t, 3>>& faces() const { return faces_; }

  std::size_t number_of_vertices() const { return vertices_.size(); }
  std::size_t number_of_edges() const { return edges_.size(); }
  std::size_t number_of_faces() const { return faces_.size(); }

  // Utility methods
  void clear() {
    vertices_.clear();
    edges_.clear();
    faces_.clear();
  }

private:
  std::vector<Sphere_3> vertices_; // Each vertex is a complete medial sphere
  std::vector<std::pair<std::size_t, std::size_t>> edges_;
  std::vector<std::array<std::size_t, 3>> faces_;
};

template <typename TriangleMesh_,
          typename GeomTraits_,
          typename VertexPointMap_ = Default>
class Variational_medial_axis
{
public:


  using VPM =
      typename Default::Get<VertexPointMap_,
                                    typename boost::property_map<TriangleMesh_, vertex_point_t>::const_type>::type;

  using GT = GeomTraits_;
  using FT = typename GT::FT;
  using Point_3 = typename GT::Point_3;
  using Vector_3 = typename GT::Vector_3;
  using Sphere_3 = typename GT::Sphere_3;
  using MSMesh = typename Internal::MedialSphereMesh<TriangleMesh_, GT>;
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

  Variational_medial_axis(
      const TriangleMesh_& tmesh,
      VPM vpm, const GT& gt = GT())
      : tmesh_(tmesh)
      , traits_(gt)
      , vpm_(vpm) {
    init();
  }

    Variational_medial_axis(const TriangleMesh_& tmesh, const GT& gt = GT())
      : Variational_medial_axis(tmesh,get(vertex_point, tmesh),  gt) {

    }

  void init() {
    namespace PMP = CGAL::Polygon_mesh_processing;


      // Build AABB-tree
    tree_ = std::make_unique<Tree>(faces(tmesh_).begin(), faces(tmesh_).end(), tmesh_, vpm_);
    tree_->accelerate_distance_queries();


    // Create property maps
    vertex_normal_map_ = get(Vertex_normal_tag(), tmesh_, Vector_3(0., 0., 0.));
    vertex_area_map_ = get(Vertex_area_tag(), tmesh_, 0.);
    vertex_error_map_ = get(Vertex_error_tag(), tmesh_, 0.);
    vertex_cluster_sphere_map_ =
        get(Vertex_cluster_sphere_tag(), tmesh_, MSMesh::INVALID_SPHERE_ID);
    vertex_medial_sphere_pos_map_ =
        get(Vertex_medial_sphere_pos_tag(), tmesh_, Point_3(0., 0., 0.));
    vertex_medial_sphere_radius_map_ =
        get(Vertex_medial_sphere_radius_tag(), tmesh_, FT(0.));
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
  }

  template <typename NamedParameters>
  void compute(const NamedParameters& np = parameters::default_values()) {
    auto result_sphere_mesh = Internal::variational_medial_axis_sampling(
        tmesh_,
        traits_,
        vpm_,
        *tree_,
        vertex_normal_map_,
        vertex_area_map_,
        face_normal_map_,
        face_area_map_,
        vertex_error_map_,
        vertex_cluster_sphere_map_,
        vertex_medial_sphere_pos_map_,
        vertex_medial_sphere_radius_map_,
        np);

    *sphere_mesh_ = std::move(result_sphere_mesh);

    /*Internal::compute_shrinking_balls_and_save_result<GT>
        (tmesh_,
            vpm_,
            *tree_,
            face_normal_map_,
            vertex_normal_map_,
            *sphere_mesh_,
            vertex_medial_sphere_pos_map_,
            vertex_medial_sphere_radius_map_);*/
  }

  void add_sphere(Sphere_ID sphere_id) {
    Internal::insert_sphere<TriangleMesh_, GT>(
        *sphere_mesh_,
        sphere_id,
        vertex_medial_sphere_pos_map_,
        vertex_medial_sphere_radius_map_);
  }

  void remove_sphere(Sphere_ID sphere_id) {
      Internal::remove_sphere<TriangleMesh_, GT>(*sphere_mesh_, sphere_id);
  }

  Skeleton<TriangleMesh_, GT> export_skeleton() const {
    Skeleton<TriangleMesh_, GT> skeleton;
    skeleton.build_skeleton_from_medial_sphere_mesh(*sphere_mesh_);
    return skeleton;
  }

private:
  const TriangleMesh_& tmesh_;


  // Computed data members
  GT traits_;
  VPM vpm_;

  std::unique_ptr<Tree> tree_;
  std::unique_ptr<MSMesh> sphere_mesh_;

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
