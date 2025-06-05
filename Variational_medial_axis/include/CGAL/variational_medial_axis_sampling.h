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
//#ifdef CGAL_EIGEN3_ENABLED
#include <Eigen/Dense>
//#endif
namespace CGAL
{

template <typename TriangleMesh, typename GT>
class MedialSphereMesh
{
public:
  using FT = typename GT::FT;
  using Point_3 = typename GT::Point_3;
  using Sphere_3 = typename GT::Sphere_3;
  using Sphere_ID = std::size_t;
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;

  struct MedialSphere
  {
    Sphere_ID id;
    Sphere_3 sphere;
    FT error;
    FT cluster_area;
    std::unordered_set<Sphere_ID> neighbors;
    std::vector<vertex_descriptor> cluster_vertices;

    MedialSphere(const Sphere_3& s, Sphere_ID i)
        : id(i)
        , sphere(s)
        , error(FT(0))
        , cluster_area(FT(0)) {}

    void reset() {
      error = FT(0);
      cluster_area = FT(0);
      neighbors.clear();
    }
    FT get_radius() const { return CGAL::sqrt(sphere.squared_radius()); }
    Point_3 get_center() const { return sphere.center(); }
    void set_center(const Point_3& p) { sphere = Sphere_3(p, get_radius()); }
    void set_radius(FT r) { sphere = Sphere_3(get_center(), r*r); }
    Sphere_ID get_id() const { return id; }
    std::vector<vertex_descriptor>& get_cluster_vertices() { return cluster_vertices; }
  };

private:
  std::vector<std::unique_ptr<MedialSphere>> spheres_;
  std::unordered_map<Sphere_ID, std::size_t> id_to_index_;
  Sphere_ID next_id_ = 0;

public:
  Sphere_ID add_sphere(const Sphere_3& s) {
    Sphere_ID id = next_id_++;
    spheres_.push_back(std::make_unique<MedialSphere>(s, id));
    id_to_index_[id] = spheres_.size() - 1;
    return id;
  }
  void remove(Sphere_ID id) {
    auto it = id_to_index_.find(id);
    if(it == id_to_index_.end())
      return;
    spheres_[it->second].reset();
    id_to_index_.erase(it);
  }
  void reset() {
    for(auto& sphere : spheres_) {
      sphere.reset();
    }
  }
  std::size_t nb_spheres() { return id_to_index_.size(); }
  MedialSphere& get(Sphere_ID id) { return *spheres_[id_to_index_.at(id)]; }
  std::vector<std::unique_ptr<MedialSphere>>& spheres() { return spheres_; }
};

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
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     must be available in `TriangleMesh`.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 */
template <class TriangleMesh, class NamedParameters = parameters::Default_named_parameters >
void variational_medial_axis_sampling(const TriangleMesh& tmesh,
                                 const NamedParameters& np= parameters::default_values())
{
  namespace PMP = CGAL::Polygon_mesh_processing;
  using parameters::choose_parameter;
  using parameters::get_parameter;

  using GT                         = typename GetGeomTraits<TriangleMesh, NamedParameters>::type;
  using FT                         = GT::FT;
  using Point_3                    = GT::Point_3;
  using Vector_3                   = GT::Vector_3;
  using Sphere_3                   = GT::Sphere_3;
  using VPM                        = typename CGAL::GetVertexPointMap<TriangleMesh, NamedParameters>::type;
  using Tree                       = AABB_tree<AABB_traits_3<GT, AABB_face_graph_triangle_primitive<TriangleMesh, VPM>>>;
                                   
  using face_descriptor            = typename boost::graph_traits<TriangleMesh>::face_descriptor;
  using halfedge_descriptor        = typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;
  using vertex_descriptor          = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;

  
#ifdef CGAL_EIGEN3_ENABLED
  using EMat                       = Eigen::Matrix<FT, Eigen::Dynamic, Eigen::Dynamic>;
  using EVec                       = Eigen::Matrix<FT, Eigen::Dynamic, 1>;
  using EVec3                      = Eigen::Matrix<FT, 3, 1>;
  using EVec4                      = Eigen::Matrix<FT, 4, 1>;
  using LDLTSolver                 = Eigen::LDLT<EMat>;
#endif
  GT traits                        = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));
  using Sphere_ID                  = MedialSphereMesh<TriangleMesh, GT>::Sphere_ID;
  using Vertex_normal_tag          = CGAL::dynamic_vertex_property_t<Vector_3>;
  using Vertex_normal_map          = typename boost::property_map<TriangleMesh, Vertex_normal_tag>::const_type;
  using Vertex_area_tag            = CGAL::dynamic_vertex_property_t<FT>;
  using Vertex_area_map            = typename boost::property_map<TriangleMesh, Vertex_area_tag>::const_type;
  using Face_normal_tag            = CGAL::dynamic_face_property_t<Vector_3>;
  using Face_normal_map            = typename boost::property_map<TriangleMesh, Face_normal_tag>::const_type;
  using Face_area_tag              = CGAL::dynamic_face_property_t<FT>;
  using Face_area_map              = typename boost::property_map<TriangleMesh, Face_area_tag>::const_type;
  using Vertex_error_tag           = CGAL::dynamic_vertex_property_t<FT>;
  using Vertex_error_map           = typename boost::property_map<TriangleMesh, Vertex_error_tag>::const_type;
  using Vertex_cluster_sphere_tag  = CGAL::dynamic_vertex_property_t<Sphere_ID>;
  using Vertex_cluster_sphere_map  = typename boost::property_map<TriangleMesh, Vertex_cluster_sphere_tag>::const_type;
  
  typedef typename internal_np::Lookup_named_param_def <
    internal_np::concurrency_tag_t,
    NamedParameters,
    Sequential_tag
  > ::type Concurrency_tag;

  constexpr bool parallel_execution = std::is_same_v<Parallel_tag, Concurrency_tag>;

  
  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, tmesh));

  // build an AABB-tree of faces
  Tree tree(faces(tmesh).begin(), faces(tmesh).end(), tmesh, vpm);
  tree.accelerate_distance_queries();

  // create property maps
  Vertex_normal_map vertex_normal_map = get(Vertex_normal_tag(), tmesh, Vector_3(0., 0., 0.));
  Vertex_area_map vertex_area_map = get(Vertex_area_tag(), tmesh, 0.);
  Vertex_error_map vertex_error_map = get(Vertex_error_tag(), tmesh, 0.);
  Vertex_cluster_sphere_map vertex_cluster_sphere_map = get(Vertex_cluster_sphere_tag(), tmesh, Sphere_ID(0));
  Face_normal_map face_normal_map = get(Face_normal_tag(), tmesh, Vector_3(0., 0., 0.));
  Face_area_map face_area_map = get(Face_area_tag(), tmesh, 0.);
  // compute vertex normal
  PMP::compute_vertex_normals(tmesh, vertex_normal_map, parameters::geom_traits(traits).vertex_point_map(vpm));
  PMP::compute_face_normals(tmesh, face_normal_map, parameters::geom_traits(traits).vertex_point_map(vpm));
  // compute vertex areas
  for(face_descriptor f : faces(tmesh)) {
    double area = PMP::face_area(f, tmesh, parameters::vertex_point_map(vpm));
    put(face_area_map, f, area);
    for(vertex_descriptor v : vertices_around_face(halfedge(f, tmesh), tmesh)) {
      put(vertex_area_map, v, get(vertex_area_map, v) + area / 3.0);
    }
  }

  
  MedialSphereMesh <TriangleMesh, GT> sphere_mesh;
  int desired_number_of_spheres = 100; // TO DO: pass desired number of sphere by parameter
  double lambda = 0.2;                 // TO DO: pass lambda by parameter
  sphere_mesh.add_sphere(Sphere_3(Point_3(0., 0., 0.), FT(1.0)));

  while(sphere_mesh.nb_spheres() < desired_number_of_spheres) {
    //Clean data
    sphere_mesh.reset();
    //------------------------------------------------------------
    //- Compute the cluster sphere for each vertex               -
    //------------------------------------------------------------
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

        //compute euclidean distance
        FT dist_eucl = CGAL::sqrt((p - center).squared_length()) - radius;
        dist_eucl *= dist_eucl;
        dist_eucl *= area; // weight by area

        //compute sqem distance
        FT dist_sqem = CGAL::scalar_product(p - center, normal) - radius;
        dist_sqem *= dist_sqem;

        FT distance = dist_sqem + lambda * dist_eucl;
        if(distance < min_distance) {
          min_distance = distance;
          closest_sphere_id = sphere->get_id();
        }

      }
      // Update the closest sphere
      put(vertex_cluster_sphere_map, v, closest_sphere_id);
      sphere_mesh.get(closest_sphere_id).cluster_area += area;
    }
      
    //------------------------------------------------------------
    //- Update the sphere by optimizing the combined metric      -
    //------------------------------------------------------------
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
          //compute sqem energy
          EVec4 lhs = EVec4::Zero();
          FT rhs = 0.0;
          for(face_descriptor f : faces_around_target(halfedge(v, tmesh), tmesh)) {
            Vector_3 normal = get(face_normal_map, f);
            EVec3 normal_eigen(normal.x(), normal.y(), normal.z());
            EVec4 n4(normal_eigen(0), normal_eigen(1), normal_eigen(2), 1.0);
            FT area = CGAL::sqrt(get(face_area_map, f) / 3.0);
            lhs += -n4 * area;
            rhs += -1.0 * ((pos - EVec3(s(0), s(1), s(2))).dot(normal_eigen) - s(3)) *
                   area; // using Eigen dot??? or CGAL::scalar_product
          }
          J.row(idx) = lhs;
          b(idx) = rhs;
          ++idx;

          // compute euclidean energy
          EVec3 d = pos - EVec3(s(0), s(1), s(2));
          FT l = d.norm();
          FT area = CGAL::sqrt(get(vertex_area_map, v));
          J.row(idx) = EVec4(-(d[0] / l), -(d[1] / l), -(d[2] / l), -1.0) * area * lambda;
          b(idx) = -(l - s(3)) * area * lambda;
          ++idx;
        }
        LDLTSolver solver(J.transpose() * J);
        EVec4 delta_s = solver.solve(J.transpose() * b);
        s += delta_s;
        if(delta_s.norm() < 1e-6) {
          break; // convergence
        }
      }
      //------------------------------------------------------------
      //- Correct Sphere using Shrinking Ball algorithm            -
      //------------------------------------------------------------
    
      auto& [c, r] = shrinking_ball_algorithm(tmesh, Point_3(s(0), s(1), s(2)), tree);
      sphere->set_center(c);
      sphere->set_radius(r);

    }
    //------------------------------------------------------------
    //- Compute error of each sphere                             -
    //------------------------------------------------------------
    
    //split sphere
    // update neighbors
  }

  
}

template <typename GT> 
inline typename GT::FT cosine_angle(const typename GT::Vector_3& v1, const typename GT::Vector_3& v2) {
  return CGAL::scalar_product(v1, v2) / (CGAL::sqrt(v1.squared_length()) * CGAL::sqrt(v2.squared_length()));
}

template <typename GT>
inline typename GT::FT
compute_radius(const typename GT::Vector_3& p, const typename GT::Vector_3& n, const typename GT::Vector_3& q) {
  using Vector_3 = typename GT::Vector_3;
  using FT = typename GT::FT;
  Vector_3 pq = q - p;
  FT d = CGAL::sqrt(pq.squared_length());
  FT cos_angle = cosine_angle(pq, n);
  return d / (2 * cos_angle);
}

template <typename TriangleMesh, typename GT, typename Tree>
std::pair<typename GT::Point_3,typename GT::FT> 
shrinking_ball_algorithm(const TriangleMesh& tmesh,
                                    const typename GT::Point_3& p, // original center
                                    const Tree& tree
                                    ) {
  using Vector_3 = typename GT::Vector_3;
  using FT = typename GT::FT;

  const FT denoise_preserve = FT(20.0) * CGAL::to_double(CGAL::PI) / FT(180.0);
  const FT delta_convergence = FT(1e-5);
  const int iteration_limit = 30;

  Point_3 closest_point = tree.closest_point(p);
  Vector_3 n = (p - closest_point)/CGAL::approximate_sqrt((p-closest_point).squared_length());
  int j = 0;
  FT r = 0.5; // initial radius
  Point_3 c = p - (r * n);
  Point_3 q = p - (2 * r * n);

  while(true) {
    Point_3 q_next = tree.closest_point(c);
    FT dist = CGAL::sqrt((c - q_next).squared_length());
    if(std::abs(d - r) <= delta_convergence || (CGAL::sqrt((q - q_next).squared_length()) < delta_convergence)) {
      // convergence
      break;
    }
    FT r_next = compute_radius(p, n, q_next);
    Point3 c_next = p - (r_next * n);
    FT seperation_angle = CGAL::approximate_angle(q - c_next, q_next - c_next);
    if(j > 0 && seperation_angle < denoise_preserve)
      break;
    c = c_next;
    r = r_next;
    q = q_next;
    j++;
    if(j > iteration_limit)
      break;

  }
  return {c, r};

} 
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
#endif // CGAL_VARIATIONAL_MEDIAL_AXIS_SAMPLING_H
