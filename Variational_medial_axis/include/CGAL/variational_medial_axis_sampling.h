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

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <Eigen/Dense>
#include <tuple>

#ifdef CGAL_LINKED_WITH_TBB
#include <functional>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#endif

namespace CGAL {

#ifndef DOXYGEN_RUNNING
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
  bool do_not_split = false; // flag to indicate if the sphere should not be split
  std::unordered_set<Sphere_ID> neighbors;
  std::vector<vertex_descriptor> cluster_vertices;
};

template <typename TriangleMesh, typename GT> class Medial_Sphere_Mesh
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
 */
template <typename TriangleMesh, typename GT> class Medial_Skeleton
{
  using Sphere_3 = typename GT::Sphere_3;
  using Point_3 = typename GT::Point_3;
  using FT = typename GT::FT;
  using Sphere_ID = std::size_t;
  using MSMesh = Medial_Sphere_Mesh<TriangleMesh, GT>;

public:
  /**
   * Write the medial skeleton to a PLY file.
   *
   * @param filepath The name of the file to write to.
   *
   * This function writes the medial skeleton to a PLY file in ASCII format.
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
   * ...        //N vertices
   * xN yN zN rN
   * v1 v2
   * ...        //M edges
   * vM vN
   * 3 v1 v2 v3
   * ...        //K faces
   * 3 vX vY vZ
   * ```
   */
  void write_to_ply_file(const std::string& filepath) const {
    std::ofstream ofs(filepath);
    if(!ofs) {
      std::cerr << "Error opening file: " << filepath << std::endl;
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
#ifndef DOXYGEN_RUNNING
  void build_skeleton_from_medial_sphere_mesh(const MSMesh& sphere_mesh) {
    clear();

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
#endif // DoXYGEN_RUNNING

  /// \name Accessor methods
  /// @{
  /**
   * Returns the container of vertices, where each vertex is a medial sphere (`Sphere_3`).
   */
  const std::vector<Sphere_3>& vertices() const { return vertices_; }
  /**
   * Returns the container of edges, where each edge is represented as a pair of indices in the vertices vector.
   */
  const std::vector<std::pair<std::size_t, std::size_t>>& edges() const { return edges_; }
  /**
   * Returns the container of faces, where each face is represented as an array of three indices in the vertices vector.
   */
  const std::vector<std::array<std::size_t, 3>>& faces() const { return faces_; }
  /**
   * Returns the number of vertices in the medial skeleton.
   */
  std::size_t number_of_vertices() const { return vertices_.size(); }
  /**
   * Returns the number of edges in the medial skeleton.
   */
  std::size_t number_of_edges() const { return edges_.size(); }
  /**
   * Returns the number of faces in the medial skeleton.
   */
  std::size_t number_of_faces() const { return faces_.size(); }

  /**
   * Clears the data for the medial skeleton.
   */
  void clear() {
    vertices_.clear();
    edges_.clear();
    faces_.clear();
  }
  /**
   * Sets the data for the medial skeleton.
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

/// \ingroup PkgVMASRef
/// \brief Function object for extracting a variational medial skeleton from a triangulated surface mesh.
///
/// This algorithm takes as input a triangulated surface mesh and iteratively samples the medial axis
/// by optimizing the placement of medial spheres. The result is a non-manifold triangle mesh that
/// captures the medial structure of the shape, consisting of both curve segments (1D) and face patches (2D),
/// derived from the adjacency of clusters associated with each medial sphere.
///
/// The process terminates when either the desired number of spheres is reached, or a maximum number
/// of iterations is exceeded.
///
/// \note This method is designed to generate coarse approximations of the medial axis. The number of
/// spheres that can be reliably generated depends on the density of the input surface sampling. For best
/// results, we recommend keeping the total number of medial spheres under nb_vertices/100.
///
/// @tparam TriangleMesh_
///         a model of `FaceListGraph`
///
/// @tparam GomTraits_
///         a model of `VMASTraits`<br>
///         <b>%Default:</b>
/// \code
///     CGAL::Kernel_traits<
///       boost::property_traits<
///          boost::property_map<TriangleMesh, CGAL::vertex_point_t>::type
///        >::value_type
///      >::Kernel
/// \endcode
///
/// @tparam VertexPointMap_
///         a model of `ReadWritePropertyMap`
///         with `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key and
///         `Traits::Point_3` as value type.<br>
///         <b>%Default:</b>
/// \code
///   boost::property_map<TriangleMesh, CGAL::vertex_point_t>::const_type.
/// \endcode
///
///

template <typename TriangleMesh_, typename GeomTraits_ = Default, typename VertexPointMap_ = Default>
class Variational_medial_axis
{
private:
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
  /* The constructor of a vmas object.
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
   *     \cgalParamExtra{This number should generally not exceed 300, as the method is designed to produce coarse
   * skeletons.}
   *   \cgalParamNEnd
   *   \cgalParamNBegin{lambda}
   *     \cgalParamDescription{A weight balancing the two energy terms (SQEM and Euclidean). Smaller values tend to
   * produce skeletons that follow local features more closely.}
   *     \cgalParamType{FT}
   *     \cgalParamDefault{FT(0.2)}
   *     \cgalParamExtra{This parameter must be strictly positive; setting it to zero may prevent correct skeleton
   * connectivity construction.}
   *   \cgalParamNEnd
   *   \cgalParamNBegin{concurrency_tag}
   *     \cgalParamDescription{A tag indicating whether the algorithm should run sequentially or in parallel.}
   *     \cgalParamType{Either `CGAL::Sequential_tag`, `CGAL::Parallel_tag`, or `CGAL::Parallel_if_available_tag`}
   *     \cgalParamDefault{`CGAL::Sequential_tag`}
   *   \cgalParamNEnd
   * \cgalNamedParamsEnd
   **/
  template <class NamedParameters = parameters::Default_named_parameters>
  Variational_medial_axis(const TriangleMesh_& tmesh,
                          VPM vpm,
                          const GT& gt = GT(),
                          const NamedParameters& np = parameters::default_values())
      : tmesh_(tmesh)
      , traits_(gt)
      , vpm_(vpm) {
    using parameters::choose_parameter;
    using parameters::get_parameter;
    desired_number_of_spheres_ = choose_parameter(get_parameter(np, internal_np::number_of_spheres), 100);
    lambda_ = choose_parameter(get_parameter(np, internal_np::lambda), FT(0.2));
    max_iteration_ = choose_parameter(get_parameter(np, internal_np::number_of_iterations), 1000);
    typedef typename internal_np::Lookup_named_param_def<internal_np::concurrency_tag_t, NamedParameters,
                                                         Sequential_tag>::type Concurrency_tag;
    parallel_execution_ = std::is_same_v<Parallel_tag, Concurrency_tag>;

#ifndef CGAL_LINKED_WITH_TBB
    if(parallel_execution_) {
      std::cerr << "Warning: Parallel execution requested but TBB is not available. Using sequential execution."
                << std::endl;
      parallel_execution_ = false;
    }
#endif
  }
  ///@}
  template <class NamedParameters = parameters::Default_named_parameters>
  Variational_medial_axis(const TriangleMesh_& tmesh,
                          const GT& gt = GT(),
                          const NamedParameters& np = parameters::default_values())
      : Variational_medial_axis(tmesh, get(vertex_point, tmesh), gt, np) {}

  /// Initialization that compute some global variable for the algorithm.
  void init() {
    namespace PMP = CGAL::Polygon_mesh_processing;

    // Build AABB-tree
    tree_ = std::make_unique<Tree>(faces(tmesh_).begin(), faces(tmesh_).end(), tmesh_, vpm_);
    tree_->accelerate_distance_queries(vertices(tmesh_).begin(), vertices(tmesh_).end(), vpm_);
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

  /**
   *
   * \brief Computes a static skeleton based on the Variational Medial Axis Sampling method.
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
   *     \cgalParamExtra{This number should generally not exceed 300, as the method is designed to produce coarse
   * skeletons.}
   *   \cgalParamNEnd
   *   \cgalParamNBegin{lambda}
   *     \cgalParamDescription{A weight balancing the two energy terms (SQEM and Euclidean). Smaller values tend to
   * produce skeletons that follow local features more closely.}
   *     \cgalParamType{FT}
   *     \cgalParamDefault{FT(0.2)}
   *     \cgalParamExtra{This parameter must be strictly positive; setting it to zero may prevent correct skeleton
   * connectivity construction.}
   *   \cgalParamNEnd
   *   \cgalParamNBegin{concurrency_tag}
   *     \cgalParamDescription{A tag indicating whether the algorithm should run sequentially or in parallel.}
   *     \cgalParamType{Either `CGAL::Sequential_tag`, `CGAL::Parallel_tag`, or `CGAL::Parallel_if_available_tag`}
   *     \cgalParamDefault{`CGAL::Sequential_tag`}
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

    typedef typename internal_np::Lookup_named_param_def<internal_np::concurrency_tag_t, NamedParameters,
                                                         Sequential_tag>::type Concurrency_tag;

    if constexpr(!std::is_same_v<Concurrency_tag, Sequential_tag>) {
      parallel_execution_ = std::is_same_v<Parallel_tag, Concurrency_tag>;
#ifndef CGAL_LINKED_WITH_TBB
      if(parallel_execution_) {
        std::cerr << "Warning: Parallel execution requested but TBB is not available. Using sequential execution."
                  << std::endl;
        parallel_execution_ = false;
      }
#endif
    }
    MSMesh sphere_mesh;
    bool success = false;
    // Algorithm variables
    iteration_count_ = 0;
    total_error_ = FT(0.0);
    total_error_diff_ = (std::numeric_limits<FT>::max)();
    last_total_error_ = total_error_;

    // Initialize with one sphere
    Sphere_3 init_sphere(Point_3(0., 0., 0.), FT(1.0));
    sphere_mesh_->add_sphere(init_sphere);
#if CGAL_LINKED_WITH_TBB
    if(parallel_execution_)
      // Compute shrinking ball of each vertex
      compute_shrinking_balls_parallel();
    else
#endif
      // Compute shrinking ball of each vertex
      compute_shrinking_balls();
    // TODO: add a paratmeter to control the output verbosity
    std::cout << "Starting variational medial axis computation..." << std::endl;
    std::cout << "Target number of spheres: " << desired_number_of_spheres_ << std::endl;
    std::cout << "Lambda: " << lambda_ << std::endl;
    std::cout << "Max iterations: " << max_iteration_ << std::endl;

    // Main algorithm loop
    while(iteration_count_ < max_iteration_) {
      bool converged = update();
      if(converged) {
        success = true;
        break;
      }
    }

    // Final neighbor update
    update_sphere_neighbors();
    // TODO: add a paratmeter to control the output verbosity
    if(success) {
      std::cout << "Algorithm completed after " << iteration_count_ << " iterations" << std::endl;
      std::cout << "Final number of spheres: " << sphere_mesh_->nb_spheres() << std::endl;
      std::cout << "Final total error: " << total_error_ << std::endl;
    } else {
      std::cout << "Algorithm did not converge after " << max_iteration_ << " iterations" << std::endl;
      std::cout << "Final number of spheres: " << sphere_mesh_->nb_spheres() << std::endl;
      std::cout << "Final total error: " << total_error_ << std::endl;
      std::cout << "Consider decreasing the target number of spheres." << std::endl;
    }
    return success;
  }
  bool update() {
    // Clean data
    sphere_mesh_->reset();
#if CGAL_LINKED_WITH_TBB
    if(parallel_execution_) {
      // Compute the cluster sphere for each vertex
      assign_vertices_to_clusters_parallel();

      // Update the sphere by optimizing the combined metric
      optimize_sphere_positions_parallel(true);
      total_error_ = compute_sphere_errors_parallel();
    } else {
#endif
      // Compute the cluster sphere for each vertex
      assign_vertices_to_clusters();
      // Update the sphere by optimizing the combined metric
      optimize_sphere_positions(true);

      // Compute error of each sphere
      total_error_ = compute_sphere_errors();
#if CGAL_LINKED_WITH_TBB
    }
#endif
    total_error_diff_ = std::abs(total_error_ - last_total_error_);
    last_total_error_ = total_error_;
    // TODO: add a paratmeter to control the output verbosity
    std::cout << "Iteration " << iteration_count_ << ": spheres=" << sphere_mesh_->nb_spheres()
              << ", error=" << total_error_ << ", error_diff=" << total_error_diff_ << std::endl;

    // Check convergence
    if((sphere_mesh_->nb_spheres() >= desired_number_of_spheres_ && total_error_diff_ < converged_threshold_)) {
      // TODO: add a paratmeter to control the output verbosity
      std::cout << "Converged: reached target number of spheres with low error change" << std::endl;
      return true;
    }

    // Split spheres periodically or when converged
    if(total_error_diff_ < converged_threshold_ || iteration_count_ % 10 == 0) {
      update_sphere_neighbors();
      split_spheres();
    }
    iteration_count_++;
    return false;
  }
  /** Add a new sphere by splitting sphere with the id sphere_id.
   * @param sphere_id
   *    The ID of the sphere to split.
   */
  void add_sphere(Sphere_ID sphere_id) {
    auto& sphere = sphere_mesh_->get_sphere(sphere_id);
    vertex_descriptor split_vertex = sphere->get_split_vertex();
    Point_3 center = get(vertex_medial_sphere_pos_map_, split_vertex);
    FT radius = get(vertex_medial_sphere_radius_map_, split_vertex);
    sphere_mesh_->add_sphere(Sphere_3(center, radius * radius));
    // update the spheres
    int max_local_iterations = 10;
    for(int i = 0; i < max_local_iterations; i++) {
      bool converged = update();
      if(converged) {
        std::cout << "Local convergence achieved after adding sphere" << std::endl;
        break;
      }
    }
  }
  /** Remove a sphere by its sphere_id.
   * @param sphere_id
   *    The ID of the sphere to remove.
   */
  void remove_sphere(Sphere_ID sphere_id) {
    sphere_mesh_->remove(sphere_id);
    // update the spheres
    int max_local_iterations = 10;
    for(int i = 0; i < max_local_iterations; i++) {
      bool converged = update();
      if(converged) {
        std::cout << "Local convergence achieved after removing sphere" << std::endl;
        break;
      }
    }
  }

  /** Export the medial skeleton as a `Medial_Skeleton` object.
   *
   * This function builds a `Medial_Skeleton` from the current state of the medial sphere mesh.
   * It extracts the vertices, edges, and faces from the medial sphere mesh and constructs
   * the medial skeleton accordingly.
   *
   * @return
   *     A `Medial_Skeleton` object containing the medial skeleton data.
   */
  Medial_Skeleton<TriangleMesh_, GT> export_skeleton() const {
    Medial_Skeleton<TriangleMesh_, GT> skeleton;
    skeleton.build_skeleton_from_medial_sphere_mesh(*sphere_mesh_);
    return skeleton;
  }
  /**Load a medial skeleton from a PLY file.
   *
   * @param filepath
   *     Filepath to the PLY file containing the medial skeleton data.
   * @return
   *     A `Medial_Skeleton` object containing the loaded skeleton data.
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
  Medial_Skeleton<TriangleMesh_, GT> read_skeleton_from_ply(std::string& filepath) const {
    std::ifstream ifs(filepath);
    if(!ifs) {
      std::cerr << "Error opening file: " << filepath << std::endl;
      return Medial_Skeleton<TriangleMesh_, GT>{};
    }

    Medial_Skeleton<TriangleMesh_, GT> skeleton;
    std::vector<Sphere_3> vertices;
    std::vector<std::pair<std::size_t, std::size_t>> edges;
    std::vector<std::array<std::size_t, 3>> faces;

    std::string line;
    std::size_t num_vertices = 0, num_edges = 0, num_faces = 0;
    bool in_header = true;
    bool is_ascii = false;

    while(std::getline(ifs, line) && in_header) {
      std::istringstream iss(line);
      std::string token;
      iss >> token;

      if(token == "ply") {
        continue;
      } else if(token == "format") {
        std::string format_type;
        iss >> format_type;
        if(format_type == "ascii") {
          is_ascii = true;
        } else {
          std::cerr << "Error: Only ASCII PLY format is supported" << std::endl;
          return Medial_Skeleton<TriangleMesh_, GT>{};
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

    vertices.reserve(num_vertices);
    for(std::size_t i = 0; i < num_vertices; ++i) {
      if(!std::getline(ifs, line)) {
        std::cerr << "Error: Unexpected end of file while reading vertices" << std::endl;
        return Medial_Skeleton<TriangleMesh_, GT>{};
      }

      std::istringstream iss(line);
      double x, y, z, radius;
      if(!(iss >> x >> y >> z >> radius)) {
        std::cerr << "Error: Invalid vertex data at line " << i + 1 << std::endl;
        return Medial_Skeleton<TriangleMesh_, GT>{};
      }

      Point_3 center(x, y, z);
      FT squared_radius = FT(radius * radius);
      vertices.emplace_back(center, squared_radius);
    }

    edges.reserve(num_edges);
    for(std::size_t i = 0; i < num_edges; ++i) {
      if(!std::getline(ifs, line)) {
        std::cerr << "Error: Unexpected end of file while reading edges" << std::endl;
        return Medial_Skeleton<TriangleMesh_, GT>{};
      }

      std::istringstream iss(line);
      std::size_t v1, v2;
      if(!(iss >> v1 >> v2)) {
        std::cerr << "Error: Invalid edge data at line " << i + 1 << std::endl;
        return Medial_Skeleton<TriangleMesh_, GT>{};
      }

      if(v1 >= num_vertices || v2 >= num_vertices) {
        std::cerr << "Error: Edge references invalid vertex indices" << std::endl;
        return Medial_Skeleton<TriangleMesh_, GT>{};
      }

      edges.emplace_back(v1, v2);
    }

    faces.reserve(num_faces);
    for(std::size_t i = 0; i < num_faces; ++i) {
      if(!std::getline(ifs, line)) {
        std::cerr << "Error: Unexpected end of file while reading faces" << std::endl;
        return Medial_Skeleton<TriangleMesh_, GT>{};
      }

      std::istringstream iss(line);
      std::size_t vertex_count;
      if(!(iss >> vertex_count)) {
        std::cerr << "Error: Invalid face data at line " << i + 1 << std::endl;
        return Medial_Skeleton<TriangleMesh_, GT>{};
      }

      if(vertex_count != 3) {
        std::cerr << "Error: Only triangular faces are supported" << std::endl;
        return Medial_Skeleton<TriangleMesh_, GT>{};
      }

      std::size_t v1, v2, v3;
      if(!(iss >> v1 >> v2 >> v3)) {
        std::cerr << "Error: Invalid face vertex indices" << std::endl;
        return Medial_Skeleton<TriangleMesh_, GT>{};
      }

      if(v1 >= num_vertices || v2 >= num_vertices || v3 >= num_vertices) {
        std::cerr << "Error: Face references invalid vertex indices" << std::endl;
        return Medial_Skeleton<TriangleMesh_, GT>{};
      }

      faces.push_back({v1, v2, v3});
    }

    ifs.close();

    skeleton.set_data(std::move(vertices), std::move(edges), std::move(faces));

    std::cout << "Successfully loaded skeleton from " << filepath << std::endl;
    std::cout << "Vertices: " << num_vertices << ", Edges: " << num_edges << ", Faces: " << num_faces << std::endl;

    return skeleton;
  }

  /// \name Parameters
  /// @{
  /** Lambda parameter for the algorithm.
   *
   * This parameter controls the balance between the SQEM and Euclidean energy terms.
   * Smaller values tend to produce skeletons that follow local features more closely.
   */
  FT lambda_param() const { return lambda_; }

  /** set function for `lambda_param()`.
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
   * Get the desired number of spheres for the algorithm.
   */
  int number_of_spheres() const { return desired_number_of_spheres_; }
  /**
   * set function for `number_of_spheres()`.
   */
  void set_number_of_spheres(int num) { desired_number_of_spheres_ = num; }

  /**
   * The maximum number of iterations for the algorithm.
   * In case the algorithm does not converge, it will stop after this number of iterations.
   */
  int max_iteration() const { return max_iteration_; }

  /**
   * set function for `max_iteration()`.
   */
  void set_max_iteration(int max_iter) { max_iteration_ = max_iter; }
  ///@}
private:
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
  shrinking_ball_algorithm(const Point_3& p,               // point on the surface
                           const Vector_3& n,              // inverse of search direction
                                                  FT denoise_radius = 1e-6, // model has to be normalized in [0, 1]^3
                                                  FT delta_convergence = FT(1e-7)) {
    using face_descriptor = typename Tree::Primitive_id;

    denoise_radius *= scale_;
    delta_convergence *= scale_;
    const FT denoise_preserve = FT(30.0); // in degree
    const int iteration_limit = 30;
    int j = 0;
    FT r = FT(1.0) * scale_; // initial radius
    Point_3 c = p - (r * n);
    Point_3 q = p - (2 * r * n);
    face_descriptor last_face;
    while(true) {
      // get the hint directly from the kd-tree only to illustrate the use of the internal kd-tree
      auto hint = tree_->kd_tree().closest_point(c);
      auto [q_next, closest_face] = tree_->closest_point_and_primitive(c, hint);

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
      if(j > 0 && closest_face == last_face) {
        // std::cout << "No change in closest face at iteration " << j << std::endl;
        break; // no change in closest face
      }
      FT r_next = compute_radius(p, n, q_next);
      if(!CGAL::is_finite(r_next) || r_next <= FT(0)) {
        // std::cerr << "Invalid radius at iteration " << j << ": " << r_next << std::endl;
        break;
      }
      if((r_next - delta_convergence) * (r_next - delta_convergence) < denoise_radius) {
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
  std::pair<Point_3, FT> shrinking_ball_algorithm_kdt(const Point_3& p,         // point on the surface
                                                  const Vector_3& n,        // inverse of search direction
                                                  FT delta_convergence = FT(1e-5)) {
    using face_descriptor = typename Tree::Primitive_id;

    delta_convergence *= scale_;
    const FT denoise_preserve = FT(20.0); // in degree
    const int iteration_limit = 30;
    int j = 0;
    FT r = FT(0.25) * scale_; // initial radius
    Point_3 c = p - (r * n);
    Point_3 q = p - (2 * r * n);
    while(true) {
      // get the hint directly from the kd-tree only to illustrate the use of the internal kd-tree
      auto hint = tree_->kd_tree().closest_point(c);
      
      Point_3 q_next = hint.first;

      FT squared_dist = (q_next - c).squared_length();
      if(squared_dist >= (r - delta_convergence) * (r - delta_convergence) || p==q_next) {
        //std::cout << "Convergence achieved at iteration " << j << ": " << CGAL::to_double(squared_dist) << std::endl;
        break; // convergence
      }
      FT r_next = compute_radius(p, n, q_next);
      if(!CGAL::is_finite(r_next) || r_next <= FT(0)) {
        std::cerr << "Invalid radius at iteration " << j << ": " << r_next << std::endl;
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
  void compute_one_vertex_shrinking_ball(vertex_descriptor v) {

    Vector_3 normal = get(vertex_normal_map_, v);
    Point_3 p = get(vpm_, v);
    auto [center, radius] = shrinking_ball_algorithm_kdt(p, normal);
    put(vertex_medial_sphere_pos_map_, v, center);
    put(vertex_medial_sphere_radius_map_, v, radius);
  }

  void compute_shrinking_balls() {
    for(auto v : vertices(tmesh_)) {
      compute_one_vertex_shrinking_ball(v);
    }
  }

  // TODO: This function will be deleted. It is only used for debugging purposes.
  void compute_shrinking_balls_and_save_result() {
    std::cout << "Compute shrinking ball for each vertex" << std::endl;
    for(vertex_descriptor v : vertices(tmesh_)) {
      Vector_3 normal = get(vertex_normal_map_, v);
      Point_3 p = get(vpm_, v);
      auto [center, radius] = shrinking_ball_algorithm_kdt(p, normal);
      put(vertex_medial_sphere_pos_map_, v, center);
      put(vertex_medial_sphere_radius_map_, v, radius);
      sphere_mesh_->add_sphere(Sphere_3(center, radius * radius));
    }
    std::string filename = "medial_sphere_mesh.ply";
    sphere_mesh_->write_to_ply_file(filename);
    std::cout << "Medial sphere mesh written to " << filename << "\n";
  }

  void assign_vertices_to_clusters() {
    for(vertex_descriptor v : vertices(tmesh_)) {

      Point_3 p = get(vpm_, v);
      FT min_distance = (std::numeric_limits<FT>::max)();
      Vector_3 normal = get(vertex_normal_map_, v);
      Sphere_ID closest_sphere_id = 0;

      // Find the sphere with smallest distance to the vertex
      for(auto& sphere : sphere_mesh_->spheres()) {
        Point_3 center = sphere->get_center();
        FT radius = sphere->get_radius();

        // compute euclidean distance
        FT dist_eucl = CGAL::approximate_sqrt((p - center).squared_length()) - radius;
        dist_eucl *= dist_eucl;

        // compute sqem distance
        FT dist_sqem = CGAL::scalar_product(p - center, normal) - radius;
        dist_sqem *= dist_sqem;

        FT distance = dist_sqem + lambda_ * dist_eucl;
        if(distance < min_distance) {
          min_distance = distance;
          closest_sphere_id = sphere->get_id();
        }
      }
      FT area = get(vertex_area_map_, v);
      // Update the closest sphere
      put(vertex_cluster_sphere_map_, v, closest_sphere_id);
      sphere_mesh_->get_sphere(closest_sphere_id)->accumulate_cluster_area(area);
      sphere_mesh_->get_sphere(closest_sphere_id)->add_cluster_vertex(v);
    }
    std::vector<Sphere_ID> sphere_ids_to_remove;
    for(auto& sphere : sphere_mesh_->spheres()) {
      auto& cluster_vertices = sphere->get_cluster_vertices();
      if(cluster_vertices.size() <= 4) {
        for(vertex_descriptor v : cluster_vertices) {
          put(vertex_cluster_sphere_map_, v, MSMesh::INVALID_SPHERE_ID);
        }
        sphere_ids_to_remove.push_back(sphere->get_id());
      }
    }
    for(Sphere_ID id : sphere_ids_to_remove) {
      // TODO: add a paratmeter to control the output verbosity
      std::cout << "Removing sphere with ID: " << id << " due to small cluster size." << std::endl;
      sphere_mesh_->remove(id); // remove spheres with small clusters
    }
  }

  void correct_sphere(std::shared_ptr<MSphere> sphere, const Eigen::Matrix<FT, 4, 1>& optimized_sphere_params) {

    Side_of_triangle_mesh<TriangleMesh_, GT, VPM, Tree> side_of(*tree_, traits_);
    Point_3 optimal_center(optimized_sphere_params(0), optimized_sphere_params(1), optimized_sphere_params(2));
    Point_3 cp = tree_->closest_point(optimal_center);
    FT len = (optimal_center - cp).squared_length();

    if(len < 1e-12) {
      std::cerr << "Warning: optimal center is too close to the closest point, skipping shrinking ball optimization.\n";
      return;
    }

    Vector_3 normal = (cp - optimal_center) / CGAL::approximate_sqrt(len);

    if((side_of(optimal_center) == CGAL::ON_UNBOUNDED_SIDE))
      normal = -normal; // if the center is outside, flip the normal

    auto [c, r] = shrinking_ball_algorithm_kdt(cp, normal);
    sphere->set_center(c);
    sphere->set_radius(r);

    /* std::cout << "Optimized center: " << optimized_sphere_params(0) << ", " << optimized_sphere_params(1) << ", "
              << optimized_sphere_params(2) << ", radius: " << optimized_sphere_params(3)
              << "\nCorrected center: " << c.x() << ", " << c.y() << ", " << c.z() << ", radius: " << r << std::endl;
  */
  }

  void optimize_single_sphere(std::shared_ptr<MSphere>& sphere,

                              bool use_shrinking_ball_correction = false) {
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
      sphere->set_do_not_split(false); // reset the split flag
    }
    if(use_shrinking_ball_correction) {
      correct_sphere(sphere, s);
    } else {
      // Update sphere with optimized parameters
      sphere->set_center(Point_3(s(0), s(1), s(2)));
      sphere->set_radius(s(3));
    }
  }

  FT compute_single_sphere_error(std::shared_ptr<MSphere>& sphere) {
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
        sphere->set_split_vertex(v); // update the split vertex
      }
    }
    error = error / sphere->get_area(); // average error
    sphere->set_error(error);

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
                          FT min_distance = std::numeric_limits<FT>::max();
                          Vector_3 normal = get(vertex_normal_map_, v);
                          Sphere_ID closest_sphere_id;

                          // Find the sphere with smallest distance to the vertex
                          for(const auto& sphere : sphere_mesh_->spheres()) {
                            Point_3 center = sphere->get_center();
                            FT radius = sphere->get_radius();

                            // compute euclidean distance
                            FT dist_eucl = CGAL::approximate_sqrt((p - center).squared_length()) - radius;
                            dist_eucl *= dist_eucl;

                            // compute sqem distance
                            FT dist_sqem = CGAL::scalar_product(p - center, normal) - radius;
                            dist_sqem *= dist_sqem;

                            FT distance = area * (dist_sqem + lambda_ * dist_eucl);
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
        put(vertex_cluster_sphere_map_, v, closest_sphere_id);
        sphere_mesh_->get_sphere(closest_sphere_id)->accumulate_cluster_area(area);
        sphere_mesh_->get_sphere(closest_sphere_id)->add_cluster_vertex(v);
      }
    }
    std::vector<Sphere_ID> sphere_ids_to_remove;
    for(auto& sphere : sphere_mesh_->spheres()) {
      auto& cluster_vertices = sphere->get_cluster_vertices();
      if(cluster_vertices.size() <= 4) {
        for(vertex_descriptor v : cluster_vertices) {
          put(vertex_cluster_sphere_map_, v, MSMesh::INVALID_SPHERE_ID);
        }
        sphere_ids_to_remove.push_back(sphere->get_id());
      }
    }
    for(Sphere_ID id : sphere_ids_to_remove) {
      // TODO: add a paratmeter to control the output verbosity
      std::cout << "Removing sphere with ID: " << id << " due to small cluster size." << std::endl;
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

    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, spheres.size()),
                      [&](const tbb::blocked_range<std::size_t>& range) {
                        for(std::size_t i = range.begin(); i != range.end(); i++) {
                          optimize_single_sphere(spheres[i], use_shrinking_ball_correction);
                        }
                      });
  }

  FT compute_sphere_errors_parallel() {
    auto& spheres = sphere_mesh_->spheres();
    FT total_error = tbb::parallel_reduce(
        tbb::blocked_range<std::size_t>(0, spheres.size()), FT(0.0),
        [&](const tbb::blocked_range<std::size_t>& range, FT local_sum) -> FT {
          for(std::size_t i = range.begin(); i != range.end(); ++i) {
            FT sphere_error = compute_single_sphere_error(spheres[i]);
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
        auto sphere1 = sphere_mesh_->get_sphere(s1);
        auto sphere2 = sphere_mesh_->get_sphere(s2);

        if(sphere1->get_neighbors().find(s2) == sphere1->get_neighbors().end()) {
          sphere1->add_neighbor(s2);
          sphere2->add_neighbor(s1);
        }
      }
    }
  }

  void split_spheres() {
    if(sphere_mesh_->nb_spheres() >= desired_number_of_spheres_) {
      return;
    }
    // TODO: add a paratmeter to control the output verbosity
    std::cout << "Start Split spheres" << std::endl;
    std::vector<std::shared_ptr<MSphere>> sorted_sphere = sphere_mesh_->spheres();
    std::sort(sorted_sphere.begin(), sorted_sphere.end(),
              [](const std::shared_ptr<MSphere>& a, const std::shared_ptr<MSphere>& b) {
                return a->get_error() > b->get_error(); // sort by error
              });

    int to_split_max = std::min(int(std::ceil(sphere_mesh_->nb_spheres() * 0.2)), 10);
    for(auto& sphere : sorted_sphere) {
      if(sphere_mesh_->nb_spheres() >= desired_number_of_spheres_)
        break;
      if(to_split_max > 0 && sphere->can_split()) {
        for(Sphere_ID neighbor_id : sphere->get_neighbors()) {
          sphere_mesh_->get_sphere(neighbor_id)->set_do_not_split(true);
        }
        vertex_descriptor split_vertex = sphere->get_split_vertex();
        Point_3 center = get(vertex_medial_sphere_pos_map_, split_vertex);
        FT radius = get(vertex_medial_sphere_radius_map_, split_vertex);
        sphere_mesh_->add_sphere(Sphere_3(center, radius * radius));
        to_split_max--;
      }
    }
  }

private:
  const TriangleMesh_& tmesh_;

  GT traits_;
  VPM vpm_;
  FT lambda_;
  bool parallel_execution_;
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
