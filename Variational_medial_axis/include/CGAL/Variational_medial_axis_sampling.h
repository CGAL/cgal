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
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Fast_winding_number.h>
#include <CGAL/point_generators_3.h>
#include <Eigen/Dense>
#include <algorithm>
#include <iterator>
#include <tuple>
#include <unordered_set>
#include <boost/iterator/zip_iterator.hpp>


#ifdef CGAL_LINKED_WITH_TBB
#include <functional>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#endif

namespace CGAL {

#ifndef DOXYGEN_RUNNING

template <typename TriangleMesh, typename GT> class Medial_sphere_mesh;

template <typename TriangleMesh, typename GT>
class Medial_sphere
{
public:
  using FT = typename GT::FT;
  using Point_3 = typename GT::Point_3;
  using Point_Index = typename Point_set_3<Point_3>::Index;
  using Sphere_3 = typename GT::Sphere_3;
  using MSMesh = Medial_sphere_mesh<TriangleMesh, GT>;
  using Sphere_ID = typename MSMesh::Sphere_ID;

  // Specalization for Compact_container_with_index
  std::size_t for_compact_container() const { return static_cast<std::size_t>(id_); }
  void for_compact_container(std::size_t idx) { id_ = Sphere_ID(idx); }

  Medial_sphere(const Sphere_3& s)
      : sphere_(s)
      , split_point_idx_(Point_Index())
      , error_(FT(0))
      , cluster_area_(FT(0))
      , do_not_split_(false)
      , id_(MSMesh::INVALID_SPHERE_ID) {}

  Medial_sphere(Sphere_3&& s)
      : sphere_(std::move(s))
      , split_point_idx_(Point_Index())
      , error_(FT(0))
      , cluster_area_(FT(0))
      , do_not_split_(false)
      , id_(MSMesh::INVALID_SPHERE_ID) {}

  void reset() {
    error_ = FT(0);
    split_point_idx_ = Point_Index();
    cluster_area_ = FT(0);
    neighbors_.clear();
    cluster_points_idx_.clear();
  }

  Sphere_3 get_sphere() const { return sphere_; }
  FT get_radius() const { return CGAL::approximate_sqrt(sphere_.squared_radius()); }
  Point_3 get_center() const { return sphere_.center(); }
  FT get_area() const { return cluster_area_; }
  Sphere_ID get_id() const { return id_; }
  const std::unordered_set<Sphere_ID>& get_neighbors() const { return neighbors_; }
  Point_Index get_split_point_idx() const { return split_point_idx_; }
  FT get_error() const { return error_; }

  void set_center(const Point_3& p) { sphere_ = Sphere_3(p, get_radius() * get_radius()); }
  void set_radius(FT r) { sphere_ = Sphere_3(get_center(), r * r); }
  void set_cluster_area(FT area) { cluster_area_ = area; }
  void set_error(FT e) { error_ = e; }
  void set_split_point(Point_Index idx) { split_point_idx_ = idx; }
  void set_id(Sphere_ID id) { id_ = id; }
  void set_do_not_split(bool value) { do_not_split_ = value; }
  std::vector<Point_Index>& get_cluster_points_idx() { return cluster_points_idx_; }
  void add_cluster_point(Point_Index idx) { cluster_points_idx_.push_back(idx); }
  void accumulate_cluster_area(FT area) { cluster_area_ += area; }
  void add_neighbor(Sphere_ID id) { neighbors_.insert(id); }
  bool can_split() const { return !do_not_split_ && split_point_idx_ != Point_Index(); }

private:
  Sphere_3 sphere_;
  Point_Index split_point_idx_; // point chosen for split
  FT error_;
  FT cluster_area_;
  bool do_not_split_ = false; // flag to indicate if the sphere should not be split
  std::unordered_set<Sphere_ID> neighbors_;
  std::vector<Point_Index> cluster_points_idx_;
  Sphere_ID id_;
};

template <typename TriangleMesh, typename GT>
class Medial_sphere_mesh
{
public:
  using FT = typename GT::FT;
  using Point_3 = typename GT::Point_3;
  using Sphere_3 = typename GT::Sphere_3;
  using MSphere = Medial_sphere<TriangleMesh, GT>;

  using Sphere_container =
      CGAL::Compact_container_with_index<MSphere, CGAL_ALLOCATOR(MSphere), Multiply_by_two_policy_for_cc_with_size<64>>;

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
 * @tparam GeomTraits The geometric traits class used for geometric computations.<br>
 *         <b>%Default:</b>
 * \code
 *     CGAL::Kernel_traits<
 *       boost::property_traits<
 *          boost::property_map<TriangleMesh, CGAL::vertex_point_t>::type
 *        >::value_type
 *      >::Kernel
 * \endcode
 */
template <typename TriangleMesh, typename GeomTraits = Default>
class Medial_skeleton
{
  using GT = typename Default::Get<
      GeomTraits,
      typename Kernel_traits<typename boost::property_traits<
          typename boost::property_map<TriangleMesh, vertex_point_t>::type>::value_type>::Kernel>::type;
  using Sphere_3 = typename GT::Sphere_3;
  using Point_3 = typename GT::Point_3;
  using FT = typename GT::FT;
  using MSMesh = Medial_sphere_mesh<TriangleMesh, GT>;
  using Sphere_ID = typename MSMesh::Sphere_ID;

public:

  bool read_PLY(const std::string& filepath) {
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

  MSMesh build_Medial_sphere_mesh_from_skeleton() const {
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

  /// \name Accessor Methods
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
   * @param vertices A list of `Sphere_3` representing the vertices (medial spheres).
   * @param edges A list of pairs representing the edges, where each pair contains indices of vertices.
   * @param faces A list of arrays representing the faces, where each array contains three indices of vertices.
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
/// results, we recommend keeping the total number of medial spheres under nb_samples/100.
/// In addition, the user can explicitly request cluster merging or splitting operations to locally simplify or refine
/// the skeleton structure.
///
/// @tparam TriangleMesh
///         a model of `FaceListGraph`
///
/// @tparam ConcurrencyTag
///         a tag indicating whether the algorithm should run sequentially or in parallel.
///         This determines the execution mode at compile time.<br>
///         <b>%Default:</b> `CGAL::Sequential_tag`<br>
///         <b>%Valid values:</b> `CGAL::Sequential_tag`, `CGAL::Parallel_tag`, `CGAL::Parallel_if_available_tag`
///
///  @tparam AccelerationType
///         a tag indicating whether the algorithm should use Kd-tree or BVH as acceleration structure.<br>
///         <b>%Default:</b> `CGAL::KD_tree_tag`<br>
///         <b>%Valid values:</b> `CGAL::KD_tree_tag`, `CGAL::BVH_tag`,
///
/// @tparam GeomTraits
///         a model of `Kernel`<br>
///         <b>%Default:</b>
/// \code
///     CGAL::Kernel_traits<
///       boost::property_traits<
///          boost::property_map<TriangleMesh, CGAL::vertex_point_t>::type
///        >::value_type
///      >::Kernel
/// \endcode
///
/// @tparam VertexPointMap
///         a model of `ReadWritePropertyMap`
///         with `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key and
///         `GeomTraits::Point_3` as value type.<br>
///         <b>%Default:</b>
/// \code
///   boost::property_map<TriangleMesh, CGAL::vertex_point_t>::const_type.
/// \endcode
///
///

template <typename TriangleMesh,
          typename ConcurrencyTag = Sequential_tag,
          typename AccelerationType = KD_tree_tag,
          typename GeomTraits = Default,
          typename VertexPointMap = Default>
class Variational_medial_axis_sampling
{
  using VPM = typename Default::Get<VertexPointMap,
                                    typename boost::property_map<TriangleMesh, vertex_point_t>::const_type>::type;
  using GT = typename Default::Get<
      GeomTraits,
      typename Kernel_traits<typename boost::property_traits<
          typename boost::property_map<TriangleMesh, vertex_point_t>::type>::value_type>::Kernel>::type;

  using Point_3 = typename GT::Point_3;
  using Vector_3 = typename GT::Vector_3;
  using MSMesh = Medial_sphere_mesh<TriangleMesh, GT>;
  using MSphere = typename MSMesh::MSphere;
  using Sphere_ID = typename MSMesh::Sphere_ID;
  using Point_set = Point_set_3<Point_3>;
  using Point_Index = typename Point_set::Index;
  using KD_tree_search_traits = CGAL::Search_traits_3<GT>;

  using Tree = AABB_tree<AABB_traits_3<GT, AABB_face_graph_triangle_primitive<TriangleMesh, VPM>>>;

  using face_descriptor = typename boost::graph_traits<TriangleMesh>::face_descriptor;
  using edge_descriptor = typename boost::graph_traits<TriangleMesh>::edge_descriptor;
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;

  // Property map types
  using Vertex_normal_tag = CGAL::dynamic_vertex_property_t<Vector_3>;
  using Vertex_normal_map = typename boost::property_map<TriangleMesh, Vertex_normal_tag>::const_type;
  using Face_normal_tag = CGAL::dynamic_face_property_t<Vector_3>;
  using Face_normal_map = typename boost::property_map<TriangleMesh, Face_normal_tag>::const_type;
  using Face_nb_samples_tag = CGAL::dynamic_face_property_t<std::size_t>;
  using Face_nb_samples_map = typename boost::property_map<TriangleMesh, Face_nb_samples_tag>::const_type;
  using Face_area_tag = CGAL::dynamic_face_property_t<typename GT::FT>;
  using Face_area_map = typename boost::property_map<TriangleMesh, Face_area_tag>::const_type;
  using Face_centroid_tag = CGAL::dynamic_face_property_t<Point_3>;
  using Face_centroid_map = typename boost::property_map<TriangleMesh, Face_centroid_tag>::const_type;

  using FWN = CGAL::Fast_winding_number<TriangleMesh, Face_normal_map, Face_area_map, Face_centroid_map, Tree>;

public:

  /// sphere type
  using Sphere_3 = typename GT::Sphere_3;

  /// number type
  using FT = typename GT::FT;

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
  ///  \cgalParamNEnd
  ///   \cgalParamNBegin{number_of_samples}
  ///     \cgalParamDescription{The number of samples on the surface mesh to use for the optimization process.}
  ///     \cgalParamType{unsigned int}
  ///     \cgalParamDefault{max(20000, number_of_spheres * 100)}
  ///     \cgalParamExtra{The number of samples should be significantly larger than the number of spheres.(x100 at least)}
  ///  \cgalParamNEnd
  /// \cgalParamNBegin{number_of_iterations}
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
  ///   \cgalParamNBegin{random_seed}
  ///     \cgalParamDescription{The random seed to sample points on the triangle mesh surface.}
  ///     \cgalParamType{unsigned int}
  ///     \cgalParamExtra{Fix the random seed so that the result can be reproduced}
  ///   \cgalParamNEnd
  ///   \cgalParamNBegin{vertex_point_map}
  ///     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
  ///     \cgalParamType{a class model of `ReadablePropertyMap` with
  ///     `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
  ///                    as key type and `%Point_3` as value type}
  ///     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
  ///   \cgalParamNEnd
  ///   \cgalParamNBegin{verbose}
  ///     \cgalParamDescription{If `true`, the algorithm will print detailed information about its progress.}
  ///     \cgalParamType{bool}
  ///     \cgalParamDefault{false}
  ///   \cgalParamNEnd
  /// \cgalNamedParamsEnd
  ///
  template <class NamedParameters = parameters::Default_named_parameters>
  Variational_medial_axis_sampling(const TriangleMesh& tmesh,
                          const NamedParameters& np = parameters::default_values(),
                          const GT& gt = GT())
      : tmesh_(tmesh)
      , traits_(gt) {

    using parameters::choose_parameter;
    using parameters::get_parameter;

    desired_number_of_spheres_ = choose_parameter(get_parameter(np, internal_np::number_of_spheres), 100);
    nb_samples_ = choose_parameter(get_parameter(np, internal_np::number_of_samples), 100 * desired_number_of_spheres_);
    nb_samples_ = (std::max)(nb_samples_, std::size_t(20000));
    lambda_ = choose_parameter(get_parameter(np, internal_np::lambda), FT(0.2));
    max_iteration_ = choose_parameter(get_parameter(np, internal_np::number_of_iterations), 1000);
    verbose_ = choose_parameter(get_parameter(np, internal_np::verbose), false);
    vpm_ = choose_parameter(get_parameter(np, internal_np::vertex_point),
                            get_const_property_map(CGAL::vertex_point, tmesh));
    seed_ = choose_parameter(get_parameter(np, internal_np::random_seed), 42);
#ifndef CGAL_LINKED_WITH_TBB
    static_assert(!std::is_same_v<ConcurrencyTag, Parallel_tag>, "Parallel_tag is enabled but TBB is unavailable.");
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
   *     \cgalParamNBegin{number_of_samples}
   *      \cgalParamDescription{The number of samples on the surface mesh to use for the optimization process.}
   *      \cgalParamType{unsigned int}
   *      \cgalParamDefault{max(20000, number_of_spheres * 100)}
   *      \cgalParamExtra{The number of samples should be significantly larger than the number of spheres.(x100 at least)}
   *   \cgalParamNEnd
   *  \cgalParamNBegin{number_of_iterations}
   *     \cgalParamDescription{The maximum number of iterations for the optimization process.}
   *     \cgalParamType{int}
   *     \cgalParamDefault{1000}
   *   \cgalParamNEnd
   *   \cgalParamNBegin{lambda}
   *     \cgalParamDescription{A weight balancing the two energy terms {SQEM(Spherical Quadrics Error Metric) and Euclidean}. Smaller values encourage
   * the skeleton to extend deeper into local geometric features of the shape.}
   *     \cgalParamType{FT}
   *     \cgalParamDefault{FT(0.2)}
   *     \cgalParamExtra{The range of this parameter is (0,1].}
   *   \cgalParamNEnd
   *   \cgalParamNBegin{random_seed}
   *      \cgalParamDescription{The random seed to sample points on the triangle mesh surface.}
   *      \cgalParamType{unsigned int}
   *      \cgalParamExtra{Fix the random seed so that the result can be reproduced}
   *    \cgalParamNEnd
   *   \cgalParamNBegin{verbose}
   *     \cgalParamDescription{If `true`, the algorithm will print detailed information about its progress on the
   * standard output.}
   *     \cgalParamType{bool}
   *     \cgalParamDefault{false}
   *   \cgalParamNEnd
   * \cgalNamedParamsEnd
   *
   *
   */
  template <class NamedParameters = parameters::Default_named_parameters>
  bool sample(const NamedParameters& np = parameters::default_values()) {

    using parameters::choose_parameter;
    using parameters::get_parameter;
    // Extract algorithm parameters from named parameters
    desired_number_of_spheres_ =
        choose_parameter(get_parameter(np, internal_np::number_of_spheres), desired_number_of_spheres_);
    lambda_ = choose_parameter(get_parameter(np, internal_np::lambda), lambda_);
    max_iteration_ = choose_parameter(get_parameter(np, internal_np::number_of_iterations), max_iteration_);
    nb_samples_ = choose_parameter(get_parameter(np, internal_np::number_of_samples), nb_samples_);
    verbose_ = choose_parameter(get_parameter(np, internal_np::verbose), verbose_);
    seed_ = choose_parameter(get_parameter(np, internal_np::random_seed), seed_);
    bool success = false;
    reset_algorithm_state();
    sphere_mesh_->spheres().reserve(desired_number_of_spheres_);
    if(nb_samples_ < desired_number_of_spheres_ * 100) {
      // Resample the surface with more points
      nb_samples_ = desired_number_of_spheres_ * 100;
      // Sample the surface mesh
      sample_surface_mesh();
      // Build AABB-tree
      if constexpr(std::is_same_v<AccelerationType, KD_tree_tag>) {
        tree_->accelerate_distance_queries(tpoints_.begin(), tpoints_.end(), tpoints_.point_map());
      }
    }

    if constexpr(std::is_same_v<ConcurrencyTag, Parallel_tag>) {
#if CGAL_LINKED_WITH_TBB
      // Compute the shrinking balls in parallel
      compute_shrinking_balls_parallel();
#endif
    } else {
      compute_shrinking_balls();
    }
    auto it = tpoints_.begin();
    Point_3 init_center = point_medial_sphere_pos_map_[*it];
    FT init_radius = point_medial_sphere_radius_map_[*it];
    Sphere_3 init_sphere(init_center, init_radius * init_radius);
    sphere_mesh_->add_sphere(init_sphere);
    if(verbose_) {
      std::cout << "Starting variational medial axis computation..." << std::endl;
      std::cout << "Target number of spheres: " << desired_number_of_spheres_ << std::endl;
      std::cout << "Lambda: " << lambda_ << std::endl;
      std::cout << "Max iterations: " << max_iteration_ << std::endl;
      std::cout << "Sampled points: " << tpoints_.size() << std::endl;
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
   * \return `true` if the algorithm has converged, `false` otherwise.
   */

  bool update_single_step(bool enable_split = false) {
    // Clean data
    sphere_mesh_->reset();

    if constexpr(std::is_same_v<ConcurrencyTag, Parallel_tag>) {
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
   * \return `true` if spheres were added successfully, `false` otherwise.
   */
  bool add_spheres(int nb_sphere) {
    if(nb_sphere == 0) {
      return false;
    }
    iteration_count_ = 0;
    int max_iteration = (std::min)(1000, (std::max)(100, 10 * nb_sphere));
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
    Point_Index split_vertex = sphere.get_split_point_idx();
    Point_3 center = point_medial_sphere_pos_map_[split_vertex];
    FT radius = point_medial_sphere_radius_map_[split_vertex];
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
   * \brief returns the medial skeleton as a `Medial_skeleton` object.
   *
   * This function builds a medial skeleton from the current state of the medial sphere mesh.
   * It extracts the vertices, edges, and faces from the medial sphere mesh and constructs
   * the medial skeleton accordingly.
   *
   * \return
   *     A `Medial_skeleton` object containing the medial skeleton data.
   */
  Medial_skeleton<TriangleMesh, GeomTraits> skeleton() const {
    Medial_skeleton<TriangleMesh, GeomTraits> skeleton;
    skeleton.build_skeleton_from_medial_sphere_mesh(*sphere_mesh_);
    return skeleton;
  }

  /// \name Parameters
  /// @{
  /**
   * \brief returns the lambda parameter for the algorithm.
   *
   * This parameter controls the balance between the SQEM and Euclidean energy terms.
   * Smaller values encourage the skeleton to extend deeper into local geometric features of the shape.
   */
  FT lambda() const { return lambda_; }

  /**
   * \brief sets function for `lambda()`.
   * Note: The value of lambda must be strictly positive; if set to zero, it will default to 0.2.
   */
  void set_lambda(FT lambda) {
    if(lambda <= FT(0)) {
      if(verbose_){
        std::cerr << "Warning: Lambda must be strictly positive. Setting to default value of 0.2." << std::endl;
      }
      lambda_ = FT(0.2);
    } else {
      lambda_ = lambda;
    }
  }
  /**
   * \brief returns the desired number of spheres for the algorithm.
   */
  int number_of_spheres() const { return desired_number_of_spheres_; }
  /**
   * \brief sets the desired number of spheres.
   */
  void set_number_of_spheres(int num) { desired_number_of_spheres_ = num; }

  /**
   * \brief returns the maximum number of iterations for the algorithm.
   *
   * In case the algorithm does not converge, it will stop after this number of iterations.
   */
  int number_of_iterations() const { return max_iteration_; }

  /**
   * \brief sets the maximum number of iterations.
   */
  void set_number_of_iterations(int max_iter) { max_iteration_ = max_iter; }
  ///@}
private:
  //Sample points on the surface mesh and compute k-nearest neighbors for each point so that we can construct the connectivity of skeleton
  void sample_surface_mesh() {

    // Reset point set
    tpoints_.clear();

    // Reset face sample counts
    for(face_descriptor f : faces(tmesh_)) {
      put(face_nb_samples_map_, f, 0);
    }
    bool success = false;
    // Create point set property maps
    std::tie(point_from_face_map_, success) = tpoints_.template add_property_map<face_descriptor>(
        "points_from_face", boost::graph_traits<TriangleMesh>::null_face());
    CGAL_assertion(success);
    std::tie(point_normal_map_, success) =
        tpoints_.template add_property_map<Vector_3>("point_normal", Vector_3(0., 0., 0.));
    CGAL_assertion(success);
    std::tie(point_area_map_, success) = tpoints_.template add_property_map<FT>("point_area", FT(0.));
    CGAL_assertion(success);
    std::tie(point_medial_sphere_pos_map_, success) =
        tpoints_.template add_property_map<Point_3>("point_medial_sphere_pos", Point_3(0., 0., 0.));
    CGAL_assertion(success);
    std::tie(point_medial_sphere_radius_map_, success) =
        tpoints_.template add_property_map<FT>("point_medial_sphere_radius", FT(0.));
    CGAL_assertion(success);
    std::tie(point_error_map_, success) = tpoints_.template add_property_map<FT>("point_error", FT(0.));
    CGAL_assertion(success);
    std::tie(point_knn_map_, success) =
        tpoints_.template add_property_map<std::vector<Point_Index>>("point_knn", std::vector<Point_Index>());
    CGAL_assertion(success);
    std::tie(point_cluster_sphere_map_, success) =
        tpoints_.template add_property_map<Sphere_ID>("point_cluster_sphere", MSMesh::INVALID_SPHERE_ID);
    CGAL_assertion(success);
    // Sample points on the surface mesh

    CGAL::Random rng(seed_);

    CGAL::Random_points_in_triangle_mesh_3<TriangleMesh, VPM> g(tmesh_, vpm_, rng);
    for(std::size_t i = 0; i < nb_samples_; ++i) {
      Point_3 p = *g;
      face_descriptor f = g.last_item_picked();
      put(face_nb_samples_map_, f, get(face_nb_samples_map_, f) + 1);
      auto it = tpoints_.insert(p);
      point_from_face_map_[*it] = f;
      point_normal_map_[*it] = get(face_normal_map_, f);
      ++g;
    }

    for(auto it = tpoints_.begin(); it != tpoints_.end(); it++) {
      face_descriptor f = point_from_face_map_[*it];
      FT area = get(face_area_map_, f);
      std::size_t nb_sample = get(face_nb_samples_map_, f);
      point_area_map_[*it] = area / FT(nb_sample);
    }


    // Compute k-nearest neighbors for each point
    // The kd-tree is built temporarily here only for k-nearest neighbor search and discarded after this function
    using Point_with_index = boost::tuple<Point_3, Point_Index>;
    using BaseTraits = CGAL::Search_traits_3<GT>;
    using Traits =
        CGAL::Search_traits_adapter<Point_with_index, CGAL::Nth_of_tuple_property_map<0, Point_with_index>, BaseTraits>;
    using KNN = CGAL::Orthogonal_k_neighbor_search<Traits>;
    using Kd_tree = typename KNN::Tree;
    std::vector<Point_3> points;
    std::vector<Point_Index> indices;
    points.reserve(tpoints_.size());
    indices.reserve(tpoints_.size());
    for(auto it = tpoints_.begin(); it != tpoints_.end(); ++it) {
      points.push_back(tpoints_.point(*it));
      indices.push_back(*it);
    }

     Kd_tree tree(boost::make_zip_iterator(boost::make_tuple(points.begin(), indices.begin())),
                 boost::make_zip_iterator(boost::make_tuple(points.end(), indices.end())));

    for(auto it = tpoints_.begin(); it != tpoints_.end(); ++it) {
       Point_Index idx = *it;
       const Point_3& query_point = tpoints_.point(idx);
       KNN knn(tree, query_point, k_);
       auto& neighbors = point_knn_map_[idx];
       neighbors.clear();
       neighbors.reserve(k_);

       for(auto knn_it = knn.begin(); knn_it != knn.end(); ++knn_it) {
         Point_Index neighbor_idx = boost::get<1>(knn_it->first);
         neighbors.push_back(neighbor_idx);
       }
     }
  }

  /// Initialization that compute some global variable for the algorithm.
  void init() {
    namespace PMP = CGAL::Polygon_mesh_processing;

    // Create Surface_mesh property maps
    vertex_normal_map_ = get(Vertex_normal_tag(), tmesh_, Vector_3(0., 0., 0.));
    face_nb_samples_map_ = get(Face_nb_samples_tag(), tmesh_, 0);
    face_normal_map_ = get(Face_normal_tag(), tmesh_, Vector_3(0., 0., 0.));
    face_area_map_ = get(Face_area_tag(), tmesh_, 0.);
    face_centroid_map_ = get(Face_centroid_tag(), tmesh_, Point_3(0., 0., 0.));
    // Compute normals
    PMP::compute_vertex_normals(tmesh_, vertex_normal_map_);
    PMP::compute_face_normals(tmesh_, face_normal_map_);

    // Compute vertex areas
    for(face_descriptor f : faces(tmesh_)) {
      double area = PMP::face_area(f, tmesh_);
      put(face_area_map_, f, area);
      FT sum = FT(0), x = FT(0), y = FT(0), z = FT(0);
      for(vertex_descriptor v : vertices_around_face(halfedge(f, tmesh_), tmesh_)) {
        const Point_3& point = get(vpm_, v);
        x += point.x();
        y += point.y();
        z += point.z();
        sum += FT(1);
      }
      CGAL_precondition(sum > FT(0));
      x /= sum;
      y /= sum;
      z /= sum;
      put(face_centroid_map_, f, Point_3(x, y, z));
      put(face_area_map_, f, area);
    }
    // Sample the surface mesh
    tree_ = std::make_unique<Tree>(faces(tmesh_).begin(), faces(tmesh_).end(), tmesh_, vpm_);
    // get bounding box of the mesh
    sample_surface_mesh();
    // Build AABB-tree
    if constexpr(std::is_same_v<AccelerationType, KD_tree_tag>) {
      tree_->accelerate_distance_queries(tpoints_.begin(), tpoints_.end(), tpoints_.point_map());
    } else {
      tree_->accelerate_distance_queries();
    }

    auto bbox = tree_->bbox();
    scale_ = (std::max)(bbox.xmax() - bbox.xmin(), (std::max)(bbox.ymax() - bbox.ymin(), bbox.zmax() - bbox.zmin()));
    converged_threshold_ = scale_ * scale_ * 1e-5;

    // Initialize the fast winding number function
    fast_winding_number_ = std::make_unique<FWN>(tmesh_, face_normal_map_, face_area_map_, face_centroid_map_, *tree_);
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
  shrinking_ball_algorithm_bvh(const std::vector<face_descriptor>& incident_faces,
                               const Point_3& p,                  // point on the surface
                               const Vector_3& n,                 // inverse of search direction
                               FT delta_convergence = FT(1e-5)) {
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
      CGAL_assertion(CGAL::is_finite(r_next) && r_next > FT(0));
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
    FT r = FT(1.0) * scale_; // initial radius
    Point_3 c = p - (r * n);

    while(true) {

      auto hint = tree_->kd_tree().closest_point(c);
      Point_3 q_next = hint.first;

      FT squared_dist = (q_next - c).squared_length();
      if(squared_dist >= (r - delta_convergence) * (r - delta_convergence) || p == q_next) {
        // std::cout << "Convergence achieved at iteration " << j << ": " << CGAL::to_double(squared_dist) <<
        // std::endl;

        break; // convergence
      }
      FT r_next = compute_radius(p, n, q_next);
      CGAL_assertion(CGAL::is_finite(r_next) && r_next > FT(0));
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

  void compute_one_vertex_shrinking_ball(Point_Index idx) {

    Vector_3 normal = point_normal_map_[idx];
    Point_3 p = tpoints_.point(idx);
    auto [center, radius] = compute_shrinking_ball_impl(p, normal, AccelerationType{});
    point_medial_sphere_pos_map_[idx] = center;
    point_medial_sphere_radius_map_[idx] = radius;
  }

  void compute_shrinking_balls() {
    for(auto it = tpoints_.begin(); it!= tpoints_.end(); ++it) {
      Point_Index idx = *it;
      compute_one_vertex_shrinking_ball(idx);
    }
  }

  void assign_vertices_to_clusters() {
    for(auto it = tpoints_.begin(); it != tpoints_.end(); ++it) {
      Point_Index idx = *it;

      Point_3 p = tpoints_.point(idx);
      FT min_distance = (std::numeric_limits<FT>::max)();
      Vector_3 normal = point_normal_map_[idx];
      Sphere_ID closest_sphere_id = MSMesh::INVALID_SPHERE_ID;

      // Find the sphere with smallest distance to the vertex
      for(auto& sphere : sphere_mesh_->spheres()) {
        Point_3 center = sphere.get_center();
        FT radius = sphere.get_radius();

        // compute Euclidean distance
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
      FT area = point_area_map_[idx];
      // Update the closest sphere
      point_cluster_sphere_map_[idx] = closest_sphere_id;
      if(closest_sphere_id != MSMesh::INVALID_SPHERE_ID) {
        sphere_mesh_->get_sphere(closest_sphere_id).accumulate_cluster_area(area);
        sphere_mesh_->get_sphere(closest_sphere_id).add_cluster_point(idx);
      }
    }
    std::vector<Sphere_ID> sphere_ids_to_remove;
    for(auto& sphere : sphere_mesh_->spheres()) {
      auto& cluster_points = sphere.get_cluster_points_idx();
      if(cluster_points.size() <= 4) {
        for(Point_Index idx : cluster_points) {
          point_cluster_sphere_map_[idx] = MSMesh::INVALID_SPHERE_ID;
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
    Point_3 optimal_center(optimized_sphere_params(0), optimized_sphere_params(1), optimized_sphere_params(2));
    auto [cp, closest_face] = tree_->closest_point_and_primitive(optimal_center);
    FT len = (optimal_center - cp).squared_length();

    Vector_3 normal = (cp - optimal_center) / CGAL::approximate_sqrt(len);

    if(!fast_winding_number_->is_inside(optimal_center))
      normal = -normal; // if the center is outside, flip the normal

    auto [c, r] = compute_shrinking_ball_impl(cp, normal, AccelerationType{});

    sphere.set_center(c);
    sphere.set_radius(r);
  }

  void optimize_single_sphere(MSphere& sphere, bool use_shrinking_ball_correction = false) {
    using EMat = Eigen::Matrix<FT, Eigen::Dynamic, Eigen::Dynamic>;
    using EVec = Eigen::Matrix<FT, Eigen::Dynamic, 1>;
    using EVec3 = Eigen::Matrix<FT, 3, 1>;
    using EVec4 = Eigen::Matrix<FT, 4, 1>;
    using LDLTSolver = Eigen::LDLT<EMat>;

    auto& cluster_points = sphere.get_cluster_points_idx();
    Point_3 center = sphere.get_center();
    FT radius = sphere.get_radius();

    auto nrows = 2 * cluster_points.size();
    EMat J(nrows, 4);
    J.setZero();
    EVec b(nrows);
    b.setZero();
    int idx = 0;
    EVec4 s;

    for(int i = 0; i < 10; i++) {
      idx = 0;

      for(Point_Index id : cluster_points) {
        Point_3 p = tpoints_.point(id);
        EVec3 pos(p.x(), p.y(), p.z());

        // compute sqem energy
        EVec4 lhs = EVec4::Zero();
        FT rhs = 0.0;

        Vector_3 normal = point_normal_map_[id];
        EVec3 normal_eigen(normal.x(), normal.y(), normal.z());
        EVec4 n4(normal_eigen(0), normal_eigen(1), normal_eigen(2), 1.0);
        FT a = CGAL::approximate_sqrt(point_area_map_[id]) / (k_ + FT(1.0));
        lhs += -n4 * a;
        rhs += -1.0 * ((pos - EVec3(s(0), s(1), s(2))).dot(normal_eigen) - s(3)) * a;
        for(Point_Index neighbor_idx : point_knn_map_[id]) {
          Point_3 pos_q = tpoints_.point(neighbor_idx);
          EVec3 pos_q_eigen(pos_q.x(), pos_q.y(), pos_q.z());
          Vector_3 n_q = point_normal_map_[neighbor_idx];
          EVec3 n_q_eigen(n_q.x(), n_q.y(), n_q.z());
          EVec4 n4_q(n_q_eigen(0), n_q_eigen(1), n_q_eigen(2), 1.0);
          FT a_q = CGAL::approximate_sqrt(point_area_map_[neighbor_idx]) / (k_ + FT(1.0));
          lhs += -n4_q * a_q;
          rhs += -1.0 * ((pos_q_eigen - EVec3(s(0), s(1), s(2))).dot(n_q_eigen) - s(3)) * a_q;
        }
        J.row(idx) = lhs;
        b(idx) = rhs;
        ++idx;

        // compute Euclidean energy
        EVec3 d = pos - EVec3(s(0), s(1), s(2));
        FT l = d.norm();

        FT area = CGAL::approximate_sqrt(point_area_map_[id]);
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
    auto& cluster_points = sphere.get_cluster_points_idx();
    Point_3 center = sphere.get_center();
    FT radius = sphere.get_radius();
    FT max_dist = (std::numeric_limits<FT>::min)();
    FT error = 0.0;
    FT area = sphere.get_area();
    if(area <= FT(0.0)) {
      if(verbose_){
        std::cerr << "Warning: sphere with zero area encountered, skipping error computation." << std::endl;
      }
      return FT(0.0);
    }

    for(Point_Index idx : cluster_points) {
      Point_3 p = tpoints_.point(idx);
      FT area = point_area_map_[idx];
      FT dist_sqem = CGAL::scalar_product(p - center, point_normal_map_[idx]) - radius;
      dist_sqem *= dist_sqem; // square the distance

      FT dist_eucl = CGAL::approximate_sqrt((p - center).squared_length()) - radius;
      dist_eucl *= dist_eucl; // square the distance

      FT distance = area * (dist_sqem + lambda_ * dist_eucl);
      error += distance;

      point_error_map_[idx] = distance;
      if(distance > max_dist) {
        max_dist = distance;
        sphere.set_split_point(idx); // update the split point
      }
    }
    error = error / sphere.get_area(); // average error
    sphere.set_error(error);

    return error;
  }

#ifdef CGAL_LINKED_WITH_TBB

  void assign_vertices_to_clusters_parallel() {

    const size_t num_points = tpoints_.end() - tpoints_.begin();
    std::vector<std::pair<Sphere_ID, FT>> point_assignments(num_points);
    auto it = tpoints_.begin();
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, num_points),
                      [&](const tbb::blocked_range<std::size_t>& range) {
                        for(std::size_t i = range.begin(); i != range.end(); ++i) {

                          Point_Index idx = *(it + i);
                          FT area = point_area_map_[idx];
                          Point_3 p = tpoints_.point(idx);
                          FT min_distance = (std::numeric_limits<FT>::max)();
                          Vector_3 normal = point_normal_map_[idx];
                          Sphere_ID closest_sphere_id;

                          // Find the sphere with smallest distance to the vertex
                          for(const auto& sphere : sphere_mesh_->spheres()) {
                            Point_3 center = sphere.get_center();
                            FT radius = sphere.get_radius();

                            // compute Euclidean distance
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
                          point_assignments[i] = {closest_sphere_id, area};
                        }
                      });
    for(std::size_t i = 0; i < num_points; ++i) {
      Point_Index idx = *(it + i);
      auto& [closest_sphere_id, area] = point_assignments[i];

      if(closest_sphere_id != MSMesh::INVALID_SPHERE_ID) {
        point_cluster_sphere_map_[idx] = closest_sphere_id;
        sphere_mesh_->get_sphere(closest_sphere_id).accumulate_cluster_area(area);
        sphere_mesh_->get_sphere(closest_sphere_id).add_cluster_point(idx);
      }
    }
    std::vector<Sphere_ID> sphere_ids_to_remove;
    for(auto& sphere : sphere_mesh_->spheres()) {
      auto& cluster_points = sphere.get_cluster_points_idx();
      if(cluster_points.size() <= 4) {
        for(Point_Index idx : cluster_points) {
          point_cluster_sphere_map_[idx] = MSMesh::INVALID_SPHERE_ID;
        }
        sphere_ids_to_remove.push_back(sphere.get_id());
      }
    }
    for(Sphere_ID id : sphere_ids_to_remove) {
      if(verbose_) {
        std::cout << "Removing sphere with ID: " << id << " due to small cluster size." << std::endl;
      }
      sphere_mesh_->remove(id);
    }
  }

  void compute_shrinking_balls_parallel() {
    const size_t num_points = tpoints_.end() - tpoints_.begin();
    auto it = tpoints_.begin();
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, num_points),
                      [&](const tbb::blocked_range<std::size_t>& range) {
                        for(std::size_t i = range.begin(); i != range.end(); i++) {
                          compute_one_vertex_shrinking_ball(*(it + i));
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

    for(auto it = tpoints_.begin(); it != tpoints_.end(); ++it) {
      Point_Index idx = *it;
      Sphere_ID s1 = point_cluster_sphere_map_[idx];
      if(s1 == MSMesh::INVALID_SPHERE_ID) continue;

      const auto& knn = point_knn_map_[idx]; // std::vector<Point_Index>
      for(const Point_Index& jdx : knn) {
        Sphere_ID s2 = point_cluster_sphere_map_[jdx];
        if(s2 == MSMesh::INVALID_SPHERE_ID || s1 == s2) continue;

        auto& sphere1 = sphere_mesh_->get_sphere(s1);
        auto& sphere2 = sphere_mesh_->get_sphere(s2);
        sphere1.add_neighbor(s2);
        sphere2.add_neighbor(s1);
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

    int to_split_max = (std::min)(int(std::ceil(sphere_mesh_->nb_spheres() * 0.2)), 10);
    for(auto& sphere_id : sorted_sphere_ids) {
      if(sphere_mesh_->nb_spheres() >= desired_number_of_spheres_)
        break;
      auto& sphere = sphere_mesh_->get_sphere(sphere_id);
      if(to_split_max > 0 && sphere.can_split()) {
        for(Sphere_ID neighbor_id : sphere.get_neighbors()) {
          sphere_mesh_->get_sphere(neighbor_id).set_do_not_split(true);
        }
        Point_Index split_vertex = sphere.get_split_point_idx();
        Point_3 center = point_medial_sphere_pos_map_[split_vertex];
        FT radius = point_medial_sphere_radius_map_[split_vertex];
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
  const TriangleMesh& tmesh_;
  Point_set tpoints_;
  GT traits_;
  VPM vpm_;
  FT lambda_;
  std::size_t desired_number_of_spheres_;
  std::size_t nb_samples_;
  unsigned int seed_;
  int max_iteration_;
  int iteration_count_;
  FT total_error_diff_;
  FT last_total_error_;
  FT total_error_;
  FT converged_threshold_;
  int k_ = 10; // number of nearest neighbors
  std::unique_ptr<Tree> tree_;
  std::unique_ptr<MSMesh> sphere_mesh_;
  std::unique_ptr<FWN> fast_winding_number_;
  FT scale_;
  bool verbose_;
  // Property maps
  typename Point_set::template Property_map<face_descriptor> point_from_face_map_;
  typename Point_set::template Property_map<Vector_3> point_normal_map_;
  typename Point_set::template Property_map<FT> point_area_map_;
  typename Point_set::template Property_map<std::vector<Point_Index>> point_knn_map_;
  typename Point_set::template Property_map<Point_3> point_medial_sphere_pos_map_;
  typename Point_set::template Property_map<FT> point_medial_sphere_radius_map_;
  typename Point_set::template Property_map<Sphere_ID> point_cluster_sphere_map_;
  typename Point_set::template Property_map<FT> point_error_map_;
  Vertex_normal_map vertex_normal_map_;
  Face_nb_samples_map face_nb_samples_map_;
  Face_normal_map face_normal_map_;
  Face_area_map face_area_map_;
  Face_centroid_map face_centroid_map_;
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
bool write_PLY(const Medial_skeleton<TriangleMesh, GeomTraits>& skeleton,
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


  /**
   * @ingroup PkgVMASRef
   * loads a medial skeleton from a PLY file.
   *
   * @param skeleton the skeleton
   * @param filepath
   *     Filepath to the PLY file containing the medial skeleton data.
   * @return
   *     `true` if the skeleton was successfully loaded, `false` otherwise.
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

template <typename TriangleMesh, typename GeomTraits>
bool read_PLY(Medial_skeleton<TriangleMesh, GeomTraits>& skeleton,
              const std::string& filepath) {
    return skeleton.read_PLY(filepath);
}

} // namespace IO

} // namespace CGAL
#endif
