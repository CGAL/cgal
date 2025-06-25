#pragma once

#include <CGAL/Graphics_scene.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_operations.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_2.h>
#include <CGAL/Point_container.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/Union_find.h>
#include <CGAL/aff_transformation_tags.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Linear_cell_complex_min_items.h>
#include <functional>
#include <array>
#include <cassert>
#include <unordered_map>
#include <variant>
#include "utils.h"


#include "../query_replace/cmap_query_replace.h"


#include <boost/container/static_vector.hpp>
#include <boost/range/join.hpp>


using uint = uint;

/**
 * @brief A generic 3D point class with basic arithmetic operations
 * @tparam T The numeric type used for coordinates (e.g., int, char)
 * 
 * This class represents a point in 3D space with x, y, and z coordinates.
 * It provides basic arithmetic operations such as addition, subtraction,
 * and division. The class is designed to be used with different numeric
 * types through template parameter T.
 */
template <typename T>
struct GenericPoint {
  T x, y, z;  ///< The x, y, z coordinates of the point

  /// @brief Default constructor, initializes coordinates to 0
  constexpr GenericPoint(): x(0), y(0), z(0) {}

  /// @brief Constructor with explicit coordinates
  /// @param x The x coordinate
  /// @param y The y coordinate
  /// @param z The z coordinate
  constexpr GenericPoint(T x, T y, T z): x(x), y(y), z(z) {}

  /// @brief Constructor from array of 3 integers
  /// @param l Array containing the coordinates [x,y,z]
  constexpr GenericPoint(std::array<int,3> l): x(l[0]), y(l[1]), z(l[2]) {}
  GenericPoint operator+(const GenericPoint& p) const { return { x + p.x, y + p.y, z + p.z }; }
  GenericPoint operator/(const GenericPoint& p) const { return { x / p.x, y / p.y, z / p.z }; }
  GenericPoint operator-(const GenericPoint& p) const { return { x - p.x, y - p.y, z - p.z }; }
  GenericPoint operator-() const { return {static_cast<T>(-x), static_cast<T>(-y), static_cast<T>(-z)}; }
  GenericPoint& operator+=(const GenericPoint& p) { *this = *this + p; return *this; }
  GenericPoint& operator-=(const GenericPoint& p) { *this = *this - p; return *this; }
  GenericPoint& operator/=(const GenericPoint& p) { *this = *this / p; return *this; }
  T& operator[](int i) { assert(i >= 0 && i <= 2); return i == 0 ? x : i == 1 ? y : i == 2 ? z : z; };
  const T& operator[](int i) const { assert(i >= 0 && i <= 2); return i == 0 ? x : i == 1 ? y : i == 2 ? z : z; };
  bool operator==(const GenericPoint& p) const { return p.x == x && p.y == y && p.z == z; }
  bool operator!=(const GenericPoint& p) const { return p.x != x or p.y != y or p.z != z; }

  /// @brief Type conversion operator
  /// @tparam K Target numeric type
  /// @return Point with coordinates converted to type K
  template <typename K>
  operator GenericPoint<K>() const { return {static_cast<K>(x), static_cast<K>(y), static_cast<K>(z)}; }
};

using PointInt = GenericPoint<int>;
using PointChar = GenericPoint<char>;

// x,y,z represents plane normals, with value between -1, 0 or 1, 0 = unset;
// If a volume is owned, multiple constraints/areas can overlap
using AreaId = PointChar;

// Basic CGAL kernel types
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT = Kernel::FT;
using Point = Kernel::Point_3;
using Vector = Kernel::Vector_3;
using Triangle = Kernel::Triangle_3;
using Polyhedron = CGAL::Polyhedron_3<Kernel>;
using Primitive = CGAL::AABB_face_graph_triangle_primitive<Polyhedron>;

using AABB_Traits = CGAL::AABB_traits<Kernel, Primitive>;
using Tree = CGAL::AABB_tree<AABB_Traits>;

/**
 * @brief Enumeration representing the type of volume in the hexahedral mesh
 * 
 * This enum class defines the different states a volume can be in during
 * the hexahedral mesh generation and refinement process.
 */
enum class VolumeType {
  NONE,         // Newly created vol_attribute
  REFINEMENT,   // Previously NONE volumes that was selected for refinement
  ID_EXPANSION, // Expansion of the Identified layer of volumes
  IDENTIFIED    // Volumes identified for refinement
};

/**
 * @brief Thread-safe producer-consumer queue implementation
 * @tparam T The type of elements stored in the queue
 * 
 * This class implements a thread-safe queue using the producer-consumer pattern.
 * It provides synchronized access to a queue where multiple threads can safely
 * produce (push) and consume (wait and pop) items.
 */
template <typename T>
class ProdCons {
public:
  /**
   * @brief Checks if the queue has any items
   * @return true if the queue is not empty, false otherwise
   */
  bool hasNext() const {
    return !queue.empty();
  }

  /**
   * @brief Waits for and retrieves the next item from the queue
   * @return The next item in the queue
   * 
   * If the queue is empty, this function blocks until an item becomes available.
   * When an item is available, it is removed from the queue and returned.
   */
  T waitNextItem(){
    std::unique_lock lock(m);

    while (!hasNext()){
      awake_signal.wait(lock);
    }

    auto item = std::move(queue.front());
    queue.pop();

    return item;
  }

  /**
   * @brief Pushes a new item into the queue
   * @param item The item to push (rvalue reference)
   * 
   * This overload pushes the item into the queue and notifies one waiting consumer.
   */
  void push(T&& item){
    std::unique_lock lock(m);
    queue.push(item);
    awake_signal.notify_one();
  }

  /**
   * @brief Pushes a new item into the queue
   * @param item The item to push (const reference)
   * 
   * This overload copies the item into the queue and notifies one waiting consumer.
   */
  void push(const T& item){
    std::unique_lock lock(m);
    queue.push(item);
    awake_signal.notify_one();
  }

private:
  std::queue<T> queue;           ///< The underlying queue container
  std::mutex m;                  ///< Mutex for thread synchronization
  std::condition_variable awake_signal;  ///< Condition variable for thread notification
};


class LCCItems
{
public:
  /**
   * @brief A wrapper class for Linear Cell Complex dart attributes
   * @tparam Storage The storage type for the Linear Cell Complex
   * 
   * This class defines the attributes associated with vertices, faces, and volumes
   * in the hexahedral mesh.
   */
  template < class Storage >
  struct Dart_wrapper
  {
    /**
     * @brief Attribute class for vertices with point information
     * 
     * Extends CGAL::Cell_attribute_with_point to store vertex coordinates and IDs.
     * Each vertex can be identified by a unique ID, which can be a simple number
     * or a composite of up to three numbers.
     */
    struct VertexAttr : public CGAL::Cell_attribute_with_point<Storage> {
      using Base = CGAL::Cell_attribute_with_point<Storage>;
      VertexAttr():
        Base(),
        id(max_id)
      {}

      VertexAttr(const Point& p):
        Base(p),
        id(max_id)
      {}
      // The id can be a simple number, or a pair up to three numbers
      // Pair identifiers are placeholders until they can be identified with a single size_t number

      static constexpr size_t max_id = std::numeric_limits<size_t>::max();
      size_t id = max_id;
    };

    /**
     * @brief Attribute values for volumes in the hexahedral mesh
     * 
     * Stores information about the state and properties of a volume cell,
     * including its iteration state, type, area ID, ownership, and connected component ID.
     */
    struct VolumeAttrValue {
      static constexpr uint max_cc_id = std::numeric_limits<uint>::max();

      char iteration = -1;
      VolumeType type = VolumeType::NONE;

      AreaId area_id = {0,0,0}; // Used for ghost cells;
      bool owned = true; // Disable this when creating others constraint area
      size_t cc_id = max_cc_id;
    };

    /**
     * @brief Attribute values for faces in the hexahedral mesh
     * 
     * Stores information about face properties including template ID,
     * plane orientation, plane ID, and connected component ID.
     */
    struct FaceAttrValue {
      static constexpr uint max_plane_id = std::numeric_limits<uint>::max();
      static constexpr uint max_cc_id = std::numeric_limits<uint>::max();

      char template_id = 0;
      std::bitset<3> plane;
      size_t plane_id = max_plane_id;
      size_t cc_id = max_cc_id;
    };

    typedef CGAL::Cell_attribute<Storage, FaceAttrValue> FaceAttr;
    typedef CGAL::Cell_attribute<Storage, VolumeAttrValue> VolumeAttr;
    typedef std::tuple<VertexAttr, void, FaceAttr, VolumeAttr> Attributes;
  };
};

using LCCTraits = CGAL::Linear_cell_complex_traits<3,Kernel>;
using LCC = CGAL::Linear_cell_complex_for_combinatorial_map<3,3, LCCTraits, LCCItems>;
using Dart_handle = typename LCC::Dart_handle;
using Vertex_handle = typename LCC::Vertex_attribute_handle;
using size_type = typename LCC::size_type;

using DartInfo = LCCItems::Dart_wrapper<LCC::Storage>;

// Thread local variables are for multithreading, but they are here since they appears in the sequential version
thread_local size_type l_debug_mark;
thread_local size_type l_debug_mark_2;
thread_local int l_thread_id;
thread_local std::unordered_set<Dart_handle> l_handles; // This is left unused

// Do not use this, only to satisfy lcc.is_isomorphic_to
bool operator==(const DartInfo::VolumeAttrValue& first, const DartInfo::VolumeAttrValue& second){
  return true;
}

// Do no use this, only to statisfy lcc.is_isomorphic_to
bool operator==(const DartInfo::FaceAttrValue& first, const DartInfo::FaceAttrValue& second){
  return true;
}

const size_t CONST_SIZE_T_MAX = std::numeric_limits<size_t>::max();

namespace CGAL::HexRefinement::TwoRefinement {
  // 未使用
  template<typename ... Ts>                                                 // (7)
  struct Overload : Ts ... {
    using Ts::operator() ...;
  };
  template<class... Ts> Overload(Ts...) -> Overload<Ts...>;

  /**
  * @brief Enumeration representing the normal direction of a plane in 3D space
  * 
  * This enum is used to specify the orientation of planes in the hexahedral mesh.
  * Each value represents a normal direction along one of the principal axes,
  * with an additional NONE value for unspecified or invalid cases.
  */
  enum PlaneNormal { 
    X,      ///< Normal direction along the X axis
    Y,      ///< Normal direction along the Y axis
    Z,      ///< Normal direction along the Z axis
    NONE = -1  ///< Represents no specific direction or invalid plane
  };

  /**
   * @brief A structure representing a 3D grid for hexahedral mesh generation
   * 
   * This structure defines a regular grid in 3D space, specified by its position,
   * cell size, and dimensions. It provides methods for creating both centered and
   * non-centered grids, as well as specialized cube grid creation.
   */
  struct Grid {
    Point pos;      ///< Starting position (origin) of the grid
    Point size;     ///< Size of each cell in the grid
    PointInt dims;  ///< Number of cells in each dimension (x, y, z)
    PointInt dims_id_start;  ///< Starting indices for each dimension used only for parallelization

    /// @brief Default constructor
    Grid() {}

    /**
     * @brief Constructs a grid with specified parameters
     * @param from Starting position of the grid
     * @param cell_size Size of each cell
     * @param dims Number of cells in each dimension
     */
    Grid(Point from, Point cell_size, PointInt dims)
    : pos(from), size(cell_size), dims(dims) {}

    /**
     * @brief Creates a grid centered around a specified point
     * @param center Center point of the grid
     * @param cell_size Size of each cell
     * @param dims Number of cells in each dimension
     * @return A Grid object centered at the specified point
     */
    static Grid make_centered_grid(Point center, Point cell_size, PointInt dims) {
      Kernel::Vector_3 offset = {
        dims.x / 2 * cell_size.x(),
        dims.y / 2 * cell_size.y(),
        dims.z / 2 * cell_size.z(),
      };

      Point from = center - offset;
      return Grid(from, cell_size, dims);
    }

    /**
     * @brief Creates a grid starting from a specified point
     * @param from Starting position of the grid
     * @param cell_size Size of each cell
     * @param dims Number of cells in each dimension
     * @return A Grid object starting at the specified point
     */
    static Grid make_grid(Point from, Point cell_size, PointInt dims){
      return Grid(from, cell_size, dims);
    }

    /**
     * @brief Creates a cubic grid centered around a point
     * @param center Center point of the cube
     * @param cell_size Size of each cubic cell
     * @param dim Number of cells in each dimension (same for x, y, z)
     * @return A Grid object representing a cubic grid
     */
    static Grid make_centered_cube(Point center, double cell_size, int dim){
      Kernel::Vector_3 offset = {
        dim / 2 * cell_size,
        dim / 2 * cell_size,
        dim / 2 * cell_size,
      };

      Point from = center - offset;
      auto& s = cell_size;
      return Grid(from, {s,s,s}, {dim,dim,dim});
    }

    /**
     * @brief Creates a cubic grid starting from a point
     * @param from Starting position of the cube
     * @param cell_size Size of each cubic cell
     * @param dim Number of cells in each dimension (same for x, y, z)
     * @return A Grid object representing a cubic grid
     */
    static Grid make_cube(Point from, double cell_size, int dim){
      auto& s = cell_size;
      return Grid(from, {s,s,s}, {dim,dim,dim});
    }
  };

  // Identifies which 3-cell should be refined
  using TrimmingFunction = std::function<bool(LCC&, Dart_handle)>;
  using MarkingFunction = std::function<bool(LCC&, Dart_handle)>;

  /**
   * @brief Container for external pattern substitution resources used in hexahedral meshing
   * 
   * This structure holds pattern substituters for both regular and partial templates.
   * These templates are used during the mesh refinement process to replace existing
   * mesh patterns with more refined ones.
   */
  struct ExternalRessources {
    Pattern_substituer<LCC> regular_templates;   ///< Pattern substituter for regular hexahedral templates
    Pattern_substituer<LCC> partial_templates;   ///< Pattern substituter for partial hexahedral templates
  };

  using PlaneCC = std::vector<Dart_handle>; // One dart per face connected components
  using PlaneSet = std::vector<PlaneCC>; // A set of planes 


  /**
   * @brief Core data structure for hexahedral mesh generation and refinement
   * 
   * This structure maintains the state of the hexahedral meshing process,
   * including the Linear Cell Complex (LCC), grid configuration, and various
   * markers used during mesh generation and refinement.
   */
  struct HexMeshingData {
    // Required initialization
    Grid grid;                  ///< Grid configuration defining the mesh structure
    ExternalRessources* ext;    ///< Pointer to external resources for pattern substitution

    // Initialized by the algorithm
    LCC lcc;                    ///< Linear Cell Complex representing the mesh
    size_type identified_mark;  ///< Mark for identified cells
    size_type template_mark;    ///< Mark for template cells
    size_type propagation_face_mark;  ///< Mark for faces during propagation
    size_type debug;           ///< Debug marker
    size_type debug2;          ///< Secondary debug marker
    int level = 0;             ///< Current refinement level
    std::array<PlaneSet, 3> first_face_of_planes;  ///< First faces of each plane set (X, Y, Z)

    /// @brief Default constructor
    HexMeshingData() {}

    /**
     * @brief Initializes the mesh data with external resources and grid configuration
     * @param ext Pointer to external resources
     * @param grid Grid configuration for the mesh
     */
    void init(ExternalRessources* ext, Grid grid) {
      this->ext = ext;
      this->grid = grid;
    }

    /**
     * @brief Fixes invalid dart handles after refinement
     * 
     * After refinement, some dart handles may become invalid. This function
     * iterates over all attributes to replace invalid handles with valid ones.
     * It ensures the consistency of the mesh structure after modifications.
     */
    void fix_dart_storage() {
      // Data that we'll use to know which plane/cc we are looking for, while iterating on attributes
      using DartRef = std::reference_wrapper<Dart_handle>;
      using CCIdToFace = std::unordered_map<size_t, DartRef>;
      using InvalidPlaneSet = std::unordered_map<size_t, CCIdToFace>;
      std::array<InvalidPlaneSet, 3> plane_non_valid_faces;

      const auto valid_face = [&](Dart_handle& face, char p, size_t cc_id){
        auto face_attr = lcc.attribute<2>(face);
        auto& face_info = face_attr->info();
        return lcc.is_dart_used(face) && face_attr != nullptr && face_attr->is_valid() && face_info.plane[p] == true && face_info.cc_id == cc_id;
      };

      bool should_fix = false;
      for (size_t p = 0; p < 3; p++){
        PlaneSet& plane_set = first_face_of_planes[p];
        for (size_t plane_id = 0; plane_id < plane_set.size(); plane_id++){
          PlaneCC& plane_cc = plane_set[plane_id];
          CCIdToFace non_valid_face;

          for (int cc_id = 0; cc_id < plane_cc.size(); cc_id++){
            Dart_handle& dart = plane_cc[cc_id];

            if (!valid_face(dart, p, cc_id)){
              non_valid_face.emplace(cc_id, std::ref(dart));
              should_fix = true;
            }
          }

          if (non_valid_face.size() > 0)
            plane_non_valid_faces[p][plane_id] = std::move(non_valid_face);
        }
      }

      if (!should_fix) return;

      // We don't care if a face has multiple planes, there is atleast one face for each cc that is not a merged face
      // Eitherway we cannot treat a face belonging to multiple planes, because the cc_id is undefined
      const auto get_plane_normal = [&](DartInfo::FaceAttrValue& face) -> std::optional<char>{
        std::bitset<3>& plane = face.plane;
        char sum = plane[0] + plane[1] + plane[2];
        if (sum > 1 or sum == 0) return {};

        for (char i = 0; i < 3; i++){
          if (plane[i]) return i;
        }

        CGAL_assertion_msg(false, "get_plane_normal: Unreachable code accessed");
        return {};
      };

      const auto fix_face_orientation = [&](Dart_handle face, char plane_normal) -> Dart_handle {
        auto n0 = lcc.attribute<0>(face)->point();
        auto n1 = lcc.attribute<0>(lcc.beta(face, 1))->point();
        auto n2 = lcc.attribute<0>(lcc.beta(face, 1, 1))->point();
        // auto n3 = lcc.attribute<0>(lcc.beta(face, 1, 1, 1));

        Vector normal = CGAL::cross_product(n1 - n0, n2 - n1);
        Vector plane_up = [&]() -> Vector {
          switch (plane_normal){
            default: return {1,0,0};
            case 1: return {0,1,0};
            case 2: return {0,0,1};
          }
        }();

        if (plane_up * normal < 0) {
          assert(!lcc.is_free<3>(face));
          return lcc.beta(face, 3);
        }

        return face;
      };

      // Invalidated handles are rare, we can afford scanning all 2 attributes
      auto& attributes = lcc.attributes<2>();
      for (auto it = attributes.begin(); it != attributes.end(); it++){
        auto& info = it->info();
        std::optional<char> normal = get_plane_normal(info);
        if (!normal) continue;

        auto& nvf_plane_set = plane_non_valid_faces[*normal];

        // Get plane
        auto plane_it = nvf_plane_set.find(info.plane_id);
        if (plane_it == nvf_plane_set.end()) continue;
        CCIdToFace& non_valid_faces = plane_it->second;

        // Get cc_id that was invalidated
        auto face_it = non_valid_faces.find(info.cc_id);
        if (face_it == non_valid_faces.end()) continue;

        // If both plane and cc_id match, reassign the handle
        DartRef handle_ptr = face_it->second;
        non_valid_faces.erase(face_it);

        handle_ptr.get() = fix_face_orientation(it->dart(), *normal);

        if (non_valid_faces.size() == 0){
          nvf_plane_set.erase(plane_it);
        }
      }

      auto all_valid = [&](){
        for (int p = 0; p < 3; p++){
        for (auto& non_valid : plane_non_valid_faces[p]){
          CGAL_assertion_msg(non_valid.second.size() == 0, "HexMeshingData::fix_planes_sets, not all stored planes were fixed");
        }}

        for (int p = 0; p < 3; p++){
          auto& plane_set = first_face_of_planes[p];
          for (int plane_id = 0; plane_id < plane_set.size(); plane_id++){
            auto& plane_cc = plane_set[plane_id];
            for (int cc_id = 0; cc_id < plane_cc.size(); cc_id++){
              auto& face = plane_cc[cc_id];
              auto face_attr = lcc.attribute<2>(face);
              CGAL_assertion_msg(valid_face(face, p, cc_id), "HexMeshingData::fix_planes_sets, face was not valid after fix");
            }
          }
        }
      };

      CGAL_assertion_code(all_valid());
    }
  };

  /**
   * @brief Data structure for managing mesh refinement operations
   * 
   * This structure contains the necessary data for performing mesh refinement,
   * including iteration direction, marked elements, and collections of elements
   * that need to be refined.
   */
  struct RefinementData {
    PlaneNormal iteration;  ///< Current iteration direction (X, Y, Z) for refinement

    std::vector<Dart_handle> marked_nodes;  ///< Collection of nodes that have been marked during refinement
    std::vector<Dart_handle> faces_of_plane;  ///< Faces lying on the current refinement plane
    std::vector<Dart_handle> additionnal_volumes_found;  ///< Additional volumes discovered during refinement

    std::vector<Dart_handle> volumes_to_refine;  ///< Collection of volumes that need to be refined
    std::vector<Dart_handle> faces_to_refine;  ///< Collection of faces that need to be refined
    std::vector<Dart_handle> partial_templates_to_refine;  ///< Collection of partial templates that need refinement
  };

  // Debug .. to remove
  ProdCons<int> debug_stream, debug_stream2;

  /**
   * @brief Marks all k-cells that are part of a given i-cell
   * 
   * This　function iterates over all darts of an i-dimensional cell and marks
   * all k-dimensional cells that are part of it.
   * 
   * @tparam i The dimension of the source cell
   * @tparam k The dimension of the cells to mark
   * @param lcc The Linear Cell Complex
   * @param dart A dart handle representing the i-cell
   * @param mark The mark to apply to the k-cells
   */
  template <uint i, uint k>
  void mark_k_cells_of_i_cell(LCC& lcc, Dart_handle dart, size_type mark){
    auto iterator = lcc.darts_of_cell<i, 0>(dart);
    for (auto dit = iterator.begin(), dend = iterator.end(); dit != dend; dit++){
      if (!lcc.is_marked(dit, mark)) {
        lcc.mark_cell<k>(dit, mark);
      }
    }
  }

  /**
   * @brief Internal function to iterate over faces and their edges in a plane
   * 
   * This function performs a breadth-first traversal of faces in a plane, ensuring each face
   * and edge is processed exactly once. It uses marks to track visited elements and applies
   * the provided operations to faces and edges.
   * 
   * @tparam FaceOp Type of the face operation functor
   * @tparam EdgeOp Type of the edge operation functor
   * @param lcc The Linear Cell Complex
   * @param to_explore Queue of dart handles to explore
   * @param face_operation Operation to apply to each face
   * @param edge_operation Operation to apply to each edge
   */
  template <typename FaceOp, typename EdgeOp>
  inline void __plane_for_each_face(LCC& lcc, std::queue<Dart_handle>& to_explore,
                                const FaceOp&& face_operation,
                                const EdgeOp&& edge_operation)
  {
    size_type face_mark = lcc.get_new_mark();
    size_type edge_mark = lcc.get_new_mark();

    while (!to_explore.empty()){
      Dart_handle face = to_explore.front();
      to_explore.pop();

      if (lcc.is_whole_cell_marked<2>(face, face_mark)) continue;
      lcc.mark_cell<2>(face, face_mark);

      auto edges = lcc.darts_of_cell<2,1>(face);
      face_operation(face, edges);

      for (auto it = edges.begin(), end = edges.end(); it != end; it++){
        if (lcc.is_whole_cell_marked<1>(it, edge_mark)) continue;
        lcc.mark_cell<1>(it, edge_mark);
        Dart_handle edge = edge_operation(it);
        if (edge != lcc.null_dart_descriptor)
          to_explore.push(edge);
      }
    }

    lcc.free_mark(face_mark);
    lcc.free_mark(edge_mark);
  }

  /**
   * @brief Retrieves all faces of a hexahedral cell
   * 
   * This function returns an array containing dart handles for all six faces of a hexahedron.
   * The faces are returned in a specific order:
   * 1. The input face (vol)
   * 2-5. The faces adjacent to the edges of the input face
   * 6. The opposite face
   * 
   * @param lcc The Linear Cell Complex
   * @param vol A dart handle representing the initial face of the hexahedron
   * @return std::array<Dart_handle, 6> Array of dart handles for all faces
   */
  std::array<Dart_handle, 6> faces_of_hex(LCC& lcc, Dart_handle vol){
    std::array<Dart_handle, 6> arr;
    int i = 0;

    arr[i++] = vol;

    auto edges = lcc.darts_of_cell<2,1>(vol);
    for (auto it = edges.begin(), end = edges.end(); it != end; it++){
        arr[i++] = lcc.beta(it, 2);
    }

    arr[i++] = lcc.beta(vol, 2, 1, 1, 2);

    return arr;
  }

  /**
   * @brief Iterates over all faces in a plane and applies specified operations
   * 
   * This function performs a breadth-first traversal of faces in a plane, starting from
   * a given face. For each face and its edges, it applies the provided operation functors.
   * The function ensures that each face and edge is processed exactly once by using marks.
   * 
   * @tparam FaceOp Type of the face operation functor
   * @tparam EdgeOp Type of the edge operation functor
   * @param lcc The Linear Cell Complex
   * @param start The starting face's dart handle
   * @param face_operation Operation to apply to each face
   * @param edge_operation Operation to apply to each edge
   */
  template <typename FaceOp, typename EdgeOp>
  void plane_for_each_face(LCC& lcc, std::vector<Dart_handle> starts,
                            const FaceOp&& face_operation,
                            const EdgeOp&& edge_operation)
  {
    std::queue<Dart_handle> to_explore;
    for (Dart_handle start : starts){
      to_explore.push(start);
    }

    __plane_for_each_face(lcc, to_explore,
      std::forward<const FaceOp>(face_operation),
      std::forward<const EdgeOp>(edge_operation));
  }

  /**
   * @brief Gets or creates an attribute for a cell in the Linear Cell Complex
   * 
   * This template function retrieves an existing attribute for a cell, or creates a new one
   * if it doesn't exist. The attribute type and cell dimension are specified through template
   * parameters.
   * 
   * @tparam i The dimension of the cell (0 for vertex, 1 for edge, 2 for face, 3 for volume)
   * @param lcc The Linear Cell Complex
   * @param dart A dart handle representing the cell
   * @return The attribute descriptor for the cell
   */
  template <uint i>
  typename LCC::Attribute_descriptor<i>::type get_or_create_attr(LCC& lcc, Dart_handle dart){
    auto attr = lcc.attribute<i>(dart);

    if (attr == nullptr){
      attr = lcc.create_attribute<i>();
      lcc.set_attribute<i>(dart, attr);
    }
    return attr;
  }

  /**
   * @brief Gets or creates a volume attribute for refinement operations
   * 
   * This function retrieves an existing volume attribute or creates a new one if it doesn't exist.
   * The created volume attribute is initialized with default values suitable for refinement:
   * - type is set to VolumeType::REFINEMENT
   * - iteration is set to -1
   * - cc_id is set to the maximum possible value
   * 
   * @param lcc The Linear Cell Complex
   * @param dart A dart handle representing the volume
   * @return The volume attribute descriptor
   */
  LCC::Attribute_descriptor<3>::type get_or_create_refinement_volume(LCC& lcc, Dart_handle dart){
    auto attr = get_or_create_attr<3>(lcc, dart);

    // Previously NONE volumes are tagged as refined
    if (attr->info().type == VolumeType::NONE)
      attr->info().type = VolumeType::REFINEMENT;

    return attr;

  }

  /**
   * @brief Retrieves all 26-connected neighboring cells of a given cell
   * 
   * This function collects all cells that share a vertex, edge, or face with the given cell.
   * In a 3D grid, a cell can have up to 26 neighbors.
   * 
   * @param lcc The Linear Cell Complex
   * @param dart A dart handle representing the central cell
   * @param include_self_vol Whether to include the input cell in the result (default: false)
   * @return A static vector containing dart handles for all neighboring cells (max size: 27)
   */
  boost::container::static_vector<Dart_handle, 27> cells_26_connectivity(LCC& lcc, Dart_handle dart, bool include_self_vol = false) {
    boost::container::static_vector<Dart_handle, 27> array;
    Dart_handle layers[3] = {
      dart,
      lcc.beta(dart, 3),
      lcc.beta(dart, 2, 1, 1, 2, 3)
    };

    for (int i = 0; i < 3; i++){
      Dart_handle mid_face_dart = layers[i];

      if (mid_face_dart == lcc.null_dart_descriptor) continue;

      if (include_self_vol || i != 0) array.push_back(mid_face_dart);

      auto edges = lcc.darts_of_cell<2, 1>(mid_face_dart);
      for (auto edge_it = edges.begin(), end = edges.end(); edge_it != end; edge_it++){

        if (!lcc.is_free<3>(lcc.beta(edge_it, 2))){
          Dart_handle other_face = lcc.beta(edge_it, 2, 3, 2);
          if (other_face != lcc.null_dart_descriptor)
            array.push_back(other_face);

          // TODO If nesting ..
          if (!lcc.is_free<3>(lcc.beta(other_face, 1, 2))){
            Dart_handle side_face = lcc.beta(other_face, 1, 2, 3, 2);
            if (side_face != lcc.null_dart_descriptor)
              array.push_back(side_face);
          }
        }

      }
    }

    return array;
  }

  /**
   * @brief Finds the origin dart of a 3-template pattern
   * 
   * This function searches for the origin dart of a 3-template pattern by traversing
   * through the marked faces.
   * 
   * @param lcc The Linear Cell Complex
   * @param marked_face A dart handle representing the marked face to start the search from
   * @param template_mark The mark used to identify template faces
   * @return A dart handle pointing to the origin of the 3-template pattern, or null if not found
   */
  Dart_handle find_3_template_origin(LCC& lcc, Dart_handle marked_face, size_type template_mark) {

    auto _edges_count = lcc.darts_of_cell<2,1>(marked_face).size();
    assert(_edges_count == 6);

    Dart_handle dart = marked_face;

    // Get the origin dart : Find the two unmarked node on the face
    // since the 3 template is created by adjacent two 2-templates
    bool found = false;
    for (int i = 0; i < 6; i++)
    {
      if (!lcc.is_marked(dart, template_mark)
        && !lcc.is_marked(lcc.beta(dart, 1), template_mark)
        && !lcc.is_marked(lcc.beta(dart, 1, 1), template_mark))
      {
        assert(lcc.is_marked(lcc.beta(dart, 1, 1, 1), template_mark));
        assert(lcc.is_marked(lcc.beta(dart, 1, 1, 1, 1), template_mark));
        assert(lcc.is_marked(lcc.beta(dart, 1, 1, 1, 1, 1), template_mark));
        assert(lcc.beta(dart, 1, 1, 1, 1, 1, 1) == dart);
        found = true;
        break;
      }

      dart = lcc.beta(dart, 1);
    }

    assert(found);

    return dart;
  }

  /**
   * @brief Marks a face without checking its current mark status
   * 
   * This function marks all darts of a 2-dimensional cell (face) with the specified mark.
   * Unlike safe marking functions, this function does not verify if the face is already marked.
   * 
   * @param lcc The Linear Cell Complex
   * @param dart A dart handle representing the face to mark
   * @param mark The mark to apply to the face
   */
  void mark_face_unchecked(LCC& lcc, Dart_handle dart, size_type mark){
    auto iterator = lcc.darts_of_cell<2, 1>(dart);
    for (auto dit = iterator.begin(), dend = iterator.end(); dit != dend; dit++) {
      lcc.mark(dit, mark);

      if (!lcc.is_free<3>(dit))
        lcc.mark(lcc.beta<3>(dit), mark);
    }

    assert(lcc.is_whole_cell_marked<2>(dart, mark));
  }

  /**
   * @brief Checks if a half-face is marked with the specified mark
   * 
   * This function checks if all darts of a half-face are marked with the given mark.
   * A half-face consists of the 4 darts that can be reached from the specified dart
   * without using beta<3> operations, out of the total 8 darts that include beta<3>.
   * 
   * @param lcc The Linear Cell Complex
   * @param dart A dart handle representing the half-face to check
   * @param mark The mark to check for
   * @return true if all darts of the half-face are marked, false otherwise
   */
  bool is_half_face_marked(const LCC& lcc, LCC::Dart_const_handle dart, size_type mark ){
    auto iterator = lcc.darts_of_cell<2, 1>(dart);
    for (auto dit = iterator.begin(), dend = iterator.end(); dit != dend; dit++) {
      if (!lcc.is_marked(dit, mark)) return false;
    }

    return true;
  }

  /**
   * @brief Marks a half-face without checking its current mark status
   * 
   * This function marks all darts of a half-face with the specified mark.
   * A half-face consists of the 4 darts that can be reached from the specified dart
   * without using beta<3> operations, out of the total 8 darts that include beta<3>.
   * Unlike safe marking functions, this function does not verify if the half-face is already marked.
   * 
   * @param lcc The Linear Cell Complex
   * @param dart A dart handle representing the half-face to mark
   * @param mark The mark to apply to the half-face
   */
  void mark_half_face_unchecked(LCC& lcc, Dart_handle dart, size_type mark){
    auto iterator = lcc.darts_of_cell<2, 1>(dart);

    for (auto dit = iterator.begin(), dend = iterator.end(); dit != dend; dit++) {
      lcc.mark(dit, mark);
    }

    assert(is_half_face_marked(lcc, dart, mark));
  }

  /**
   * @brief Internal function to find the adjacent face on a given plane
   * 
   * This function searches for a face adjacent to the given edge that lies on the specified plane.
   * It traverses through the mesh using beta operations to find faces that share the edge
   * and checks if they belong to the target plane. The function may fail if the adjacent face
   * is not directly accessible from beta3 to the current face. To complete the search, you should
   * call this function with both edge and lcc.beta<3>(edge)
   * 
   * @param lcc The Linear Cell Complex
   * @param plane The plane normal direction to search for
   * @param edge A dart handle representing the edge to find adjacent faces for
   * @return A dart handle to the adjacent face on the plane, or null if not found
   */
  Dart_handle __adjacent_face_on_plane(LCC& lcc, PlaneNormal plane, Dart_handle edge){
    Dart_handle other_face = lcc.beta(edge, 2);
    auto this_face_handle = lcc.attribute<2>(edge);
    auto other_face_handle = lcc.attribute<2>(other_face);

    assert(this_face_handle != nullptr);

    // Early exit if the beta(2) face is already on the plane
    if (other_face_handle != nullptr
      && other_face_handle->info().plane[plane])
      return other_face;

    bool found = false;
    for (int i = 0; i < 5 && !found; i++){

      if (lcc.is_free<3>(other_face)) {
        // We hit the grid border, there are no adjacent face possible to that edge
        break;
      }

      other_face = lcc.beta(other_face, 3, 2);
      other_face_handle = lcc.attribute<2>(other_face);
      auto other_vol_handle = lcc.attribute<3>(other_face);

      // Exit if we fall back on the same face or outside of the domain
      if (this_face_handle == other_face_handle)
        return lcc.null_dart_descriptor;

      found = other_face_handle != nullptr
        && other_face_handle->info().plane[plane];
    }

    return found ? other_face : lcc.null_dart_descriptor;
  }

  /**
   * @brief Internal function to find the adjacent face on a given plane with additional volume collection
   * 
   * This overload of __adjacent_face_on_plane additionally collects volumes that couldn't be reached
   * directly from the plane during the search. It searches for a face adjacent to the given edge
   * that lies on the specified plane, and stores any encountered volumes that are inside the
   * refinement domain but not directly accessible from the plane.
   * 
   * @param lcc The Linear Cell Complex
   * @param plane The plane normal direction to search for
   * @param edge A dart handle representing the edge to find adjacent faces for
   * @param additional_volumes Vector to store volumes encountered during the search
   * @return A dart handle to the adjacent face on the plane, or null if not found
   */
  Dart_handle __adjacent_face_on_plane(LCC& lcc, PlaneNormal plane, Dart_handle edge, std::vector<Dart_handle>& additional_volumes){
    Dart_handle other_face = lcc.beta(edge, 2);
    auto this_face_handle = lcc.attribute<2>(edge);
    auto other_face_handle = lcc.attribute<2>(other_face);
    boost::container::static_vector<Dart_handle, 5> __additional_volumes;

    assert(this_face_handle != nullptr);

    // Early exit if the beta(2) face is already on the plane
    if (other_face_handle != nullptr && other_face_handle->info().plane[plane])
      return other_face;

    bool found = false;
    for (int i = 0; i < 5 && !found; i++){

      if (lcc.is_free<3>(other_face)) {
        // We hit the grid border, there are no adjacent face possible to that edge
        break;
      }

      other_face = lcc.beta(other_face, 3, 2);
      other_face_handle = lcc.attribute<2>(other_face);
      auto other_vol_handle = lcc.attribute<3>(other_face);

      // Exit if we fall back on the same face or outside of the domain
      if (this_face_handle == other_face_handle)
        break;

      found = other_face_handle != nullptr && other_face_handle->info().plane[plane];

      // Add the incident volume only if it is inside the refinement domain
      if (!found && other_vol_handle != nullptr)
        __additional_volumes.push_back(other_face);
    }

    additional_volumes.insert(additional_volumes.end(), __additional_volumes.begin(), __additional_volumes.end());

    return found ? other_face : lcc.null_dart_descriptor;
  }

  /**
   * @brief Finds the adjacent face on a given plane
   * 
   * This function searches for a face adjacent to the given edge that lies on the specified plane.
   * It uses the internal __adjacent_face_on_plane function and handles cases where the adjacent face
   * might not be directly accessible from the current orientation. If the first search fails,
   * it tries again with the opposite orientation (beta<3> of the edge).
   * 
   * @param lcc The Linear Cell Complex
   * @param plane The plane normal direction to search for
   * @param edge A dart handle representing the edge to find adjacent faces for
   * @return A dart handle to the adjacent face on the plane, or null if not found
   */
  Dart_handle adjacent_face_on_plane(LCC& lcc, PlaneNormal plane, Dart_handle edge){
    Dart_handle other_face = __adjacent_face_on_plane(lcc, plane, edge);

    // This is needed to iterate on 3 templates that are on a grid border, missing a connected volume
    // And preventing iteration the first way
    if (other_face == lcc.null_dart_descriptor && !lcc.is_free<3>(edge))
      other_face = __adjacent_face_on_plane(lcc, plane, lcc.beta<3>(edge));

    return other_face;
  }

  /**
   * @brief Finds the adjacent face on a given plane with additional volume collection
   * 
   * This overload of adjacent_face_on_plane additionally collects volumes that couldn't be reached
   * directly from the plane during the search. It uses the internal __adjacent_face_on_plane function
   * with volume collection and handles cases where the adjacent face might not be directly accessible
   * from the current orientation. If the first search fails, it tries again with the opposite orientation.
   * 
   * @param lcc The Linear Cell Complex
   * @param plane The plane normal direction to search for
   * @param edge A dart handle representing the edge to find adjacent faces for
   * @param additionnal_volumes Pointer to vector to store volumes encountered during the search
   * @return A dart handle to the adjacent face on the plane, or null if not found
   */
  Dart_handle adjacent_face_on_plane(LCC& lcc, PlaneNormal plane, Dart_handle edge, std::vector<Dart_handle>* additionnal_volumes){
    Dart_handle other_face = __adjacent_face_on_plane(lcc, plane, edge, *additionnal_volumes);

    if (other_face == lcc.null_dart_descriptor && !lcc.is_free<3>(edge))
      other_face = __adjacent_face_on_plane(lcc, plane, lcc.beta<3>(edge), *additionnal_volumes);

    return other_face;
  }

  /**
   * @brief Finds all faces around a node that lie on a specific plane
   * 
   * This function collects all faces that are adjacent to a given node and lie on the specified plane.
   * There may be up to 8 faces adjacent to a vertex, where each dart handle belongs to a unique edge
   * and a unique face on the plane around the selected node. The function performs a traversal around
   * the node to gather all connected faces on the plane.
   * 
   * @param lcc The Linear Cell Complex
   * @param rdata Refinement data containing the current iteration direction
   * @param node A dart handle representing the node to find faces around
   * @return A static vector containing dart handles for all faces around the node on the plane (max size: 8)
   */
  boost::container::static_vector<Dart_handle, 8>
  plane_faces_around_node(LCC& lcc,
                                    RefinementData &rdata,
                                    Dart_handle node){
    boost::container::static_vector<Dart_handle, 8> arr;

    // Add initial face
    arr.push_back(node);
    auto this_face_attr = lcc.attribute<2>(node);

    auto turn_around_edge = [&](Dart_handle edge){
      Dart_handle adjacent_face = adjacent_face_on_plane(lcc, rdata.iteration, edge);
      auto other_face_attr = lcc.attribute<2>(adjacent_face);

      while (adjacent_face != lcc.null_dart_descriptor && this_face_attr != other_face_attr){
        arr.push_back(adjacent_face);

        int f = lcc.belong_to_same_cell<0>(adjacent_face, node) ? 0 : 1;
        adjacent_face = adjacent_face_on_plane(lcc, rdata.iteration, lcc.beta(adjacent_face, f));
        if (adjacent_face != lcc.null_dart_descriptor){
          // Check if the dart is expectedly turning around the node
          assert(lcc.belong_to_same_cell<0>(adjacent_face, node) or lcc.belong_to_same_cell<0>(lcc.other_extremity(adjacent_face), node));
        }
        other_face_attr = lcc.attribute<2>(adjacent_face);
      } ;

      return adjacent_face;
    };

    // Turn around the first edge, gathering faces encountered
    Dart_handle end = turn_around_edge(node);

    // If the iteration doesn't loop back to the start, do it again on the opposite edge to that node
    if (lcc.attribute<2>(end) != this_face_attr){
      end = turn_around_edge(lcc.beta(node, 0));
    }

    // Check if all found faces are unique, if they are not, this means that the LCC was not properly refined
    auto assertion = [&](){
      std::unordered_set<DartInfo::FaceAttrValue*> faces_set;
      for (Dart_handle face : arr){
        auto f = lcc.attribute<2>(face);
        CGAL_postcondition_msg(f != nullptr, "plane_faces_around_node: returned array contains nullptr, Is the refinement correctly done ?");
        CGAL_postcondition_msg(true, "plane_faces_around_node: array contains a duplicate face, Is the refinement correctly done ?");
        assert(faces_set.count(&f->info()) == 0);
        faces_set.insert(&f->info());
      }
    };

    CGAL_postcondition_code(assertion());

    return arr;
  }

  /**
   * @brief Asserts that all dart attributes in the given arrays are unique
   * 
   * This template function checks that all dart attributes of dimension i in the provided
   * arrays are unique (i.e., no duplicate attributes exist across the arrays). It uses
   * a hash map to track encountered attributes and raises assertions if duplicates are found.
   * This function is typically used for debugging to ensure data integrity during mesh operations.
   * 
   * @tparam i The dimension of the cell attributes to check (0 for vertex, 1 for edge, 2 for face, 3 for volume)
   * @tparam DartArray Variadic template parameter for the dart handle arrays
   * @param lcc The Linear Cell Complex
   * @param array Variadic parameter containing arrays of dart handles to check
   */
  template <uint i, typename... DartArray>
  void assert_dart_attr_are_unique(LCC& lcc, DartArray... array){
    auto assertion = [&](){
      std::unordered_map<void*, int> attributes;
      ([&]{
        int index = 0;
        for (auto dart : array){
          auto attr = lcc.attribute<i>(dart);
          CGAL_assertion_msg(attr != nullptr, "assert_dart_attr_are_unique, nullptr found in arrays ");
          auto ptr = &attr->info();
          CGAL_assertion_msg(attributes.count(ptr) == 0, "assert_dart_attr_are_unique, duplicate attribute found in arrays");
          attributes[ptr] = index;
          index++;
        }
      }(),...);
    };

    CGAL_assertion_code(assertion());
  }

  /**
   * @brief Asserts that all faces of the plane are valid
   * 
   * This function validates that all faces in the current refinement plane have consistent
   * template IDs. It checks that each face's template_id matches the number of marked nodes
   * on that face. This validation ensures that the mesh refinement process has correctly
   * assigned template IDs to faces based on their marked node patterns.
   * 
   * @param hdata Hexahedral meshing data containing the mesh and marks
   * @param rdata Refinement data containing the faces of the current plane
   */
  void assert_faces_of_plane_valid(HexMeshingData& hdata, RefinementData& rdata){
    LCC& lcc = hdata.lcc;
    assert_dart_attr_are_unique<2>(hdata.lcc, rdata.faces_of_plane);

    auto assertion = [&](){
      for (Dart_handle face : rdata.faces_of_plane){
        auto& attr = lcc.attribute<2>(face)->info();
        auto nodes = lcc.darts_of_cell<2, 0>(face);
        int marked = 0;
        for (auto it = nodes.begin(), end = nodes.end(); it != end; it++){
          if (lcc.is_marked(it, hdata.template_mark)) marked++;
        }
        CGAL_assertion_msg(attr.template_id == marked, "assert_faces_of_plane_valid: Found an invalid face");
      }
    };

    CGAL_assertion_code(assertion());
  }

  // hdata は未使用
  /**
   * @brief Asserts that all faces in the Linear Cell Complex are quadrilateral
   * 
   * This function validates that every face in the mesh has exactly 4 edges,
   * ensuring the mesh consists only of quadrilateral faces. It iterates through
   * all faces in the Linear Cell Complex and checks that each face has exactly
   * 4 edges. If any face is found with a different number of edges, an assertion
   * is triggered.
   * 
   * @param lcc The Linear Cell Complex to validate
   * @param hdata Hexahedral meshing data
   */
  void assert_all_faces_are_quadrilateral(LCC &lcc, HexMeshingData& hdata){
    auto assertion = [&](){
      auto iterator = lcc.one_dart_per_cell<2>();
      for (auto it = iterator.begin(); it != iterator.end(); it++){
        auto edges = lcc.darts_of_cell<2, 0>(it);
        CGAL_assertion_msg(edges.size() == 4, "assert_all_faces_are_quadrilateral: Found a non quadrilateral face");
      }
    };

    CGAL_assertion_code(assertion());
  }

  /**
   * @brief Asserts that all volumes in the Linear Cell Complex are hexahedral
   * 
   * This function validates that every volume in the mesh has exactly 24 edges
   * (4 edges per face × 6 faces), ensuring the mesh consists only of hexahedral
   * volumes. It iterates through all volumes in the Linear Cell Complex and
   * checks that each volume has exactly 24 edges. If any volume is found with
   * a different number of edges, an assertion is triggered.
   * 
   * @param lcc The Linear Cell Complex to validate
   */
  void assert_all_volumes_are_hexes(LCC &lcc){
    auto assertion = [&](){
      auto iterator = lcc.one_dart_per_cell<3>();
      for (auto it = iterator.begin(); it != iterator.end(); it++){
        auto edges = lcc.darts_of_cell<3, 0>(it);
        CGAL_assertion_msg(edges.size() == (4 * 6), "assert_all_volumes_are_hexes: Found a non hexahedral volume");
      }
    };

    CGAL_assertion_code(assertion());
  }

  template <typename HexData>
  void thread_number_vertex_in_edge(HexData& hdata,
    Dart_handle node, Dart_handle extremity0, Dart_handle extremity1){}

  template <typename HexData>
  void thread_number_vertex_in_1t_face(HexData& hdata, Dart_handle node) {}

  template <typename HexData>
  void thread_number_vertex_in_1t_vol(HexData& hdata, Dart_handle v_signature_start) {}

  template <typename HexData>
  void thread_join_3_template_vertex__pair(HexData& hdata, Dart_handle edge) {}

  template <typename HexData>
  void thread_join_3_template_vertex__pairpair(HexData& hdata, Dart_handle edge) {}

  /**
   * @brief Expands the marking of a 0-cell (node) to affect surrounding faces
   * 
   * This function handles the propagation of node marking by finding all faces
   * around the marked node and updating their template IDs. When a node is marked,
   * all faces that share this node need to have their template_id incremented.
   * Faces that transition from template_id 0 to 1 are added to the refinement list,
   * and faces with template_id 2 or 3 are added to the faces_to_check queue for
   * further processing.
   * 
   * @param lcc The Linear Cell Complex
   * @param rdata Refinement data containing the current iteration direction and collections
   * @param faces_to_check Queue of faces that need to be checked for template conflicts
   * @param edge A dart handle representing the marked node (0-cell)
   */
  void __expand_0_cell_marking(LCC &lcc, RefinementData &rdata, std::queue<Dart_handle>& faces_to_check, Dart_handle &edge) {
    auto faces = plane_faces_around_node(lcc, rdata, edge);
    int s = faces.size();

    for (Dart_handle face : faces){
      assert( lcc.attribute<2>(face) != nullptr);
      auto& face_attr =  lcc.attribute<2>(face)->info();

      // If the face didn't have any template before, it will have one, so add it in faces to refine
      if (face_attr.template_id == 0) {
        rdata.faces_to_refine.push_back(face);
      }

      face_attr.template_id++;
      assert(face_attr.template_id <= 4);

      if (face_attr.template_id == 2 || face_attr.template_id == 3)
        faces_to_check.push(face);
    }
  }

  /**
   * @brief Fixes adjacent 3-templates by transforming them into 4-templates when they share an unmarked node
   * 
   * This function handles the case where two adjacent faces both have template_id 3 and share
   * an unmarked node. In such cases, the refinement pattern becomes impossible to resolve,
   * so the function marks the shared unmarked node to force both faces to become 4-templates.
   * 
   * The function finds the unmarked edge on the given face, checks if the adjacent face
   * on the same plane also has template_id 3, and if so, marks the unmarked node to resolve
   * the conflict.
   * 
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and marks
   * @param rdata Refinement data containing the current iteration direction
   * @param faces_to_check Queue of faces that need to be checked for template conflicts
   * @param face A dart handle representing the face to check for adjacent 3-template conflicts
   * @return true if a fix was applied (node was marked), false otherwise
   */
  bool fix_adjacent_3_templates(HexMeshingData &hdata, RefinementData &rdata, std::queue<Dart_handle>& faces_to_check, Dart_handle &face)
  {
    // Transform nearby 3 templates into two 4 templates
    // if they both share the unmarked node

    LCC& lcc = hdata.lcc;

    auto &face_attr = lcc.attribute<2>(face)->info();

    auto edges = lcc.darts_of_cell<2, 1>(face);

    // find the unmarked edge;
    bool found = false;
    Dart_handle edges_unmarked[2];

    for (auto it = edges.begin(); it != edges.end() && !found; it++)
    {
      if (!lcc.is_marked(it, hdata.template_mark)){
        found = true;
        edges_unmarked[0] = it;
      }
    }

    assert(found);
    edges_unmarked[1] = lcc.beta(edges_unmarked[0], 0);

    bool is_impossible = false;
    for (int i = 0; i < 2; i++){
      Dart_handle edge = edges_unmarked[i];
      Dart_handle other_face = adjacent_face_on_plane(lcc, rdata.iteration, edge);

      if (other_face == lcc.null_dart_descriptor) continue;

      auto other_face_handle = lcc.attribute<2>(other_face);

      assert(other_face_handle != nullptr);

      if (other_face_handle->info().template_id == 3){
        is_impossible = true;
        break;
      }

    }

    if (is_impossible) {
      lcc.mark_cell<0>(edges_unmarked[0], hdata.template_mark);
      rdata.marked_nodes.push_back(edges_unmarked[0]);
      __expand_0_cell_marking(lcc, rdata, faces_to_check, edges_unmarked[0]);

      return true;
    }

    return false;
  }

  /**
   * @brief Fixes diagonal marking patterns on faces with template_id=2 by converting them to template_id=4
   * 
   * This function handles the specific case where a face has template_id=2 with marks placed
   * diagonally (e.g., nodes 1 and 3 are marked but node 2 is not). In such cases, the
   * refinement pattern becomes impossible to resolve, so the function marks all unmarked
   * nodes on the face to create a complete 4-template pattern.
   * 
   * The function examines the face's edges in sequence and detects diagonal marking patterns
   * where marked nodes are not consecutively chained. If such a pattern is found, it marks
   * all remaining unmarked nodes and sets the face's template_id to 4.
   * 
   * @pre The given face must have template_id=2
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and marks
   * @param rdata Refinement data containing the current iteration direction
   * @param faces_to_check Queue of faces that need to be checked for template conflicts
   * @param face A dart handle representing the face with template_id=2 to check for diagonal marking
   * @return true if a fix was applied (diagonal pattern was converted to template_id=4), false if the marks were already consecutive
   */
  bool fix_mark_connectivity(HexMeshingData &hdata, RefinementData &rdata, std::queue<Dart_handle>& faces_to_check, Dart_handle face){
    LCC& lcc = hdata.lcc;

    auto edges = lcc.darts_of_cell<2, 1>(face);

    bool connected = true;
    Dart_handle edge1 = face;
    if (lcc.is_marked(edge1, hdata.template_mark)
      && !lcc.is_marked(lcc.beta(edge1, 1), hdata.template_mark)
      &&  lcc.is_marked(lcc.beta(edge1, 1, 1), hdata.template_mark))
    { connected = false;}

    Dart_handle edge2 = lcc.beta(face, 1);
    if (lcc.is_marked(edge2, hdata.template_mark)
      && !lcc.is_marked(lcc.beta(edge2, 1), hdata.template_mark)
      &&  lcc.is_marked(lcc.beta(edge2, 1, 1), hdata.template_mark))
    { connected = false;}

    if (connected) return false;

    // mark the face
    for (auto edge = edges.begin(); edge != edges.end(); edge++) {
      if (lcc.is_marked(edge, hdata.template_mark)) continue;

      lcc.mark_cell<0>(edge, hdata.template_mark);
      rdata.marked_nodes.push_back(edge);
      __expand_0_cell_marking(lcc, rdata, faces_to_check, edge);
    }

    auto &face_attr = lcc.attribute<2>(face)->info();
    face_attr.template_id = 4;

    return true;

  }

  /**
   * @brief Fixes impossible template cases by resolving conflicts in 2-templates and 3-templates
   * 
   * This function identifies and resolves impossible refinement patterns that cannot be
   * processed by the standard template substitution. It handles two main cases:
   * 1. Faces with template_id=2 that have diagonal marking patterns (non-consecutive marks)
   * 2. Faces with template_id=3 that are adjacent to other 3-templates sharing unmarked nodes
   * 
   * The function performs an iterative process: it first examines all existing faces in the
   * current plane, then processes any new faces that are added to the faces_to_check queue
   * during the fixing process. This continues until no more faces need to be checked.
   * 
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and marks
   * @param rdata Refinement data containing the faces of the current plane and collections
   * @return The total number of fixes applied (sum of diagonal 2-templates and adjacent 3-templates fixed)
   */
  int fix_impossible_cases(HexMeshingData &hdata, RefinementData &rdata){
    int fix_c_count = 0, fix_3_count = 0;

    // Condition for checking faces : 2 or 3 templates.
    std::queue<Dart_handle> faces_to_check;

    // First iteration

    // faces_of_plane will grow, we only want to examine only faces that existed right now
    int faces_end = rdata.faces_of_plane.size();

    for (int i = 0; i < faces_end; i++){
      Dart_handle face = rdata.faces_of_plane[i];
      auto& face_attr = hdata.lcc.attribute<2>(face)->info();

      if (face_attr.template_id == 2 && fix_mark_connectivity(hdata, rdata, faces_to_check, face))
        fix_c_count++;

      if (face_attr.template_id == 3 && fix_adjacent_3_templates(hdata, rdata, faces_to_check, face))
        fix_3_count++;
    }

    // Repeat until there are no more faces to check
    while (!faces_to_check.empty()){
      Dart_handle front_face = faces_to_check.front();
      faces_to_check.pop();

      auto& face_attr = hdata.lcc.attribute<2>(front_face)->info();

      if (face_attr.template_id == 2 && fix_mark_connectivity(hdata, rdata, faces_to_check, front_face))
        fix_c_count++;

      if (face_attr.template_id == 3 && fix_adjacent_3_templates(hdata, rdata, faces_to_check, front_face))
        fix_3_count++;
    }

    std::cout << "Diagonal 2 templates fixed: " << fix_c_count << std::endl;
    std::cout << "Neighboring 3 templates repaired: " << fix_3_count << std::endl;

    return fix_c_count + fix_3_count;
  }

  /**
   * @brief Cleans up a 3-template by removing temporary elements and merging face attributes
   * 
   * This function performs the cleanup phase after a 3-template refinement operation.
   * It removes temporary elements created during the refinement process and merges
   * the face attributes of the two neighboring faces that were split during the
   * 3-template creation.
   * 
   * The function performs the following operations:
   * 1. Preserves face attributes before they are potentially destroyed
   * 2. Inserts barycenters in the lower edge and contracts them
   * 3. Contracts the two remaining 2-dart faces that were created during refinement
   * 4. Removes the temporary volume created by the partial template substitution
   * 5. Sews the two neighboring volumes together if both faces are valid
   * 6. Merges the plane information from the two original faces into a single face attribute
   * 
   * @tparam HexData Type of the hexahedral meshing data structure
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex
   * @param origin_dart The origin dart of the 3-template volume to be cleaned up
   * @param upper_edge The upper edge of the 3-template
   * @param lower_edge The lower edge of the 3-template
   * @param face1 Reference to the first neighboring face (may be modified)
   * @param face2 Reference to the second neighboring face (may be modified)
   */
  template <typename HexData>
  void clean_up_3_template(HexData &hdata, const Dart_handle &origin_dart, const Dart_handle &upper_edge, const Dart_handle lower_edge, Dart_handle &face1, Dart_handle &face2)
  {
    LCC& lcc = hdata.lcc;

    // Ne pas maintenir des refs sur les attributs, ils vont disparaitre
    DartInfo::FaceAttrValue face1_attr, face2_attr;
    if (lcc.attribute<2>(face1) != nullptr) face1_attr = lcc.attribute<2>(face1)->info();
    if (lcc.attribute<2>(face2) != nullptr) face2_attr = lcc.attribute<2>(face2)->info();

    Dart_handle lower_mid_1 = lcc.insert_barycenter_in_cell<1>(lower_edge);
    Dart_handle lower_mid_2 = lcc.beta(lower_mid_1, 2, 1);

    thread_join_3_template_vertex__pairpair(hdata, lower_edge);

    lcc.contract_cell<1>(lower_mid_1);
    lcc.contract_cell<1>(lower_mid_2);

    // Contract the two remaining 2-darts faces
    Dart_handle face_to_remove_1 = lcc.beta(upper_edge, 2, 1);
    Dart_handle face_to_remove_2 = lcc.beta(upper_edge, 2, 1, 1, 2, 1, 2);

    assert(lcc.darts_of_orbit<1>(face_to_remove_1).size() == 2);
    assert(lcc.darts_of_orbit<1>(face_to_remove_2).size() == 2);

    // Contract the two remaining 2-darts faces
    lcc.contract_cell<2>(face_to_remove_1);
    lcc.contract_cell<2>(face_to_remove_2);

    // Remove the created volume from the partial query replace and sew the two neighboring volumes
    lcc.remove_cell<3>(origin_dart);

    if (face1 != lcc.null_dart_descriptor && face2 != lcc.null_dart_descriptor){
      lcc.sew<3>(face1, face2);
      assert(lcc.attribute<2>(face1) == lcc.attribute<2>(face2));
    }

    // Requires at least one face for merging/adding 2-attributes
    std::bitset<3> merged_planes = face1_attr.plane | face2_attr.plane;
    auto merge_face = face1 != lcc.null_dart_descriptor ? face1 : face2;

    assert(merge_face != lcc.null_dart_descriptor);

    // Merge the two previous face attributes bitsets
    auto& merge_face_attr = get_or_create_attr<2>(lcc, merge_face)->info();
    // TODO don't entirely reset...
    merge_face_attr = DartInfo::FaceAttrValue();
    merge_face_attr.plane = merged_planes;
  }

  /**
   * @brief Refines all 3-template patterns in the current refinement data
   *
   * This function processes all partial templates (3-templates) collected in the
   * refinement data and applies the necessary topological and combinatorial operations
   * to convert each 3-template into two regular hexahedral volumes. For each 3-template:
   *   - The function identifies the origin dart and its adjacent volume
   *   - It locates the relevant edges and faces for both volumes
   *   - It inserts barycenters into the upper edge and contracts the resulting edges
   *   - It calls clean_up_3_template to remove temporary elements and merge face attributes
   *   - The process is repeated for both the original and adjacent volumes
   *
   * @tparam HexData Type of the hexahedral meshing data structure
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex
   * @param rdata Refinement data containing the list of partial templates to refine
   * @return The total number of volumic substitutions performed (should be twice the number of 3-templates)
   */
  template <typename HexData>
  int refine_3_template(HexData &hdata, RefinementData& rdata)
  {
    int nb_sub = 0;
    for (auto origin_dart : rdata.partial_templates_to_refine){
      LCC& lcc = hdata.lcc;

      int nb_edges = lcc.darts_of_cell<2,1>(origin_dart).size();
      assert(nb_edges == 3);

      Dart_handle vol2_origin_dart = lcc.beta(origin_dart, 3);

      Dart_handle upper_edge = lcc.beta(origin_dart, 0); //?
      Dart_handle vol2_upper_edge = lcc.beta(upper_edge, 3);

      // Assert the origin dart have different directions on the adjacent volume
      // So assert that their opposite are equal
      assert(lcc.beta(origin_dart, 1) == lcc.beta(vol2_origin_dart, 0, 3));

      // Face of the two neighboring volumes to the created volume
      Dart_handle face1 = lcc.beta(origin_dart, 2, 3);
      Dart_handle face2 = lcc.beta(origin_dart, 1, 2, 3);
      Dart_handle lower_edge = lcc.beta(upper_edge, 2, 1, 1);

      Dart_handle vol2_face1 = lcc.beta(vol2_origin_dart, 2, 3);
      Dart_handle vol2_face2 = lcc.beta(vol2_origin_dart, 0, 2, 3); // 0 because opposite directions
      Dart_handle vol2_lower_edge = lcc.beta(vol2_upper_edge, 2, 1, 1);

      // Contract upper and lower edge into its barycenter

      Dart_handle upper_mid_1 = lcc.insert_barycenter_in_cell<1>(upper_edge);
      Dart_handle upper_mid_2 = lcc.beta(upper_mid_1, 2, 1);

      thread_join_3_template_vertex__pair(hdata, upper_edge);

      // contract the shared edge between the two volumes
      lcc.contract_cell<1>(upper_mid_1);
      lcc.contract_cell<1>(upper_mid_2);

      clean_up_3_template(hdata, origin_dart, upper_edge, lower_edge, face1, face2);
      clean_up_3_template(hdata, vol2_origin_dart, vol2_upper_edge, vol2_lower_edge, vol2_face1, vol2_face2);
      nb_sub += 2;
    }

    return nb_sub;
  }

  /**
   * @brief Refines all marked faces using regular template substitution
   *
   * This function processes all faces that need refinement (from both the current plane
   * and the faces_to_refine collection) and applies regular template substitution to each
   * face. For each face, it:
   *   - Computes the face signature based on the marked nodes pattern
   *   - Applies the appropriate regular template substitution
   *   - Handles propagation marking for faces that need to propagate to adjacent volumes
   *   - Marks half-faces for propagation when necessary
   *
   * The function uses the regular_templates pattern substituter to replace each face
   * with a refined pattern based on its signature. It also handles the propagation
   * of refinement marks to adjacent faces and volumes.
   *
   * @tparam HexData Type of the hexahedral meshing data structure
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and external resources
   * @param rdata Refinement data containing the faces to refine and propagation information
   */
  template <typename HexData>
  void refine_marked_faces(HexData& hdata, RefinementData& rdata){
    LCC& lcc = hdata.lcc;
    int nbsub = 0;
    for (Dart_handle& dart : boost::join(rdata.faces_of_plane, rdata.faces_to_refine))
    {
      //Utile uniquement si les faces marqués ne sont pas 100% templatés
      auto edges_range = lcc.darts_of_cell<2, 1>(dart);
      int propagation = 0, beta3_propagation = 0;
      bool has_beta3 = !lcc.is_free<3>(dart);

      std::vector<Dart_handle> edges_vec;
      for (auto it = edges_range.begin(), end = edges_range.end(); it != end; it++){
        edges_vec.push_back(it);
        if (lcc.is_marked(it, hdata.propagation_face_mark)) propagation++;
        if (lcc.is_marked(lcc.beta<3>(it), hdata.propagation_face_mark)) beta3_propagation++;
      }

      auto& substituer = hdata.ext->regular_templates;
      Signature signature;
      Dart_handle f_signature_start = fsignature_of_face(lcc, dart, hdata.template_mark, signature);

      size_t temp_id = hdata.ext->regular_templates.replace_one_face_from_signature
        (lcc, dart, signature, f_signature_start);

      assert(temp_id < CONST_SIZE_T_MAX);
      nbsub++;

      if (temp_id == 0) {
        // lcc.mark_cell<1>(lcc.beta(f_signature_start, 1, 1), l_debug_mark_2);
        // lcc.mark_cell<1>(lcc.beta(f_signature_start, 1, 1, 2, 1), l_debug_mark);
        // debug_stream.push(l_thread_id);
        // assert(false);
        thread_number_vertex_in_1t_face(hdata, f_signature_start);
      }

      if (propagation >= 1){
        for (auto dart : edges_vec){
          mark_half_face_unchecked(lcc, dart, hdata.propagation_face_mark);
        }
      }

      if (beta3_propagation >= 1){
        for (auto dart : edges_vec){
          mark_half_face_unchecked(lcc, lcc.beta<3>(dart), hdata.propagation_face_mark);
        }
      }
    }

    // Cannot easily assert if all faces has been correctly treated, because some faces don't have attr
    // and we don't refine 3/4 template faces.

    std::cout << nbsub << " face substitution was made" << std::endl;
  }

  /**
   * @brief Refines all marked volumes using regular template substitution
   *
   * This function processes all volumes that need refinement (from the volumes_to_refine collection)
   * and applies regular template substitution to each volume. For each volume, it:
   *   - Computes the volume signature based on the marked nodes pattern
   *   - Applies the appropriate regular template substitution using the regular_templates substituter
   *   - Handles any special cases for 1-template volumes (e.g., vertex numbering)
   *
   * The function returns the total number of successful volume substitutions performed.
   *
   * @tparam HexData Type of the hexahedral meshing data structure
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and external resources
   * @param rdata Refinement data containing the volumes to refine
   * @return The total number of successful volume substitutions
   */
  template <typename HexData>
  int refine_regular_templates(HexData& hdata, RefinementData& rdata){

    LCC& lcc = hdata.lcc;
    int nbsub = 0;

    nbsub = 0;
    for (auto& dart : rdata.volumes_to_refine)
    {
      Pattern_substituer<LCC>& substituer = hdata.ext->regular_templates;
      Signature signature;
      Dart_handle v_signature_start = vsignature_of_volume(lcc, dart, hdata.template_mark, signature);
      size_type temp_id = substituer.replace_one_volume_from_signature(lcc, dart, signature, v_signature_start);
      if (temp_id < CONST_SIZE_T_MAX)
        nbsub++;

      if (temp_id == 0){
        thread_number_vertex_in_1t_vol(hdata, v_signature_start);
      }

      // if (temp > 0) nbsub++;
      // if (temp == 2) temp = 3;
      // if (temp < 10) temp++;
      // assert(vol_attr.template_id >= 0 && vol_attr.template_id <= 4 && vol_attr.template_id != 3 || vol_attr.template_id == 8);
      // assert(vol_attr.template_id == temp || temp == CONST_SIZE_T_MAX && (vol_attr.template_id == 0 || vol_attr.template_id == 8));
    }

    return nbsub;
  }

  /**
   * @brief Refines all partial 3-template patterns using partial template substitution
   *
   * This function processes all faces marked as partial templates (3-templates) in the
   * refinement data and applies partial template substitution to each. For each such face:
   *   - The origin dart of the 3-template is found and used for substitution
   *   - The missing edge is created to enable the query_replace operation
   *   - The partial template substituter is used to replace the volume(s) with the correct pattern
   *   - Both the original and the adjacent 3-template volumes are processed
   *
   * This operation is necessary to handle cases where a regular template cannot be applied
   * due to the presence of a 3-template configuration.
   *
   * @tparam HexData Type of the hexahedral meshing data structure
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and external resources
   * @param rdata Refinement data containing the list of partial templates to refine
   */
  template <typename HexData>
  void refine_partial_templates(HexData& hdata, RefinementData& rdata)
  {
    LCC& lcc = hdata.lcc;

    // Partially refine the 3 templates
    for (auto& marked_face : rdata.partial_templates_to_refine){
      assert(lcc.attribute<2>(marked_face) != nullptr);
      assert(lcc.attribute<2>(marked_face)->info().template_id == 3);

      // Query replace with the partial 3-template, making it into two volumes
      Dart_handle origin_dart = find_3_template_origin(lcc, marked_face, hdata.template_mark);
      marked_face = origin_dart;

      Dart_handle upper_d1 = origin_dart;
      Dart_handle upper_d2 = lcc.beta(origin_dart, 1, 1);

      // Create the missing edge to perform the query_replace
      Dart_handle upper_edge = lcc.insert_cell_1_in_cell_2(upper_d1, upper_d2);

      size_type p = hdata.ext->partial_templates.query_replace_one_volume(lcc, marked_face, hdata.template_mark);
      assert(p == 0);

      // Also replace the other connected volume that is 3 template
      p = hdata.ext->partial_templates.query_replace_one_volume(lcc, lcc.beta(marked_face, 3), hdata.template_mark);
      assert(p == 0);
    }
  }

  template <typename HexData>
  void thread_communicate_marked_nodes(HexData&, RefinementData&, size_type) {}

  template <typename HexData>
  void thread_communicate_cells_id_and_3t(HexData&, RefinementData&){}

  template <typename HexData>
  void thread_remove_ghosts(HexData& hdata) {}

  /**
   * @brief Explores a face in a plane to identify marked nodes and gather adjacent faces
   *
   * This function performs a breadth-first exploration of faces in a plane, starting from
   * a given face. It identifies marked nodes (0-cells), gathers adjacent faces (2-cells),
   * and collects incident volumes (3-cells) that need refinement. The function:
   *   - Marks the current face as explored to avoid revisiting
   *   - Initializes the face's template_id to 0
   *   - Iterates through all edges of the face
   *   - Marks identified nodes that haven't been marked yet
   *   - Increments the face's template_id for each identified node
   *   - Explores unvisited edges and finds adjacent faces on the same plane
   *   - Adds adjacent faces to the exploration queue
   *   - Collects additional volumes that are incident to the explored edges
   *   - Adds the face to the plane's face list if it has any identified nodes
   *
   * This function is part of the plane exploration process that identifies which faces
   * and volumes need to be refined based on the identified nodes.
   *
   * @tparam HexData Type of the hexahedral meshing data structure
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and marks
   * @param rdata Refinement data to store marked nodes, faces, and additional volumes
   * @param queue Queue of faces to explore (breadth-first traversal)
   * @param face The face to explore
   * @param explored_edge Mark to track explored edges
   * @param explored_face Mark to track explored faces
   */
  template <typename HexData>
  void explore_face_of_plane(HexData& hdata, RefinementData& rdata, std::queue<Dart_handle>& queue,
        Dart_handle face, size_type explored_edge, size_type explored_face) {
    LCC& lcc = hdata.lcc;

    auto& face_attr = lcc.attribute<2>(face)->info();

    if (!lcc.is_whole_cell_unmarked<2>(face, explored_face)) return ;
    lcc.mark_cell<2>(face, explored_face);

    face_attr.template_id = 0;

    // In multi threaded code, check if both volumes are owned, otherwise skip
    bool is_markable = true;

    // Might not be needed,

    // if constexpr (std::is_same_v<HexData, ProcessData>){
    //   auto& front_vol = lcc.attribute<3>(face)->info();
    //   auto back_vol_attr = lcc.attribute<3>(lcc.beta<3>(face));
    //   is_markable = front_vol.owned or (back_vol_attr != nullptr && back_vol_attr->info().owned);
    // }

    auto edges = lcc.darts_of_cell<2,1>(face);
    assert(edges.size() == 4);
    // Add neighboring faces
    for (auto dit = edges.begin(), dend = edges.end(); dit != dend; dit++){
      bool edge_explored = lcc.is_whole_cell_marked<1>(dit, explored_edge);
      bool marked = lcc.is_marked(dit, hdata.template_mark);
      bool identified = lcc.is_marked(dit, hdata.identified_mark);

      if (is_markable && !marked && identified){
        lcc.mark_cell<0>(dit, hdata.template_mark);
        rdata.marked_nodes.push_back(dit);
      }

      if (identified){
        face_attr.template_id++;
      }

      if (!edge_explored) {
        lcc.mark_cell<1>(dit, explored_edge);

        // Also add incident volumes while iterating
        Dart_handle other_face = __adjacent_face_on_plane(lcc, rdata.iteration, dit, rdata.additionnal_volumes_found);
        // Also do the other way around to catch volumes in the opposite dir
        if (!lcc.is_free<3>(dit)){
          Dart_handle other_face2 = __adjacent_face_on_plane(lcc, rdata.iteration, lcc.beta(dit, 3), rdata.additionnal_volumes_found);
          if (other_face == lcc.null_dart_descriptor) other_face  = other_face2;
        }

        if (other_face != lcc.null_dart_descriptor)
          queue.push(other_face);
      }
    }


    if (face_attr.template_id > 0) {
      rdata.faces_of_plane.push_back(face);
    }
  }

  /**
   * @brief Marks faces for propagation based on the template pattern of the given face
   * 
   * This function identifies which faces need to be marked for propagation based on the
   * template_id of the given face. It handles two cases:
   * 1. For template_id=1: Finds the single marked node and marks three adjacent faces
   * 2. For template_id=2: Finds the edge with both endpoints marked and marks two adjacent faces
   * 
   * The function locates the marked edge(s) on the face and then marks the appropriate
   * half-faces for propagation using the propagation_face_mark. This ensures that
   * refinement patterns can properly propagate to adjacent volumes.
   * 
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and marks
   * @param face A dart handle representing the face to analyze for propagation marking
   * @param face_attr The face attribute containing the template_id information
   */
  void mark_template_for_propagation(HexMeshingData &hdata, Dart_handle face, DartInfo::FaceAttrValue &face_attr)
  {
    LCC& lcc = hdata.lcc;
    Dart_handle marked_edge;
    auto nodes = lcc.darts_of_cell<2, 0>(face);

    bool found = false;

    if (face_attr.template_id == 1)
    {
      for (auto it = nodes.begin(); it != nodes.end() && !found; it++)
      {
        if (lcc.is_marked(it, hdata.template_mark))
        {
          found = true;
          marked_edge = it;
        }
      }
    }

    else if (face_attr.template_id == 2)
    {
      for (auto it = nodes.begin(); it != nodes.end() && !found; it++)
      {
        if (lcc.is_marked(it, hdata.template_mark) && lcc.is_marked(lcc.other_extremity(it), hdata.template_mark))
        {
          found = true;
          marked_edge = it;
        }
      }
    }

    assert(found);

    Dart_handle face1 = lcc.beta(marked_edge, 2, 1, 1, 2);
    Dart_handle face2 = lcc.beta(face1, 0, 2);
    Dart_handle face3 = lcc.beta(face1, 0, 0, 2);

    mark_half_face_unchecked(lcc, face1, hdata.propagation_face_mark);
    mark_half_face_unchecked(lcc, face3, hdata.propagation_face_mark);

    if (face_attr.template_id == 1)
    {
      mark_half_face_unchecked(lcc, face2, hdata.propagation_face_mark);
    }
  }

  /**
   * @brief Propagates refinement from a marked face to adjacent volumes and faces
   * 
   * This function handles the propagation of refinement patterns from a marked face
   * to its adjacent volumes and faces. It performs the following operations:
   * 1. Marks the current volume for refinement and sets its iteration
   * 2. Identifies the back face and its associated volume for refinement
   * 3. Adds the back volume to the refinement list if not already processed
   * 4. Marks the back face for refinement if not already explored
   * 5. Iterates through all edges of the back volume face
   * 6. Marks unmarked nodes and adds incident faces to the refinement list
   * 
   * This function ensures that refinement patterns properly propagate through the
   * mesh structure, maintaining consistency across adjacent volumes and faces.
   * 
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and marks
   * @param rdata Refinement data containing the current iteration and collections
   * @param face A dart handle representing the face to propagate from
   * @param explored_face_mark Mark used to track already explored faces
   */
  void propagate_face(HexMeshingData &hdata, RefinementData &rdata, const Dart_handle &face, size_type explored_face_mark)
  {
    LCC& lcc = hdata.lcc;


    assert(lcc.attribute<3>(face) != nullptr);
    auto &vol_attr = get_or_create_refinement_volume(lcc, face)->info();
    vol_attr.iteration = rdata.iteration;

    // /!\ ALSO: Mark add the hex attached to the back_face for refinement.
    // It is garanted that faces and volumes added will be unique in our array
    // Because that back_hex is only accessible within the template itself, and not
    // accessible from the plane.
    Dart_handle back_face = lcc.beta(face, 2, 1, 1, 2);
    Dart_handle back_volume_face = lcc.beta(back_face, 3);
    assert(back_volume_face != nullptr);
    auto &back_vol_attr = get_or_create_refinement_volume(lcc, back_volume_face)->info();

    if (back_vol_attr.iteration != rdata.iteration){
      rdata.volumes_to_refine.push_back(back_volume_face);
      back_vol_attr.iteration = rdata.iteration;
    }

    // Also add the neighboring faces to that 4-template face forrefinement
    if (!lcc.is_marked(back_volume_face, explored_face_mark)){
      rdata.faces_to_refine.push_back(back_volume_face);
      mark_face_unchecked(lcc, back_volume_face, explored_face_mark);
    }

    auto edges = lcc.darts_of_cell<2, 1>(back_volume_face);

    for (auto it = edges.begin(); it != edges.end(); it++)
    {
      auto top_face = lcc.beta(it, 2);

      if (!lcc.is_marked(it, hdata.template_mark)){
        lcc.mark_cell<0>(it, hdata.template_mark);
        rdata.marked_nodes.push_back(it);
      }

      if (!lcc.is_whole_cell_marked<2>(top_face, explored_face_mark)){
        rdata.faces_to_refine.push_back(top_face);
        mark_face_unchecked(lcc, top_face, explored_face_mark);
      }
    }
  }

  /**
   * @brief Executes the propagation stage of the refinement algorithm
   * 
   * This function handles the propagation of refinement patterns across the mesh.
   * It consists of two main phases:
   * 
   * Phase 1 (for iterations > 0): Propagates from 4-template faces
   * - Iterates through all faces in the current plane
   * - For faces with template_id=4, checks if they are marked for propagation
   * - Calls propagate_face for both the face and its opposite face (if it exists)
   * 
   * Phase 2 (for iterations < 2): Marks faces for future propagation
   * - Iterates through all faces in the current plane
   * - For faces with template_id=1 or 2, calls mark_template_for_propagation
   * - Handles both the face and its opposite face (if it exists)
   * 
   * The function skips propagation entirely for iteration 0 and stops marking
   * for propagation after iteration 2, as the refinement patterns are fully
   * established by that point.
   * 
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and marks
   * @param rdata Refinement data containing the current iteration and faces of the plane
   * @param explored_face_mark Mark used to track already explored faces
   */
  void propagation_stage(HexMeshingData &hdata, RefinementData &rdata, size_type explored_face_mark)
  {
    LCC& lcc = hdata.lcc;

    bool skip_propagation = rdata.iteration == 0;
    int propagated_count = 0;
    int marked_for_prop_count = 0;

    if (!skip_propagation) for (Dart_handle face : rdata.faces_of_plane){
      auto& face_attr = lcc.attribute<2>(face)->info();

      if (face_attr.template_id == 4){
        if (is_half_face_marked(lcc, face, hdata.propagation_face_mark)){
          propagate_face(hdata, rdata, face, explored_face_mark);
          propagated_count++;
        }

        Dart_handle other_face = lcc.beta(face, 3);
        if (!lcc.is_free<3>(face) && is_half_face_marked(lcc, other_face, hdata.propagation_face_mark)) {
          propagate_face(hdata, rdata, other_face, explored_face_mark);
          propagated_count++;
        }
      }
    }

    if (rdata.iteration >= 2) return;

    for (Dart_handle face : rdata.faces_of_plane){
      auto& face_attr = lcc.attribute<2>(face)->info();
      if (face_attr.template_id < 1 or face_attr.template_id > 2) continue;

      mark_template_for_propagation(hdata, face, face_attr);
      marked_for_prop_count++;

      if (!lcc.is_free<3>(face)){
        mark_template_for_propagation(hdata, lcc.beta(face, 3), face_attr);
        marked_for_prop_count++;
      }
    }

    std::cout << "Number of faces propagated : " << propagated_count << std::endl;
    std::cout << "Number of faces marked for propagation : " << marked_for_prop_count << std::endl;
  }

  /**
   * @brief Collects cells that need refinement from the current plane
   * 
   * This function processes all faces in the current refinement plane and identifies
   * which cells (faces and volumes) need to be refined. It performs the following operations:
   * 
   * 1. For each face in the plane, iterates through all its edges
   * 2. For edges with marked nodes, finds incident faces normal to the plane
   * 3. Adds unmarked incident faces to the faces_to_refine collection
   * 4. Processes the volumes associated with each face:
   *    - Creates or updates volume attributes with the current iteration
   *    - Adds volumes to volumes_to_refine or partial_templates_to_refine based on template_id
   *    - Handles both the current face's volume and its opposite volume (if it exists)
   * 
   * The function ensures that all cells affected by the refinement patterns are properly
   * identified and categorized for subsequent refinement operations.
   * 
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and marks
   * @param rdata Refinement data containing the faces of the plane and collections to populate
   * @param explored_faces Mark used to track already explored faces to avoid duplicates
   */
  void get_cells_to_refine_from_plane(HexMeshingData &hdata, RefinementData &rdata, size_type explored_faces)
  {
    LCC& lcc = hdata.lcc;

    // Iterate over all faces of even planes
    for (Dart_handle face : rdata.faces_of_plane)
    {
      auto edges = lcc.darts_of_cell<2, 0>(face);

      for (auto edge = edges.begin(); edge != edges.end(); edge++)
      {
        // get incident faces to marked node to be refined
        // Incident faces normal to the plane
        // We also need to prevent adding twice the faces by marking them

        if (!lcc.is_marked(edge, hdata.template_mark) && !lcc.is_marked(lcc.other_extremity(edge), hdata.template_mark))
          continue;

        auto top_face_1 = lcc.beta(edge, 2);
        auto top_face_2 = lcc.beta(edge, 3, 2);

        if (!lcc.is_whole_cell_marked<2>(top_face_1, explored_faces))
        {
          rdata.faces_to_refine.push_back(top_face_1);
          mark_face_unchecked(lcc, top_face_1, explored_faces);
        }

        if (top_face_2 != lcc.null_dart_descriptor && !lcc.is_whole_cell_marked<2>(top_face_2, explored_faces))
        {
          rdata.faces_to_refine.push_back(top_face_2);
          mark_face_unchecked(lcc, top_face_2, explored_faces);
        }
      }
      // Also add the adjacent volumes if there is atleast one marked node
      // Because we are on odd/even layers, we can't accidently add twice a volume

      auto &face_attr = lcc.attribute<2>(face)->info();
      auto &vol_attr = get_or_create_refinement_volume(lcc, face)->info();

      bool vol_processed = vol_attr.iteration == rdata.iteration;

      // Two faces point to the same volume. We are on a merged cell (caused by 3 templates)
      // if (vol_processed){
      //   mark_all_0_cells<3>(lcc, face, hdata.template_mark);
      // }

      // Create first volume attr and push back the volume
      if (!vol_processed) {
        vol_attr.iteration = rdata.iteration;

        if (face_attr.template_id != 3){
          rdata.volumes_to_refine.push_back(face);
        }
        else
          rdata.partial_templates_to_refine.push_back(face); // Both volumes are treated if it is 3 template
      }

      // Create second volume attr and push the face from the other volume
      if (!lcc.is_free<3>(face))
      {
        Dart_handle other_face = lcc.beta(face, 3);
        auto &other_vol_attr = get_or_create_refinement_volume(lcc, other_face)->info();
        bool other_vol_processed = other_vol_attr.iteration == rdata.iteration;

        // if (other_vol_processed)
          // mark_all_0_cells<3>(lcc, other_face, hdata.template_mark);
        if (!other_vol_processed) {
          other_vol_attr.iteration = rdata.iteration;

          if (face_attr.template_id != 3){
            rdata.volumes_to_refine.push_back(other_face);
          }
        }
      }
    }
  }

  /**
   * @brief Collects cells that need refinement from additional volumes discovered during plane exploration
   * 
   * This function processes additional volumes that were discovered during the plane exploration
   * phase but are not directly accessible from the current plane. These volumes typically
   * represent cells that are inside the refinement domain but not directly connected to
   * the current refinement plane.
   * 
   * For each additional volume, the function:
   * 1. Checks if either endpoint of the volume's edge is marked
   * 2. If marked nodes are found, adds the volume to the refinement list
   * 3. Identifies the four adjacent faces of the volume
   * 4. Adds unmarked adjacent faces to the faces_to_refine collection
   * 5. Handles face marking based on which specific nodes are marked
   * 
   * The function ensures that all volumes and faces affected by the refinement patterns
   * are properly identified, even if they are not directly accessible from the main plane.
   * 
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and marks
   * @param rdata Refinement data containing the additional volumes and collections to populate
   * @param explored_face Mark used to track already explored faces to avoid duplicates
   */
  void get_cells_to_refine_from_additionnal_volumes(HexMeshingData &hdata, RefinementData &rdata, size_type explored_face)
  {
    LCC& lcc = hdata.lcc;

    // No additionnal volumes should be found on the first iteration
    assert(rdata.iteration != 0 || rdata.iteration == 0 && rdata.additionnal_volumes_found.size() == 0);

    for (Dart_handle initial_edge : rdata.additionnal_volumes_found)
    {
      auto &vol_attr = lcc.attribute<3>(initial_edge)->info();

      bool node_1_marked = lcc.is_marked(initial_edge, hdata.template_mark);
      bool node_2_marked = lcc.is_marked(lcc.other_extremity(initial_edge), hdata.template_mark);

      Dart_handle adjacent_faces[] = {
        initial_edge,
        lcc.beta(initial_edge, 2),

        lcc.beta(initial_edge, 0, 2),
        lcc.beta(initial_edge, 1, 2)
      };

      if (node_1_marked || node_2_marked)
      {
        if (vol_attr.iteration != rdata.iteration)
          rdata.volumes_to_refine.push_back(initial_edge);

        vol_attr.iteration = rdata.iteration;

        if (!lcc.is_whole_cell_marked<2>(adjacent_faces[0], explored_face))
        {
          mark_face_unchecked(lcc, adjacent_faces[0], explored_face);
          rdata.faces_to_refine.push_back(adjacent_faces[0]);
        }

        if (!lcc.is_whole_cell_marked<2>(adjacent_faces[1], explored_face))
        {
          mark_face_unchecked(lcc, adjacent_faces[1], explored_face);
          rdata.faces_to_refine.push_back(adjacent_faces[1]);
        }
      }

      if (node_1_marked && !lcc.is_whole_cell_marked<2>(adjacent_faces[2], explored_face))
      {
        mark_face_unchecked(lcc, adjacent_faces[2], explored_face);
        rdata.faces_to_refine.push_back(adjacent_faces[2]);
      }

      if (node_2_marked && !lcc.is_whole_cell_marked<2>(adjacent_faces[3], explored_face))
      {
        mark_face_unchecked(lcc, adjacent_faces[3], explored_face);
        rdata.faces_to_refine.push_back(adjacent_faces[3]);
      }
    }
  }

  /**
   * @brief Main function to identify and collect all cells that need refinement for a given plane direction
   * 
   * This function orchestrates the complete process of identifying cells that need refinement
   * for a specific plane direction (X, Y, or Z). It performs the following sequence of operations:
   * 
   * 1. **Plane Exploration**: Explores all even planes in the specified direction using breadth-first search
   * 2. **Validation**: Validates that all faces in the plane have consistent template IDs
   * 3. **Conflict Resolution**: Fixes impossible template cases (diagonal 2-templates and adjacent 3-templates)
   * 4. **Node Communication**: Communicates marked nodes between threads (in parallel version)
   * 5. **Propagation**: Executes the propagation stage to mark faces for future propagation
   * 6. **Cell Collection**: Collects cells to refine from both the plane and additional volumes
   * 
   * The function manages marks for tracking explored edges and faces, and ensures proper cleanup
   * of these marks at the end. This is the central coordination function for the refinement process.
   * 
   * @tparam HexData Type of the hexahedral meshing data structure
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and plane information
   * @param rdata Refinement data to be populated with cells that need refinement
   * @param iterationPlane The plane direction (X, Y, or Z) for the current refinement iteration
   */
  template <typename HexData>
  void get_cells_to_refine(HexData &hdata, RefinementData &rdata, PlaneNormal iterationPlane)
  {
    LCC& lcc = hdata.lcc;
    rdata.iteration = iterationPlane;

    size_type explored_edge = lcc.get_new_mark();
    size_type explored_face = lcc.get_new_mark();

    PlaneSet& plane_set = hdata.first_face_of_planes[iterationPlane];

    // Explore all even planes
    for (int i = 1; i < plane_set.size(); i += 2) {
      std::queue<Dart_handle> to_explore;

      for (auto start : plane_set[i])
        to_explore.push(start);

      while (!to_explore.empty()) {
        Dart_handle front = to_explore.front();
        to_explore.pop();
        explore_face_of_plane(hdata, rdata, to_explore, front, explored_edge, explored_face);
      }
    }

    assert_faces_of_plane_valid(hdata, rdata);

    // Communicate found marked nodes,
    // fix procedure and propagation procedure will produce on two threads the same result

    fix_impossible_cases(hdata, rdata);
    assert_faces_of_plane_valid(hdata, rdata);

    thread_communicate_marked_nodes(hdata, rdata, explored_edge);
    propagation_stage(hdata, rdata, explored_face);

    get_cells_to_refine_from_plane(hdata, rdata,  explored_face);
    get_cells_to_refine_from_additionnal_volumes(hdata, rdata, explored_face);

    lcc.free_mark(explored_edge);
    lcc.free_mark(explored_face);

  }

  /**
   * @brief Cleans up mesh attributes and reevaluates cell identification status after refinement
   * 
   * This function performs cleanup operations after a refinement level is completed.
   * It handles the removal of attributes for cells outside the refinement domain and
   * reevaluates the identification status of cells within the domain. The function
   * performs the following operations:
   * 
   * 1. **Face Attribute Cleanup**: Removes face attributes for faces that are not
   *    associated with volumes inside the refinement domain or unowned ghost zones
   * 2. **Volume Attribute Cleanup**: Removes volume attributes for volumes outside
   *    the refinement domain that are owned by the current process
   * 3. **Volume Reset**: Resets all volume attributes within the refinement domain
   *    to their default state
   * 4. **Identification Reevaluation**: Re-evaluates the identification status of
   *    previously identified volumes using the provided cellIdentifier function
   * 5. **Validation**: Ensures that the cleanup operations correctly removed the
   *    expected number of attributes
   * 
   * Since the refinement process subdivides each original cell into 8 new cells
   * uniformly, this function only needs to reevaluate cells that were previously
   * identified to determine if they should remain identified in the new refined mesh.
   * 
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex
   * @param cellIdentifier Function that determines whether a cell should be identified
   *                       for refinement based on its position and properties
   */
  void clean_up_and_reevaluate_attributes(HexMeshingData& hdata, MarkingFunction& cellIdentifier){
    // Erase all volume / faces attributes outside the domain,
    // Reset all volumes inside the domain, and reevaluate their identification status

    // Since the refinement is uniform, we subdivided each original cell in 8 new cells.
    // This means we only need to reevaluate identified cells, if they are still identified or not.
    LCC& lcc = hdata.lcc;

    auto& face_attributes = lcc.attributes<2>();
    auto& vol_attributes = lcc.attributes<3>();

    std::vector<LCC::Attribute_descriptor<2>::type> faces_to_delete;
    std::vector<LCC::Attribute_descriptor<3>::type> volumes_to_delete;

    for (auto it = face_attributes.begin(), end = face_attributes.end(); it != end; it++){
      auto vol_handle = lcc.attribute<3>(it->dart());
      auto other_vol_handle = lcc.attribute<3>(lcc.beta(it->dart(), 3));

      if ((vol_handle == nullptr || vol_handle->info().type <= VolumeType::REFINEMENT)
        && (other_vol_handle == nullptr || other_vol_handle->info().type <= VolumeType::REFINEMENT))
        faces_to_delete.push_back(it);
    }

    for (auto it = vol_attributes.begin(), end = vol_attributes.end(); it != end; it++){
      DartInfo::VolumeAttrValue old_info = it->info();
      DartInfo::VolumeAttrValue& current_info = it->info();

      // Delete volumes outside of refinement domain, and outside of unowned ghosted zone
      if (old_info.type <= VolumeType::REFINEMENT && old_info.owned){
        volumes_to_delete.push_back(it);
        continue;
      }

      // Reset all volumes inside the refinement domain
      current_info = DartInfo::VolumeAttrValue();

      // Reevaluate the identification status of those who were identified
      if (old_info.type == VolumeType::IDENTIFIED && cellIdentifier(lcc, it->dart()))
        current_info.type = VolumeType::IDENTIFIED;
    }

    size_type size_face_before = lcc.attributes<2>().size();
    size_type size_vol_before = lcc.attributes<3>().size();

    for (auto attr_desc : faces_to_delete){
      assert(lcc.is_attribute_used<2>(attr_desc));
      assert(lcc.is_valid_attribute<2>(attr_desc));
      lcc.set_attribute<2>(attr_desc->dart(), nullptr);
    }

    for (auto attr_desc : volumes_to_delete){
      assert(lcc.is_attribute_used<3>(attr_desc));
      assert(lcc.is_valid_attribute<3>(attr_desc));
      lcc.set_attribute<3>(attr_desc->dart(), nullptr);
    }

    assert(size_face_before - faces_to_delete.size() == lcc.attributes<2>().size());
    assert(size_vol_before - volumes_to_delete.size() == lcc.attributes<3>().size());
  }

  /**
   * @brief Expands the set of identified cells by propagating identification to neighboring cells
   * 
   * This function expands the identification of cells that need refinement by propagating
   * the identification status to neighboring cells in the mesh. The expansion is performed
   * multiple times based on the current refinement level and total number of levels.
   * 
   * The function calculates the number of expansion iterations needed using a height-based
   * formula that depends on the current level and total number of refinement levels:
   * - For levels < nb_levels-2: 3 iterations
   * - For level nb_levels-2: 2 iterations  
   * - For level nb_levels-1: 0 iterations
   * 
   * During each expansion iteration, the function:
   * 1. Collects all volumes that have type >= VolumeType::ID_EXPANSION
   * 2. For each identified volume, finds all 26-connected neighboring volumes
   * 3. Sets the type of neighboring volumes to VolumeType::ID_EXPANSION if they
   *    currently have a lower type
   * 
   * This expansion ensures that the refinement region has sufficient padding around
   * the originally identified cells to maintain mesh quality and prevent artifacts
   * at the boundaries of the refinement domain.
   * 
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex
   * @param current_lvl The current refinement level (0-based)
   * @param nb_levels The total number of refinement levels to be performed
   */
  void expand_identified_cells(HexMeshingData& hdata, int current_lvl, int nb_levels){
    LCC& lcc = hdata.lcc;

    assert(nb_levels >= 1 && current_lvl >= 0);

    // TODO Reformulate, information not up to date

    // Calculate the totaling cells needed per level (height), with i ranging from 0 (lowest level) to n (highest level)
    // Very simple sequence : 4, 6, 6, 6, 6, ....               (height of n cells, cell sized from current i+1-subdivision)
    // Divide this by 2 : height of n cells, cell sized from current i-subdivision
    auto height_of_refinement_level = [](int i, int nb_levels) -> int {
      return i < nb_levels - 2 ? 3
          : i < nb_levels - 1 ? 2
          : 0
          ;
    };

    // Expand the set of identified cells n times
    for (int i = 0 ; i < height_of_refinement_level(current_lvl, nb_levels); i++){
      std::vector<Dart_handle> to_propagate;
      auto& vol_attributes = lcc.attributes<3>();

      for(auto it = vol_attributes.begin(), end = vol_attributes.end(); it != end; it++){
        if (it->info().type >= VolumeType::ID_EXPANSION)
          to_propagate.push_back(it->dart());
      }

      for (auto dart : to_propagate){
        for (auto vol : cells_26_connectivity(lcc, dart)){
          auto& vol_attr = get_or_create_attr<3>(lcc, vol)->info();

          if (vol_attr.type <= VolumeType::ID_EXPANSION)
            vol_attr.type = VolumeType::ID_EXPANSION;
        }
      }
    }

  }

  /**
   * @brief Sets up the initial plane structure for the hexahedral mesh refinement algorithm
   * 
   * This function initializes the plane structure that will be used throughout the refinement
   * process. It creates and organizes faces into planes along the X, Y, and Z axes, establishing
   * the foundation for the refinement algorithm's plane-based approach.
   * 
   * The function performs two main operations:
   * 
   * 1. **Plane Extraction**: For each axis (X, Y, Z), it extracts the first faces of each plane
   *    by traversing the mesh structure using beta operations. It uses helper functions:
   *    - `__first_plane`: Finds the starting face for each plane direction
   *    - `__next_plane`: Navigates to the next plane in the sequence
   *    The extraction iterates through all planes in each direction based on the grid dimensions.
   * 
   * 2. **Attribute Creation**: For all faces in all planes, it creates and initializes face
   *    attributes with the following properties:
   *    - Sets the appropriate plane bit in the plane bitset
   *    - Assigns a plane_id based on the plane's position in the sequence
   *    - Sets cc_id (connected component ID) to 0 for all faces
   * 
   * The function uses `plane_for_each_face` to traverse all faces in each plane and ensure
   * proper attribute initialization. This setup is crucial for the refinement algorithm to
   * properly identify and process faces during the refinement stages.
   * 
   * Note: The last odd plane in each direction is not considered in the current implementation,
   * as indicated by the commented code at the end of the function.
   * 
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and grid information
   */
  void setup_initial_planes(HexMeshingData& hdata){
    LCC& lcc = hdata.lcc;

    auto __first_plane = [](LCC& lcc, PlaneNormal plane){
      switch (plane){
        case X: return lcc.beta(lcc.first_dart(), 0, 2);
        case Y: return lcc.beta(lcc.first_dart(), 2, 1);
        case Z: return lcc.first_dart();
      }

      CGAL_assertion_msg(false, "Unexpected value for plane");
      return lcc.null_dart_descriptor;
    };

    auto __next_plane = [](LCC& lcc, Dart_handle start_plane, PlaneNormal plane){
      return lcc.beta(start_plane, 1, 2, 3, 2, 1);
    };

    // Extract the first faces of each plane
    // Note: The last odd plane is not considered
    for (int p = 0; p < 3; p++){
      PlaneNormal plane = (PlaneNormal)p;

      auto& plane_set = hdata.first_face_of_planes[p];

      Dart_handle start_plane = __first_plane(hdata.lcc, plane);

      for (int z = 0; z < hdata.grid.dims[p]; z++){
        std::vector<Dart_handle> starts;
        starts.push_back(lcc.beta(start_plane, 0, 2));

        start_plane = __next_plane(lcc, start_plane, plane);

        plane_set.push_back(std::move(starts));
      }
    }

    // Create attributes of all faces of all planes
    // To be able to iterate properly at later stages
    for (int p = 0; p < 3; p++){
      auto& plane_set = hdata.first_face_of_planes[p];
      for (int i = 0; i < plane_set.size(); i++) {
        plane_for_each_face(lcc, plane_set[i],
          [&](Dart_handle face, auto& edges){
            auto& face_attr = get_or_create_attr<2>(lcc, face)->info();
            face_attr.plane[p] = true;
            face_attr.plane_id = i;
            face_attr.cc_id = 0;
          },
          [&](Dart_handle edge){
            return lcc.beta(edge, 2, 3, 2);
          });
      }

       // Last odd plane
       // int size = hdata.first_face_of_planes[p].size();
       // Dart_handle start = hdata.first_face_of_planes[p][size-1][0];
       // start = lcc.beta(start, 2, 1, 1, 2);

       // plane_for_each_face(lcc, start,
       //  [&](Dart_handle face){
       //    auto& face_attr = get_or_create_attr<2>(lcc, face)->info();
       //    face_attr.plane[p] = true;
       //    face_attr.cc_id = 0;
       //  },
       //  [&](Dart_handle edge){
       //    return lcc.beta(edge, 2, 3, 2);
       //  });
    }
  }

  /**
   * @brief Performs the initial setup for the hexahedral mesh refinement algorithm
   * 
   * This function performs the essential initialization steps required before starting
   * the refinement process. It establishes the foundational structure that the refinement
   * algorithm will use throughout its execution.
   * 
   * The function performs two main operations:
   * 
   * 1. **Plane Setup**: Calls `setup_initial_planes` to initialize the plane structure
   *    along the X, Y, and Z axes. This creates the organizational framework for
   *    the refinement algorithm's plane-based approach.
   * 
   * 2. **Cell Identification**: Iterates through all volumes in the Linear Cell Complex
   *    and identifies which cells should be marked for refinement. For each volume:
   *    - Creates a volume attribute if it doesn't exist
   *    - Uses the provided `cellIdentifier` function to determine if the volume
   *      should be marked as identified for refinement
   *    - Sets the volume type to `VolumeType::IDENTIFIED` if the cell should be refined
   * 
   * This initialization is crucial because it establishes which cells will be the
   * starting points for the refinement process. The identified cells will be used
   * to determine the refinement patterns and guide the subsequent refinement operations.
   * 
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and grid information
   * @param cellIdentifier Function that determines whether a cell should be identified
   *                       for refinement based on its position and properties
   */
  void initial_setup(HexMeshingData& hdata, MarkingFunction& cellIdentifier){
    LCC& lcc = hdata.lcc;
    setup_initial_planes(hdata);

    // Mark initial identified cells
    auto volumes = lcc.one_dart_per_cell<3>();
    for (auto dart = volumes.begin(), end = volumes.end(); dart != end; dart++){
      // Create a 3-attr for all 3-cells in the LCC
      auto& vol_attr = get_or_create_attr<3>(lcc, dart)->info();

      // Mark those who are identified
      if (cellIdentifier(lcc, dart))
        vol_attr.type = VolumeType::IDENTIFIED;
    }
  }

  /**
   * @brief Sets up face attributes and union-find structures for the next refinement level
   * 
   * This function processes a single face during the setup phase for the next refinement level.
   * It handles the creation and management of union-find structures for both odd and even planes,
   * ensuring proper connectivity tracking for faces that belong to identified volumes.
   * 
   * The function performs the following operations:
   * 
   * 1. **Volume Identification Check**: Determines if the current face and its opposite face
   *    belong to volumes that are identified for refinement (type >= VolumeType::ID_EXPANSION)
   * 
   * 2. **Face Marking**: Marks the current face to avoid reprocessing and adds unvisited
   *    adjacent faces to the exploration queue
   * 
   * 3. **Edge Processing**: For each edge of the face:
   *    - Finds adjacent faces on the same plane using `__adjacent_face_on_plane`
   *    - Adds unvisited adjacent faces to the exploration queue
   *    - Performs union-find operations for odd planes if either volume is identified
   *    - Performs union-find operations for even planes if the current volume is identified
   * 
   * 4. **Union-Find Management**: 
   *    - For odd planes: Unifies connected components when both volumes are identified
   *    - For even planes: Creates and manages connected components for the back face
   *    - Handles the creation of new union-find sets when needed
   * 
   * 5. **Attribute Assignment**: Assigns the appropriate connected component IDs to face
   *    attributes based on the union-find results
   * 
   * This function is part of the plane reconstruction process that occurs after each
   * refinement level to maintain proper plane connectivity for the next refinement iteration.
   * 
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex
   * @param odd_face_to_handle Map from face attributes to odd plane union-find handles
   * @param even_face_to_handle Map from face attributes to even plane union-find handles
   * @param to_explore Queue of faces to be explored in breadth-first traversal
   * @param odd_union_find Union-find structure for odd plane connected components
   * @param even_union_find Union-find structure for even plane connected components
   * @param face The face to be processed
   * @param planeIteration The current plane direction (X, Y, or Z)
   * @param plane_id The ID of the current plane
   * @param edge_mark Mark used to track visited edges
   * @param face_mark Mark used to track visited faces
   */
  void setup_next_level_face(HexMeshingData& hdata,
                            std::unordered_map<LCC::Attribute_handle<2>::type, Union_find<Dart_handle>::handle>& odd_face_to_handle,
                            std::unordered_map<LCC::Attribute_handle<2>::type, Union_find<Dart_handle>::handle>& even_face_to_handle,
                            std::queue<Dart_handle>& to_explore,
                            Union_find<Dart_handle>& odd_union_find,
                            Union_find<Dart_handle>& even_union_find,
                            Dart_handle face,
                            PlaneNormal planeIteration,
                            size_t plane_id,
                            int edge_mark, int face_mark){
    LCC& lcc = hdata.lcc;
    auto vol_handle = lcc.attribute<3>(face);
    auto beta3_vol_handle = lcc.attribute<3>(lcc.beta<3>(face));

    bool identified = vol_handle != nullptr && vol_handle->info().type >= VolumeType::ID_EXPANSION;
    bool beta3_identified = beta3_vol_handle != nullptr && beta3_vol_handle->info().type >= VolumeType::ID_EXPANSION;

    if (lcc.is_whole_cell_marked<2>(face, face_mark)) return;
    lcc.mark_cell<2>(face, face_mark);

    auto edges = lcc.darts_of_cell<2, 1>(face);
    Union_find<Dart_handle>::handle odd_cc_id, even_cc_id;

    Dart_handle back_face = lcc.beta(face, 2, 1, 1, 2);

    for (auto it = edges.begin(), end = edges.end(); it != end; it++){
      // Iterate
      assert(lcc.belong_to_same_cell<3>(it, face));

      // Here we must ensure that the face is always oriented the same way by using __adjacent_face_on_plane
      // We don't care if we miss a grid-border 3-template, because if we didn't found it using this function,
      // it has a missing volume (in the same beta3 orientation of our iteration).
      // This means that we can discard this face because there is nothing left to check.
      Dart_handle adjacent_face = __adjacent_face_on_plane(lcc, planeIteration, it);

      if (adjacent_face == lcc.null_dart_descriptor) continue;

      if (lcc.is_whole_cell_unmarked<1>(it, edge_mark)){
        to_explore.push(adjacent_face);
        lcc.mark_cell<1>(it, edge_mark);
      }

      // Union find, only if a volume attached to the face is identified
      if (identified || beta3_identified){
        auto adj_face_attr = lcc.attribute<2>(adjacent_face);
        auto adj_cc_id = odd_face_to_handle.find(adj_face_attr);

        if (odd_cc_id != nullptr
          && adj_cc_id != odd_face_to_handle.end()
          && !odd_union_find.same_set(odd_cc_id, adj_cc_id->second)){
          odd_union_find.unify_sets(odd_cc_id, adj_cc_id->second);
        }

        if (odd_cc_id == nullptr
          && adj_cc_id != odd_face_to_handle.end())
          odd_cc_id = odd_union_find.find(adj_cc_id->second);
      }

      // Even plane union find, only if the studied volume is identified
      if (identified){
        auto back_face_attr = lcc.attribute<2>(lcc.beta(adjacent_face, 2, 1, 1, 2));
        // Skip back face with no attributes
        if (back_face_attr == nullptr) continue;
        auto back_cc_id = even_face_to_handle.find(back_face_attr);

        if (even_cc_id != nullptr
          && back_cc_id != even_face_to_handle.end()
          && !even_union_find.same_set(even_cc_id, back_cc_id->second)){
          even_union_find.unify_sets(even_cc_id, back_cc_id->second);
        }

        if (even_cc_id == nullptr
          && back_cc_id != even_face_to_handle.end())
          even_cc_id = even_union_find.find(back_cc_id->second);
      }
    }

    if (identified or beta3_identified){
      auto face_attr = lcc.attribute<2>(face);
      auto& face_info = face_attr->info();
      auto cc_id = odd_cc_id != nullptr ? odd_cc_id : odd_union_find.make_set(face);
      odd_face_to_handle[face_attr] = cc_id;
    }

    if (identified) {
      if (lcc.attribute<2>(back_face) != nullptr){
        lcc.mark_cell<2>(back_face, hdata.debug);
        lcc.mark_cell<2>(face, hdata.debug2);
        debug_stream.push(l_thread_id);
        std::this_thread::sleep_for(std::chrono::hours(1));
      }
      assert(lcc.attribute<2>(back_face) == nullptr);
      assert(!lcc.is_free<3>(back_face));

      auto back_face_attr = get_or_create_attr<2>(lcc, back_face);
      auto& back_face_info = back_face_attr->info();
      auto cc_id = even_cc_id != nullptr ? even_cc_id : even_union_find.make_set(back_face);
      even_face_to_handle[back_face_attr] = cc_id;
    }
  }

  /**
   * @brief Reconstructs the plane structure for the next refinement level
   * 
   * This function reconstructs the plane structure after a refinement level has been completed.
   * It processes the existing planes and creates new plane sets that reflect the refined mesh
   * structure. The function handles the transition from the old plane structure to the new one
   * that will be used in the next refinement iteration.
   * 
   * The function performs the following operations:
   * 
   * 1. **Plane Processing**: For each axis (X, Y, Z), processes all existing planes:
   *    - Iterates through each plane in the current plane set
   *    - Uses breadth-first traversal to explore all faces in each plane
   *    - Calls `setup_next_level_face` for each face to build union-find structures
   * 
   * 2. **Union-Find Analysis**: For each plane, creates separate union-find structures:
   *    - `odd_union_find` and `even_union_find` for tracking connected components
   *    - Maps face attributes to union-find handles for both odd and even planes
   *    - Uses `get_partitions` to extract connected component partitions
   * 
   * 3. **New Plane Creation**: Creates new plane sets based on the union-find results:
   *    - Each old plane is split into two new planes (odd and even)
   *    - Uses the `create_plane` lambda function to build new plane structures
   *    - Assigns new plane IDs and connected component IDs to face attributes
   * 
   * 4. **Attribute Management**: Updates face attributes with new plane information:
   *    - Sets the appropriate plane bit in the plane bitset
   *    - Assigns new plane_id values (pid * 2 for odd, pid * 2 + 1 for even)
   *    - Sets cc_id based on the connected component partition
   * 
   * 5. **Cleanup**: Manages marks and updates the global plane structure:
   *    - Frees edge and face marks after processing each axis
   *    - Replaces the old plane structure with the new one
   * 
   * This reconstruction is essential because after refinement, the mesh topology has changed,
   * and the plane structure needs to be updated to reflect the new connectivity patterns
   * for the next refinement iteration.
   * 
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and plane information
   */
  void setup_next_level_plane(HexMeshingData& hdata){
    LCC& lcc = hdata.lcc;
    std::array<PlaneSet, 3> new_planes;

    using FaceToHandle = std::unordered_map<LCC::Attribute_handle<2>::type, Union_find<Dart_handle>::handle>;
    using Union_find = Union_find<Dart_handle>;
    using UF_Partition = std::vector<Union_find::handle>;

    size_type edge_mark = lcc.get_new_mark();
    size_type face_mark = lcc.get_new_mark();

    const auto find_next_plane_id = [&](Dart_handle face, char axis) -> std::optional<size_t>{
      if (face == lcc.null_dart_descriptor) return {};
      face = lcc.beta(face, 2, 1, 1, 2);
      auto attr = lcc.attribute<2>(face);

      if (attr != nullptr){
        assert(attr->info().plane[axis]);
        return attr->info().plane_id;
      }

      face = lcc.beta(face, 3);
      if (face == lcc.null_dart_descriptor) return {};

      // Checks the second volume, if we are on a 1/2 template, we might have to check
      // all faces of the hex to find the next plane of 'axis'

      for (Dart_handle f : faces_of_hex(lcc, face)){
        attr = lcc.attribute<2>(f);
        if (attr != nullptr && attr->info().plane[axis]){
          return attr->info().plane_id;
        }
      }

      return {};
    };

    const auto create_plane = [&](char plane_normal, size_t plane_id, const UF_Partition& partition,
        const Union_find& uf, const FaceToHandle& face_to_handle) -> std::vector<Dart_handle>
    {
      std::vector<Dart_handle> plane_cc;
      std::unordered_map<Union_find::pointer, size_t> ptr_to_id;

      for (int i = 0; i < partition.size(); i++){
        ptr_to_id[partition[i].ptr()] = i;
      }

      for (auto uf_handle : partition){
        Dart_handle face_cc = *uf_handle;
        plane_cc.push_back(face_cc);
      }

      for (auto& [face, handle] : face_to_handle){
        auto id_it = ptr_to_id.find(uf.find(handle).ptr());
        if (id_it != ptr_to_id.end()){
          face->info() = DartInfo::FaceAttrValue();
          face->info().plane[plane_normal] = true;
          face->info().plane_id = plane_id;
          face->info().cc_id = id_it->second;
        }
      }

      return plane_cc;
    };

    for (int p = 0; p < 3; p++){
      PlaneSet& old_plane_set = hdata.first_face_of_planes[p];
      PlaneSet new_plane_set;

      for (int pid = 0; pid < old_plane_set.size(); pid++){
        std::queue<Dart_handle> to_explore;
        Union_find odd_union_find, even_union_find;
        FaceToHandle odd_face_to_handle, even_face_to_handle;


        for (auto start : old_plane_set[pid]){
          // Check the orientation of the face, facing up or down relative to its plane axis

          // bool is_facing_backwards = false;
          // bool found = false;

          // std::optional<size_t> id;
          // if ((id = find_next_plane_id(start, p))) {
          //   is_facing_backwards = *id < pid;
          //   found = true;
          // }
          // else if ((id = find_next_plane_id(lcc.beta(start, 3), p))){
          //   is_facing_backwards = *id > pid;
          //   found = true;
          // }

          // assert(found);
          // to_explore.push(is_facing_backwards ?  lcc.beta(start, 3) : start );
          to_explore.push( start );
        }

        while (!to_explore.empty()){
          Dart_handle face = to_explore.front();
          to_explore.pop();

          // TODO write this in a better way, there is so much args / lines
          setup_next_level_face(hdata, odd_face_to_handle, even_face_to_handle, to_explore, odd_union_find,
              even_union_find, face, (PlaneNormal)p, pid,
              edge_mark, face_mark);
        }

        PlaneCC odd_plane_cc, even_plane_cc;
        odd_plane_cc.reserve(odd_union_find.number_of_sets());
        even_plane_cc.reserve(even_union_find.number_of_sets());

        // Attribute new planes cc_id
        auto odd_partitions = get_partitions(odd_union_find);
        auto even_partitions = get_partitions(even_union_find);

        int plane_id = pid * 2;

        new_plane_set.push_back(create_plane(p, plane_id, odd_partitions, odd_union_find, odd_face_to_handle));
        new_plane_set.push_back(create_plane(p, plane_id+1, even_partitions, even_union_find, even_face_to_handle));
      }

      new_planes[p] = std::move(new_plane_set);

      lcc.unmark_all(edge_mark);
      lcc.unmark_all(face_mark);
    }

    hdata.first_face_of_planes = std::move(new_planes);

    lcc.free_mark(edge_mark);
    lcc.free_mark(face_mark);
  }

  /**
   * @brief Sets up the mesh structure for the next refinement level
   * 
   * This function performs the necessary setup operations to prepare the mesh
   * for the next refinement level. It handles the transition from the current
   * refinement level to the next one by updating the plane structure and
   * reevaluating cell identification status.
   * 
   * The function performs three main operations in sequence:
   * 
   * 1. **Plane Reconstruction**: Calls `setup_next_level_plane` to reconstruct
   *    the plane structure after the previous refinement level. This updates
   *    the plane organization to reflect the new mesh topology created by
   *    the refinement operations.
   * 
   * 2. **Attribute Cleanup and Reevaluation**: Calls `clean_up_and_reevaluate_attributes`
   *    to remove attributes for cells outside the refinement domain and
   *    reevaluate the identification status of cells within the domain.
   *    This ensures that only relevant cells are tracked for the next level.
   * 
   * 3. **Level Increment**: Increments the refinement level counter in the
   *    hexahedral meshing data structure to track the current refinement progress.
   * 
   * This function is called at the beginning of each refinement level (except
   * the first level, which uses `initial_setup`). It ensures that the mesh
   * structure is properly prepared for the refinement operations that will
   * follow in the current level.
   * 
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and plane information
   * @param cellIdentifier Function that determines whether a cell should be identified
   *                       for refinement based on its position and properties
   */
  void setup_next_level(HexMeshingData& hdata, MarkingFunction& cellIdentifier){
    setup_next_level_plane(hdata);
    clean_up_and_reevaluate_attributes(hdata, cellIdentifier);
    hdata.level++;
  }

  /**
   * @brief Marks all nodes (0-cells) that belong to identified volumes for refinement
   * 
   * This function iterates through all volume attributes in the Linear Cell Complex
   * and marks all nodes (vertices) that belong to volumes that have been identified
   * for refinement. The marking is done using the identified_mark to track which
   * nodes should be considered during the refinement process.
   * 
   * The function performs the following operations:
   * 
   * 1. **Volume Attribute Iteration**: Iterates through all volume attributes (3-cells)
   *    in the Linear Cell Complex using `lcc.attributes<3>()`
   * 
   * 2. **Identification Check**: For each volume, checks if:
   *    - The volume's type is greater than `VolumeType::NONE` (i.e., it has been
   *      identified for refinement, is in expansion, or is marked for refinement)
   *    - The volume is owned by the current process (`owned` flag is true)
   * 
   * 3. **Node Marking**: For volumes that meet the criteria, calls `mark_k_cells_of_i_cell<3, 0>`
   *    to mark all nodes (0-cells) that belong to the volume with the `identified_mark`
   * 
   * This function is called at the beginning of each plane iteration in the refinement
   * algorithm to ensure that all nodes belonging to identified volumes are properly
   * marked before the refinement process begins. This marking is essential for the
   * subsequent template identification and refinement operations.
   * 
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and marks
   */
  void mark_identified_cells_from_3_attrs(HexMeshingData& hdata) {
    LCC& lcc = hdata.lcc;

    auto& attributes = lcc.attributes<3>();

    for (auto it = attributes.begin(), end = attributes.end(); it != end; it++){
      if (it->info().type > VolumeType::NONE && it->info().owned){
        mark_k_cells_of_i_cell<3, 0>(lcc, it->dart(), hdata.identified_mark);
      }
    }
  }

  /**
   * @brief Generates a regular hexahedral grid in the Linear Cell Complex
   * 
   * This function creates a regular 3D grid of hexahedral cells based on the provided
   * grid configuration. It generates the initial mesh structure that will be used
   * as the starting point for the hexahedral refinement algorithm.
   * 
   * The function performs the following operations:
   * 
   * 1. **Grid Iteration**: Iterates through all grid positions in three dimensions
   *    (x, y, z) based on the grid dimensions specified in `grid.dims`
   * 
   * 2. **Coordinate Calculation**: For each grid position, calculates the coordinates
   *    of the eight vertices that define a hexahedral cell:
   *    - Bottom face vertices: (x1,y1,z1), (x2,y1,z1), (x2,y2,z1), (x1,y2,z1)
   *    - Top face vertices: (x1,y1,z2), (x2,y1,z2), (x2,y2,z2), (x1,y2,z2)
   *    Where x1, y1, z1 are the coordinates of the current grid position, and
   *    x2, y2, z2 are the coordinates of the next grid position
   * 
   * 3. **Hexahedron Creation**: Creates a hexahedral cell at each grid position using
   *    `lcc.make_hexahedron()` with the calculated vertex coordinates
   * 
   * 4. **Face Sewing**: After creating all hexahedra, calls `lcc.sew3_same_facets()`
   *    to properly connect adjacent faces and establish the correct topological
   *    relationships between neighboring cells
   * 
   * The resulting mesh is a regular grid where each cell is a hexahedron with
   * faces properly connected to its neighbors. This provides the foundation for
   * the subsequent refinement operations.
   * 
   * @param lcc The Linear Cell Complex where the grid will be created
   * @param grid Grid configuration containing position, size, and dimensions
   */
  void generate_grid(LCC& lcc, Grid& grid) {
    for (int x = 0; x < grid.dims.x; x++) {
      for (int y = 0; y < grid.dims.y; y++) {
        for (int z = 0; z < grid.dims.z; z++) {
          double x1 = grid.pos.x() + x * grid.size.x();
          double y1 = grid.pos.y() + y * grid.size.y();
          double z1 = grid.pos.z() + z * grid.size.z();

          double x2 = grid.pos.x() + (x+1)*grid.size.x();
          double y2 = grid.pos.y() + (y+1)*grid.size.y();
          double z2 = grid.pos.z() + (z+1)*grid.size.z();

          lcc.make_hexahedron(Point(x1,y1,z1), Point(x2,y1,z1),
                              Point(x2,y2,z1), Point(x1,y2,z1),
                              Point(x1,y2,z2), Point(x1,y1,z2),
                              Point(x2,y1,z2), Point(x2,y2,z2));
        }
      }
    }

    lcc.sew3_same_facets();
  }

  /**
   * @brief Creates new vertices for template refinement by subdividing edges
   * 
   * This function creates new vertices that are needed for the template refinement
   * process. It identifies edges that need to be subdivided based on the marked
   * nodes and creates new vertices at the barycenter of these edges.
   * 
   * The function follows a specific rule for vertex creation:
   * - Two adjacent marked nodes do not produce a new vertex
   * - One marked node adjacent to an unmarked node produces a new vertex
   * 
   * The function performs the following operations:
   * 
   * 1. **Edge Identification**: Iterates through all marked nodes and examines
   *    their incident edges using `lcc.one_dart_per_incident_cell<1, 0>()`
   * 
   * 2. **Edge Filtering**: For each edge incident to a marked node:
   *    - Skips edges that have already been processed (marked with `arrete_done`)
   *    - Skips edges where both endpoints are marked (no vertex creation needed)
   *    - Adds edges to the subdivision list where only one endpoint is marked
   * 
   * 3. **Vertex Creation**: For each edge that needs subdivision:
   *    - Inserts a barycenter in the edge using `lcc.insert_barycenter_in_cell<1>()`
   *    - Calls `thread_number_vertex_in_edge()` to handle any thread-specific
   *      vertex numbering operations
   * 
   * 4. **Cleanup**: Frees the temporary mark used for tracking processed edges
   * 
   * This function is called during the refinement stage to prepare the mesh
   * for template substitution by creating the necessary intermediate vertices.
   * 
   * @tparam HexData Type of the hexahedral meshing data structure
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and marks
   * @param rdata Refinement data containing the marked nodes for vertex creation
   */
  template <typename HexData>
  void create_vertices_for_templates(HexData& hdata, RefinementData& rdata)
  {
    // 2 noeuds marqué l'un à coté de l'autre ne produit pas de sommet
    // 1 noeud marqué a coté d'un noeud non marqué produit un sommet

    std::vector<Dart_handle> edges_to_subdivide;
    LCC& lcc = hdata.lcc;

    auto arrete_done = lcc.get_new_mark();

    int vertices_created = 0;
    for (auto dart : rdata.marked_nodes)
    {
      for (auto nit = lcc.one_dart_per_incident_cell<1, 0>(dart).begin(),
                nend = lcc.one_dart_per_incident_cell<1, 0>(dart).end();
          nit != nend;
          nit++)
      {
        if (lcc.is_marked(nit, arrete_done))
          continue;

        // If the node is next to an other marked node, we don't have to create vertices
        if (lcc.is_marked(lcc.beta<1>(nit), hdata.template_mark)){
          lcc.mark_cell<1>(nit, arrete_done);
          continue;
        }

        vertices_created++;
        edges_to_subdivide.push_back(nit);
        lcc.mark_cell<1>(nit, arrete_done);
      }
    }

    for (Dart_handle dart : edges_to_subdivide)
    {
      Dart_handle other_ext = lcc.other_extremity(dart);
      Dart_handle new_node = lcc.insert_barycenter_in_cell<1>(dart);
      thread_number_vertex_in_edge(hdata, new_node, dart, other_ext);
    }

    std::cout << "Vertices created: " << vertices_created << std::endl;

    lcc.free_mark(arrete_done);
  }

  /**
   * @brief Main algorithm for hexahedral mesh refinement
   * 
   * This function implements the complete two-refinement algorithm for hexahedral
   * mesh generation. It performs multiple levels of refinement, each consisting
   * of plane-based refinement operations along the X, Y, and Z axes.
   * 
   * The algorithm performs the following operations:
   * 
   * 1. **Initialization**: Sets up marks, generates the initial grid, and
   *    performs initial setup for the first refinement level
   * 
   * 2. **Refinement Levels**: For each refinement level:
   *    - Sets up the next level structure (except for the first level)
   *    - Expands identified cells to include neighboring cells
   *    - For each plane direction (X, Y, Z):
   *      - Marks identified cells from volume attributes
   *      - Identifies cells that need refinement
   *      - Creates vertices for templates
   *      - Refines marked faces and volumes
   *      - Validates mesh quality
   *      - Cleans up marks
   * 
   * 3. **Cleanup**: Removes ghost cells and performs final cleanup
   * 
   * @tparam HexData Type of the hexahedral meshing data structure
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and configuration
   * @param cellIdentifier Function that determines which cells should be refined
   * @param nb_levels Number of refinement levels to perform
   * @param thread_id Thread identifier for parallel processing (default: 0)
   */
  template <typename HexData>
  void two_refinement_algorithm(HexData& hdata, MarkingFunction& cellIdentifier, int nb_levels, int thread_id = 0){
    // static_assert(std::is_same_v<HexData, HexMeshingData> or std::is_same_v<HexData, ProcessData>);
    static_assert(std::is_base_of_v<HexMeshingData, HexData>);

    LCC& lcc = hdata.lcc;

    hdata.debug = lcc.get_new_mark();
    hdata.debug2 = lcc.get_new_mark();
    hdata.identified_mark = lcc.get_new_mark();
    hdata.template_mark = lcc.get_new_mark();
    hdata.propagation_face_mark = lcc.get_new_mark();
    l_debug_mark = hdata.debug;
    l_debug_mark_2 = hdata.debug2;
    l_thread_id = thread_id;

    // if constexpr (std::is_same_v<HexData, ProcessData>){
    if constexpr (!std::is_same_v<HexData, HexMeshingData>) {
      hdata.three_template_node_mark = lcc.get_new_mark();
      hdata.reset_temp_vertex_ids();
    }

    generate_grid(lcc, hdata.grid);

    // Levels of refinement
    for (int r = 0; r < nb_levels; r++){
      if (r == 0) initial_setup(hdata, cellIdentifier);
      else setup_next_level(hdata, cellIdentifier);

      // debug_stream.push(l_thread_id);
      // std::this_thread::sleep_for(std::chrono::hours(1));

      expand_identified_cells(hdata, r, nb_levels);

      // For each plane
      for (int p = 0; p < 3; p++) {
        RefinementData rdata;

        mark_identified_cells_from_3_attrs(hdata);

        get_cells_to_refine(hdata, rdata, (PlaneNormal)p);

        assert_dart_attr_are_unique<3>(lcc, rdata.volumes_to_refine, rdata.partial_templates_to_refine);

        // Refinement stage
        int total_sub = 0;
        create_vertices_for_templates(hdata, rdata);
        refine_marked_faces(hdata, rdata);
        total_sub += refine_regular_templates(hdata, rdata);
        refine_partial_templates(hdata, rdata);
        total_sub += refine_3_template(hdata, rdata);
        std::cout << total_sub << " volumic substitution was made" << std::endl;

        assert_all_faces_are_quadrilateral(lcc, hdata);
        assert_all_volumes_are_hexes(lcc);

        hdata.fix_dart_storage();

        thread_communicate_cells_id_and_3t(hdata, rdata);

        lcc.unmark_all(hdata.identified_mark);
        lcc.unmark_all(hdata.template_mark);
      }


      lcc.unmark_all(hdata.propagation_face_mark);
    }

    thread_remove_ghosts(hdata);
  }

  void mark_1template_face(LCC& lcc, size_type mark){
    auto dart = lcc.one_dart_per_cell<3>().begin();
    lcc.mark_cell<0>(lcc.beta(dart, 0, 2, 0, 0), mark);
  }

  void mark_2template_face(LCC& lcc, size_type mark){
    auto dart = lcc.one_dart_per_cell<3>().begin();
    lcc.mark_cell<0>(dart, mark);
    lcc.mark_cell<0>(lcc.beta(dart, 1), mark);
  }

  void mark_1template_volume(LCC& lcc, size_type mark){
    auto first = lcc.first_dart();
    lcc.mark_cell<0>(lcc.beta(first, 1, 2, 3, 1, 2, 0), mark);
  }

  void mark_2template_volume(LCC& lcc, size_type mark){
    auto first = lcc.first_dart();

    auto node1 = lcc.beta(first, 0, 2, 3, 1, 1, 2, 0);
    auto node2 = lcc.beta(node1, 0);

    lcc.mark_cell<0>(node1, mark);
    lcc.mark_cell<0>(node2, mark);
  }

  void mark_3template_partial_volume(LCC& lcc, size_type mark){
    auto first = lcc.first_dart();

    auto m2 = lcc.get_new_mark();

    auto node1 = lcc.beta(first, 0, 0, 2, 1, 1, 2);
    auto node2 = lcc.beta(node1, 1);
    auto node3 = lcc.beta(node2, 1);

    lcc.mark_cell<0>(node1, mark);
    lcc.mark_cell<0>(node2, mark);
    lcc.mark_cell<0>(node3, mark);
  }

  void mark_4template_volume(LCC& lcc, size_type mark){
    auto first = lcc.first_dart();

    auto node1 = first;
    auto node2 = lcc.beta(node1, 0);
    auto node3 = lcc.beta(node2, 0);
    auto node4 = lcc.beta(node3, 0);

    lcc.mark_cell<0>(node1, mark);
    lcc.mark_cell<0>(node2, mark);
    lcc.mark_cell<0>(node3, mark);
    lcc.mark_cell<0>(node4, mark);
  }

  /**
   * @brief Loads pattern substitution templates for hexahedral mesh refinement
   * 
   * This function loads the pattern substitution templates that are used during
   * the hexahedral mesh refinement process. It loads both regular templates for
   * complete refinement patterns and partial templates for 3-template cases.
   * 
   * The function performs the following operations:
   * 
   * 1. **Memory Reservation**: Reserves memory for pattern storage to avoid
   *    reallocation issues during loading (temporary fix for missing move operators)
   * 
   * 2. **Face Pattern Loading**: Loads face patterns for regular template substitution:
   *    - pattern1-face.moka: 1-template face pattern
   *    - pattern2-face.moka: 2-template face pattern
   * 
   * 3. **Volume Pattern Loading**: Loads volume patterns for regular template substitution:
   *    - pattern1.moka: 1-template volume pattern
   *    - pattern2.moka: 2-template volume pattern  
   *    - pattern4.moka: 4-template volume pattern
   * 
   * 4. **Partial Template Loading**: Loads partial volume patterns for 3-template cases:
   *    - pattern3.moka: 3-template partial volume pattern
   * 
   * Each pattern is loaded with its corresponding marking function that identifies
   * the nodes that should be marked for the specific template type.
   * 
   * @param regular_templates Pattern substituter for regular hexahedral templates
   * @param partial_3_template Pattern substituter for partial 3-template patterns
   */
  void load_patterns(Pattern_substituer<LCC> &regular_templates, Pattern_substituer<LCC>& partial_3_template) {
    // BUG Manque des opérateurs de déplacements pour Pattern, le reserve est un fix temporaire
    // Pour pouvoir charger les patterns correctement sans réallocation
    regular_templates.m_fpatterns.reserve(10);
    regular_templates.m_vpatterns.reserve(10);

    regular_templates.load_additional_fpattern(CGAL::data_file_path("hexmeshing/fpattern/pattern1-face.moka"), mark_1template_face);
    regular_templates.load_additional_fpattern(CGAL::data_file_path("hexmeshing/fpattern/pattern2-face.moka"), mark_2template_face);

    regular_templates.load_additional_vpattern(CGAL::data_file_path("hexmeshing/vpattern/complete/pattern1.moka"), mark_1template_volume);
    regular_templates.load_additional_vpattern(CGAL::data_file_path("hexmeshing/vpattern/complete/pattern2.moka"), mark_2template_volume);
    regular_templates.load_additional_vpattern(CGAL::data_file_path("hexmeshing/vpattern/complete/pattern4.moka"), mark_4template_volume);

    // TODO: Chargé séparément pour le moment, voir si je laisse comme ça
    partial_3_template.load_additional_vpattern(CGAL::data_file_path("hexmeshing/vpattern/partial/pattern3.moka"), mark_3template_partial_volume);
  }

  bool is_intersect(double x1, double y1, double z1,
                  double x2, double y2, double z2,
                  double x3, double y3, double z3,
                  double x4, double y4, double z4,
                  const Tree& t)
  {
    Kernel::Point_3 p1(x1,y1,z1);
    Kernel::Point_3 p2(x2,y2,z2);
    Kernel::Point_3 p3(x3,y3,z3);
    Kernel::Point_3 p4(x4,y4,z4);

    // And compute the two triangles
    Triangle t1(p1, p2, p3);
    if(t.do_intersect(t1))
    { return true; }

    t1=Triangle(p1, p3, p4);
    if(t.do_intersect(t1))
    { return true; }

    return false;
  }

  bool is_intersect(double x1, double y1, double z1,
                  double x2, double y2, double z2,
                  const Tree& t)
  {
    return
        is_intersect(x1,y1,z1, x2,y1,z1, x2,y1,z2, x1,y1,z2, t) || // f1 y1
        is_intersect(x2,y2,z1, x2,y1,z1, x2,y1,z2, x2,y2,z2, t) || // f2 x2
        is_intersect(x1,y2,z1, x2,y2,z1, x2,y2,z2, x1,y2,z2, t) || // f3 y2
        is_intersect(x1,y1,z1, x1,y1,z2, x1,y2,z2, x1,y2,z1, t) || // f4 x1
        is_intersect(x1,y1,z1, x2,y1,z1, x2,y2,z1, x1,y2,z1, t) || // f5 z1
        is_intersect(x1,y1,z2, x2,y1,z2, x2,y2,z2, x1,y2,z2, t);   // f6 z2
  }

  bool is_intersect(const LCC& lcc, LCC::Dart_const_handle dh, const Tree& t)
  {
    CGAL::Bbox_3 bbox=lcc.point(dh).bbox();
    // For each vertex of the volume
    for(auto it=lcc.one_dart_per_incident_cell<0,3>(dh).begin(),
        itend=lcc.one_dart_per_incident_cell<0,3>(dh).end(); it!=itend; ++it)
    { bbox+=lcc.point(it).bbox(); }

    return is_intersect(bbox.xmin(), bbox.ymin(), bbox.zmin(),
                        bbox.xmax(), bbox.ymax(), bbox.zmax(), t);
  }

  auto is_volume_intersecting_poly(Tree& tree) {
    return [&](LCC& lcc, Dart_handle dart){
      return is_intersect(lcc, dart, tree);
    };
  }

  /**
   * @brief Removes volumes from the Linear Cell Complex based on a trimming function
   * 
   * This function iterates through all volumes in the Linear Cell Complex and removes
   * those that do not satisfy the criteria specified by the trimming function.
   * 
   * The function performs the following operations:
   * 
   * 1. **Volume Iteration**: Iterates through all volumes (3-cells) in the Linear Cell Complex
   *    using `lcc.one_dart_per_cell<3>()`
   * 
   * 2. **Volume Evaluation**: For each volume, calls the provided trimming function
   *    `func(lcc, it)` to determine if the volume should be kept
   * 
   * 3. **Volume Removal**: If the trimming function returns false (indicating the volume
   *    should be removed), calls `lcc.remove_cell<3>(it)` to remove the volume from the mesh
   * 
   * This function is typically used to remove volumes that are outside the desired domain
   * or do not meet certain geometric or topological criteria. It is commonly called after
   * the refinement process to clean up the final mesh.
   * 
   * @param lcc The Linear Cell Complex containing the mesh
   * @param func Trimming function that determines whether a volume should be kept
   *             (returns true to keep, false to remove)
   */
  void trim_excedent_volumes(LCC& lcc, TrimmingFunction func){
    auto volumes = lcc.one_dart_per_cell<3>();
    for (auto it = volumes.begin(); it != volumes.end(); it++){
      if (func(lcc, it)) continue;
      lcc.remove_cell<3>(it);
    }
  }
}

namespace CGAL::HexRefinement {
  /**
   * @brief Renders the result of the two-refinement algorithm as a graphics scene
   * 
   * This function creates a visual representation of the hexahedral mesh generated
   * by the two-refinement algorithm. It sets up rendering options and displays
   * the mesh with optional trimming based on intersection with an AABB tree.
   * 
   * The function performs the following operations:
   * 
   * 1. **Scene Options Setup**: Creates `LCCSceneOptions` with custom rendering behavior:
   *    - Enables colored volume rendering for all volumes
   *    - Assigns random colors to volumes using `rand_color_from_dart`
   *    - Configures volume visibility based on intersection with AABB tree when trimming is enabled
   * 
   * 2. **Graphics Scene Creation**: Creates a `CGAL::Graphics_scene` buffer and adds
   *    the Linear Cell Complex to it using the configured scene options
   * 
   * 3. **Scene Display**: Calls `CGAL::draw_graphics_scene` to render the mesh
   * 
   * When trimming is enabled (trim = true), only volumes that intersect with the
   * provided AABB tree are displayed. This is useful for visualizing only the
   * relevant portions of the mesh that intersect with a specific geometric domain.
   * 
   * @param lcc The Linear Cell Complex containing the refined hexahedral mesh
   * @param aabb AABB tree used for intersection testing when trimming is enabled
   * @param trim Whether to trim the display based on AABB intersection (default: true)
   * @param title Title for the graphics window (default: "TwoRefinement Result")
   */
  void render_two_refinement_result(const LCC& lcc, Tree& aabb, bool trim = true, const char* title = "TwoRefinement Result"){
    LCCSceneOptions<LCC> gso;

    gso.colored_volume = [&](const LCC& lcc, LCC::Dart_const_handle dart){ return true; };
    gso.volume_color = [&](const LCC& lcc, LCC::Dart_const_handle dart){
      return rand_color_from_dart(lcc, dart);
    };
    gso.draw_volume = [&](const LCC& lcc, LCC::Dart_const_handle dart){
      return !trim or TwoRefinement::is_intersect(lcc, dart, aabb);
    };

    CGAL::Graphics_scene buffer;
    add_to_graphics_scene(lcc, buffer, gso);
    CGAL::draw_graphics_scene(buffer);
  }

  /**
   * @brief Main entry point for the two-refinement hexahedral mesh generation algorithm
   * 
   * This function provides a high-level interface for generating hexahedral meshes
   * using the two-refinement algorithm. It orchestrates the complete process from
   * initialization to final mesh generation, including pattern loading, algorithm
   * execution, and optional trimming.
   * 
   * The function performs the following operations:
   * 
   * 1. **Resource Initialization**: Creates `HexMeshingData` and `ExternalRessources`
   *    structures to hold the mesh data and pattern substitution resources
   * 
   * 2. **Pattern Loading**: Calls `load_patterns` to load all necessary template
   *    patterns for both regular and partial template substitution
   * 
   * 3. **Data Initialization**: Initializes the hexahedral meshing data with the
   *    provided grid configuration and external resources
   * 
   * 4. **Algorithm Execution**: Calls `two_refinement_algorithm` to perform the
   *    complete refinement process with the specified number of levels
   * 
   * 5. **Optional Trimming**: Currently commented out, but can be enabled to remove
   *    excess volumes using the provided trimming function
   * 
   * The function returns a Linear Cell Complex containing the refined hexahedral mesh.
   * The mesh quality and refinement patterns are determined by the provided grid
   * configuration and cell identification function.
   * 
   * @param grid Grid configuration defining the initial mesh structure and dimensions
   * @param cellIdentifier Function that determines which cells should be refined
   *                       based on their position and properties
   * @param trimmingFunction Function that determines which volumes should be kept
   *                         in the final mesh (currently not used)
   * @param nb_levels Number of refinement levels to perform (default: 1)
   * @return Linear Cell Complex containing the refined hexahedral mesh
   */
  LCC two_refinement(
      TwoRefinement::Grid grid,
      TwoRefinement::MarkingFunction cellIdentifier,
      TwoRefinement::TrimmingFunction trimmingFunction,
      int nb_levels = 1,
      bool trim = false)
  {
    using namespace TwoRefinement;

    HexMeshingData hdata;
    ExternalRessources res;

    load_patterns(res.regular_templates, res.partial_templates);
    hdata.init(&res, grid);

    two_refinement_algorithm(hdata, cellIdentifier, nb_levels);
    if(trim)
    // trim_excedent_volumes(hdata.lcc, trimmingFunction);

    // debug_render(hdata);

    return hdata.lcc;
  }
}
///////////////////////
/////     栞     //////
///////////////////////

