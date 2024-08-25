#pragma once

#include <CGAL/Graphics_scene.h>
#include <CGAL/IO/Color.h>
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
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <algorithm>
#include <bitset>
#include <boost/container_hash/hash_fwd.hpp>
#include <chrono>
#include <cstddef>
#include <exception>
#include <functional>
#include "utils.h"
#include <CGAL/Simple_Cartesian.h>

#include "../query_replace/cmap_query_replace.h"
#include "../query_replace/init_to_preserve_for_query_replace.h"

#include <boost/container/static_vector.hpp>
#include <boost/range/join.hpp>
#include <limits>
#include <memory>
#include <mutex>
#include <optional>
#include <qnamespace.h>
#include <qpointingdevice.h>
#include <thread>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>
#include <winnt.h>
#include <CGAL/Simple_cartesian.h>

template <typename T>
struct GenericPoint {
  T x, y, z;
  constexpr GenericPoint(): x(0), y(0), z(0) {}
  constexpr GenericPoint(T x, T y, T z): x(x), y(y), z(z) {}
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

  template <typename K>
  operator GenericPoint<K>() const { return {static_cast<K>(x), static_cast<K>(y), static_cast<K>(z)}; }
};

using PointInt = GenericPoint<int>;
using PointChar = GenericPoint<char>;

// x,y,z represents plane normals, with value between -1, 0 or 1, 0 = unset;
// If a volume is owned, multiple constraints/areas can overlap
using AreaId = PointChar;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT       FT;
typedef Kernel::Point_3  Point;
typedef Kernel::Vector_3 Vector;

typedef typename Kernel::Triangle_3 Triangle;

typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> AABB_Traits;
typedef CGAL::AABB_tree<AABB_Traits> Tree;

enum class VolumeType {
  NONE,         // Newly created vol_attribute
  REFINEMENT,   // Previously NONE volumes that was selected for refinement
  ID_EXPANSION, // Expansion of the Identified layer of volumes
  IDENTIFIED    // Volumes identified for refinement
};

struct AreaIdHasher
{
  std::size_t operator()(const AreaId& area) const noexcept
  {
    // Perfect hashing
    assert(sizeof(size_t) > 3); // should be always true
    int h = 8 * sizeof(char); // should be 8;
    uchar x = static_cast<uchar>(area.x);
    uchar y = static_cast<uchar>(area.y);
    uchar z = static_cast<uchar>(area.z);
    size_t r = z;
    return (((r << 8) + y) << 8) + x;
  }
};

template <typename T>
using AreaIDMap = std::unordered_map<AreaId, T, AreaIdHasher>;

// Parameters order are important
// Checks if two areas are equivalent
bool are_areas_equivalent(const AreaId& area, const AreaId& sub_area, bool owned){
  // An area contains the other if all non zero axes from the contained area are equal to those in area
  if (!owned) return area == sub_area;

  for (int i = 0; i < 3; i++){
    if (sub_area[i] != 0 && area[i] != sub_area[i]) return false;
  }

  return true;
}

template <typename T>
class ProdCons {
public:
  bool hasNext() const {
    return !queue.empty();
  }

  T waitNextItem(){
    std::unique_lock lock(m);

    while (!hasNext()){
      awake_signal.wait(lock);
    }

    auto item = std::move(queue.front());
    queue.pop();

    return item;
  }

  void push(T&& item){
    std::unique_lock lock(m);
    queue.push(item);
    awake_signal.notify_one();
  }

  void push(const T& item){
    std::unique_lock lock(m);
    queue.push(item);
    awake_signal.notify_one();
  }

private:
  std::queue<T> queue;
  std::mutex m;
  std::condition_variable awake_signal;
};

class LCCItems
{
public:
  template < class Storage >
  struct Dart_wrapper
  {
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

    struct VolumeAttrValue {
      char iteration = -1;
      VolumeType type = VolumeType::NONE;

      AreaId area_id = {0,0,0}; // Used for ghost cells;
      bool owned = true; // Disable this when creating others constraint area
    };

    struct FaceAttrValue {
      char template_id = 0;
      std::bitset<3> plane;
      typename CGAL::Union_find<typename Storage::Dart_handle>::handle cc_id;
      bool cc_wave = false;
    };

    typedef CGAL::Cell_attribute<Storage, FaceAttrValue> FaceAttr;
    typedef CGAL::Cell_attribute<Storage, VolumeAttrValue> VolumeAttr;
    typedef std::tuple<VertexAttr, void, FaceAttr, VolumeAttr> Attributes;
  };
};


typedef CGAL::Linear_cell_complex_traits<3,Kernel> LCCTraits;
typedef CGAL::Linear_cell_complex_for_combinatorial_map<3,3, LCCTraits, LCCItems> LCC;
typedef typename LCC::Dart_handle Dart_handle;
typedef typename LCC::Vertex_attribute_handle Vertex_handle;
typedef typename LCC::size_type  size_type;

using DartInfo = LCCItems::Dart_wrapper<LCC::Storage>;

thread_local size_type l_debug_mark;
thread_local size_type l_debug_mark_2;
thread_local int l_thread_id;
thread_local std::unordered_set<Dart_handle> l_handles;

// Do not use this, only to satisfy lcc.is_isomorphic_to
bool operator==(const DartInfo::VolumeAttrValue& first, const DartInfo::VolumeAttrValue& second){
  return true;
}

// Do no use this, only to statisfy lcc.is_isomorphic_to
bool operator==(const DartInfo::FaceAttrValue& first, const DartInfo::FaceAttrValue& second){
  return true;
}

const size_t SIZE_T_MAX = std::numeric_limits<size_t>::max();

namespace CGAL::HexRefinement::TwoRefinement {

  template<typename ... Ts>                                                 // (7)
  struct Overload : Ts ... {
      using Ts::operator() ...;
  };
  template<class... Ts> Overload(Ts...) -> Overload<Ts...>;

  enum PlaneNormal { X, Y, Z, NONE = -1};
  enum class NumberableCell {None, False, True};

  struct Grid {
    Point pos;
    Point size;
    PointInt dims;
    PointInt dims_id_start;

    // For each plane normals, the starting index of those planes
    // PointInt dims_id_start;

    Grid() {}
    Grid(Point from, Point cell_size, PointInt dims)
    : pos(from), size(cell_size), dims(dims) {}

    static Grid make_centered_grid(Point center, Point cell_size, PointInt dims) {
      Kernel::Vector_3 offset = {
        dims.x / 2 * cell_size.x(),
        dims.y / 2 * cell_size.y(),
        dims.z / 2 * cell_size.z(),
      };

      Point from = center - offset;
      return Grid(from, cell_size, dims);
    }

    static Grid make_grid(Point from, Point cell_size, PointInt dims){
      return Grid(from, cell_size, dims);
    }

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

    static Grid make_cube(Point from, double cell_size, int dim){
      auto& s = cell_size;
      return Grid(from, {s,s,s}, {dim,dim,dim});
    }
  };

  // Identifies which 3-cell should be refined
  using TrimmingFunction = std::function<bool(LCC&, Dart_handle)>;
  using MarkingFunction = std::function<bool(LCC&, Dart_handle)>;

  struct ExternalRessources {
    Pattern_substituer<LCC> regular_templates, partial_templates;
  };

  using PlaneCC = std::vector<Dart_handle>; // One dart per face connected components
  using PlaneSet = std::vector<PlaneCC>; // A set of planes

  struct HexMeshingData {
    // Must initialize
    Grid grid;
    ExternalRessources* ext;

    // Initialized by the algorithm
    LCC lcc;
    size_type identified_mark, template_mark, propagation_face_mark, debug, debug2;
    int level = 0;
    std::array<PlaneSet, 3> first_face_of_planes;

    HexMeshingData() {}
    void init(ExternalRessources* ext, Grid grid) {
      this->ext = ext;
      this->grid = grid;
    }
  };

  struct RefinementData {
    PlaneNormal iteration;

    std::vector<Dart_handle> marked_nodes;
    std::vector<Dart_handle> faces_of_plane;
    std::vector<Dart_handle> additionnal_volumes_found;

    std::vector<Dart_handle> volumes_to_refine;
    std::vector<Dart_handle> faces_to_refine;
    std::vector<Dart_handle> partial_templates_to_refine;
  };

  using AltVertexIdSingle = std::pair<size_t, size_t>;
  using AltVertexIdPair = std::pair<AltVertexIdSingle, AltVertexIdSingle>;
  using AltVertexIdPairPair = std::pair<AltVertexIdPair, AltVertexIdPair>; // Maximum level of pairing possible with a 3 template
  using AltVertexIdTripletPair = std::tuple<AltVertexIdPair, AltVertexIdPair, AltVertexIdPair>;

  // Total size = 64 Bytes + 1
  using AltVertexIdVariants = std::variant<
    std::monostate,
    AltVertexIdSingle,
    AltVertexIdPair,
    AltVertexIdPairPair,
    AltVertexIdTripletPair
  >;

  struct AltVertexIdHasher {
    std::size_t operator()(const size_t s) const noexcept {
      return std::hash<size_t>{}(s);
    }

    template <typename T, typename K>
    std::size_t operator()(const std::pair<T,K>& pair) const noexcept {
      size_t seed = 0;
      boost::hash_combine(seed, AltVertexIdHasher{}(pair.first));
      boost::hash_combine(seed, AltVertexIdHasher{}(pair.second));
      return seed;
    }

    template <typename T, typename K, typename U>
    std::size_t operator()(const std::tuple<T,K,U>& tuple) const noexcept {
      auto& [first, second, third] = tuple;
      size_t seed = 0;
      boost::hash_combine(seed, AltVertexIdHasher{}(first));
      boost::hash_combine(seed, AltVertexIdHasher{}(second));
      boost::hash_combine(seed, AltVertexIdHasher{}(third));
      return seed;
    }
  };

  // using AltIdsTuple = std::tuple<AltVertexIdSingle, AltVertexIdPair, AltVertexIdPairPair, AltVertexIdTripletPair>;

  // TODO find better names
  template <typename T>
  struct AltVertexIdMap {
    std::unordered_map<AltVertexIdSingle, T, AltVertexIdHasher> vertexId;
    std::unordered_map<AltVertexIdPair, T, AltVertexIdHasher> vertexIdPair;
    std::unordered_map<AltVertexIdTripletPair, T, AltVertexIdHasher> vertexIdTripletPair;
    std::unordered_map<AltVertexIdPairPair, T, AltVertexIdHasher> vertexIdPairPair;
  };

  struct TempVertexIdStorage {
    std::unordered_map<size_t, AltVertexIdSingle> vertexId;
    std::unordered_map<size_t, AltVertexIdPair> vertexIdPair;
    std::unordered_map<size_t, AltVertexIdTripletPair> vertexIdTripletPair;
    std::unordered_map<size_t, AltVertexIdPairPair> vertexIdPairPair;
  };

  // Attributes the id of cells after a refinement and their positions
  struct ThreadCellsIdMsg {
    AltVertexIdMap<size_t> ids;
    std::unordered_map<size_t, Point> positions;
  };
  using ThreadMarkedCellsMsg = std::unordered_set<size_t>;
  using ThreadMarkedCellsPtr = std::shared_ptr<ThreadMarkedCellsMsg>;
  using ThreadCellsIdPtr = std::shared_ptr<ThreadCellsIdMsg>;

  using ThreadMsg = std::variant<
    std::monostate,
    ThreadMarkedCellsPtr,
    ThreadCellsIdPtr
  >;

  template <typename T>
  struct IOThreadStream {
    // If both threads calls receive() and no data are inside the streams
    // => Deadlock, however, this won't happen because the thread always push data first, then receive
    using MsgStream = ProdCons<T>;

    IOThreadStream(){};
    IOThreadStream(MsgStream* in, MsgStream* out): in(in), out(out){}
    bool is_valid() const { return in && out; }
    void send(const ThreadMsg& msg){ out->push(msg);}
    void send(ThreadMsg&& msg){ out->push(msg);}
    ThreadMsg receive(){ return in->waitNextItem(); }

    operator bool() const { return is_valid(); }
    IOThreadStream& operator<<(const ThreadMsg& msg) { send(msg); return *this;}
    IOThreadStream& operator<<(ThreadMsg&& msg) { send(msg); return *this;}
    IOThreadStream& operator>>(ThreadMsg& msg) { msg = receive(); return *this; }

    MsgStream *in, *out;
  };

  // Debug .. to remove
  ProdCons<int> debug_stream, debug_stream2;
  using IOMsgStream = IOThreadStream<ThreadMsg>;

  struct ProcessData : public HexMeshingData {
    struct Neighbor {
      Neighbor(){};
      Neighbor(int thread_id, IOMsgStream stream) : thread_id(thread_id), stream(stream){};
      int thread_id;
      IOMsgStream stream;

    };
    AreaIDMap<Neighbor> neighboring_threads;

    TempVertexIdStorage vertex_temp_ids;
    AltVertexIdMap<size_t> vertices_to_share;

    // TODO : Handles may be deleted and point to an other dart
    // Find a way to keep them on any valid volume / face (for HexMeshingData::first_face_of_planes)
    Dart_handle owned_ghosts_begin, unowned_ghosts_begin;

    std::array<bool, 3> positive_axes = {false, false, false};
    std::array<bool, 3> negative_axes = {false, false, false};;

    PointInt pos;
    PointInt full_grid;
    PointInt max_position; // TODO, we don't need to store that much, refactor the
    PointInt cell_dimension;

    size_t v_domain_current, v_domain_end, v_temp_id_begin;
    int thread_id;
    size_type three_template_node_mark;

    LCC::Attribute_handle<1>::type numberable_edge_attr, unnumberable_edge_attr;

    size_t get_new_vertex_id(){
      assert(v_domain_current + 1 < v_domain_end);
      return v_domain_current++;
    }

    void reset_temp_vertex_ids(){
      v_temp_id_begin = v_domain_end - 1;
    }

    size_t get_temp_vertex_id(){
      assert(v_temp_id_begin - 1 > v_domain_current);
      return v_temp_id_begin--;
    }

    bool is_temp_id(size_t id){
      return id >= v_temp_id_begin && id < v_domain_end;
    }
  };


  template <unsigned int i, unsigned int k>
  void mark_k_cells_of_i_cell(LCC& lcc, Dart_handle dart, size_type mark){
    auto iterator = lcc.darts_of_cell<i, 0>(dart);
    for (auto dit = iterator.begin(), dend = iterator.end(); dit != dend; dit++){
      if (!lcc.is_marked(dit, mark)) {
        lcc.mark_cell<k>(dit, mark);
      }
    }
  }

  template <typename FaceOp, typename EdgeOp>
  // Iterates on all faces and its edges
  // Doesn't iterate twice on the same edge
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
      face_operation(face);

      auto edges = lcc.darts_of_cell<2,1>(face);
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

  std::array<Dart_handle, 8> nodes_of_hex(LCC& lcc, Dart_handle vol){
    std::array<Dart_handle, 8> arr;
    int x = 0;

    Dart_handle it = vol;
    for (int i = 0; i < 4; it = lcc.beta(it, 1), i++){
      arr[x++] = it;
    }
    assert(it == vol);

    Dart_handle start = lcc.beta(vol, 2, 1, 1, 2);
    it = start;
    for (int i = 0; i < 4; it = lcc.beta(it, 1), i++){
      arr[x++] = it;
    }

    assert(it == start);
    return arr;
  }
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

  template <typename VolumeOp, typename VolumeSelector>
  inline void for_each_hex(LCC& lcc, Dart_handle start,
                            const VolumeOp&& volume_operation,
                            const VolumeSelector&& volume_selector)
  {
    std::queue<Dart_handle> to_explore;
    to_explore.push(start);
    size_type volume_mark = lcc.get_new_mark();
    size_type face_mark = lcc.get_new_mark();
    lcc.mark_cell<3>(start, volume_mark);

    while (!to_explore.empty()){
      Dart_handle volume = to_explore.front();
      to_explore.pop();

      auto hex_faces = faces_of_hex(lcc, volume);
      volume_operation(volume, hex_faces);

      for (auto face : hex_faces){
        Dart_handle other_vol = lcc.beta(face,3);
        if (other_vol != lcc.null_dart_descriptor && !lcc.is_marked(other_vol, volume_mark)){
          lcc.mark_cell<3>(other_vol, volume_mark);
          if (volume_selector(other_vol)) to_explore.push(other_vol);
        }
      }
    }

    lcc.free_mark(volume_mark);
  }


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

  template <typename FaceOp, typename EdgeOp>
  void plane_for_each_face(LCC& lcc, Dart_handle start,
                            const FaceOp&& face_operation,
                            const EdgeOp&& edge_operation)
  {
    std::queue<Dart_handle> to_explore;
    to_explore.push(start);
    __plane_for_each_face(lcc, to_explore,
      std::forward<const FaceOp>(face_operation),
      std::forward<const EdgeOp>(edge_operation));
  }

  PointInt grid_cell_from_index_2x2xN(int i){
    // 3 === 2
    // /     /
    // 0 === 1
    return {  i%4 == 3 or i%4 == 0 ? 0 : 1,
              i%4 < 2 ? 0 : 1,
              i/4};
  }

  PointInt grid_cell_from_index_2x1xN(int i){
    // 0 === 1
    return { i%2, 0, i/2};
  }

  void link_threads(std::vector<ProcessData>& processes, std::vector<ProdCons<ThreadMsg>>& streams, int a, int b, int& nextStream){
    auto& p1 = processes[a];
    auto& p2 = processes[b];

    PointChar rel_pos = p2.pos - p1.pos;

    p1.neighboring_threads[rel_pos] = ProcessData::Neighbor(
      b, IOMsgStream(&streams[nextStream], &streams[nextStream+1])
    );

    p2.neighboring_threads[-rel_pos] = ProcessData::Neighbor(
      a, IOMsgStream(&streams[nextStream+1], &streams[nextStream])
    );

    assert(nextStream < streams.size());
    nextStream += 2;
  }

  // Creates threads data, in multiples of 4
  void create_threads_2x2xN(std::vector<ProcessData>& procs, std::vector<ProdCons<ThreadMsg>>& streams, int nb_threads){
    assert(nb_threads % 4 == 0);

    streams = std::vector<ProdCons<ThreadMsg>>(nb_threads * 6 * 2 + (nb_threads - 1) * 14 * 2);

    int nextStream = 0;

    for (int i = 0; i < nb_threads; i++){
      procs[i].thread_id = i;
      procs[i].pos = grid_cell_from_index_2x2xN(i);
    }

    for (int i = 0; i < nb_threads / 4; i++){
      int x = i * 4;
      int prev_x = (i-1) * 4;

      // Full mesh the 4 cells
      link_threads(procs, streams, x + 0, x + 1, nextStream);
      link_threads(procs, streams, x + 1, x + 2, nextStream);
      link_threads(procs, streams, x + 2, x + 3, nextStream);
      link_threads(procs, streams, x + 3, x + 0, nextStream);

      link_threads(procs, streams, x + 0, x + 2, nextStream);
      link_threads(procs, streams, x + 1, x + 3, nextStream);

      if (x == 0) continue;

      // // Full mesh with previous layer
      // Link one to one with prev layer
      link_threads(procs, streams, x + 0, prev_x + 0, nextStream);
      link_threads(procs, streams, x + 1, prev_x + 1, nextStream);
      link_threads(procs, streams, x + 2, prev_x + 2, nextStream);
      link_threads(procs, streams, x + 3, prev_x + 3, nextStream);

      // Inner diagonals
      link_threads(procs, streams, x + 0, prev_x + 2, nextStream);
      link_threads(procs, streams, x + 2, prev_x + 0, nextStream);
      link_threads(procs, streams, x + 1, prev_x + 3, nextStream);
      link_threads(procs, streams, x + 3, prev_x + 1, nextStream);

      // Bottom layer diagonals
      link_threads(procs, streams, x + 1, prev_x + 0, nextStream);
      link_threads(procs, streams, x + 0, prev_x + 1, nextStream);

      // Right side diagonals
      link_threads(procs, streams, x + 0, prev_x + 3, nextStream);
      link_threads(procs, streams, x + 3, prev_x + 0, nextStream);

      // Top layer diagonals
      link_threads(procs, streams, x + 2, prev_x + 3, nextStream);
      link_threads(procs, streams, x + 3, prev_x + 2, nextStream);

      // Left side diagonals
      link_threads(procs, streams, x + 1, prev_x + 2, nextStream);
      link_threads(procs, streams, x + 2, prev_x + 1, nextStream);
    }
  }

  // TODO test this function
  void create_threads_2x1xN(std::vector<ProcessData>& procs, std::vector<ProdCons<ThreadMsg>>& streams, int nb_threads){
    assert(nb_threads % 2 == 0);

    streams = std::vector<ProdCons<ThreadMsg>>(nb_threads * 2 + (nb_threads - 1) * 4 * 2);

    int nextStream = 0;

    for (int i = 0; i< nb_threads; i++){
      procs[i].thread_id = i;
      procs[i].pos = grid_cell_from_index_2x1xN(i);
    }

    for (int i = 0; i < nb_threads / 2; i++){
      int x = i * 2;
      int prev_x = (i-1)*2;

      // Link the 2 procs
      link_threads(procs, streams, x + 0, x + 1, nextStream);

      if (x == 0) continue;

      // Link with previous procs
      link_threads(procs, streams, x + 0, prev_x + 0, nextStream);
      link_threads(procs, streams, x + 1, - + 1, nextStream);
      link_threads(procs, streams, x + 0, prev_x + 1, nextStream);
      link_threads(procs, streams, x + 1, prev_x + 0, nextStream);
    }
  }

  void create_threads_1x1xN(std::vector<ProcessData>& procs, std::vector<ProdCons<ThreadMsg>>& streams, int nb_threads){
    for (int i = 0; i < nb_threads; i++){
      procs[i].thread_id = i;
      procs[i].pos = {0, 0, i};
    }

    streams = std::vector<ProdCons<ThreadMsg>>(nb_threads * 2);

    int nextStream = 0;
    for (int i = 1; i < nb_threads; i++){
      link_threads(procs, streams, i-1, i, nextStream);
    }
  }

  template <unsigned int i>
  typename LCC::Attribute_descriptor<i>::type get_or_create_attr(LCC& lcc, Dart_handle dart){
    auto attr = lcc.attribute<i>(dart);

    if (attr == nullptr){
      attr = lcc.create_attribute<i>();
      lcc.set_attribute<i>(dart, attr);
    }
    return attr;
  }

  LCC::Attribute_descriptor<3>::type get_or_create_refinement_volume(LCC& lcc, Dart_handle dart){
    auto attr = get_or_create_attr<3>(lcc, dart);

    // Previously NONE volumes are tagged as refined
    if (attr->info().type == VolumeType::NONE)
      attr->info().type = VolumeType::REFINEMENT;

    return attr;

  }

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

  void mark_face_unchecked(LCC& lcc, Dart_handle dart, size_type mark){
    auto iterator = lcc.darts_of_cell<2, 1>(dart);
    for (auto dit = iterator.begin(), dend = iterator.end(); dit != dend; dit++) {
      lcc.mark(dit, mark);

      if (!lcc.is_free<3>(dit))
        lcc.mark(lcc.beta<3>(dit), mark);
    }

    assert(lcc.is_whole_cell_marked<2>(dart, mark));
  }

  bool is_half_face_marked(const LCC& lcc, LCC::Dart_const_handle dart, size_type mark ){
    auto iterator = lcc.darts_of_cell<2, 1>(dart);
    for (auto dit = iterator.begin(), dend = iterator.end(); dit != dend; dit++) {
      if (!lcc.is_marked(dit, mark)) return false;
    }

    return true;
  }

  void mark_half_face_unchecked(LCC& lcc, Dart_handle dart, size_type mark){
    auto iterator = lcc.darts_of_cell<2, 1>(dart);

    for (auto dit = iterator.begin(), dend = iterator.end(); dit != dend; dit++) {
      lcc.mark(dit, mark);
    }

    assert(is_half_face_marked(lcc, dart, mark));
  }

  // Gets the adjacent face to a edge on a given plane.
  // This method can fail if the adjacent face is not directly accessible from the current beta3 of the face.
  // Otherwise, you should compute this function with both edge and lcc.beta<3>(edge)
  // TODO : Give a better name
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

  // This overload adds encountered volume that couldn't be reached directly from the plane
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

  // Gets the adjacent face to a edge on a given plane.
  // Beware that the orientation of the face (beta3) can change.
  // This can be verified if the result and the edge are oriented the same way (TODO Verify)
  Dart_handle adjacent_face_on_plane(LCC& lcc, PlaneNormal plane, Dart_handle edge){
    Dart_handle other_face = __adjacent_face_on_plane(lcc, plane, edge);

  // This is needed to iterate on 3 templates that are on a grid border, missing a connected volume
    // And preventing iteration the first way
    if (other_face == lcc.null_dart_descriptor && !lcc.is_free<3>(edge))
      other_face = __adjacent_face_on_plane(lcc, plane, lcc.beta<3>(edge));

    return other_face;
  }

  Dart_handle adjacent_face_on_plane(LCC& lcc, PlaneNormal plane, Dart_handle edge, std::vector<Dart_handle>* additionnal_volumes){
    Dart_handle other_face = __adjacent_face_on_plane(lcc, plane, edge, *additionnal_volumes);

    if (other_face == lcc.null_dart_descriptor && !lcc.is_free<3>(edge))
      other_face = __adjacent_face_on_plane(lcc, plane, lcc.beta<3>(edge), *additionnal_volumes);

    return other_face;
  }


  // There may be up to 8 faces adjacent to a vertex.
  // Each dart_handle belongs to a unique edge, and an unique face on the plane and around the selected node
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
          assert(lcc.belong_to_same_cell<0>(adjacent_face, node) or (!lcc.is_free<3>(adjacent_face) && lcc.belong_to_same_cell<0>(lcc.beta(adjacent_face, 3), node)));
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

  template <unsigned int i, typename... DartArray>
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

  void thread_assert_no_temporary_ids(ProcessData& proc){
    LCC& lcc = proc.lcc;
    auto assertion = [&](){
      auto iterator = lcc.one_dart_per_cell<0>();
      for (auto it = iterator.begin(); it != iterator.end(); it++){
        bool is_temp_id = proc.is_temp_id(lcc.attribute<0>(it)->id);
        CGAL_assertion_msg(!is_temp_id, "thread_assert_no_temporary_ids: A vertex has not been assigned a true ID after communication");
      }
    };

    CGAL_assertion_code(assertion());
  }

  template <typename T, typename K>
  T get_temp_vertex_id(K& left, K& right){
    auto [id_min, id_max] = std::minmax(left, right);
    return {id_min, id_max};
  }

  template <typename T, typename K>
  T get_temp_vertex_id(K& a, K& b, K& c){
    std::array<K,3> arr = {a,b,c};
    std::sort(arr.begin(), arr.end());
    return {arr[0], arr[1], arr[2]};
  }

  NumberableCell is_edge_numberable(ProcessData& proc, LCC::Dart_const_handle edge) {
    LCC& lcc = proc.lcc;
    auto face_orbit = lcc.darts_of_orbit<2,3>(edge);

    bool sharable_owned = false;
    bool sharable_unowned = false;
    bool lowest_id = true;

    for (auto it = face_orbit.begin(); it != face_orbit.end(); it++){
      auto vol = lcc.attribute<3>(it)->info();

      if (!vol.owned)
        sharable_unowned = true;

      if (vol.owned && vol.area_id != AreaId{0,0,0})
        sharable_owned = true;

      if (vol.owned)
        continue;

      assert(proc.neighboring_threads.count(vol.area_id));
      auto other_id = proc.neighboring_threads[vol.area_id].thread_id;

      if (other_id < proc.thread_id)
        lowest_id = false;
    }

    if (!sharable_owned && sharable_unowned)
      return NumberableCell::False;

    if (sharable_owned && !sharable_unowned)
      return NumberableCell::True;

    if (sharable_owned && sharable_unowned)
      return lowest_id ? NumberableCell::True : NumberableCell::False;

    return NumberableCell::None;
  };

  NumberableCell is_face_numberable(ProcessData& proc, LCC::Dart_const_handle face){
    LCC& lcc = proc.lcc;
    auto& front_vol = lcc.attribute<3>(face)->info();
    auto back_vol = lcc.attribute<3>(lcc.beta(face,3));

    constexpr auto NoneArea = AreaId{0,0,0};

    // Don't number faces inside non shared volumes
    if (front_vol.area_id == NoneArea && (back_vol == nullptr or back_vol->info().area_id == NoneArea))
      return NumberableCell::None;

    // === For faces inside sharable volumes

    // Number them if they are both owned volumes
    if (front_vol.owned && (back_vol == nullptr or back_vol->info().owned))
      return NumberableCell::True;

    // Let other threads number them if both volumes are unowned
    if (!front_vol.owned && (back_vol == nullptr or !back_vol->info().owned))
      return NumberableCell::False;

    // If this is a boundary, check which processor has the lowest id
    assert(back_vol != nullptr && back_vol->info().owned != front_vol.owned);
    auto& unowned_vol = !front_vol.owned ? front_vol : back_vol->info();
    assert(proc.neighboring_threads.count(unowned_vol.area_id));

    return proc.thread_id < proc.neighboring_threads[unowned_vol.area_id].thread_id
      ? NumberableCell::True
      : NumberableCell::False;
  }

  NumberableCell is_volume_numberable(ProcessData& proc, LCC::Dart_const_handle vol){
    constexpr auto NoneArea = AreaId{0,0,0};
    LCC& lcc = proc.lcc;
    auto vol_attr = lcc.attribute<3>(vol);
    if (vol_attr->info().area_id == NoneArea)
      return NumberableCell::None;

    return vol_attr->info().owned ? NumberableCell::True : NumberableCell::False;
  }

  template <typename HexData>
  void thread_number_vertex_in_edge(HexData& hdata,
    Dart_handle node, Dart_handle extremity0, Dart_handle extremity1){}

  template <>
  void thread_number_vertex_in_edge(ProcessData& proc,
    Dart_handle node, Dart_handle extremity0, Dart_handle extremity1)
  {
    LCC& lcc = proc.lcc;
    auto vertex_attr = lcc.attribute<0>(node);
    auto e_numberable = is_edge_numberable(proc, extremity0);

    // No need to assign ids for edges that wont ever be communicated
    if (e_numberable == NumberableCell::None) return;

    size_t id0 = lcc.attribute<0>(extremity0)->id;
    size_t id1 = lcc.attribute<0>(extremity1)->id;
    assert(id0 != DartInfo::VertexAttr::max_id);
    assert(id1 != DartInfo::VertexAttr::max_id);
    auto tid = get_temp_vertex_id<AltVertexIdSingle>(id0, id1);

    if (e_numberable == NumberableCell::True) {
      vertex_attr->id = proc.get_new_vertex_id();
      proc.vertices_to_share.vertexId[tid] = vertex_attr->id;
    } else {
      vertex_attr->id = proc.get_temp_vertex_id();
    }

    proc.vertex_temp_ids.vertexId[vertex_attr->id] = tid;
  }


  template <typename HexData>
  void thread_number_vertex_in_1t_face(HexData& hdata, Dart_handle node) {}

  template<>
  void thread_number_vertex_in_1t_face(ProcessData& proc, Dart_handle f_signature_start) {
    LCC& lcc = proc.lcc;

    NumberableCell f_numberable = is_face_numberable(proc, f_signature_start);
    // Assign temporary or fixed ids for the new vertices.
    // No need to assign ids for faces that wont ever be communicated
    if (f_numberable == NumberableCell::None)
      return;

    // Only one templates can create a vertex within a face
    Dart_handle extremity0 = lcc.beta(f_signature_start, 1);
    Dart_handle node = lcc.beta(extremity0, 1);
    Dart_handle extremity1 = lcc.beta(node, 1);
    // Dart_handle outgoing_edge = lcc.beta(f_signature_start, 1, 1, 2, 1);

    auto vertex_attr = lcc.attribute<0>(node);
    assert(vertex_attr->id == vertex_attr->max_id);

    size_t id0 = lcc.attribute<0>(extremity0)->id;
    size_t id1 = lcc.attribute<0>(extremity1)->id;

    constexpr size_t max_v_id = DartInfo::VertexAttr::max_id;
    assert(id0 != max_v_id && proc.vertex_temp_ids.vertexId.count(id0));
    assert(id1 != max_v_id && proc.vertex_temp_ids.vertexId.count(id1));
    auto tid0 = proc.vertex_temp_ids.vertexId[id0];
    auto tid1 = proc.vertex_temp_ids.vertexId[id1];

    auto tid = get_temp_vertex_id<AltVertexIdPair>(tid0, tid1);

    if (f_numberable == NumberableCell::True){
      vertex_attr->id = proc.get_new_vertex_id();
      proc.vertices_to_share.vertexIdPair[tid] = vertex_attr->id;
    } else {
      vertex_attr->id = proc.get_temp_vertex_id();
    }
    proc.vertex_temp_ids.vertexIdPair[vertex_attr->id] = tid;
  }

  template <typename HexData>
  void thread_number_vertex_in_1t_vol(HexData& hdata, Dart_handle v_signature_start) {}

  template<>
  void thread_number_vertex_in_1t_vol(ProcessData& proc, Dart_handle v_signature_start) {
    LCC& lcc = proc.lcc;

    // lcc.mark_cell<1>(v_signature_start, l_debug_mark);
    //v_signature_start, 1, 2 => ext0
    //v_signature_start, 1, 2, 0 => node
    //v_signature_start, 1, 2, 0, 0 => ext1
    //1, 2, 0, 2, 1, +1 => ext3
    // outgoing cell :

    // lcc.mark_cell<1>(lcc.beta(v_signature_start, 1, 2, 0, 2, 1), l_debug_mark_2);
    // debug_stream.push(l_thread_id);
    // assert(false);

    NumberableCell v_numberable = is_volume_numberable(proc, v_signature_start);
    if (v_numberable == NumberableCell::None) return;

    Dart_handle extremity0 = lcc.beta(v_signature_start, 1, 2);
    Dart_handle node = lcc.beta(extremity0, 0);
    Dart_handle extremity1 = lcc.beta(node, 0);
    Dart_handle extremity2 = lcc.beta(extremity1, 2, 0);

    auto vertex_attr = lcc.attribute<0>(node);
    assert(vertex_attr->id == vertex_attr->max_id);

    size_t id0 = lcc.attribute<0>(extremity0)->id;
    size_t id1 = lcc.attribute<0>(extremity1)->id;
    size_t id2 = lcc.attribute<0>(extremity2)->id;

    constexpr size_t max_v_id = DartInfo::VertexAttr::max_id;
    assert(id0 != max_v_id && proc.vertex_temp_ids.vertexIdPair.count(id1));
    assert(id1 != max_v_id && proc.vertex_temp_ids.vertexIdPair.count(id0));
    assert(id2 != max_v_id && proc.vertex_temp_ids.vertexIdPair.count(id2));
    auto tid0 = proc.vertex_temp_ids.vertexIdPair[id0];
    auto tid1 = proc.vertex_temp_ids.vertexIdPair[id1];
    auto tid2 = proc.vertex_temp_ids.vertexIdPair[id2];

    auto tid = get_temp_vertex_id<AltVertexIdTripletPair>(tid0, tid1, tid2);

    if (v_numberable == NumberableCell::True){
      vertex_attr->id = proc.get_new_vertex_id();
      proc.vertices_to_share.vertexIdTripletPair[tid] = vertex_attr->id;
    } else {
      vertex_attr->id = proc.get_temp_vertex_id();
    }

    proc.vertex_temp_ids.vertexIdTripletPair[vertex_attr->id] = tid;
  }

  template <typename HexData>
  void thread_join_3_template_vertex__pair(HexData& hdata, Dart_handle edge) {}

  void thread_join_3_template_vertex__pair(ProcessData& proc, Dart_handle edge) {
    LCC& lcc = proc.lcc;

    NumberableCell v_numberable = is_face_numberable(proc, edge);
    if (v_numberable == NumberableCell::None)
      return;

    Dart_handle ext0 = edge;
    Dart_handle node = lcc.beta(edge, 1);
    Dart_handle ext1 = lcc.beta(node, 1);

    auto vertex_attr = lcc.attribute<0>(node);
    size_t id0 = lcc.attribute<0>(ext0)->id;
    size_t id1 = lcc.attribute<0>(ext1)->id;

    constexpr size_t max_v_id = DartInfo::VertexAttr::max_id;
    assert(id0 != max_v_id && proc.vertex_temp_ids.vertexId.count(id0));
    assert(id1 != max_v_id && proc.vertex_temp_ids.vertexId.count(id1));
    auto tid0 = proc.vertex_temp_ids.vertexId[id0];
    auto tid1 = proc.vertex_temp_ids.vertexId[id1];
    auto tid = get_temp_vertex_id<AltVertexIdPair>(tid0, tid1);

    if (v_numberable == NumberableCell::True){
      vertex_attr->id = proc.get_new_vertex_id();
      proc.vertices_to_share.vertexId[tid0] = vertex_attr->id;
      proc.vertices_to_share.vertexId[tid1] = vertex_attr->id;
      proc.vertices_to_share.vertexIdPair[tid] = vertex_attr->id;
    } else {
      vertex_attr->id = proc.get_temp_vertex_id();
    }

    proc.vertex_temp_ids.vertexIdPair[vertex_attr->id] = tid;

  }

  template <typename HexData>
  void thread_join_3_template_vertex__pairpair(HexData& hdata, Dart_handle edge) {}

  void thread_join_3_template_vertex__pairpair(ProcessData& proc, Dart_handle edge) {
    LCC& lcc = proc.lcc;

    NumberableCell v_numberable = is_volume_numberable(proc, edge);
    if (v_numberable == NumberableCell::None)
      return;

    Dart_handle ext0 = edge;
    Dart_handle node = lcc.beta(edge, 1);
    Dart_handle ext1 = lcc.beta(edge, 1, 1);

    auto vertex_attr = lcc.attribute<0>(node);
    size_t id0 = lcc.attribute<0>(ext0)->id;
    size_t id1 = lcc.attribute<0>(ext1)->id;

    constexpr size_t max_v_id = DartInfo::VertexAttr::max_id;
    assert(id0 != max_v_id && proc.vertex_temp_ids.vertexIdPair.count(id0));
    assert(id1 != max_v_id && proc.vertex_temp_ids.vertexIdPair.count(id1));
    auto tid0 = proc.vertex_temp_ids.vertexIdPair[id0];
    auto tid1 = proc.vertex_temp_ids.vertexIdPair[id1];
    proc.vertex_temp_ids.vertexIdPair.erase(id0);
    proc.vertex_temp_ids.vertexIdPair.erase(id1);

    auto tid = get_temp_vertex_id<AltVertexIdPairPair>(tid0, tid1);

    if (v_numberable == NumberableCell::True){
      vertex_attr->id = proc.get_new_vertex_id();
      proc.vertices_to_share.vertexIdPair[tid0] = vertex_attr->id;
      proc.vertices_to_share.vertexIdPair[tid1] = vertex_attr->id;
      proc.vertices_to_share.vertexIdPairPair[tid] = vertex_attr->id;
    } else {
      vertex_attr->id = proc.get_temp_vertex_id();
    }

    proc.vertex_temp_ids.vertexIdPairPair[vertex_attr->id] = tid;

  }

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

      assert(temp_id < SIZE_T_MAX);
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
      if (temp_id < SIZE_T_MAX)
        nbsub++;

      if (temp_id == 0){
        thread_number_vertex_in_1t_vol(hdata, v_signature_start);
      }

      // if (temp > 0) nbsub++;
      // if (temp == 2) temp = 3;
      // if (temp < 10) temp++;
      // assert(vol_attr.template_id >= 0 && vol_attr.template_id <= 4 && vol_attr.template_id != 3 || vol_attr.template_id == 8);
      // assert(vol_attr.template_id == temp || temp == SIZE_T_MAX && (vol_attr.template_id == 0 || vol_attr.template_id == 8));
    }

    return nbsub;
  }

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

  template <>
  void thread_communicate_marked_nodes(ProcessData& proc, RefinementData& rdata, size_type edge_mark){
    LCC& lcc = proc.lcc;

    size_type node_mark = lcc.get_new_mark();

    ThreadMarkedCellsPtr marked_cells = std::make_shared<ThreadMarkedCellsMsg>();
    std::vector<ThreadMarkedCellsPtr> others_msg;

    others_msg.reserve(proc.neighboring_threads.size());

    auto retrieve_nodes_in_face = [&](Dart_handle face){
      auto nodes =  lcc.darts_of_cell<2,0>(face);
      for (auto it = nodes.begin(); it != nodes.end(); it++){
        if (!lcc.is_marked(it, proc.template_mark)) continue;

        size_t id = lcc.attribute<0>(it)->id;
        assert(id != DartInfo::VertexAttr::max_id);
        marked_cells->insert(id);
      }
    };

    auto mark_nodes_in_face = [&](Dart_handle face){
      auto nodes =  lcc.darts_of_cell<2,0>(face);
      for (auto it = nodes.begin(); it != nodes.end(); it++){
        if (lcc.is_marked(it, node_mark)) continue;
        lcc.mark_cell<0>(it, node_mark);

        size_t id = lcc.attribute<0>(it)->id;

        bool marked = false;
        for (ThreadMarkedCellsPtr& set : others_msg){
          if (set->count(id)) {
            marked = true;
            break;
          }
        }

        if (!marked or lcc.is_marked(it, proc.template_mark))
          continue;

        // Mark the node and update neighboring faces
        lcc.mark_cell<0>(it, proc.template_mark);
        rdata.marked_nodes.push_back(it);

        // Find faces associated with this node
        for (Dart_handle face : plane_faces_around_node(lcc, rdata, it)){
          auto& face_attr =  lcc.attribute<2>(face)->info();

          if (face_attr.template_id == 0)
            rdata.faces_of_plane.push_back(face);

          face_attr.template_id++;
          assert(face_attr.template_id <= 4);

          // Also gather additionnal volumes if the edge is not marked
          // faces_around_node also gives all unique edges on the plane
          if (lcc.is_marked(face, edge_mark)) continue;
          lcc.mark_cell<1>(face, edge_mark);

          // TODO This could be used with plane_faces_around_node code instead of doing twice the work
          __adjacent_face_on_plane(lcc, rdata.iteration, face, rdata.additionnal_volumes_found);
          // Also do the other way around to catch volumes in the opposite dir
          if (!lcc.is_free<3>(face)){
            __adjacent_face_on_plane(lcc, rdata.iteration, lcc.beta(face, 3), rdata.additionnal_volumes_found);
          }
        }
      }
    };

    auto& plane_set = proc.first_face_of_planes[rdata.iteration];
    for (int i = 1; i < plane_set.size(); i += 2){
      plane_for_each_face(lcc, plane_set[i],
        [&](Dart_handle face){
          auto back_vol = lcc.attribute<3>(lcc.beta(face, 3));
          auto front_vol = lcc.attribute<3>(face);

          bool owned_zone = front_vol->info().owned or (back_vol != nullptr && back_vol->info().owned);

          if (owned_zone){
            retrieve_nodes_in_face(face);
          }
        },
        [&](Dart_handle edge){
          return adjacent_face_on_plane(lcc, rdata.iteration, edge);
        }
      );
    }

    for (auto& [area, thread] : proc.neighboring_threads){
      // Use a shared pointer to avoid copy
      thread.stream << marked_cells;
    }

    for (auto& [area, thread] : proc.neighboring_threads){
      ThreadMsg msg;
      thread.stream >> msg;

      if (!std::holds_alternative<ThreadMarkedCellsPtr>(msg)){
        assert(false);
        exit(0);
      }

      others_msg.push_back(std::get<ThreadMarkedCellsPtr>(msg));
    }

    for (int i = 1; i < plane_set.size(); i += 2){
      plane_for_each_face(lcc, plane_set[i],
        [&](Dart_handle face){
          auto back_vol = lcc.attribute<3>(lcc.beta(face, 3));
          auto front_vol = lcc.attribute<3>(face);
          auto face_attr = lcc.attribute<2>(face);

          bool unowned_zone = !front_vol->info().owned or (back_vol != nullptr && !back_vol->info().owned);

          if (unowned_zone){
            mark_nodes_in_face(face);
          }
        },
        [&](Dart_handle edge){
          return adjacent_face_on_plane(lcc, rdata.iteration, edge);
        });
    }

    lcc.free_mark(node_mark);
  }

  // In assertion mode, we check if no duplicates exists
  template <bool assertion = false>
  AltVertexIdVariants get_node_tid(const ProcessData& proc, const size_t& id) {
    AltVertexIdVariants ret;

    if (proc.vertex_temp_ids.vertexId.count(id)){
      if (assertion) CGAL_assertion_msg(ret.index() == 0, "Duplicate temp id found");
      ret = proc.vertex_temp_ids.vertexId.at(id);
      if (!assertion) return ret;
    }
    else if (proc.vertex_temp_ids.vertexIdPair.count(id)){
      if (assertion) CGAL_assertion_msg(ret.index() == 0, "Duplicate temp id found");
      ret = proc.vertex_temp_ids.vertexIdPair.at(id);
      if (!assertion) return ret;
    }
    else if (proc.vertex_temp_ids.vertexIdPairPair.count(id)){
      if (assertion) CGAL_assertion_msg(ret.index() == 0, "Duplicate temp id found");
      ret = proc.vertex_temp_ids.vertexIdPairPair.at(id);
      if (!assertion) return ret;
    }
    else if (proc.vertex_temp_ids.vertexIdTripletPair.count(id)){
      if (assertion) CGAL_assertion_msg(ret.index() == 0, "Duplicate temp id found");
      ret = proc.vertex_temp_ids.vertexIdTripletPair.at(id);
      if (!assertion) return ret;
    }

    CGAL_assertion_msg(ret.index() != 0, "Temporary node ID not found");
    return ret;
  };

  template <bool assertion = false>
  size_t get_node_rid_from_tid(const std::vector<ThreadCellsIdPtr>& others_msg, const AltVertexIdVariants& vtid) {
    size_t ret = DartInfo::VertexAttr::max_id;

    for (const ThreadCellsIdPtr& ptr : others_msg){
      if (std::holds_alternative<AltVertexIdSingle>(vtid)){
        const AltVertexIdSingle& tid = std::get<AltVertexIdSingle>(vtid);
        if (ptr->ids.vertexId.count(tid)){
          if (assertion) CGAL_assertion_msg(ret == DartInfo::VertexAttr::max_id, "Duplicate real id found");
          ret = ptr->ids.vertexId[tid];
          if (!assertion) return ret;
        }
      }

      if (std::holds_alternative<AltVertexIdPair>(vtid)){
        const AltVertexIdPair& tid = std::get<AltVertexIdPair>(vtid);
        if (ptr->ids.vertexIdPair.count(tid)){
          if (assertion) CGAL_assertion_msg(ret == DartInfo::VertexAttr::max_id, "Duplicate real id found");
          ret = ptr->ids.vertexIdPair[tid];
          if (!assertion) return ret;
        }
      }

      if (std::holds_alternative<AltVertexIdPairPair>(vtid)){
        const AltVertexIdPairPair& tid = std::get<AltVertexIdPairPair>(vtid);
        if (ptr->ids.vertexIdPairPair.count(tid)){
          if (assertion) CGAL_assertion_msg(ret == DartInfo::VertexAttr::max_id, "Duplicate real id found");
          ret = ptr->ids.vertexIdPairPair[tid];
          if (!assertion) return ret;
        }
      }

      if (std::holds_alternative<AltVertexIdTripletPair>(vtid)){
        const AltVertexIdTripletPair& tid = std::get<AltVertexIdTripletPair>(vtid);
        if (ptr->ids.vertexIdTripletPair.count(tid)){
          if (assertion) CGAL_assertion_msg(ret == DartInfo::VertexAttr::max_id, "Duplicate real id found");
          ret = ptr->ids.vertexIdTripletPair[tid];
          if (!assertion) return ret;
        }
      }
    }

    CGAL_assertion_msg(ret != DartInfo::VertexAttr::max_id, "Real id not found");
    return ret;
  };

  template <typename HexData>
  void thread_communicate_cells_id_and_3t(HexData&, RefinementData&){}

  template <>
  void thread_communicate_cells_id_and_3t(ProcessData& proc, RefinementData& rdata){
    LCC& lcc = proc.lcc;

    ThreadCellsIdPtr threadMsg = std::make_shared<ThreadCellsIdMsg>();
    threadMsg->ids = std::move(proc.vertices_to_share);

    std::vector<ThreadCellsIdPtr> others_msg;

    for (auto [aid, thread] : proc.neighboring_threads){
      thread.stream << threadMsg;
    }

    for (auto& [area, thread] : proc.neighboring_threads){
      ThreadMsg msg;
      thread.stream >> msg;

      if (!std::holds_alternative<ThreadCellsIdPtr>(msg)){
        assert(false);
        exit(0);
      }

      others_msg.push_back(std::get<ThreadCellsIdPtr>(msg));
    }

    const auto node_operation = [&](Dart_handle node){
      size_t& v_id = lcc.attribute<0>(node)->id;

      if (!proc.is_temp_id(v_id))
        return;

      // Attribute node id
      size_t r_id = [&](){
        #ifdef CGAL_ASSERTIONS_ENABLED
          return get_node_rid_from_tid<true>(others_msg, get_node_tid<true>(proc, v_id));
        #else
          return get_node_rid_from_tid(others_msg, get_node_tid(proc, v_id));
        #endif
      }();

      v_id = r_id;
    };

    // TODO carefull, make a version that takes a vector to put {owned, unowned}, one of them may be null;
    // IF both are null then we have nothing to process.
    // Explore all assignable volumes
    for_each_hex(lcc, proc.owned_ghosts_begin,
      [&](Dart_handle vol, auto faces){
        for (Dart_handle node : nodes_of_hex(lcc, vol)){
          node_operation(node);
        }
      },
      [&](Dart_handle vol) -> bool{
        auto attr = lcc.attribute<3>(vol)->info();
        return attr.area_id != AreaId{0,0,0};
      });


    thread_assert_no_temporary_ids(proc);

    // lcc.unmark_all(proc.three_template_node_mark);
  }

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

  // If not all marks are chained consecutively, mark the whole face
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
   * Mark 0-cells
   * Gather 2-cells adjacent to marked 0-cells
   * Gather 3-cells adjacents to marked 0-cells
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
    if constexpr (std::is_same_v<HexData, ProcessData>){
      auto& front_vol = lcc.attribute<3>(face)->info();
      auto back_vol_attr = lcc.attribute<3>(lcc.beta<3>(face));
      is_markable = front_vol.owned or (back_vol_attr != nullptr && back_vol_attr->info().owned);
    }

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
    assert_faces_of_plane_valid(hdata, rdata);

    fix_impossible_cases(hdata, rdata);
    assert_faces_of_plane_valid(hdata, rdata);

    propagation_stage(hdata, rdata, explored_face);

    get_cells_to_refine_from_plane(hdata, rdata,  explored_face);
    get_cells_to_refine_from_additionnal_volumes(hdata, rdata, explored_face);

    lcc.free_mark(explored_edge);
    lcc.free_mark(explored_face);

  }

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

      // Delete volumes outside of refinement domain
      if (old_info.type <= VolumeType::REFINEMENT ){
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

  void expand_identified_cells(HexMeshingData& hdata, int current_lvl, int nb_levels){
    LCC& lcc = hdata.lcc;

    assert(nb_levels >= 1 && current_lvl >= 0);

    // TODO Reformulate, information not up to date

    // Calculate the totaling cells needed per level (height), with i ranging from 0 (lowest level) to n (highest level)
    // Very simple sequence : 4, 6, 6, 6, 6, ....               (height of n cells, cell sized from current i+1-subdivision)
    // Divide this by 2 : height of n cells, cell sized from current i-subdivision
    auto height_of_refinement_level = [](int i, int nb_levels) -> int{
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
       for (auto& plane_cc :  hdata.first_face_of_planes[p]) {
         plane_for_each_face(lcc, plane_cc,
           [&](Dart_handle face){
             auto& face_attr = get_or_create_attr<2>(lcc, face)->info();
             face_attr.plane[p] = true;
           },
           [&](Dart_handle edge){
             return lcc.beta(edge, 2, 3, 2);
           });
       }

       // Last odd plane
       int size = hdata.first_face_of_planes[p].size();
       Dart_handle start = hdata.first_face_of_planes[p][size-1][0];
       start = lcc.beta(start, 2, 1, 1, 2);

       plane_for_each_face(lcc, start,
         [&](Dart_handle face){
           auto& face_attr = get_or_create_attr<2>(lcc, face)->info();
           face_attr.plane[p] = true;
         },
         [&](Dart_handle edge){
           return lcc.beta(edge, 2, 3, 2);
         });
     }

  }

  // TODO, messy code but it works
  // Only works on a uniform grid
  Dart_handle find_area_on_grid(ProcessData& proc, Dart_handle from, PlaneNormal fromPlane, const AreaId& area, bool owned){
    LCC& lcc = proc.lcc;

    const auto can_move_right = [&](Dart_handle d) { return !lcc.is_free<3>(lcc.beta(d, 1, 2)); };
    const auto right = [&](Dart_handle d){ return lcc.beta(d, 1, 2, 3, 2, 1); };
    const auto can_move_up = [&](Dart_handle d) { return !lcc.is_free<3>(lcc.beta(d, 1, 1, 2)); };
    const auto up = [&](Dart_handle d){ return lcc.beta(d, 1, 1, 2, 3, 2); };

    const int x_axis = (fromPlane + 1)%3;
    const int y_axis = (fromPlane + 2)%3;
    Dart_handle position = from;

    // Offset rules for owned areas
    // If area searched is inside a owned area, offset the starting dart to be in a owned area
    if (owned){
      if (proc.negative_axes[y_axis]){
        assert(!lcc.attribute<3>(position)->info().owned);
        position = up(up(position));
      }

      if (proc.negative_axes[x_axis]){
        assert(!lcc.attribute<3>(position)->info().owned);
        position = right(right(position));
      }

      // assert(is_volume_owned(position));
    }

    // Offset rules for unowned area
    // Because two areas cannot overlap
    else if (!owned){
      bool x_offset = proc.negative_axes[x_axis] && area[x_axis] == 0;
      bool y_offset = proc.negative_axes[y_axis] && area[y_axis] == 0;

      if (x_offset) position = right(right(position));
      if (y_offset) position = up(up(position));
    }

    const DartInfo::VolumeAttrValue* vol_attr = &lcc.attribute<3>(position)->info();

    if (area[x_axis] > 0){
      while (can_move_right(position) && (vol_attr->area_id[x_axis] <= 0 or vol_attr->owned != owned)){
        position = right(position);
        vol_attr = &lcc.attribute<3>(position)->info();
      }
    }

    if (area[y_axis] > 0){
      while (can_move_up(position) && (vol_attr->area_id[y_axis] <= 0 or vol_attr->owned != owned)){
        position = up(position);
        vol_attr = &lcc.attribute<3>(position)->info();
      }
    }

    if (!are_areas_equivalent(vol_attr->area_id, area, owned) or vol_attr->owned != owned){
      return lcc.null_dart_descriptor;
    }

    return position;
  };

  void assign_ghost_handles(ProcessData& proc){
    auto& volumes = proc.lcc.attributes<3>();
    bool found_owned = false, found_unowned = false;
    proc.unowned_ghosts_begin = nullptr;
    proc.owned_ghosts_begin = nullptr;

    for (auto& handle : volumes){
      if (!found_owned && handle.info().owned && handle.info().area_id != AreaId(0,0,0)){
        found_owned = true;
        proc.owned_ghosts_begin = handle.dart();
      }

      if (!found_unowned && !handle.info().owned){
        found_unowned = true;
        proc.unowned_ghosts_begin = handle.dart();
      }

      if (found_owned && found_unowned) break;
    }

    if (proc.level == 0){
      assert(found_owned && found_unowned);
    }
  }

  void setup_ghost_areas(ProcessData& proc){
    LCC& lcc = proc.lcc;
    // TODO Le code peut etre mieux écrit que ça ?

    // Attributes a AreaID to each volumes and sets the ownership of an area
    // If the area is owned, areas id can overlap on each other, equality is sufficient if all non zero dimensions are equals
    // If the area is unowned, no overlaps happens, two id are must be strictly equal to belong in the same area
    auto mark_ghost_area = [&](std::vector<Dart_handle> starts, int& plane, bool positive_axis, bool owned){
      plane_for_each_face(lcc, starts,
        [&](Dart_handle face){
          auto& front_vol = lcc.attribute<3>(face)->info();
          auto back_vol_attr = lcc.attribute<3>(lcc.beta(face, 3));

          // Return if we are stepping on unowned area, only when owned is set
          if (owned && (!front_vol.owned or back_vol_attr != nullptr && !back_vol_attr->info().owned))
            return;

          front_vol.area_id[plane] = positive_axis ? 1 : -1;
          front_vol.owned = owned;

          if (back_vol_attr != nullptr){
            auto& back_vol = back_vol_attr->info();
            back_vol.area_id[plane] = positive_axis ? 1 : -1;
            back_vol.owned = owned;
          }
        },
        [&](Dart_handle edge){
          return lcc.beta(edge, 2, 3, 2);
        }
        );
    };

    // First, mark unowned areas
    for (int p = 0; p < 3; p++){
      auto& plane_set = proc.first_face_of_planes[p];
      assert(plane_set.size() > 4);

      if (proc.positive_axes[p]){
        auto& last_even_plane = plane_set[plane_set.size() - 1]; // -1 because we don't store the last odd plane
        // Mark the last odd-even set of volumes outside the process frontier
        mark_ghost_area(last_even_plane, p, true, false);
      }

      if (proc.negative_axes[p]){
        auto& first_even_plane = plane_set[1];
        // Mark the last odd-even set of volumes outside the process frontier
        mark_ghost_area(first_even_plane, p, false, false);
      }
    }

    // then mark ghost owned areas
    for (int p = 0; p < 3; p++){
      auto& plane_set = proc.first_face_of_planes[p];

      assert(plane_set.size() > 4);

      if (proc.positive_axes[p]){
        auto& last_owning_plane = plane_set[plane_set.size() - 3];
        // Mark the last odd-even set of volumes inside the process frontier
        mark_ghost_area(last_owning_plane, p, true, true);
      }

      if (proc.negative_axes[p]){
        auto& first_owning_plane = plane_set[3];
        // Mark the last odd-even set of volumes inside the process frontier
        mark_ghost_area(first_owning_plane, p, false, true);
      }
    }

    assign_ghost_handles(proc);
  }

  void setup_vertex_ids(ProcessData& proc){
    // Setup ids for volumes inside a ghosted cell

    LCC& lcc = proc.lcc;
    Grid& grid = proc.grid;

    const auto can_move_up = [&](Dart_handle d) { return !lcc.is_free<3>(lcc.beta(d, 1, 1, 2)); };
    const auto up = [&](Dart_handle d){ return lcc.beta(d, 1, 1, 2, 3, 2); };

    size_t size_of_planes = proc.first_face_of_planes[0].size();

    // Plane xyz origin => yzx in global space;
    size_t line_size = proc.full_grid.y + 1;
    size_t plane_size = (proc.full_grid.y + 1) * (proc.full_grid.z + 1);

    size_t z_dim_offset = grid.dims_id_start.x * plane_size;
    size_t y_dim_offset = grid.dims_id_start.z * line_size;
    size_t x_dim_offset = grid.dims_id_start.y;

    const auto assign_plane = [&](Dart_handle start, int plane_id, int forward){
      const auto can_move_right = [&](Dart_handle it) -> bool {
        return !lcc.is_free<3>(lcc.beta(it, forward, 2));
      };

      const auto can_move_backwards = [&](Dart_handle it) -> bool {
        return !lcc.is_free<3>(lcc.beta(it, (int)!forward, 2));
      };

      const auto right = [&](Dart_handle it){
        return lcc.beta(it, forward, 2, 3, 2, forward);
      };

      const auto backwards = [&](Dart_handle it){
        return lcc.beta(it, (int)!forward, 2, 3, 2, (int)!forward);
      };

      Dart_handle it = start;
      int y = 0;
      while (true){
        size_t id = z_dim_offset + y_dim_offset + x_dim_offset
                  + plane_id * plane_size
                  + y * line_size;

        // Assign ids, left to right, bottom to up
        if (!forward) // Number the first vertex
          lcc.attribute<0>(lcc.other_extremity(it))->id = id++;

        while (can_move_right(it)){
          size_t& v_id = lcc.attribute<0>(it)->id;
          assert(v_id == DartInfo::VertexAttr::max_id);
          v_id = id++;
          it = right(it);
        }

        // Number the last dart
        size_t& v_id = lcc.attribute<0>(it)->id;
        assert(v_id == DartInfo::VertexAttr::max_id);
        v_id = id++;

        if (forward) // Number the last vertex
          lcc.attribute<0>(lcc.other_extremity(it))->id = id++;

        y++;
        if (!can_move_up(start)) break;
        start = it = up(start);
      }

      // Assign ids on the last line
      start = it = lcc.beta(start, 0, 0);
      size_t id = z_dim_offset + y_dim_offset + x_dim_offset
                + plane_id * plane_size
                + y * line_size;

      if (forward) // Number the first vertex if the dart handle is of opposite direction
        lcc.attribute<0>(lcc.other_extremity(it))->id = id++;

      while (can_move_backwards(it)){
        size_t& v_id = lcc.attribute<0>(it)->id;
        assert(v_id == DartInfo::VertexAttr::max_id);
        v_id = id++;
        it = backwards(it);
      }

      // Number the last dart
      size_t& v_id = lcc.attribute<0>(it)->id;
      assert(v_id == DartInfo::VertexAttr::max_id);
      v_id = id++;

      if (!forward) // Number the last vertex
        lcc.attribute<0>(lcc.other_extremity(it))->id = id++;
    };

    Dart_handle start, it;
    for (int x = 0; x < size_of_planes; x++){
      start = it = proc.first_face_of_planes[0][x][0];
      size_t id = z_dim_offset + y_dim_offset + x_dim_offset + plane_size * x;
      assign_plane(start, x, true);
    }

    // Assign ids on the last odd plane,
    start = proc.first_face_of_planes[0][size_of_planes - 1][0];
    start = it = lcc.beta(start, 2, 0, 0, 2);
    assert(lcc.is_free<3>(start));
    assign_plane(start, size_of_planes, false);
  }

  // Sets up the plane between two processors boundaries.
  // The thread with the lowest id owns the interface
  void setup_volume_boundaries(ProcessData& proc){
    // LCC& lcc = proc.lcc;
    // const auto can_move = [&](Dart_handle it, int f) { return !lcc.is_free<3>(lcc.beta(it, f, 2)); };
    // const auto move = [&](Dart_handle it, int f) { return lcc.beta(it, f, 2, 3, 2, f); };

    // // Since we only use the attr to determine whether a edge is owned or not.
    // // We just need create only one attribute for all darts
    // proc.numberable_edge_attr = lcc.create_attribute<1>(true);
    // proc.unnumberable_edge_attr = lcc.create_attribute<1>(false);

    // const auto is_edge_numberable = [&](Dart_handle edge) -> bool {
    //   auto face_orbit = lcc.darts_of_orbit<2,3>(edge);

    //   for (auto it = face_orbit.begin(); it != face_orbit.end(); it++){
    //     auto vol = lcc.attribute<3>(it)->info();
    //     if (vol.owned or vol.area_id == AreaId{0,0,0}) {
    //       continue;
    //     }

    //     assert(proc.neighboring_threads.count(vol.area_id));
    //     auto other_id = proc.neighboring_threads[vol.area_id].thread_id;
    //     if (other_id < proc.thread_id)
    //       return false;
    //   }

    //   return true;
    // };

    // const auto assign_edge_boundaries = [&](Dart_handle vol, std::array<Dart_handle,6> faces){
    //   auto edges = lcc.darts_of_cell<3,1>(vol);
    //   for (auto it = edges.begin(); it != edges.end(); it++){
    //     auto edge_attr = lcc.attribute<1>(it);
    //     if (edge_attr != nullptr) continue;

    //     auto& attr = is_edge_numberable(it) ? proc.numberable_edge_attr : proc.unnumberable_edge_attr;
    //     lcc.set_attribute<1>(it, attr);
    //   }
    // };

    // const auto assign_edge_unowned_boundaries = [&](Dart_handle vol, std::array<Dart_handle,6> faces){
    //   auto edges = lcc.darts_of_cell<3,1>(vol);
    //   for (auto it = edges.begin(); it != edges.end(); it++){
    //     auto edge_attr = lcc.attribute<1>(it);
    //     if (edge_attr != nullptr) continue;

    //     lcc.set_attribute<1>(it, proc.unnumberable_edge_attr);
    //   }
    // };


    // size_type face_mark = lcc.get_new_mark();

    // for_each_hex(lcc, proc.owned_ghosts_begin,
    //   [&](Dart_handle vol, std::array<Dart_handle,6> faces) {
    //     assign_edge_boundaries(vol, faces);
    //   },
    //   [&](Dart_handle vol){
    //     auto attr = lcc.attribute<3>(vol);
    //     return attr != nullptr && attr->info().owned && attr->info().area_id != AreaId{0,0,0};
    //   });

    // for_each_hex(lcc, proc.unowned_ghosts_begin,
    //   [&](Dart_handle vol, std::array<Dart_handle,6> faces) {
    //     assign_edge_unowned_boundaries(vol, faces);
    //   },
    //   [&](Dart_handle vol){
    //     auto attr = lcc.attribute<3>(vol);
    //     return attr != nullptr && !attr->info().owned;
    //   });

    // lcc.free_mark(face_mark);

    // debug_stream.push(l_thread_id);
    // assert(false);
  }

  // Gather the faces first, if we don't we will have issues later manually
  // Iterating to the next plane after a first refinement stage
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

  void initial_setup(ProcessData& proc, MarkingFunction& cellIdentifier){
    LCC& lcc = proc.lcc;

    auto volumes = lcc.one_dart_per_cell<3>();
    // Mark initial identified cells
    for (auto dart = volumes.begin(), end = volumes.end(); dart != end; dart++){
      // Create a 3-attr for all 3-cells in the LCC
      auto& vol_attr = get_or_create_attr<3>(lcc, dart)->info();
    }

    setup_initial_planes(proc);
    setup_ghost_areas(proc);
    setup_vertex_ids(proc);
    setup_volume_boundaries(proc);

    // Mark initial identified cells
    for (auto dart = volumes.begin(), end = volumes.end(); dart != end; dart++){
      auto& vol_attr = lcc.attribute<3>(dart)->info();

      // Mark those who are identified
      if (vol_attr.owned && cellIdentifier(lcc, dart))
        vol_attr.type = VolumeType::IDENTIFIED;
    }
  }


  void setup_next_level_face(HexMeshingData& hdata,
                            std::queue<Dart_handle>& to_explore,
                            Union_find<Dart_handle>& odd_union_find,
                            Union_find<Dart_handle>& even_union_find,
                            Dart_handle face,
                            PlaneNormal planeIteration,
                            int edge_mark, int face_mark){
    LCC& lcc = hdata.lcc;
    auto vol_handle = lcc.attribute<3>(face);
    auto beta3_vol_handle = lcc.attribute<3>(lcc.beta<3>(face));

    bool identified = vol_handle != nullptr && vol_handle->info().type >= VolumeType::ID_EXPANSION;
    bool beta3_identified = beta3_vol_handle != nullptr && beta3_vol_handle->info().type >= VolumeType::ID_EXPANSION;
    bool cc_wave = (hdata.level % 2) == 0;  // Cycling between true and false each level
                                            // To prevent doing one full iteration to remove all cc_ids

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
        auto& adj_face_attr = lcc.attribute<2>(adjacent_face)->info();

        if (odd_cc_id != nullptr
          && adj_face_attr.cc_wave == cc_wave
          && adj_face_attr.cc_id != nullptr
          && !odd_union_find.same_set(odd_cc_id, adj_face_attr.cc_id)){
          odd_union_find.unify_sets(odd_cc_id, adj_face_attr.cc_id);
        }

        if (odd_cc_id == nullptr
          && adj_face_attr.cc_wave == cc_wave
          && adj_face_attr.cc_id != nullptr)
          odd_cc_id = odd_union_find.find(adj_face_attr.cc_id);
      }

      // Even plane union find, only if the studied volume is identified
      if (identified){
        auto back_face_handle = lcc.attribute<2>(lcc.beta(adjacent_face, 2, 1, 1, 2));
        // Skip back face with no attributes
        if (back_face_handle == nullptr) continue;
        auto& back_face_attr = back_face_handle->info();

        if (even_cc_id != nullptr
          && back_face_attr.cc_wave == cc_wave
          && back_face_attr.cc_id != nullptr
          && !even_union_find.same_set(even_cc_id, back_face_attr.cc_id)){
          even_union_find.unify_sets(even_cc_id, back_face_attr.cc_id);
        }

        if (even_cc_id == nullptr
          && back_face_attr.cc_wave == cc_wave
          && back_face_attr.cc_id != nullptr)
          even_cc_id = even_union_find.find(back_face_attr.cc_id);
      }
    }

    if (identified or beta3_identified){
      auto& face_attr = lcc.attribute<2>(face)->info();
      face_attr = DartInfo::FaceAttrValue();
      face_attr.plane[planeIteration] = true;
      face_attr.cc_wave = cc_wave;
      face_attr.cc_id = odd_cc_id != nullptr ? odd_cc_id : odd_union_find.make_set(face);
    }

    if (identified) {
      assert(lcc.attribute<2>(back_face) == nullptr);
      assert(!lcc.is_free<3>(back_face));

      auto& back_face_attr = get_or_create_attr<2>(lcc, back_face)->info();
      back_face_attr.plane[planeIteration] = true;
      back_face_attr.cc_wave = cc_wave;
      back_face_attr.cc_id = even_cc_id != nullptr ? even_cc_id : even_union_find.make_set(lcc.beta<3>(back_face));
    }
  }

  void setup_next_level_plane(HexMeshingData& hdata){
    LCC& lcc = hdata.lcc;
    std::array<PlaneSet, 3> new_planes;

    size_type edge_mark = lcc.get_new_mark();
    size_type face_mark = lcc.get_new_mark();

    for (int p = 0; p < 3; p++){
      PlaneSet new_plane_set;

      int i = 0;
      for (auto& old_plane_set : hdata.first_face_of_planes[p]){
        std::queue<Dart_handle> to_explore;
        Union_find<Dart_handle> odd_union_find, even_union_find;


        for (auto start : old_plane_set)
          to_explore.push(start);

        while (!to_explore.empty()){
          Dart_handle face = to_explore.front();
          to_explore.pop();

          setup_next_level_face(hdata, to_explore, odd_union_find, even_union_find, face, (PlaneNormal)p, edge_mark, face_mark);
        }

        PlaneCC odd_plane_cc, even_plane_cc;
        odd_plane_cc.reserve(odd_union_find.number_of_sets());
        even_plane_cc.reserve(odd_union_find.number_of_sets());

        for (auto uf_handle : get_partitions(odd_union_find)){
          Dart_handle face_cc = uf_handle->value;
          odd_plane_cc.push_back(face_cc);
        }

        for (auto uf_handle : get_partitions(even_union_find)){
          Dart_handle face_cc = uf_handle->value;
          even_plane_cc.push_back(face_cc);
        }

        new_plane_set.push_back(std::move(odd_plane_cc));
        new_plane_set.push_back(std::move(even_plane_cc));
      }

      new_planes[p] = std::move(new_plane_set);

      lcc.unmark_all(edge_mark);
      lcc.unmark_all(face_mark);
    }

    hdata.first_face_of_planes = std::move(new_planes);

    lcc.free_mark(edge_mark);
    lcc.free_mark(face_mark);
  }

  void setup_next_level(HexMeshingData& hdata, MarkingFunction& cellIdentifier){
    setup_next_level_plane(hdata);
    clean_up_and_reevaluate_attributes(hdata, cellIdentifier);
    hdata.level++;
  }

  void setup_next_level(ProcessData& proc, MarkingFunction& cellIdentifier){
    setup_next_level(static_cast<HexMeshingData&>(proc), cellIdentifier);
  }

  void mark_identified_cells_from_3_attrs(HexMeshingData& hdata) {
    LCC& lcc = hdata.lcc;

    auto& attributes = lcc.attributes<3>();

    for (auto it = attributes.begin(), end = attributes.end(); it != end; it++){
      if (it->info().type > VolumeType::NONE && it->info().owned){
        mark_k_cells_of_i_cell<3, 0>(lcc, it->dart(), hdata.identified_mark);
      }
    }
  }

  void generate_grid(HexMeshingData& hdata) {
    Grid& grid = hdata.grid;
    for (int x = 0; x < grid.dims.x; x++) {
      for (int y = 0; y < grid.dims.y; y++) {
        for (int z = 0; z < grid.dims.z; z++) {
          double x1 = grid.pos.x() + x * grid.size.x();
          double y1 = grid.pos.y() + y * grid.size.y();
          double z1 = grid.pos.z() + z * grid.size.z();

          double x2 = grid.pos.x() + (x+1)*grid.size.x();
          double y2 = grid.pos.y() + (y+1)*grid.size.y();
          double z2 = grid.pos.z() + (z+1)*grid.size.z();

          hdata.lcc.make_hexahedron(Point(x1,y1,z1), Point(x2,y1,z1),
                              Point(x2,y2,z1), Point(x1,y2,z1),
                              Point(x1,y2,z2), Point(x1,y1,z2),
                              Point(x2,y1,z2), Point(x2,y2,z2));
        }
      }
    }

    hdata.lcc.sew3_same_facets();
  }

  Grid calculate_thread_grid(ProcessData& proc, const Grid& full_grid, int thread_id, int nb_threads){
    Grid grid;

    // Starting point can be displaced if neighboring cells
    // exists on -X, -Y, -Z planes

    for (auto& [rel_pos, _] : proc.neighboring_threads){
      for (int i = 0; i < 3; i++){
        if (rel_pos[i] < 0) proc.negative_axes[i] = true;
        if (rel_pos[i] > 0) proc.positive_axes[i] = true;
      }
    }

    proc.full_grid = full_grid.dims;

    // TODO : Calculate this outside of this function
    proc.cell_dimension = nb_threads % 4 == 0
      ? full_grid.dims / PointInt{2, 2, nb_threads / 4}
      : full_grid.dims / PointInt{2, 1, nb_threads / 2};

    grid.size = full_grid.size;

    grid.pos = full_grid.pos
    + Vector{
      full_grid.size.x() * proc.cell_dimension.x * proc.pos.x - proc.negative_axes[0] * 2 * full_grid.size.x(),
      full_grid.size.y() * proc.cell_dimension.y * proc.pos.y - proc.negative_axes[1] * 2 * full_grid.size.y(),
      full_grid.size.z() * proc.cell_dimension.z * proc.pos.z - proc.negative_axes[2] * 2 * full_grid.size.z(),
    };

    grid.dims = proc.cell_dimension + PointInt{
      (proc.negative_axes[0] + proc.positive_axes[0]) * 2,
      (proc.negative_axes[1] + proc.positive_axes[1]) * 2,
      (proc.negative_axes[2] + proc.positive_axes[2]) * 2
    };

    grid.dims_id_start = {
      proc.pos.x * proc.cell_dimension.x - proc.negative_axes[0] * 2,
      proc.pos.y * proc.cell_dimension.y - proc.negative_axes[1] * 2,
      proc.pos.z * proc.cell_dimension.z - proc.negative_axes[2] * 2,
    };

    return grid;
  }

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

  template <typename HexData>
  void two_refinement_algorithm(HexData& hdata, MarkingFunction& cellIdentifier, int nb_levels, int thread_id = 0){
    static_assert(std::is_same_v<HexData, HexMeshingData> or std::is_same_v<HexData, ProcessData>);

    LCC& lcc = hdata.lcc;

    hdata.debug = lcc.get_new_mark();
    hdata.debug2 = lcc.get_new_mark();
    hdata.identified_mark = lcc.get_new_mark();
    hdata.template_mark = lcc.get_new_mark();
    hdata.propagation_face_mark = lcc.get_new_mark();
    l_debug_mark = hdata.debug;
    l_debug_mark_2 = hdata.debug2;
    l_thread_id = thread_id;

    if constexpr (std::is_same_v<HexData, ProcessData>){
      hdata.three_template_node_mark = lcc.get_new_mark();
      hdata.reset_temp_vertex_ids();
    }

    generate_grid(hdata);

    // Levels of refinement
    for (int r = 0; r < nb_levels; r++){
      if (r == 0) initial_setup(hdata, cellIdentifier);
      else setup_next_level(hdata, cellIdentifier);

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

        thread_communicate_cells_id_and_3t(hdata, rdata);

        lcc.unmark_all(hdata.identified_mark);
        lcc.unmark_all(hdata.template_mark);
      }


      lcc.unmark_all(hdata.propagation_face_mark);
    }
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

  void trim_excedent_volumes(LCC& lcc, TrimmingFunction func){
    auto volumes = lcc.one_dart_per_cell<3>();
    for (auto it = volumes.begin(); it != volumes.end(); it++){
      if (func(lcc, it)) continue;
      lcc.contract_cell<3>(it);
    }
  }
}

namespace CGAL::HexRefinement {

  void render_two_refinement_result(const LCC& lcc, Tree& aabb, bool trim = true, const char* title = "TwoRefinement Result"){
    LCCSceneOptions<LCC> gso;

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

  void debug_render(TwoRefinement::HexMeshingData& hdata){
    LCCSceneOptions<LCC> gso;
    LCC& lcc = hdata.lcc;

    AreaIDMap<CGAL::IO::Color> colors;
    colors[{0,0,0}] = blue();
    int i = 0;

    gso.draw_face = [](const LCC& lcc, LCC::Dart_const_handle dart){
      return true;
      auto a = lcc.attribute<2>(dart);

      return a != nullptr && a->info().plane[1];
    };
    gso.face_color = [&](const LCC& lcc, LCC::Dart_const_handle dart){
      auto a = lcc.attribute<2>(dart);
      auto b = lcc.attribute<3>(dart);

      // return TwoRefinement::is_volume_numberable(hdata, dart) != TwoRefinement::NumberableCell::None ? red(): blue();

      // return b->info().type >= VolumeType::REFINEMENT ? red() : blue();


      // if (b == nullptr) return black();
      // if (b->info().type > VolumeType::NONE) return red();
      // return blue();

      // auto n = TwoRefinement::is_face_numberable(hdata, dart);
      // if (n != TwoRefinement::NumberableCell::None) {
      //   if (n == TwoRefinement::NumberableCell::True) return purple();
      //   if (n == TwoRefinement::NumberableCell::False) return yellow();
      // };

      // return blue();
      // if (lcc.is_whole_cell_marked<2>(dart, hdata.debug)) return white();
      // if (lcc.is_whole_cell_marked<2>(dart, hdata.debug2)) return black();
      // if (a != nullptr && a->info().template_id == 3 && lcc.is_whole_cell_marked<2>(dart, hdata.debug2)) return purple();
      // if (a != nullptr && a->info().template_id == 3) return red();
      // if (a == nullptr) return blue();
      // if (a->info().plane[0]) return yellow();
      // if (a->info().plane[1]) return green();
      // if (a != nullptr && a->info().plane[hdata.d]) return purple();
      // return blue();

      colors[{0,0,0}] = blue();

      if (colors.count(b->info().area_id) == 0){
        CGAL::Random random((i+=3));
        colors[b->info().area_id] = CGAL::get_random_color(random);
      }

      return colors[b->info().area_id];
    };
    gso.colored_face = [](const LCC& lcc, LCC::Dart_const_handle dart){ return true; };

    gso.draw_volume = [](const LCC& lcc, LCC::Dart_const_handle dart){
      return true;
      return lcc.attribute<3>(dart)->info().owned;
    };

    gso.colored_vertex = [](const LCC& lcc, LCC::Dart_const_handle dart){
      return true;
    };
    gso.draw_vertex = [](const LCC& lcc, LCC::Dart_const_handle dart){
      return true;
    };
    gso.vertex_color = [&](const LCC& lcc, LCC::Dart_const_handle dart){
      return lcc.attribute<0>(dart)->id == DartInfo::VertexAttr::max_id ? red() : blue();
      // return lcc.is_whole_cell_marked<0>(dart, hdata.template_mark) ? red() : black();
      auto a = lcc.is_whole_cell_marked<0>(dart, hdata.debug);
      auto b = lcc.is_whole_cell_marked<0>(dart, hdata.debug2);
      return a && b ? purple() : a ? red() : b ? green() : black();
      // auto a = lcc.attribute<0>(dart);
      // return a->id != a->max_id ? red() : blue();
    };

    // gso.colored_edge = [](const LCC& lcc, LCC::Dart_const_handle dart){
    //   return true;
    // };
    // gso.draw_edge = [](const LCC& lcc, LCC::Dart_const_handle dart){
    //   return true;
    // };
    // gso.edge_color = [&](const LCC& lcc, LCC::Dart_const_handle dart){
    //   // auto n = TwoRefinement::is_edge_numberable(hdata, dart);
    //   // if (n == TwoRefinement::NumberableCell::True)
    //   //   return red();
    //   // if (n == TwoRefinement::NumberableCell::False)
    //   //   return blue();

    //   // return black();

    //   return lcc.is_marked(dart, hdata.debug)
    //   ? red()
    //   : lcc.is_marked(dart, hdata.debug2)
    //     ? orange()
    //     : black();
    // };

    CGAL::Graphics_scene buffer;
    add_to_graphics_scene(lcc, buffer, gso);
    CGAL::draw_graphics_scene(buffer);
  }

  LCC two_refinement(
      TwoRefinement::Grid grid,
      TwoRefinement::MarkingFunction cellIdentifier,
      TwoRefinement::TrimmingFunction trimmingFunction,
      int nb_levels = 1)
  {
    using namespace TwoRefinement;

    HexMeshingData hdata;
    ExternalRessources res;

    load_patterns(res.regular_templates, res.partial_templates);
    hdata.init(&res, grid);

    two_refinement_algorithm(hdata, cellIdentifier, nb_levels);
    // trim_excedent_volumes(hdata.lcc, trimmingFunction);

    return hdata.lcc;
  }

  LCC two_refinement_mt(
    TwoRefinement::Grid grid,
    TwoRefinement::MarkingFunction cellIdentifier,
    int nb_threads,
    int nb_levels = 1)
  {
    using namespace TwoRefinement;
    nb_threads = 4;

    std::vector<ProcessData> proc_datas(nb_threads);
    std::vector<ProdCons<ThreadMsg>> streams;
    std::vector<std::thread> threads(nb_threads);

    ExternalRessources res;
    load_patterns(res.regular_templates, res.partial_templates);

    // TODO : Make the N dimension on the largest axis
    // TODO : Verify if the layout is possible, else try an other layout.
    // If no layout matches : cerr + return false
    if (nb_threads % 4 == 0)
      create_threads_2x2xN(proc_datas, streams, nb_threads);
    else if (nb_threads % 2 == 0)
      create_threads_2x1xN(proc_datas, streams, nb_threads);
    else {
      // not yet done
      // create_threads_1x1xN(proc_datas, streams, nb_threads);
    };

    constexpr size_t max_size_t = std::numeric_limits<size_t>::max();
    size_t max_vertices = (grid.dims.x + 1) * (grid.dims.y + 1) * (grid.dims.z + 1);

    for (int i = 0; i < proc_datas.size(); i++){
      auto& proc = proc_datas[i];
      proc.init(&res, calculate_thread_grid(proc, grid, i, nb_threads));

      size_t domain_length = (max_size_t - max_vertices) / nb_threads;

      proc.v_domain_current = max_vertices + i       * domain_length;
      proc.v_domain_end     = max_vertices + (i + 1) * domain_length;
    }

    proc_datas[proc_datas.size() - 1].v_domain_end = max_size_t;

    // two_refinement_algorithm<ProcessData>(proc_datas[1], cellIdentifier, nb_levels);

    for (int i = 0; i < nb_threads; i++){
      threads[i] = std::thread(two_refinement_algorithm<ProcessData>, std::ref(proc_datas[i]), std::ref(cellIdentifier), nb_levels, i);
    }

    while (true){
    int i = debug_stream.waitNextItem();
    debug_render(proc_datas[i]);
    debug_stream2.push(1);
    }

    for (int i = 0; i < nb_threads; i++){
      threads[i].join();
    }

    // LCC& combined = proc_datas[0].lcc;
    // for (int i = 1; i < proc_datas.size(); i++){
    //   combined.copy(proc_datas[i].lcc);
    // }

    // debug_render(proc_datas[0]);

    for (auto& proc : proc_datas){
      debug_render(proc);
    }

    // ProcessData& proc = proc_datas[2];
    // int i = 0;
    // for (auto& proc : proc_datas){
      // two_refinement_algorithm(proc, cellIdentifier, 1);
      // i++;

      // debug3 = proc.lcc.get_new_mark();
      // LCC& lcc = proc.lcc;
      // for (auto& plane : proc.owned_ghost_areas[0]){
      //   for (auto& [constraint, dartvec] : plane){
      //     for (auto dart : dartvec){
      //       // lcc.mark(dart, debug3);
      //       // lcc.mark(lcc.other_extremity(dart), debug3);
      //       if(lcc.is_whole_cell_marked<1>(dart, debug3)) continue;
      //       lcc.mark_cell<1>(dart, debug3);
      //     }
      //   }
      // }
      // debug_render(proc);
    // }
    return {};
  }
}
