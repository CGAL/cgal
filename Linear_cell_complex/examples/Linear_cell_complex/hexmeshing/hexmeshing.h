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
#include <algorithm>
#include <bitset>
#include <boost/container_hash/hash_fwd.hpp>
#include <exception>
#include <functional>
#include "utils.h"
#include <CGAL/Simple_Cartesian.h>

#include "../query_replace/cmap_query_replace.h"
#include "../query_replace/init_to_preserve_for_query_replace.h"

#include <boost/container/static_vector.hpp>
#include <boost/range/join.hpp>
#include <limits>
#include <mutex>
#include <qnamespace.h>
#include <qpointingdevice.h>
#include <thread>
#include <unordered_map>
#include <variant>
#include <vector>
#include <winnt.h>
#include <CGAL/Simple_cartesian.h>

template <typename T>
struct GenericPoint {
  T x, y, z;
  GenericPoint(): x(0), y(0), z(0) {}
  GenericPoint(T x, T y, T z): x(x), y(y), z(z) {}
  GenericPoint(std::array<int,3> l): x(l[0]), y(l[1]), z(l[2]) {}
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
    typedef std::tuple<CGAL::Cell_attribute_with_point<Storage>, void, FaceAttr, VolumeAttr> Attributes;
  };
};


typedef CGAL::Linear_cell_complex_traits<3,Kernel> LCCTraits;
typedef CGAL::Linear_cell_complex_for_combinatorial_map<3,3, LCCTraits, LCCItems> LCC;
typedef typename LCC::Dart_handle Dart_handle;
typedef typename LCC::Vertex_attribute_handle Vertex_handle;
typedef typename LCC::size_type  size_type;

using DartInfo = LCCItems::Dart_wrapper<LCC::Storage>;

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
  enum PlaneNormal { X, Y, Z, NONE = -1};

  struct Grid {
    Point pos;
    Point size;
    PointInt dims;

    // For each plane normals, the starting index of those planes
    PointInt dims_id_start;

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
  using MarkingFunction = std::function<bool(LCC&, Polyhedron&, Dart_handle)>;
  using SetupFunction = std::function<Grid(Polyhedron&)>;

  struct ExternalRessources {
    Pattern_substituer<LCC> regular_templates, partial_templates;
    Polyhedron surface;
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

  struct MarkedCellsNode {
    static constexpr size_t none_child() { return std::numeric_limits<size_t>::max(); };

    std::array<size_t, 4> childs {none_child(), none_child(), none_child(), none_child()};
    std::array<bool, 4> marked_nodes = {false, false, false, false};
  };


  struct CellPositionNode {
    static constexpr size_t none_child() { return std::numeric_limits<size_t>::max(); };
    static constexpr char none_cell_id() { return -1; };

    std::array<size_t, 4> childs {none_child(), none_child(), none_child(), none_child()};
    char cell0_id = none_cell_id();
    Point position;
  };

  using MarkedCellsTree = std::vector<MarkedCellsNode>;
  using CellsPositionTree = std::vector<CellPositionNode>;

  using MCT_CC = std::vector<MarkedCellsTree>;
  using CPT_CC = std::vector<CellsPositionTree>;

  // Mapping from plane_id => Tree
  using MCT_Set = std::unordered_map<size_t, MCT_CC>;
  using CPT_Set = std::unordered_map<size_t, CPT_CC>;

  using ThreadMsg = std::variant<
    std::monostate,
    MCT_Set,
    CPT_Set
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

  using IOMsgStream = IOThreadStream<ThreadMsg>;

  using GhostCC = std::vector<Dart_handle>;
  using AreaIDToGhostCC = AreaIDMap<GhostCC>;
  using GhostAreaSet = std::vector<AreaIDToGhostCC>;

  struct ProcessData : public HexMeshingData{
    AreaIDMap<IOMsgStream> neighboring_threads;
    std::array<GhostAreaSet, 3> owned_ghost_areas, unowned_ghost_areas;
    std::array<bool, 3> positive_axes = {false, false, false};
    std::array<bool, 3> negative_axes = {false, false, false};;
    PointInt pos;

    size_type three_template_node_mark;
    std::vector<Dart_handle> marked_three_template;
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

    p1.neighboring_threads[rel_pos] = IOMsgStream(&streams[nextStream], &streams[nextStream+1]);
    p2.neighboring_threads[-rel_pos] = IOMsgStream(&streams[nextStream+1], &streams[nextStream]);

    assert(nextStream < streams.size());
    nextStream += 2;
  }

  // Creates threads data, in multiples of 4
  void create_threads_2x2xN(std::vector<ProcessData>& procs, std::vector<ProdCons<ThreadMsg>>& streams, int nb_threads){
    assert(nb_threads % 4 == 0);

    streams = std::vector<ProdCons<ThreadMsg>>(nb_threads * 6 * 2 + (nb_threads - 1) * 14 * 2);

    int nextStream = 0;

    for (int i = 0; i < nb_threads; i++){
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
      procs[i].pos = grid_cell_from_index_2x1xN(i);
      procs[i].grid.dims /= {2, 1, nb_threads / 2 };
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

  boost::container::static_vector<Dart_handle, 4> incident_faces_to_0_cell_on_plane(LCC& lcc, RefinementData &rdata, Dart_handle edge){
    boost::container::static_vector<Dart_handle, 4> arr;

    // Add initial face
    arr.push_back(edge);
    Dart_handle d2 = lcc.beta(edge, 0);

    // Left neighbour
    Dart_handle adjacent_face = adjacent_face_on_plane(lcc, rdata.iteration, d2);
    bool left_neigh = adjacent_face != lcc.null_dart_descriptor;
    if (left_neigh) arr.push_back(adjacent_face);

    //Right neighbour
    adjacent_face = adjacent_face_on_plane(lcc, rdata.iteration, edge);
    bool right_neigh = adjacent_face != lcc.null_dart_descriptor;
    if (right_neigh) arr.push_back(adjacent_face);

    // Diagonal neighbour to 'edge' face
    if (left_neigh && right_neigh){
      // Adjacent face might change the orientation of the face (because its choosing the beta3 face instead)
      int f = lcc.belong_to_same_cell<0>(arr[1], edge) ? 0 : 1;
      Dart_handle d3 = adjacent_face_on_plane(lcc, rdata.iteration, lcc.beta(arr[1], f));
      assert(d3 != lcc.null_dart_descriptor);
      arr.push_back(d3);

      #ifndef NDEBUG
        // Assert that we fall back on the same edge after doing a final adjacent face
        f = lcc.belong_to_same_cell<0>(d3, edge) ? 0 : 1;
        Dart_handle d4 = adjacent_face_on_plane(lcc, rdata.iteration, lcc.beta(d3, f));
        f = lcc.belong_to_same_cell<0>(d4, edge) ? 0 : 1;
        assert(d4 != lcc.null_dart_descriptor);
        d4 = lcc.beta(d4, f);
        assert(lcc.belong_to_same_cell<1>(arr[2], d4));
      #endif
    }

    return arr;
  }

  void clean_up_3_template(HexMeshingData &hdata, const Dart_handle &origin_dart, const Dart_handle &upper_edge, const Dart_handle lower_edge, Dart_handle &face1, Dart_handle &face2)
  {
    LCC& lcc = hdata.lcc;

    // Ne pas maintenir des refs sur les attributs, ils vont disparaitre
    DartInfo::FaceAttrValue face1_attr, face2_attr;
    if (lcc.attribute<2>(face1) != nullptr) face1_attr = lcc.attribute<2>(face1)->info();
    if (lcc.attribute<2>(face2) != nullptr) face2_attr = lcc.attribute<2>(face2)->info();


    Dart_handle lower_mid_1 = lcc.insert_barycenter_in_cell<1>(lower_edge);
    Dart_handle lower_mid_2 = lcc.beta(lower_mid_1, 2, 1);

    // Contract the lower edge asweel, this leads to two faces being partially contracted.
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

    if (face1 != lcc.null_dart_descriptor && face2 != lcc.null_dart_descriptor)
      lcc.sew<3>(face1, face2);

    // Requires at least one face for merging/adding 2-attributes
    std::bitset<3> merged_planes = face1_attr.plane | face2_attr.plane;
    auto merge_face = face1 != lcc.null_dart_descriptor ? face1 : face2;

    if (merge_face == lcc.null_dart_descriptor) return;

    // Merge the two previous face attributes bitsets
    auto& merge_face_attr = get_or_create_attr<2>(lcc, merge_face)->info();
    merge_face_attr = DartInfo::FaceAttrValue();
    merge_face_attr.plane = merged_planes;
  }

  template <unsigned int i, typename... DartArray>
  void assert_dart_attr_are_unique(LCC& lcc, DartArray... array){
    #ifndef NDEBUG
      std::unordered_map<void*, int> attributes;

      ([&]{
        int index = 0;
        for (auto dart : array){
          auto attr = lcc.attribute<i>(dart);
          assert(attr != nullptr);
          auto ptr = &attr->info();
          assert(attributes.count(ptr) == 0);
          attributes[ptr] = index;
          index++;
        }
      }(),...);
    #endif
  }

  void assert_faces_of_plane_valid(HexMeshingData& hdata, RefinementData& rdata){
    #ifndef NDEBUG
      LCC& lcc = hdata.lcc;
      assert_dart_attr_are_unique<2>(hdata.lcc, rdata.faces_of_plane);

      for (Dart_handle face : rdata.faces_of_plane){
        auto& attr = lcc.attribute<2>(face)->info();
        auto nodes = lcc.darts_of_cell<2, 0>(face);
        int marked = 0;
        for (auto it = nodes.begin(), end = nodes.end(); it != end; it++){
          if (lcc.is_marked(it, hdata.template_mark)) marked++;
        }

        assert( attr.plane[rdata.iteration]  == true );
        assert( attr.template_id == marked );
      }
    #endif
  }

  void assert_all_faces_are_quadrilateral(LCC &lcc, HexMeshingData& hdata){
    #ifndef NDEBUG
      auto iterator = lcc.one_dart_per_cell<2>();
      for (auto it = iterator.begin(); it != iterator.end(); it++){
        auto edges = lcc.darts_of_cell<2, 0>(it);
        assert(edges.size() == 4);
      }
    #endif
  }

  void assert_all_volumes_are_hexes(LCC &lcc){
    #ifndef NDEBUG
      auto iterator = lcc.one_dart_per_cell<3>();
      for (auto it = iterator.begin(); it != iterator.end(); it++){
        auto edges = lcc.darts_of_cell<3, 0>(it);
        assert(edges.size() == (4 * 6));
      }
    #endif
  }

  template <typename HexData>
  void refine_3_template(HexData &hdata, Dart_handle marked_face)
  {
    // TODO Might be written better
    LCC& lcc = hdata.lcc;

    Dart_handle origin_dart = find_3_template_origin(lcc, marked_face, hdata.template_mark);
    Dart_handle vol2_origin_dart = lcc.beta(origin_dart, 3);

    Dart_handle upper_d1 = origin_dart;
    Dart_handle upper_d2 = lcc.beta(origin_dart, 1, 1);

    Dart_handle upper_edge = lcc.insert_cell_1_in_cell_2(upper_d1, upper_d2);
    Dart_handle vol2_upper_edge = lcc.beta(upper_edge, 3);

    if constexpr (std::is_same_v<HexData, ProcessData>){
      // TODO : Instead of adding in a seperate vec, calculate origin dart before inserting
      // the volume in partial_templates

      // Weird, normally we don't process twice the same face
      if (lcc.is_whole_cell_unmarked<0>(origin_dart, hdata.three_template_node_mark)) {
        lcc.mark_cell<0>(origin_dart, hdata.three_template_node_mark);
        hdata.marked_three_template.push_back(origin_dart);
      }
    }

    // Query replace with the partial 3-template, making it into two volumes
    size_type p = hdata.ext->partial_templates.query_replace_one_volume(lcc, origin_dart, hdata.template_mark);
    assert(p == 0);

    // Also replace the other connected volume that is 3 template
    p = hdata.ext->partial_templates.query_replace_one_volume(lcc, lcc.beta(origin_dart, 3), hdata.template_mark);
    assert(p == 0);

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

    // contract the shared edge between the two volumes
    lcc.contract_cell<1>(upper_mid_1);
    lcc.contract_cell<1>(upper_mid_2);

    clean_up_3_template(hdata, origin_dart, upper_edge, lower_edge, face1, face2);
    clean_up_3_template(hdata, vol2_origin_dart, vol2_upper_edge, vol2_lower_edge, vol2_face1, vol2_face2);
  }

  void refine_marked_faces(HexMeshingData& hdata, RefinementData& rdata){
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

      if (hdata.ext->regular_templates.query_replace_one_face(lcc, dart, hdata.template_mark) != SIZE_T_MAX) nbsub++;

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
  void refine_marked_hexes(HexData& hdata, RefinementData& rdata)
  {
    LCC& lcc = hdata.lcc;
    int nbsub = 0;
    int nb_3_tp = 0;

    nbsub = 0;
    for (auto& dart : rdata.volumes_to_refine)
    {
      size_type temp = hdata.ext->regular_templates.query_replace_one_volume(lcc, dart, hdata.template_mark);
      if (temp != SIZE_T_MAX) nbsub++;

      if (temp == 2) temp = 3;
      if (temp < 10) temp++;

      // assert(vol_attr.template_id >= 0 && vol_attr.template_id <= 4 && vol_attr.template_id != 3 || vol_attr.template_id == 8);
      // assert(vol_attr.template_id == temp || temp == SIZE_T_MAX && (vol_attr.template_id == 0 || vol_attr.template_id == 8));
    }

    // assert(nbsub == rdata.volumes_to_refine.size());

    // Refine remaining 3 patterns
    for (auto marked_face : rdata.partial_templates_to_refine){
      assert(lcc.attribute<2>(marked_face) != nullptr);
      assert(lcc.attribute<2>(marked_face)->info().template_id == 3);
      refine_3_template(hdata, marked_face);
      nbsub += 2;
    }

    std::cout << nbsub << " volumic substitution was made" << std::endl;
  }



  // TODO : All those functions are duplicated code, refactor ?
  // =============================================================================================================
  // =============================================================================================================
  MarkedCellsTree create_marked_cells_tree(ProcessData& proc, Dart_handle start, PlaneNormal plane, AreaId area_id){
    LCC& lcc = proc.lcc;

    MarkedCellsTree tree(1);
    tree.reserve(256);

    using FaceNode = std::pair<Dart_handle, size_t>;
    std::queue<FaceNode> to_explore;

    size_type face_mark = lcc.get_new_mark();
    size_type edge_mark = lcc.get_new_mark();

    to_explore.push({start, 0});
    lcc.mark_cell<2>(start, face_mark);

    // I changed the loop to only iterate on unique element, because we associate a new node element to a face prior
    // to iterate on it.
    while (!to_explore.empty()){
      auto [face, face_id] = to_explore.front();
      to_explore.pop();

      Dart_handle edge = face;
      for (int i = 0; i < 4; edge = lcc.beta(edge, 1), i++){
        if (lcc.is_whole_cell_marked<1>(edge, edge_mark)) continue;
        lcc.mark_cell<1>(edge, edge_mark);

        if (lcc.is_marked(edge, proc.template_mark))
          tree[face_id].marked_nodes[i] = true;

        Dart_handle other_face = adjacent_face_on_plane(lcc, plane, edge);

        // Verify if the face is not already marked
        if (other_face == lcc.null_dart_descriptor or lcc.is_whole_cell_marked<2>(other_face, face_mark))
          continue;

        lcc.mark_cell<2>(other_face, face_mark);
        auto front_vol = lcc.attribute<3>(other_face);
        auto back_vol = lcc.attribute<3>(other_face);

        bool same_area = front_vol != nullptr
          && front_vol->info().owned
          && are_areas_equivalent(front_vol->info().area_id, area_id, true)
        or back_vol != nullptr
          && back_vol->info().owned
          && are_areas_equivalent(back_vol->info().area_id, area_id, true);

        if (same_area){
          size_t new_child_id = tree.size();
          tree.push_back(MarkedCellsNode());
          to_explore.push({other_face, new_child_id});
          tree[face_id].childs[i] = new_child_id;
        }
      }

      assert(edge == face);
    }

    lcc.free_mark(face_mark);
    lcc.free_mark(edge_mark);

    assert(tree.size() > 1);

    return tree;
  }

  CellsPositionTree create_cells_position_tree(ProcessData& proc, Dart_handle start, PlaneNormal plane, AreaId area_id){
    LCC& lcc = proc.lcc;

    CellsPositionTree tree(1);
    tree.reserve(256);

    using FaceNode = std::pair<Dart_handle, size_t>;
    std::queue<FaceNode> to_explore;

    size_type face_mark = lcc.get_new_mark();
    size_type edge_mark = lcc.get_new_mark();

    to_explore.push({start, 0});
    lcc.mark_cell<2>(start, face_mark);

    // I changed the loop to only iterate on unique element, because we associate a new node element to a face prior
    // to iterate on it.
    while (!to_explore.empty()){
      auto [face, face_id] = to_explore.front();
      to_explore.pop();

      Dart_handle edge = face;
      bool found = false;
      for (int i = 0; i < 4; edge = lcc.beta(edge, 1), i++){
        if (lcc.is_whole_cell_marked<1>(edge, edge_mark)) continue;
        lcc.mark_cell<1>(edge, edge_mark);

        if (!found && lcc.is_marked(edge, proc.three_template_node_mark)){
          tree[face_id].cell0_id = i;
          tree[face_id].position = lcc.attribute<0>(edge)->point();
          found = true;
        }

        Dart_handle other_face = adjacent_face_on_plane(lcc, plane, edge);

        // Verify if the face is not already marked
        if (other_face == lcc.null_dart_descriptor or lcc.is_whole_cell_marked<2>(other_face, face_mark))
          continue;

        lcc.mark_cell<2>(other_face, face_mark);
        auto front_vol = lcc.attribute<3>(other_face);
        auto back_vol = lcc.attribute<3>(other_face);

        bool same_area = front_vol != nullptr
          && front_vol->info().owned
          && are_areas_equivalent(front_vol->info().area_id, area_id, true)
        or back_vol != nullptr
          && back_vol->info().owned
          && are_areas_equivalent(back_vol->info().area_id, area_id, true);

        if (same_area){
          size_t new_child_id = tree.size();
          tree.push_back(CellPositionNode());
          to_explore.push({other_face, new_child_id});
          tree[face_id].childs[i] = new_child_id;
        }
      }

      assert(edge == face);
    }

    lcc.free_mark(face_mark);
    lcc.free_mark(edge_mark);

    assert(tree.size() > 1);

    return tree;
  }

  // TODO renommer
  void mct_action(ProcessData& proc, RefinementData& rdata, Dart_handle node){
    LCC& lcc = proc.lcc;
    lcc.mark_cell<0>(node, proc.template_mark);
    rdata.marked_nodes.push_back(node);

    for (Dart_handle face : incident_faces_to_0_cell_on_plane(lcc, rdata, node)){
      auto& face_attr =  lcc.attribute<2>(face)->info();

      // If the face didn't have any template before, it will have one
      if (face_attr.template_id == 0)
        rdata.faces_of_plane.push_back(face);

      face_attr.template_id++;
      assert(face_attr.template_id <= 4);
    }
  }

  void parse_marked_cells_tree(ProcessData& proc, RefinementData& rdata, const AreaId& area_id, const MarkedCellsTree& tree, Dart_handle start){
    assert(tree.size() != 0);
    LCC& lcc = proc.lcc;
    using FaceNode = std::pair<Dart_handle, size_t>;

    std::queue<FaceNode> to_explore;
    to_explore.push({start, 0});

    while (!to_explore.empty()){
      auto [face, index] = to_explore.front();
      to_explore.pop();

      const MarkedCellsNode& node = tree[index];

      Dart_handle edge = face;
      for (int i = 0; i < 4; edge = lcc.beta(edge, 1), i++){
        size_t child = node.childs[i];
        bool marked = node.marked_nodes[i];

        if (marked && !lcc.is_whole_cell_marked<0>(edge, proc.template_mark)){
          mct_action(proc, rdata, edge);
        }

        if (child == node.none_child())
          continue;

        Dart_handle next_face = adjacent_face_on_plane(lcc, rdata.iteration, edge);
        auto back_vol_attr = lcc.attribute<3>(next_face);
        auto vol_attr = lcc.attribute<3>(next_face);

        bool is_valid_face = next_face != lcc.null_dart_descriptor && (
          vol_attr != nullptr
          && !vol_attr->info().owned
          && are_areas_equivalent(vol_attr->info().area_id, area_id, false)
          or
          back_vol_attr != nullptr
          && !back_vol_attr->info().owned
          && are_areas_equivalent(back_vol_attr->info().area_id, area_id, false)
        );

        assert(is_valid_face);

        to_explore.push({next_face, child});
      }

      assert(edge == face);
    }
  }

  void parse_cells_position_tree(ProcessData& proc, const AreaId& area_id, const CellsPositionTree& tree, Dart_handle start, PlaneNormal plane){
    assert(tree.size() != 0);
    LCC& lcc = proc.lcc;
    using FaceNode = std::pair<Dart_handle, size_t>;

    std::queue<FaceNode> to_explore;
    to_explore.push({start, 0});

    while (!to_explore.empty()){
      auto [face, index] = to_explore.front();
      to_explore.pop();

      const CellPositionNode& node = tree[index];

      auto edges = lcc.darts_of_cell<2,1>(face);

      Dart_handle edge = face;
      for (int i = 0; i < 4; edge = lcc.beta(edge, 1), i++){
        size_t child = node.childs[i];

        if (i == node.cell0_id)
          lcc.attribute<0>(edge)->point() = node.position;

        if (child == node.none_child())
          continue;

        Dart_handle next_face = adjacent_face_on_plane(lcc, plane, edge);
        auto back_vol_attr = lcc.attribute<3>(next_face);
        auto vol_attr = lcc.attribute<3>(next_face);

        bool is_valid_face = next_face != lcc.null_dart_descriptor && (
          vol_attr != nullptr
          && !vol_attr->info().owned
          && are_areas_equivalent(vol_attr->info().area_id, area_id, false)
          or
          back_vol_attr != nullptr
          && !back_vol_attr->info().owned
          && are_areas_equivalent(back_vol_attr->info().area_id, area_id, false)
        );

        assert(is_valid_face);

        to_explore.push({next_face, child});
      }

      assert(edge == face);
    }
  }

  template <typename HexData>
  void thread_communicate_marked_nodes(HexData&, RefinementData&) {}

  template <>
  void thread_communicate_marked_nodes(ProcessData& proc, RefinementData& rdata){
    LCC& lcc = proc.lcc;

    // First send node information to other threads
    // Computes the marked cell tree for a specific thread (identified by the AreaId)
    // Multiples trees are computed if the size of CC > 1
    for (auto& [area_id, stream] : proc.neighboring_threads){
      GhostAreaSet& areas = proc.owned_ghost_areas[rdata.iteration];
      MCT_Set marked_cells_set;

      for (int i = 1; i < areas.size(); i += 2){
        if (areas[i].count(area_id) == 0) continue;
        // Calculate the plane_id, because the plane index won't be the same on other threads
        int plane_id = proc.grid.dims_id_start[rdata.iteration] + i;

        GhostCC& ghost_cc = areas[i][area_id];
        MCT_CC& marked_cells_cc = marked_cells_set[plane_id];

        marked_cells_cc = MCT_CC(ghost_cc.size());

        for (int cc = 0; cc < ghost_cc.size(); cc++)
          marked_cells_cc[cc] = create_marked_cells_tree(proc, ghost_cc[cc], rdata.iteration, area_id);
      }

      stream << std::move(marked_cells_set);
    }

    // Receive other threads trees, execute the path and mark flagged nodes from the tree.
    for (auto& [area_id, stream] : proc.neighboring_threads){
      ThreadMsg msg;
      stream >> msg;

      if (!std::holds_alternative<MCT_Set>(msg)){
        // what to do ?
        return;
      }

      MCT_Set& mct_set = std::get<MCT_Set>(msg);

      for (auto& [plane_id, tree_cc] : mct_set){
        // Calculate the plane_index from the plane_id
        int plane_index = plane_id - proc.grid.dims_id_start[rdata.iteration];
        assert(plane_index >= 0 && plane_index <= proc.unowned_ghost_areas[rdata.iteration].size());

        AreaIDToGhostCC& ghosts = proc.unowned_ghost_areas[rdata.iteration][plane_index];

        assert(ghosts.count(area_id) != 0);
        auto& ghost_cc = ghosts[area_id];
        assert(ghost_cc.size() == tree_cc.size());

        for (int i = 0; i < tree_cc.size(); i++){
          parse_marked_cells_tree(proc, rdata, area_id, tree_cc[i], ghost_cc[i]);
        }
      }
    }
  }

  template <typename HexData>
  void thread_communicate_3_templates(HexData&, RefinementData&){}

  template <>
  void thread_communicate_3_templates(ProcessData& proc, RefinementData& rdata){
    LCC& lcc = proc.lcc;

    // First send node information to other threads
    // Computes the marked cell tree for a specific thread (identified by the AreaId)
    // Multiples trees are computed if the size of CC > 1
    for (auto& [area_id, stream] : proc.neighboring_threads){
      GhostAreaSet& areas = proc.owned_ghost_areas[rdata.iteration];
      CPT_Set cells_position_tree;

      for (int i = 1; i < areas.size(); i += 2){
        if (areas[i].count(area_id) == 0) continue;
        // Calculate the plane_id, because the plane index won't be the same on other threads
        int plane_id = proc.grid.dims_id_start[rdata.iteration] + i;

        GhostCC& ghost_cc = areas[i][area_id];
        CPT_CC& marked_cells_cc = cells_position_tree[plane_id];

        marked_cells_cc = CPT_CC(ghost_cc.size());

        for (int cc = 0; cc < ghost_cc.size(); cc++)
          marked_cells_cc[cc] = create_cells_position_tree(proc, ghost_cc[cc], rdata.iteration, area_id);
      }

      stream << std::move(cells_position_tree);
    }

    // Receive other threads trees, execute the path and mark flagged nodes from the tree.
    for (auto& [area_id, stream] : proc.neighboring_threads){
      ThreadMsg msg;
      stream >> msg;

      if (!std::holds_alternative<CPT_Set>(msg)){
        // Terminate thread
        exit(EXIT_FAILURE);
        return;
      }

      CPT_Set& cell_position_set = std::get<CPT_Set>(msg);

      for (auto& [plane_id, tree_cc] : cell_position_set){
        // Calculate the plane_index from the plane_id
        int plane_index = plane_id - proc.grid.dims_id_start[rdata.iteration];
        assert(plane_index >= 0 && plane_index <= proc.unowned_ghost_areas[rdata.iteration].size());

        AreaIDToGhostCC& ghosts = proc.unowned_ghost_areas[rdata.iteration][plane_index];

        assert(ghosts.count(area_id) != 0);
        auto& ghost_cc = ghosts[area_id];
        assert(ghost_cc.size() == tree_cc.size());

        for (int i = 0; i < tree_cc.size(); i++){
          parse_cells_position_tree(proc, area_id, tree_cc[i], ghost_cc[i], rdata.iteration);
        }
      }
    }
  }

  // =============================================================================================================
  // =============================================================================================================

  void __expand_0_cell_marking(LCC &lcc, RefinementData &rdata, std::queue<Dart_handle>& faces_to_check, Dart_handle &edge) {
    auto faces = incident_faces_to_0_cell_on_plane(lcc, rdata, edge);
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
      lcc.mark_cell<0>(edges_unmarked[0], hdata.debug2);
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

  void fix_impossible_cases(HexMeshingData &hdata, RefinementData &rdata){
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
  }


  /**
   * Mark 0-cells
   * Gather 2-cells adjacent to marked 0-cells
   * Gather 3-cells adjacents to marked 0-cells
   */
  template <typename HexData>
  void explore_face_of_plane(HexData& hdata, RefinementData& rdata, std::queue<Dart_handle>& queue,
        Dart_handle face, size_type explored_mark, size_type explored_edge, size_type explored_face) {
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
      is_markable = front_vol.owned && (back_vol_attr == nullptr or back_vol_attr->info().owned);
    }

    auto edges = lcc.darts_of_cell<2,1>(face);
    // Add neighboring faces
    for (auto dit = edges.begin(), dend = edges.end(); dit != dend; dit++){
      bool explored = lcc.is_marked(dit, explored_mark);
      bool edge_explored = lcc.is_whole_cell_marked<1>(dit, explored_edge);
      bool identified = lcc.is_marked(dit, hdata.identified_mark);

      if (is_markable){
        if (!explored){
          lcc.mark_cell<0>(dit, explored_mark);
        }

        if (!explored && identified){
          lcc.mark_cell<0>(dit, hdata.template_mark);
          rdata.marked_nodes.push_back(dit);
        }
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

  void propagate_face(HexMeshingData &hdata, RefinementData &rdata, const Dart_handle &face, size_type explored_node_mark, size_type explored_face_mark)
  {
    LCC& lcc = hdata.lcc;
    Dart_handle back_face = lcc.beta(face, 2, 1, 1, 2);

    mark_k_cells_of_i_cell<2, 0>(lcc, back_face, hdata.template_mark);
    assert(lcc.attribute<3>(face) != nullptr);

    auto &vol_attr = get_or_create_refinement_volume(lcc, face)->info();
    vol_attr.iteration = rdata.iteration;
    // /!\ ALSO: Mark add the hex attached to the back_face for refinement.
    // It is garanted that faces and volumes added will be unique in our array
    // Because that back_hex is only accessible within the template itself, and not
    // accessible from the plane.

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

      if (!lcc.is_marked(it, explored_node_mark)){
        rdata.marked_nodes.push_back(it);
        lcc.mark_cell<0>(it, explored_node_mark);
      }

      if (!lcc.is_whole_cell_marked<2>(top_face, explored_face_mark)){
        rdata.faces_to_refine.push_back(top_face);
        mark_face_unchecked(lcc, top_face, explored_face_mark);
      }
    }
  }

  void propagation_stage(HexMeshingData &hdata, RefinementData &rdata, size_type explored_node_mark, size_type explored_face_mark)
  {
    LCC& lcc = hdata.lcc;

    bool skip_propagation = rdata.iteration == 0;
    int propagated_count = 0;
    int marked_for_prop_count = 0;

    if (!skip_propagation) for (Dart_handle face : rdata.faces_of_plane){
      auto& face_attr = lcc.attribute<2>(face)->info();

      if (face_attr.template_id == 4){
        if (is_half_face_marked(lcc, face, hdata.propagation_face_mark)){
          propagate_face(hdata, rdata, face, explored_node_mark, explored_face_mark);
          propagated_count++;
        }

        Dart_handle other_face = lcc.beta(face, 3);
        if (!lcc.is_free<3>(face) && is_half_face_marked(lcc, other_face, hdata.propagation_face_mark)) {
          propagate_face(hdata, rdata, other_face, explored_node_mark, explored_face_mark);
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

    // TODO refactor?
    size_type explored_mark = lcc.get_new_mark();
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
        explore_face_of_plane(hdata, rdata, to_explore, front,
          explored_mark, explored_edge, explored_face);
      }
    }

    assert_faces_of_plane_valid(hdata, rdata);

    fix_impossible_cases(hdata, rdata);

    assert_faces_of_plane_valid(hdata, rdata);

    propagation_stage(hdata, rdata, explored_mark, explored_face);

    thread_communicate_marked_nodes(hdata, rdata);
    get_cells_to_refine_from_plane(hdata, rdata,  explored_face);
    get_cells_to_refine_from_additionnal_volumes(hdata, rdata, explored_face);

    lcc.free_mark(explored_mark);
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
      if (old_info.type == VolumeType::IDENTIFIED && cellIdentifier(lcc, hdata.ext->surface, it->dart()))
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

  void setup_ghost_areas(ProcessData& proc){
    LCC& lcc = proc.lcc;
    // TODO Le code peut etre mieux écrit que ça ?

    // Attributes a AreaID to each volumes and sets the ownership of an area
    // If the area is owned, areas id can overlap on each other, equality is sufficient if all non zero dimensions are equals
    // If the area is unowned, no overlaps happens, two id are must be strictly equal to belong in the same area
    auto mark_ghost_area = [&](std::vector<Dart_handle> starts, int& plane, bool positive, bool owned){
      plane_for_each_face(lcc, starts,
        [&](Dart_handle face){
          auto& front_vol = lcc.attribute<3>(face)->info();
          auto back_vol_attr = lcc.attribute<3>(lcc.beta(face, 3));

          // Return if we are stepping on unowned area, only when owned is set
          if (owned && (!front_vol.owned or back_vol_attr != nullptr && !back_vol_attr->info().owned))
            return;

          front_vol.area_id[plane] = positive ? 1 : -1;
          front_vol.owned = owned;

          if (back_vol_attr != nullptr){
            auto& back_vol = back_vol_attr->info();
            back_vol.area_id[plane] = positive ? 1 : -1;
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

    // Then, for each possible areas, find the starting dart.
    // Since the grid was just created, we can iterate like in a cartesian grid
    // Later on, this wont be the case, but we can rely on the parralel consistency of
    // creating planes connected components
    for (int p = 0; p < 3; p++){
      // TODO, looks very messy can it be better?

      PlaneSet& plane_set = proc.first_face_of_planes[p];
      GhostAreaSet unowned_ghost_set(plane_set.size());
      GhostAreaSet owned_ghost_set(plane_set.size());

      int own_start = proc.negative_axes[p] ? 2 : 0;
      int own_end = proc.negative_axes[p] ? plane_set.size() - 2 : 0;

      for (int s = 0; s < plane_set.size(); s++){
        PlaneCC& plane_cc = plane_set[s];

        // if (s >= own_start && s <= own_end){
        AreaIDToGhostCC& owned_ghost_cc = owned_ghost_set[s];
        for (auto& [area_id, _] : proc.neighboring_threads){
          Dart_handle area_pos = find_area_on_grid(proc, plane_cc[0], (PlaneNormal)p, area_id, true);
          if (area_pos != lcc.null_dart_descriptor)
            owned_ghost_cc[area_id] = {area_pos};
        }
        // }

        AreaIDToGhostCC& unowned_ghost_cc = unowned_ghost_set[s];
        for (auto& [area_id, _] : proc.neighboring_threads){
          Dart_handle area_pos = find_area_on_grid(proc, plane_cc[0], (PlaneNormal)p, area_id, false);
          if (area_pos != lcc.null_dart_descriptor)
            unowned_ghost_cc[area_id] = {area_pos};
        }
      }

      proc.owned_ghost_areas[p] = std::move(owned_ghost_set);
      proc.unowned_ghost_areas[p] = std::move(unowned_ghost_set);
    }
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
      if (cellIdentifier(lcc, hdata.ext->surface, dart))
        vol_attr.type = VolumeType::IDENTIFIED;
    }
  }

  // thread_initial_setup
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

    // Mark initial identified cells
    for (auto dart = volumes.begin(), end = volumes.end(); dart != end; dart++){
      auto& vol_attr = lcc.attribute<3>(dart)->info();

      // Mark those who are identified
      if (vol_attr.owned && cellIdentifier(lcc, proc.ext->surface, dart))
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

  void setup_next_level(ProcessData& hdata, MarkingFunction& cellIdentifier){}

  void mark_identified_cells_from_3_attrs(HexMeshingData& hdata) {
    LCC& lcc = hdata.lcc;

    auto& attributes = lcc.attributes<3>();

    for (auto it = attributes.begin(), end = attributes.end(); it != end; it++){
      if (it->info().type > VolumeType::NONE){
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

  Grid calculate_thread_grid(ProcessData& proc, const Grid& full_grid, int nb_threads){
    Grid grid;

    // Starting point can be displaced if neighboring cells
    // exists on -X, -Y, -Z planes

    for (auto& [rel_pos, _] : proc.neighboring_threads){
      for (int i = 0; i < 3; i++){
        if (rel_pos[i] < 0) proc.negative_axes[i] = true;
        if (rel_pos[i] > 0) proc.positive_axes[i] = true;
      }
    }

    // TODO : Calculate this outside of this function
    PointInt proc_cell_dimension = nb_threads % 4 == 0
      ? full_grid.dims / PointInt{2, 2, nb_threads / 4}
      : full_grid.dims / PointInt{2, 1, nb_threads / 2};

    grid.size = full_grid.size;

    grid.pos = full_grid.pos
    + Vector{
      full_grid.size.x() * proc_cell_dimension.x * proc.pos.x - proc.negative_axes[0] * 2 * full_grid.size.x(),
      full_grid.size.y() * proc_cell_dimension.y * proc.pos.y - proc.negative_axes[1] * 2 * full_grid.size.y(),
      full_grid.size.z() * proc_cell_dimension.z * proc.pos.z - proc.negative_axes[2] * 2 * full_grid.size.z(),
    };



    grid.dims = proc_cell_dimension + PointInt{
      (proc.negative_axes[0] + proc.positive_axes[0]) * 2,
      (proc.negative_axes[1] + proc.positive_axes[1]) * 2,
      (proc.negative_axes[2] + proc.positive_axes[2]) * 2
    };

    grid.dims_id_start = {
      proc.pos.x * proc_cell_dimension.x,
      proc.pos.y * proc_cell_dimension.y,
      proc.pos.z * proc_cell_dimension.z,
    };

    grid.dims_id_start -= {
      proc.negative_axes[0] ? 2 : 0,
      proc.negative_axes[1] ? 2 : 0,
      proc.negative_axes[2] ? 2 : 0
    };

    return grid;
  }

  void create_vertices_for_templates(HexMeshingData& hdata, RefinementData& rdata)
  {

    // 2 noeuds marqué l'un à coté de l'autre ne produit pas de sommet
    // 1 noeud marqué a coté d'un noeud non marqué produit un sommet

    // TODO A changer pour une itération sur les arretes ou pas selon comment l'algo va s'executer

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
      lcc.insert_barycenter_in_cell<1>(dart);
    }

    std::cout << "Vertices created: " << vertices_created << std::endl;

    lcc.free_mark(arrete_done);
  }

  template <typename HexData>
  void two_refinement_algorithm(HexData& hdata, MarkingFunction& cellIdentifier, int nb_levels){
    static_assert(std::is_same_v<HexData, HexMeshingData> or std::is_same_v<HexData, ProcessData>);

    LCC& lcc = hdata.lcc;

    hdata.debug = lcc.get_new_mark();
    hdata.debug2 = lcc.get_new_mark();
    hdata.identified_mark = lcc.get_new_mark();
    hdata.template_mark = lcc.get_new_mark();
    hdata.propagation_face_mark = lcc.get_new_mark();

    if constexpr (std::is_same_v<HexData, ProcessData>){
      hdata.three_template_node_mark = lcc.get_new_mark();
    }

    generate_grid(hdata);

    for (int i = 0; i < nb_levels; i++){
      if (i == 0) initial_setup(hdata, cellIdentifier);
      else setup_next_level(hdata, cellIdentifier);

      expand_identified_cells(hdata, i, nb_levels);

      for (int p = 0; p < 3; p++) {
        RefinementData rdata;

        mark_identified_cells_from_3_attrs(hdata);

        get_cells_to_refine(hdata, rdata, (PlaneNormal)p);

        assert_dart_attr_are_unique<3>(lcc, rdata.volumes_to_refine, rdata.partial_templates_to_refine);

        // Refinement stage
        create_vertices_for_templates(hdata, rdata);
        refine_marked_faces(hdata, rdata);
        refine_marked_hexes(hdata, rdata);

        assert_all_faces_are_quadrilateral(lcc, hdata);
        assert_all_volumes_are_hexes(lcc);

        // Empty if not multi threaded data
        thread_communicate_3_templates(hdata, rdata);
        return;


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

  void load_surface(const std::string& file, Polyhedron& out_surface) {
    std::ifstream off_file(file);
    CGAL_precondition_msg(off_file.good(), ("Input .off couldn't be read : " + file).c_str());

    Polyhedron surface;
    off_file>>out_surface;
  }

  void load_resources(const std::string& file, ExternalRessources& res){
    load_patterns(res.regular_templates, res.partial_templates);
    load_surface(file, res.surface);
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

  auto simple_voxelisation_setup(Tree& tree) {
    return [&](Polyhedron& poly){
      // Triangulate before AABB
      CGAL::Polygon_mesh_processing::triangulate_faces(poly);
      // Compute AABB tree
      tree.insert(faces(poly).first, faces(poly).second, poly);
      tree.accelerate_distance_queries();
      tree.bbox();
    };
  }

  auto mark_intersecting_volume_with_poly(Tree& tree) {
    return [&](LCC& lcc, Polyhedron& poly, Dart_handle dart){
      return is_intersect(lcc, dart, tree);
    };
  }
}

namespace CGAL::HexRefinement {

  void render_two_refinement_result(const LCC& lcc, Tree& aabb, const char* title = "TwoRefinement Result"){
    LCCSceneOptions<LCC> gso;

    gso.draw_volume = [&](const LCC& lcc, LCC::Dart_const_handle dart){
      return TwoRefinement::is_intersect(lcc, dart, aabb);
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
    };
    gso.face_color = [&](const LCC& lcc, LCC::Dart_const_handle dart){
      auto a = lcc.attribute<2>(dart);
      auto b = lcc.attribute<3>(dart);

      return lcc.is_marked(dart, hdata.debug) ? red() : blue();

      // if (a == nullptr) return blue();
      // if (a->info().plane[0]) return yellow();
      // if (a->info().plane[1]) return red();
      // if (a->info().plane[2]) return green();
      // return blue();

      // if (b == nullptr) return black();
      // if (b->info().type > VolumeType::NONE) return red();
      // return blue();

      // if (colors.count(b->info().area_id) == 0){
      //   CGAL::Random random((i+=3));
      //   colors[b->info().area_id] = CGAL::get_random_color(random);
      // }

      // return colors[b->info().area_id];

      // return lcc.is_whole_cell_marked<3>(dart, debug3)? red() : a != nullptr && a->info().type == VolumeType::IDENTIFIED ? green() : blue();
    };
    gso.colored_face = [](const LCC& lcc, LCC::Dart_const_handle dart){ return true; };

    gso.colored_vertex = [](const LCC& lcc, LCC::Dart_const_handle dart){
      return true;
    };
    gso.draw_vertex = [](const LCC& lcc, LCC::Dart_const_handle dart){
      return true;
    };
    gso.vertex_color = [&](const LCC& lcc, LCC::Dart_const_handle dart){
      return lcc.is_whole_cell_marked<0>(dart, hdata.debug) ? yellow()
      : lcc.is_whole_cell_marked<0>(dart, hdata.template_mark) ? red() : black();
    };

    // gso.colored_edge = [](const LCC& lcc, LCC::Dart_const_handle dart){
    //   return true;
    // };
    // gso.draw_edge = [](const LCC& lcc, LCC::Dart_const_handle dart){
    //   return true;
    // };
    // gso.edge_color = [&](const LCC& lcc, LCC::Dart_const_handle dart){
    //   return lcc.is_whole_cell_marked<1>(dart, debug3) ? red() : black();
    // };
    CGAL::Graphics_scene buffer;
    add_to_graphics_scene(lcc, buffer, gso);
    CGAL::draw_graphics_scene(buffer);
  }

  using R = std::tuple<LCC, Polyhedron>;

  // Pas convaincu que c'est la meilleure façon d'appeller la fonction ...
  // Les lambdas permettent de capturer des variables en dehors de ce qui est donné en argument
  // Si une fonction a besoins d'un AABB tree, dans markingfunction par exemple.
  // C'est à la fonction setup de capturer et générer l'aabb pour que markingfunction puisse l'utiliser.
  R two_refinement(
        const std::string& file,
        TwoRefinement::SetupFunction setupFunction,
        TwoRefinement::MarkingFunction cellIdentifier,
        int nb_levels = 1)
  {
    using namespace TwoRefinement;

    HexMeshingData hdata;
    ExternalRessources res;
    load_resources(file, res);

    hdata.init(&res, setupFunction(res.surface));

    two_refinement_algorithm(hdata, cellIdentifier, nb_levels);

    debug_render(hdata);

    return {hdata.lcc, res.surface};
  }

  R two_refinement_mt(
        const std::string& file,
        int nb_threads,
        TwoRefinement::SetupFunction setupFunction,
        TwoRefinement::MarkingFunction cellIdentifier,
        int nb_levels = 1)
  {
    using namespace TwoRefinement;
    nb_threads = 8;

    std::vector<ProcessData> proc_datas(nb_threads);
    std::vector<ProdCons<ThreadMsg>> streams;
    std::vector<std::thread> threads(nb_threads);

    ExternalRessources res;
    load_resources(file, res);

    Grid grid = setupFunction(res.surface);

    // TODO : Make the N dimension on the largest axis
    // TODO : Verify if the layout is possible, else try an other layout.
    // If no layout matches : cerr + return false
    if (nb_threads % 4 == 0)
      create_threads_2x2xN(proc_datas, streams, nb_threads);
    else if (nb_threads % 2 == 0)
      create_threads_2x1xN(proc_datas, streams, nb_threads);
    else {
      create_threads_1x1xN(proc_datas, streams, nb_threads);
    };

    for (ProcessData& proc : proc_datas){
      proc.init(&res, calculate_thread_grid(proc, grid, nb_threads));
    }

    // two_refinement_algorithm<ProcessData>(proc_datas[1], cellIdentifier, nb_levels);

    for (int i = 0; i < nb_threads; i++){
      threads[i] = std::thread(two_refinement_algorithm<ProcessData>, std::ref(proc_datas[i]), std::ref(cellIdentifier), nb_levels);
    }

    for (int i = 0; i < nb_threads; i++){
      threads[i].join();
    }

    LCC& combined = proc_datas[0].lcc;
    for (int i = 1; i < proc_datas.size(); i++){
      combined.copy(proc_datas[i].lcc);
    }

    debug_render(proc_datas[0]);

    // for (auto& proc : proc_datas){
    //   debug_render(proc);
    // }

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
