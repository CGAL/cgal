#pragma once

#include "hexmeshing_sequential.h"
#include <queue>
#include <mutex>
#include <thread>
#include <condition_variable>
#include <variant>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <cassert>


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

namespace CGAL::HexRefinement::TwoRefinement {
  enum class NumberableCell {None, False, True};

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

    std::vector<Dart_handle> owned_ghost_area, unowned_ghost_area;

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

    void fix_dart_storage(){
      HexMeshingData::fix_dart_storage();

      bool should_fix = false;
      using CCToHandle = std::unordered_map<size_t, std::reference_wrapper<Dart_handle>>;
      CCToHandle owned_cc_to_handle, unowned_cc_to_handle;

      const auto gather_fixable_volumes = [&](std::vector<Dart_handle>& ghost_area, CCToHandle& cc_to_handle, bool owned){
        for (int cc_id = 0; cc_id < ghost_area.size(); cc_id++){
          Dart_handle& start = ghost_area[cc_id];
          auto attr = lcc.attribute<3>(start);
          if (!lcc.is_dart_used(start) or attr == nullptr
            or !attr->is_valid()
            or attr->info().owned != owned
            or attr->info().cc_id != cc_id)
          {
            should_fix = true;
            cc_to_handle.emplace(cc_id, start);
          }
        }
      };

      const auto has_invalid_volume = [&](std::vector<Dart_handle>& ghost_area, CCToHandle& cc_to_handle, bool owned){
        for (int cc_id = 0; cc_id < ghost_area.size(); cc_id++){
          Dart_handle& start = ghost_area[cc_id];
          auto attr = lcc.attribute<3>(start);
          if (!lcc.is_dart_used(start) or attr == nullptr
            or !attr->is_valid()
            or attr->info().owned != owned
            or attr->info().cc_id != cc_id)
          {
            return true;
          }
        }

        return false;
      };

      gather_fixable_volumes(owned_ghost_area, owned_cc_to_handle, true);
      gather_fixable_volumes(unowned_ghost_area, unowned_cc_to_handle, false);

      if (!should_fix) return;

      auto& attributes = lcc.attributes<3>();
      for (auto it = attributes.begin(); it != attributes.end(); it++){
        auto& info = it->info();
        auto& cc_to_handle = info.owned
        ? owned_cc_to_handle
        : unowned_cc_to_handle;

        auto dart_it = cc_to_handle.find(info.cc_id);
        if (dart_it == cc_to_handle.end()) continue;

        auto& dart_ref = dart_it->second;
        dart_ref.get() = it->dart();

        cc_to_handle.erase(dart_it);

        if (owned_cc_to_handle.empty() && unowned_cc_to_handle.empty())
          break;
      }

      const auto all_valid = [&](){
        CGAL_postcondition_msg(owned_cc_to_handle.empty() && unowned_cc_to_handle.empty(), "");
        bool invalid = has_invalid_volume(owned_ghost_area, owned_cc_to_handle, true)
                    or has_invalid_volume(unowned_ghost_area, unowned_cc_to_handle, false);
        CGAL_postcondition_msg(invalid, "");
      };

      CGAL_postcondition_code(all_valid());
    }
  };

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

  template <typename VolumeOp, typename VolumeSelector>
  inline void for_each_hex(LCC& lcc, std::vector<Dart_handle> starts,
                            const VolumeOp&& volume_operation,
                            const VolumeSelector&& volume_selector)
  {
    size_type volume_mark = lcc.get_new_mark();
    size_type face_mark = lcc.get_new_mark();

    std::queue<Dart_handle> to_explore;
    for (Dart_handle start : starts){
      to_explore.push(start);
      lcc.mark_cell<3>(start, volume_mark);
    }

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

  template <>
  void thread_communicate_marked_nodes(ProcessData& proc, RefinementData& rdata, size_type edge_mark){
    LCC& lcc = proc.lcc;

    size_type node_mark = lcc.get_new_mark();

    ThreadMarkedCellsMsg all_sent_cells;
    ThreadMarkedCellsPtr marked_cells; // TODO: remove all_sent_cells if we avoid iterating the plane for getting nodes, see todo below
    std::vector<ThreadMarkedCellsPtr> others_msg;
    others_msg.reserve(proc.neighboring_threads.size());

    auto retrieve_nodes_in_face = [&](Dart_handle face){
      auto nodes =  lcc.darts_of_cell<2,0>(face);
      for (auto it = nodes.begin(); it != nodes.end(); it++){
        size_t id = lcc.attribute<0>(it)->id;

        if (!lcc.is_marked(it, proc.template_mark) or all_sent_cells.count(id)) continue;
        assert(id != DartInfo::VertexAttr::max_id);
        marked_cells->insert(id);
        all_sent_cells.insert(id);
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

    AreaIDMap<ProcessData::Neighbor> threads = proc.neighboring_threads;

    while (true){
      // Create a new vector, while other threads might be working on the older one
      marked_cells = std::make_shared<ThreadMarkedCellsMsg>();

      // Retrieves nodes to send
      // TODO (?): Avoid iterating again, by directly adding any vertices that is connected to a owned ghost volume when created
      auto& plane_set = proc.first_face_of_planes[rdata.iteration];
      for (int i = 1; i < plane_set.size(); i += 2){
        plane_for_each_face(lcc, plane_set[i],
          [&](Dart_handle face, auto& edges){
            auto back_vol = lcc.attribute<3>(lcc.beta(face, 3));
            auto front_vol = lcc.attribute<3>(face);

            const auto owned_sharable = [&](LCC::Attribute_handle<3>::type v){
              return v->info().owned && v->info().area_id != AreaId{0,0,0};
            };

            if (owned_sharable(front_vol) or (back_vol != nullptr && owned_sharable(back_vol))){
              retrieve_nodes_in_face(face);
            }
          },
          [&](Dart_handle edge){
            return adjacent_face_on_plane(lcc, rdata.iteration, edge);
          }
        );
      }

      // Send them
      for (auto& [area, thread] : proc.neighboring_threads){
        // Use a shared pointer to avoid copy
        thread.stream << marked_cells;
      }

      // Receive all newly marked nodes
      int empty_count = 0;
      for (auto& [area, thread] : proc.neighboring_threads){
        ThreadMsg msg;
        thread.stream >> msg;

        if (std::holds_alternative<ThreadMarkedCellsPtr>(msg)){
          auto& others_cells = std::get<ThreadMarkedCellsPtr>(msg);
          if (others_cells->empty()) empty_count++;
          others_msg.push_back(std::move(others_cells));
        }
        else {
          // TODO what now?
          assert(false);
          exit(0);
        }
      }

      for (int i = 1; i < plane_set.size(); i += 2){
        plane_for_each_face(lcc, plane_set[i],
        [&](Dart_handle face, auto& edges){
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

      // EXIT: if no more changes are received
      if (empty_count >= proc.neighboring_threads.size())
        break;

      // Clear all received nodes
      others_msg.clear();
      empty_count = 0;

      size_t nb_changes = fix_impossible_cases(proc, rdata);
      assert_faces_of_plane_valid(proc, rdata);

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

    // Explore all assignable volumes
    for_each_hex(lcc, proc.owned_ghost_area,
      [&](Dart_handle vol, auto faces){
        for (Dart_handle node : nodes_of_hex(lcc, vol)){
          node_operation(node);
        }
      },
      [&](Dart_handle vol) -> bool{
        auto attr = lcc.attribute<3>(vol);
        return attr != nullptr && attr->info().area_id != AreaId{0,0,0};
      });


    thread_assert_no_temporary_ids(proc);

    // lcc.unmark_all(proc.three_template_node_mark);
  }

  template <>
  void thread_remove_ghosts(ProcessData& proc){
    LCC& lcc = proc.lcc;
    std::vector<Dart_handle> vol_to_delete;

    Dart_handle start;
    auto& attributes = lcc.attributes<3>();
    for (auto it = attributes.begin(); it != attributes.end(); it++){
      if (!it->info().owned) {
        vol_to_delete.push_back(it->dart());
      }
    }

    for (Dart_handle dart : vol_to_delete){
      lcc.remove_cell<3>(dart);
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
    proc.unowned_ghost_area = {};
    proc.owned_ghost_area = {};

    for (auto& handle : volumes){
      if (!found_owned && handle.info().owned && handle.info().area_id != AreaId(0,0,0)){
        found_owned = true;
        proc.owned_ghost_area = { handle.dart() };
      }

      if (!found_unowned && !handle.info().owned){
        found_unowned = true;
        proc.unowned_ghost_area = { handle.dart() };
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
        [&](Dart_handle face, auto& edges){
          auto& front_vol = lcc.attribute<3>(face)->info();
          auto back_vol_attr = lcc.attribute<3>(lcc.beta(face, 3));

          // Return if we are stepping on unowned area, only when owned is set
          if (owned && (!front_vol.owned or back_vol_attr != nullptr && !back_vol_attr->info().owned))
            return;

          front_vol.area_id[plane] = positive_axis ? 1 : -1;
          front_vol.owned = owned;
          front_vol.cc_id = 0;

          if (back_vol_attr != nullptr){
            auto& back_vol = back_vol_attr->info();
            back_vol.area_id[plane] = positive_axis ? 1 : -1;
            back_vol.owned = owned;
            back_vol.cc_id = 0;
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

  void initial_setup(ProcessData& proc, MarkingFunction& cellIdentifier){
    initial_setup(static_cast<HexMeshingData&>(proc), cellIdentifier);
    setup_ghost_areas(proc);
    setup_vertex_ids(proc);
  }

  void setup_ghost_cell_cc(ProcessData& proc){
    LCC& lcc = proc.lcc;
    using Union_find = Union_find<Dart_handle>;
    using HandleToId = std::unordered_map<Union_find::pointer, size_t>;
    using VolToHandle = std::unordered_map<LCC::Attribute_handle<3>::type, Union_find::handle>;

    Union_find owned_uf, unowned_uf;
    VolToHandle owned_vol_to_handles, unowned_vol_to_handles;

    const auto process_volume = [&](Union_find& uf, VolToHandle& vol_to_handles){
      return [&](Dart_handle vol, auto& faces){
        Union_find::handle cc_id = nullptr;
        auto vol_attr = lcc.attribute<3>(vol);

        assert(vol != lcc.null_dart_descriptor);
        if (vol_attr->info().type < VolumeType::ID_EXPANSION)
          return;

        for (Dart_handle other_face : faces){
          Dart_handle other_vol = lcc.beta(other_face, 3);
          if (other_vol == lcc.null_dart_descriptor) continue;
          auto other_vol_attr = lcc.attribute<3>(other_vol);
          if (other_vol_attr == nullptr or other_vol_attr->info().type < VolumeType::ID_EXPANSION) continue;

          auto other_cc_it = owned_vol_to_handles.find(other_vol_attr);

          if (other_cc_it != owned_vol_to_handles.end()){
            Union_find::handle other_cc_id = other_cc_it->second;
            if (cc_id == nullptr)
              cc_id = other_cc_id;
            else if (!uf.same_set(cc_id, other_cc_id))
              uf.unify_sets(cc_id, other_cc_id);
          }
        }

        if (cc_id == nullptr){
          cc_id = uf.make_set(vol);
        }

        owned_vol_to_handles[vol_attr] = cc_id;
      };
    };

    const auto assign_cc_ids = [&](std::vector<Dart_handle>& old_ghost_vol_cc ,const Union_find& uf,
                                                  const std::vector<Union_find::handle>& partition, VolToHandle& vol_to_handle){
      HandleToId handle_to_ids;

      // Assign the new CC to the processdata::owned/unowned volumes cc
      old_ghost_vol_cc.clear();
      for (int i = 0; i < partition.size(); i++){
        handle_to_ids.emplace(partition[i].ptr(), i);
        old_ghost_vol_cc.push_back(*partition[i]);
      }

      for (auto& [vol, handle] : vol_to_handle){
        Union_find::handle cc_handle = uf.find(handle);
        auto cc_id = handle_to_ids.find(cc_handle.ptr());
        assert(cc_id != handle_to_ids.end());
        vol->info().cc_id = cc_id->second;
      }
    };

    // Calculate connected components of owned / unowned hexes
    for_each_hex(lcc, proc.owned_ghost_area, process_volume(owned_uf, owned_vol_to_handles),
      [&](Dart_handle vol){
        auto attr = lcc.attribute<3>(vol);
        auto& info = attr->info();
        return attr != nullptr && info.area_id != AreaId{0,0,0} && info.owned;
      }
    );

    for_each_hex(lcc, proc.unowned_ghost_area, process_volume(unowned_uf, unowned_vol_to_handles),
      [&](Dart_handle vol){
        auto attr = lcc.attribute<3>(vol);
        auto& info = attr->info();
        return attr != nullptr && !info.owned;
      }
    );

    assign_cc_ids(proc.owned_ghost_area, owned_uf, get_partitions(owned_uf), owned_vol_to_handles);
    assign_cc_ids(proc.unowned_ghost_area, unowned_uf, get_partitions(unowned_uf), unowned_vol_to_handles);
  }

  void setup_next_level(ProcessData& proc, MarkingFunction& cellIdentifier){
    setup_next_level_plane(proc);
    setup_ghost_cell_cc(proc);
    clean_up_and_reevaluate_attributes(proc, cellIdentifier);
    proc.level++;
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
}

namespace CGAL::HexRefinement {
  void debug_render(TwoRefinement::HexMeshingData& hdata){
    LCCSceneOptions<LCC> gso;
    LCC& lcc = hdata.lcc;

    AreaIDMap<CGAL::IO::Color> colors;
    colors[{0,0,0}] = blue();
    int i = 0;

    // これただ true を返すだけの関数になってるけど？
    gso.draw_face = [&](const LCC& lcc, LCC::Dart_const_handle dart){
      auto a = lcc.attribute<2>(dart);
      return true;
      return a != nullptr or lcc.is_whole_cell_marked<2>(dart, hdata.debug) or lcc.is_whole_cell_marked<2>(dart, hdata.debug2);
    };
    // これも最初の return 以降無意味だし、a は全く使われてない
    gso.face_color = [&](const LCC& lcc, LCC::Dart_const_handle dart){
      auto a = lcc.attribute<2>(dart);
      auto b = lcc.attribute<3>(dart);

      return b->info().cc_id != DartInfo::VolumeAttrValue::max_cc_id ? red() : blue();

      if (b == nullptr) return blue();

      colors[{0,0,0}] = blue();

      if (colors.count(b->info().area_id) == 0){
        CGAL::Random random((i+=3));
        colors[b->info().area_id] = CGAL::get_random_color(random);
      }

      return colors[b->info().area_id];
    };
    gso.colored_face = [](const LCC& lcc, LCC::Dart_const_handle dart){ return true; };

    // gso.draw_volume = [](const LCC& lcc, LCC::Dart_const_handle dart){
    //   return true;
    //   return lcc.attribute<3>(dart)->info().owned;
    // };

    gso.colored_vertex = [](const LCC& lcc, LCC::Dart_const_handle dart){
      return true;
    };
    gso.draw_vertex = [](const LCC& lcc, LCC::Dart_const_handle dart){
      return true;
    };
    // これも最初の行で完結しちゃってる。もしかしたら、なにかやりたかったことがあったのかも
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

    // while (true) {
    //   auto i = debug_stream.waitNextItem();
    //   debug_render(proc_datas[i]);
    // }

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
///////////////////////
/////     栞     //////
///////////////////////
