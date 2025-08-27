
#ifndef HEXMESHING_SETUP_NEXT_LEVEL_H
#define HEXMESHING_SETUP_NEXT_LEVEL_H

#include <CGAL/hexmeshing/LCC_items_for_hexmeshing.h>
#include <CGAL/hexmeshing/Hexmeshing_function_alias.h>
#include <vector>
#include <cassert>


namespace CGAL::internal::Hexmeshing {
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
   * @tparam HexData Type of the hexahedral meshing data structure
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex
   * @param cellIdentifier Function that determines whether a cell should be identified
   *                       for refinement based on its position and properties
   */
  template <typename HexData>
  void clean_up_and_reevaluate_attributes(HexData& hdata, MarkingFunction& cellIdentifier){
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
   * @tparam HexData Type of the hexahedral meshing data structure
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
  template <typename HexData>
  void setup_next_level_face(HexData& hdata,
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
        // debug_stream.push(l_thread_id);
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
   * @tparam HexData Type of the hexahedral meshing data structure
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and plane information
   */
  template <typename HexData>
  void setup_next_level_plane(HexData& hdata){
    using FaceToHandle = std::unordered_map<LCC::Attribute_handle<2>::type, Union_find<Dart_handle>::handle>;
    using Union_find = Union_find<Dart_handle>;
    using UF_Partition = std::vector<Union_find::handle>;
    using PlaneCC = std::vector<Dart_handle>;
    using PlaneSet = std::vector<PlaneCC>;

    LCC& lcc = hdata.lcc;
    std::array<PlaneSet, 3> new_planes;

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
   * @tparam HexData Type of the hexahedral meshing data structure
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and plane information
   * @param cellIdentifier Function that determines whether a cell should be identified
   *                       for refinement based on its position and properties
   */
  template <typename HexData>
  void setup_next_level(HexData& hdata, MarkingFunction& cellIdentifier){
    setup_next_level_plane(hdata);
    clean_up_and_reevaluate_attributes(hdata, cellIdentifier);
    hdata.level++;
  }
}



#endif