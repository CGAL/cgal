
#ifndef HEXMESHING_TWO_REFINEMENT_TEMPLATE_UTILS_H
#define HEXMESHING_TWO_REFINEMENT_TEMPLATE_UTILS_H

#include <CGAL/hexmeshing/LCC_items_for_hexmeshing.h>
#include <CGAL/hexmeshing/Hexmeshing_two_refinement_mark_utils.h>
#include <CGAL/query_replace/cmap_signature.h>
#include <CGAL/query_replace/cmap_query_replace.h>
#include <iostream>
#include <boost/range/join.hpp>


namespace CGAL::Hexmeshing {
  const size_t CONST_SIZE_T_MAX = std::numeric_limits<size_t>::max();
  
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
   * @tparam HexData Type of the hexahedral meshing data structure
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and marks
   * @param rdata Refinement data containing the faces of the current plane and collections
   * @return The total number of fixes applied (sum of diagonal 2-templates and adjacent 3-templates fixed)
   */
  template <typename HexData>
  int fix_impossible_cases(HexData &hdata, RefinementData &rdata){
    int fix_c_count = 0, fix_3_count = 0;

    // Condition for checking faces : 2 or 3 templates.
    std::queue<Dart_handle> faces_to_check;

    // First iteration

    // faces_of_plane will grow, we only want to examine only faces that existed right now
    int faces_end = rdata.faces_of_plane.size();
    
    for (int i = 0; i < faces_end; i++){
      Dart_handle face = rdata.faces_of_plane[i];
      auto& face_attr = hdata.lcc.template attribute<2>(face)->info();

      if (face_attr.template_id == 2 && fix_mark_connectivity(hdata.lcc, rdata, hdata.template_mark, faces_to_check, face))
        fix_c_count++;

      if (face_attr.template_id == 3 && fix_adjacent_3_templates(hdata.lcc, rdata, hdata.template_mark, faces_to_check, face))
        fix_3_count++;
    }

    // Repeat until there are no more faces to check
    while (!faces_to_check.empty()){
      Dart_handle front_face = faces_to_check.front();
      faces_to_check.pop();

      auto& face_attr = hdata.lcc.template attribute<2>(front_face)->info();

      if (face_attr.template_id == 2 && fix_mark_connectivity(hdata.lcc, rdata, hdata.template_mark, faces_to_check, front_face))
        fix_c_count++;

      if (face_attr.template_id == 3 && fix_adjacent_3_templates(hdata.lcc, rdata, hdata.template_mark, faces_to_check, front_face))
        fix_3_count++;
    }

    std::cout << "Diagonal 2 templates fixed: " << fix_c_count << std::endl;
    std::cout << "Neighboring 3 templates repaired: " << fix_3_count << std::endl;

    return fix_c_count + fix_3_count;
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
}




#endif