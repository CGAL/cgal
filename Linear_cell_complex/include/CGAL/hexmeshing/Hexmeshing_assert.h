#ifndef HEXMESHING_ASSERT_H
#define HEXMESHING_ASSERT_H

#include <CGAL/hexmeshing/LCC_items_for_hexmeshing.h>
#include <CGAL/Hexmeshing_for_linear_cell_complex_sequential.h>

namespace CGAL::Hexmeshing {
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
   * @tparam HexData Type of the hexahedral meshing data structure
   * @param hdata Hexahedral meshing data containing the mesh and marks
   * @param rdata Refinement data containing the faces of the current plane
   */
  template <typename HexData>
  void assert_faces_of_plane_valid(HexData& hdata, RefinementData& rdata){
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
   */
  void assert_all_faces_are_quadrilateral(LCC &lcc){
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




}


#endif