#ifndef HEXMESHING_TWO_REFINEMENT_ALGORITHM_H
#define HEXMESHING_TWO_REFINEMENT_ALGORITHM_H

#include <CGAL/hexmeshing/LCC_items_for_hexmeshing.h>
#include <CGAL/hexmeshing/Hexmeshing_assert.h>
#include <CGAL/hexmeshing/Hexmeshing_create_vertices_for_templates.h>
#include <CGAL/hexmeshing/Hexmeshing_expand_identified_cells.h>
#include <CGAL/hexmeshing/Hexmeshing_get_cells_to_refine.h>
#include <CGAL/hexmeshing/Hexmeshing_initial_setup.h>
#include <CGAL/hexmeshing/Hexmeshing_setup_next_level.h>
#include <CGAL/hexmeshing/Hexmeshing_two_refinement_3_template_utils.h>
#include <CGAL/hexmeshing/Hexmeshing_two_refinement_template_utils.h>
#include <CGAL/hexmeshing/Hexmeshing_function_alias.h>


namespace CGAL::Hexmeshing {
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
   * @pre HexData is HexMeshingData or ProcessData
   */
  template <typename HexData>
  void two_refinement_algorithm(HexData& hdata, MarkingFunction& cellIdentifier, int nb_levels, int thread_id = 0){
    // static_assert(std::is_same_v<HexData, HexMeshingData> or std::is_same_v<HexData, ProcessData>);
    // static_assert(std::is_base_of_v<HexMeshingData, HexData>);

    LCC& lcc = hdata.lcc;

    hdata.debug = lcc.get_new_mark();
    hdata.debug2 = lcc.get_new_mark();
    hdata.identified_mark = lcc.get_new_mark();
    hdata.template_mark = lcc.get_new_mark();
    hdata.propagation_face_mark = lcc.get_new_mark();
    // l_debug_mark = hdata.debug;
    // l_debug_mark_2 = hdata.debug2;
    // l_thread_id = thread_id;

    hdata.three_template_node_mark = lcc.get_new_mark();
    hdata.reset_temp_vertex_ids();

    Grid::generate_grid(lcc, hdata.grid);

    // Levels of refinement
    for (int r = 0; r < nb_levels; r++){
      if (r == 0) initial_setup(hdata, cellIdentifier);
      else setup_next_level(hdata, cellIdentifier);

      expand_identified_cells(hdata, r, nb_levels);

      // For each plane
      for (int p = 0; p < 3; p++) {
        RefinementData rdata;

        mark_identified_cells_from_3_attrs(hdata.lcc, hdata.identified_mark);

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

        assert_all_faces_are_quadrilateral(lcc);
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
}



#endif