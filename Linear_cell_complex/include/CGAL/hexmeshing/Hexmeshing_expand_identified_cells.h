#ifndef CGAL_HEXMESHING_EXPAND_IDENTIFIED_CELLS_H
#define CGAL_HEXMESHING_EXPAND_IDENTIFIED_CELLS_H

#include <CGAL/hexmeshing/LCC_items_for_hexmeshing.h>

namespace CGAL::Hexmeshing {
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
   * @tparam HexData Type of the hexahedral meshing data structure
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex
   * @param current_lvl The current refinement level (0-based)
   * @param nb_levels The total number of refinement levels to be performed
   */
  template <typename HexData>
  void expand_identified_cells(HexData& hdata, int current_lvl, int nb_levels){
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
}

#endif // CGAL_HEXMESHING_EXPAND_IDENTIFIED_CELLS_H
