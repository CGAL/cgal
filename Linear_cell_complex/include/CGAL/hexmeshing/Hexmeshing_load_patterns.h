#ifndef HEXMESHING_LOAD_PATTERNS_H
#define HEXMESHING_LOAD_PATTERNS_H

#include <CGAL/hexmeshing/LCC_items_for_hexmeshing.h>
#include <CGAL/query_replace/cmap_query_replace.h>


namespace CGAL::Hexmeshing {
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
}




#endif