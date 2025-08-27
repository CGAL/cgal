#ifndef HEXMESHING_POST_PROCESSING_H
#define HEXMESHING_POST_PROCESSING_H

#include <CGAL/hexmeshing/LCC_items_for_hexmeshing.h>
#include <CGAL/hexmeshing/Hexmeshing_set_attributes.h>
#include <CGAL/hexmeshing/Hexmeshing_move_points_onto_mesh.h>
#include <CGAL/hexmeshing/Hexmeshing_laplacian_smoothing.h>
#include <CGAL/hexmeshing/Hexmeshing_two_refinement_utils.h>

namespace CGAL::internal::Hexmeshing {
  void post_processing(LCC& lcc, double length_of_4_template, bool trim, MarkingFunction cellIdentifier, DecideInsideFunction decideFunc) {
    size_type move_mark = lcc.get_new_mark();
    size_type inner_mark = lcc.get_new_mark();

    set_fraction(lcc, length_of_4_template, cellIdentifier, decideFunc);

    move_points_onto_mesh_with_volume_fraction(lcc, move_mark, inner_mark);

    // smoothing
    surface_smoothing(lcc, move_mark, inner_mark);
    volume_smoothing(lcc, move_mark);

    // trimming
    if(trim)
      trim_excedent_volumes(lcc, is_marked_volume(inner_mark));

    lcc.free_mark(move_mark);
    lcc.free_mark(inner_mark);
  }
}




#endif