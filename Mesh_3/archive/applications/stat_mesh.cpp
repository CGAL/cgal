#include "weighted_types.h"

// ios
#include <iostream>
#include <fstream>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/Mesh_3/IO.h>

#include "utils.h"
#include <CGAL/min_dihedral_angle.h>

int main(int , char** argv)
{
  Tr tr;
  C2T3 c2t3(tr);

  std::ifstream ifs(argv[1]);
  if( !ifs )
  {
    std::cerr << "Usage:\n"
              << "  " << argv[0] << " FILE\n"
              << "\n"
              << "  FILE must be " << format_cgal_description
              << "\n";
    return EXIT_FAILURE;
  }

  std::cout << "  Reading " << argv[1] << std::endl;

  if( CGAL::Mesh_3::input_mesh(ifs, c2t3,
                               true,
                               &std::cerr) )
//   if( CGAL::input_pslg_from_medit(ifs,
//                                   c2t3,
//                                   true,         // debug
//                                   &std::cout) ) // debug to cout
  {
    display_faces_counts(tr, "    ", &std::cout);

    std::cout << "\n  Statistics:\n";

    std::cout << "(vertices)\n";
    display_vertices_by_surface_indices_statistics(tr, "    ", &std::cout);

    std::cout << "(facets)\n";
    display_facets_by_surface_indices_statistics(c2t3, "    ", &std::cout);

    Compute_min_angle<Tr> min_angle(tr);

    double min = 180;
    for(Tr::Finite_cells_iterator cit = tr.finite_cells_begin();
        cit != tr.finite_cells_end();
        ++cit)
      if(cit->is_in_domain())
      {
        const double angle = min_angle(cit);
        if( angle < min ) min = angle;
      }

    std::cout << "\nmin angle: " << min << "\n";

    return EXIT_SUCCESS;
  }
  else
    return EXIT_FAILURE;
}
