#include "types.h"

// ios
#include <iostream>
#include <fstream>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/IO/File_medit.h>

#include "utils.h"

int main(int , char** argv)
{
  Tr tr;
  C2T3 c2t3(tr);

  std::ifstream ifs(argv[1]);
  if( !ifs )
  {
    std::cerr << "usage:\n"
              << "  " << argv[0] << " file.mesh\n";
    return EXIT_FAILURE;
  }

  std::cout << "  Reading " << argv[1] << std::endl;
  if( CGAL::input_pslg_from_medit(ifs,
                                  c2t3,
                                  true,         // debug
                                  &std::cout) ) // debug to cout
  {
    display_faces_counts(tr, "    ", &std::cout);

    std::cout << "\n  Statistics:\n";

    std::cout << "(vertices)\n";
    display_vertices_by_surface_indices_statistics(tr, "    ", &std::cout);

    std:: cout << "(facets)\n";
    display_facets_by_surface_indices_statistics(c2t3, "    ", &std::cout);

    return EXIT_SUCCESS;
  }
  else
    return EXIT_FAILURE;
}
