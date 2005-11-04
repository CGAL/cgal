#include "weighted_types.h"

// ios
#include <iostream>
#include <fstream>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/IO/File_medit.h>
#include <CGAL/Mesh_3/IO.h>

#include "utils.h"

int main(int , char** argv)
{
  Tr tr;
  C2T3 c2t3(tr);

  std::ifstream ifs(argv[1]);
  if( !ifs )
  {
    std::cerr << "Usage:\n"
              << "  " << argv[0] << " FILE\n"
	      << "  FILE must be a file of format .mesh.cgal,\n"
              << "  produced by CGAL::Mesh_3::output_mesh(),\n"
              << "  with points CGAL::Weighted_point_with_surface_index\n"
              << "  and cells CGAL::Mesh_3::"
      "Complex_2_in_triangulation_cell_base_3.\n";
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
    return EXIT_SUCCESS;
  }
  else
    return EXIT_FAILURE;
}
