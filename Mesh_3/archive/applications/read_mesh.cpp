#include "weighted_types.h"

// ios
#include <iostream>
#include <fstream>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>

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
              << "\n"
              << "  FILE must be " << format_cgal_description
              << "\n";
    return EXIT_FAILURE;
  }

  std::cout << "  Reading " << argv[1] << std::endl;

  if( CGAL::Mesh_3::input_mesh(ifs, c2t3,
                               true,
                               &std::cerr) )
  {
    display_faces_counts(tr, "    ", &std::cout);
    return EXIT_SUCCESS;
  }
  else
    return EXIT_FAILURE;
}
