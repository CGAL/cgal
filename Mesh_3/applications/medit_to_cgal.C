#include "weighted_types.h"

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
  std::ofstream ofs(argv[2]);
  if( !ifs || !ofs )
  {
    std::cerr << "Usage:\n"
              << "  " << argv[0] << " INPUT OUPUT\n" <<
"\n"
"  INPUT must be a medit file.\n"
"\n"
"  OUTPUT will be " << format_cgal_description <<
"\n";
    return EXIT_FAILURE;
  }

  std::cout << "  Reading " << argv[1] << std::endl;

  if( CGAL::input_from_medit(ifs,
                             c2t3,
                             true,         // debug
                             &std::cout) ) // debug to cout
  {
    std::cout << "  Writing " << argv[2] << std::endl;
    CGAL::Mesh_3::output_mesh(ofs, c2t3);
    return EXIT_SUCCESS;
  }
  else
    return EXIT_FAILURE;
}
