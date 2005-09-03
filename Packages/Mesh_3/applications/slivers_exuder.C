#include "weighted_types.h"

// ios
#include <iostream>
#include <fstream>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/IO/File_medit.h>

#include <CGAL/Mesh_3/Slivers_exuder.h>

#include "utils.h"

int main(int , char** argv)
{
  Tr tr;
  C2T3 c2t3(tr);

  std::ifstream ifs(argv[1]);
  std::ofstream ofs(argv[2]);
  if( !ifs || !ofs)
  {
    std::cerr << "usage:\n"
              << "  " << argv[0] << " input.mesh output.mesh\n";
    return EXIT_FAILURE;
  }

  std::cout << "  Reading " << argv[1] << std::endl;
  if(! CGAL::input_pslg_from_medit(ifs,
                                   c2t3,
                                   true,         // debug
                                   &std::cout) ) // debug to cout
    return EXIT_FAILURE;
  ifs.close();

  CGAL::Mesh_3::Slivers_exuder<Tr> exuder(tr);

  std::cout << "  Pumping" << std::endl;
  exuder.init();
  exuder.pump_vertices();

  std::cout << "  Writing " << argv[2] << std::endl;
  CGAL::output_pslg_to_medit(ofs, c2t3);

  return EXIT_SUCCESS;
}
