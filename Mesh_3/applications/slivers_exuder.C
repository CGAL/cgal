#include "weighted_types.h"

// ios
#include <iostream>
#include <fstream>
#include <sstream>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>

#include <CGAL/Mesh_3/Slivers_exuder.h>

#include "utils.h"

int main(int argc, char** argv)
{
  Tr tr;
  C2T3 c2t3(tr);

  std::ifstream ifs(argv[1]);
  std::ofstream ofs(argv[2]);
  if( argc != 4 || !ifs || !ofs )
  {
    std::cerr <<
"Usage:\n" <<
"  " << argv[0] << " INPUT OUTPUT bound\n" <<
"\n" <<
"  INPUT must be " << format_cgal_description <<
"\n" <<
"  OUPUT is the name of the ouput file. It will be in the same format.\n" <<
"\n"
"  0<=bound<=1.0 is a bound on pumping process.\n";
    return EXIT_FAILURE;
  }

  std::stringstream str_stream(argv[3]);

  double pumping_bound;

  str_stream >> pumping_bound;

  std::cout << "  Reading " << argv[1] << std::endl;
  if(! CGAL::Mesh_3::input_mesh(ifs,
                                c2t3,
                                true,         // debug
                                &std::cerr) ) // debug to cerr
    return EXIT_FAILURE;
  ifs.close();

  CGAL::Mesh_3::Slivers_exuder<Tr> exuder(tr);

  std::cout << "  Pumping" << std::endl;
  exuder.init(pumping_bound);
  exuder.pump_vertices(pumping_bound);

  std::cout << "  Writing " << argv[2] << std::endl;
  CGAL::Mesh_3::output_mesh(ofs, c2t3);

  return EXIT_SUCCESS;
}
