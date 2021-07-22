#include "weighted_types.h"

// ios
#include <iostream>
#include <fstream>
#include <sstream>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/Timer.h>

#include <CGAL/Mesh_3/Slivers_exuder.h>

#include "utils.h"

int main(int argc, char** argv)
{
  Tr tr;
  C2T3 c2t3(tr);

  std::ifstream ifs(argv[1]);
  std::ofstream ofs(argv[2]);
  if( !(argc == 4 ||argc == 5) || !ifs || !ofs )
  {
    std::cerr <<
"Usage:\n" <<
"  " << argv[0] << " INPUT OUTPUT bound [filename.off]\n" <<
"\n" <<
"  INPUT must be " << format_cgal_description <<
"\n" <<
"  OUPUT is the name of the ouput file. It will be in the same format.\n" <<
"\n"
"  0<=bound<=1.0 is a bound on pumping process.\n"
"\n"
"  If non-nul, filename.off is the name of an output off file, where\n"
"  only slivers are displayed.\n";
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

  if(pumping_bound == 0.)
    return EXIT_SUCCESS;

  CGAL::Mesh_3::Slivers_exuder<C2T3> exuder(c2t3);

  std::cout << "  Pumping" << std::endl;
  CGAL::Timer timer;
  timer.start();
  exuder.init(pumping_bound);
  if(pumping_bound > 0.)
    // pump vertices on surfaces
    exuder.pump_vertices(pumping_bound);
  else {
    // do not pump vertices on surfaces
    pumping_bound = -pumping_bound;
    exuder.pump_vertices<false>(pumping_bound);
  }
  timer.stop();
  std::cout << "  Pumping done. CPU time: " << timer.time() << std::endl;

  CGAL::IO::set_binary_mode(ofs);
  std::cout << "  Writing " << argv[2] << std::endl;
  CGAL::Mesh_3::output_mesh(ofs, c2t3);
  ofs.close();
  if( argc >= 5)
  {
    ofs.open(argv[4]);
    if(ofs) {
      std::cout << "  Writing slivers to " << argv[4] << std::endl;
      CGAL::output_slivers_to_off(ofs, tr, pumping_bound);
    }
  }

  return EXIT_SUCCESS;
}
