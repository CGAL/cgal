#include "algo/3d/OutwardMeshOffset.h"

#include <filesystem>
#include <iostream>
#include <list>

int main(int argc, char** argv)
{
  // Can be any format read by CGAL::IO::read_polygon_mesh()
  const char* input_filename = argv[1];

  // Weight File Format:
  //   x1 vx1
  //   x2 vx2
  //   y1 vy1
  //   y2 vy2
  // optional:
  //   z1 vz1
  //   z2 vz2
  const char* weights_filename = (argc > 2) ? argv[2] : nullptr;

  // path of the directory where the output is written
  const std::filesystem::path save_path = (argc > 3) ? argv[3] : std::filesystem::current_path();

  // the offset(s) are global, but the geometric translation of the face
  // takes into account their respective speed.
  std::list<CGAL::FT> save_offsets = { -1 };

  bool success = algo::_3d::OutwardMeshOffset::run(input_filename,
                                                   weights_filename,
                                                   save_offsets,
                                                   save_path);

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
