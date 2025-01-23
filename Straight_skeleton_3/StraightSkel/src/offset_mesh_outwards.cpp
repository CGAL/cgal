#include "algo/3d/OutwardMeshOffset.h"

#include <filesystem>
#include <iostream>
#include <list>

int main(int argc, char** argv)
{
  const char* input_filename = argv[1];
  const char* weights_filename = (argc > 2) ? argv[2] : nullptr;
  const std::filesystem::path save_path = (argc > 3) ? argv[3] : std::filesystem::current_path();

  // the offset(s) are global, but the geometric translation of the face
  // takes into account their respective speed.
  std::list<CGAL::FT> save_offsets = { -1 };

  bool success = algo::_3d::OutwardMeshOffset::run(input_filename,
                                                   weights_filename,
                                                   save_offsets,
                                                   save_path);
  std::cout << "success = " << success << std::endl;

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
