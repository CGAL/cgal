#include "algo/3d/OutwardMeshOffset.h"

int main(int argc, char** argv)
{
  const char* input_filename = argv[1];
  const char* weights_filename = (argc > 2) ? argv[2] : nullptr;
  const char* output_filename = (argc > 3) ? argv[3] : "out.ply";

  bool success = algo::_3d::OutwardMeshOffset::run(input_filename, weights_filename, output_filename);
  std::cout << "success = " << success << std::endl;

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
