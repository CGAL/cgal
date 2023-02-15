#include <CGAL/Classification/ETHZ/Random_forest_classifier.h>

#include <cstdlib>
#include <fstream>

int main (int argc, char** argv)
{
  if (argc != 3)
    std::cerr << "Usage: " << argv[0] << " input.gz output.bin" << std::endl;
  else
  {
    std::ifstream ifile (argv[1], std::ios_base::binary);
    std::ofstream ofile (argv[2], std::ios_base::binary);

    CGAL::Classification::ETHZ::Random_forest_classifier::
      convert_deprecated_configuration_to_new_format(ifile, ofile);
  }

  return EXIT_SUCCESS;
}
