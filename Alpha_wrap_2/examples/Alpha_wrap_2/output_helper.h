#ifndef CGAL_ALPHA_WRAP_2_EXAMPLES_OUTPUT_HELPER_H
#define CGAL_ALPHA_WRAP_2_EXAMPLES_OUTPUT_HELPER_H

#include <fstream>
#include <string>
#include <filesystem>

std::filesystem::path generate_output_name(const std::filesystem::path& input_name,
                                           const double alpha,
                                           const double offset)
{
  std::filesystem::path output_name = input_name.stem();
  std::string suffix("_" + std::to_string(static_cast<int>(alpha))
                         + "_" + std::to_string(static_cast<int>(offset)) + "-wrap.wkt");
  output_name += suffix;

  return output_name;
}

#endif // CGAL_ALPHA_WRAP_2_EXAMPLES_OUTPUT_HELPER_H
