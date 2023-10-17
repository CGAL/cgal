#ifndef CGAL_ALPHA_WRAP_3_EXAMPLES_OUTPUT_HELPER_H
#define CGAL_ALPHA_WRAP_3_EXAMPLES_OUTPUT_HELPER_H

#include <string>

std::string generate_output_name(std::string input_name,
                                 const double alpha,
                                 const double offset)
{
  input_name = input_name.substr(input_name.find_last_of("/") + 1, input_name.length() - 1);
  input_name = input_name.substr(0, input_name.find_last_of("."));
  std::string output_name = input_name
                            + "_" + std::to_string(static_cast<int>(alpha))
                            + "_" + std::to_string(static_cast<int>(offset)) + ".off";

  return output_name;
}

#endif // CGAL_ALPHA_WRAP_3_EXAMPLES_OUTPUT_HELPER_H
