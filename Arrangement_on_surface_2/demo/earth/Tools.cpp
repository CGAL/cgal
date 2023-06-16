
#include "Tools.h"

#include <fstream>
#include <iostream>
#include <vector>


std::string read_file(const std::string& file_name)
{
  const auto flags = std::ios::in | std::ios::binary | std::ios::ate;
  std::ifstream ifs(file_name.c_str(), flags);

  if (ifs.is_open() == false)
  {
    std::cout << "could not open file: " << file_name << std::endl;
    return "";
  }

  std::ifstream::pos_type file_size = ifs.tellg();
  ifs.seekg(0, std::ios::beg);

  std::vector<char> bytes(file_size);
  ifs.read(&bytes[0], file_size);

  return std::string(&bytes[0], file_size);
}
