
#ifndef TOOLS_H
#define TOOLS_H

#include <fstream>
#include <string>
#include <vector>


std::string read_file(const std::string& file_name)
{
  std::ifstream ifs(file_name.c_str(), std::ios::in | std::ios::binary | 
                                                                 std::ios::ate);

  std::ifstream::pos_type file_size = ifs.tellg();
  ifs.seekg(0, std::ios::beg);

  std::vector<char> bytes(file_size);
  ifs.read(&bytes[0], file_size);

  return std::string(&bytes[0], file_size);
}


#endif
