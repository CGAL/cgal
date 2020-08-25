#ifndef READ_OBJECTS_H
#define READ_OBJECTS_H

#include <CGAL/basic.h>

#include <fstream>
#include <string>

// Open an input stream.
inline bool open_stream(std::ifstream& is, const char* filename)
{
  is.open(filename);
  if (is.is_open()) return true;
  std::string error_message("Failed to open ");
  error_message.append(filename).append("!");
  CGAL_error_msg(error_message.c_str());
  return false;
}

// Read objects from an input stream.
template<typename Type, typename Output_iterator>
Output_iterator
read_stream(std::ifstream& is, unsigned int n, Output_iterator oi)
{
  for (size_t i = 0; i < n; ++i) {
    Type s;
    is >> s;
    *oi++ = s;
  }
  return oi;
}

// Close an input stream.
inline void close_stream(std::ifstream& is) { is.close(); }

// Read geometric objects from a file and write them to a given output iterator.
template<typename Type, typename Output_iterator>
Output_iterator read_objects(const char* filename, Output_iterator oi)
{
  std::ifstream  is;
  if (! open_stream(is, filename)) return oi;   // open the input file
  unsigned int n;
  is >> n;                                      // read number of objects
  read_stream<Type>(is, n, oi);                 // read objects from the file
  close_stream(is);                             // close the input file
  return oi;
}

#endif
