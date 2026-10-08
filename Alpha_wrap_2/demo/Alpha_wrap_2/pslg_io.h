// Copyright (c) 2019-2020 X, The Moonshot Factory (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later
//
//
// Author(s)     : Pierre Alliez pierre.alliez@inria.fr
//               : Michael Hemmer mhsaar@gmail.com
//               : Cedric Portaneri cportaneri@gmail.com
//
#ifndef CGAL_ALPHA_WRAP_2_IO_H
#define CGAL_ALPHA_WRAP_2_IO_H

#include "pslg_2.h"

#include <CGAL/IO/io.h>

#include <fstream>
#include <iostream>
#include <string>

namespace CGAL {
namespace Alpha_wraps_2 {
namespace IO {

enum IO_exit_code
{
  VALID_INPUT_POLYLINES = 0,
  UNREADABLE_INPUT = 1,
  INPUT_IS_EMPTY = 2,
};

inline int string_to_int(const std::string &str)
{
  char *end;
  return static_cast<int>(strtol(str.c_str(), &end, 10));  // NO Large INT
}

inline double string_to_double(const std::string &str)
{
  return atof(str.c_str());
}

inline bool check_extension(std::string file_name,
                            const std::string& extension)
{
  size_t find_ext = file_name.find_last_of(".");
  if (find_ext == std::string::npos)
    return false;
  std::string file_name_ext = file_name.substr(find_ext, file_name.size() - 1);
  return (file_name_ext.compare(extension) == 0);
}

template <typename GeomTraits>
inline bool read_input_dat_file(std::ifstream &in,
                                internal::Pslg<GeomTraits>& input_polylines)
{
  using Point_2 = typename GeomTraits::Point_2;

  internal::Component<GeomTraits> component;
  std::string line;
  while(std::getline(in, line))
  {
    if(line[0] == 'n') {
      input_polylines.push_back(component);
      component.clear();
    } else {
      double x,y;
      std::istringstream iss(line);
      if(iss >> x >> y)
        component.push_back(Point_2(x,y));
      else return false;
    }
  }
  return true;
}

template <typename GeomTraits>
inline IO_exit_code read_input_polylines_file(const std::string &file_name,
                                              internal::Pslg<GeomTraits>& input_polylines)
{
  std::cout << "Read surface polylines... " << file_name << std::endl;
  std::ifstream in(file_name, std::ios::in);
  if (!in) {
    std::cout << "Unable to open the file." << std::endl;
    return UNREADABLE_INPUT;
  }

  if (!check_extension(file_name, ".dat")) {
    return UNREADABLE_INPUT;
  } else if (!read_input_dat_file(in, input_polylines)) {
    std::cout << "Unable to read stl file." << std::endl;
    return UNREADABLE_INPUT;
  }
  in.close();

  if (input_polylines.empty()) {
    std::cout << "Input is empty." << std::endl;
    return INPUT_IS_EMPTY;
  }

  return VALID_INPUT_POLYLINES;
}

template <typename GeomTraits>
inline bool write_outputput_dat_file(std::ofstream &out,
                                     const internal::Pslg<GeomTraits>& output_polylines)
{
  for(const auto& component : output_polylines) {
    for(const auto& point : component) {
      out << point.x() << " " << point.y() << "\n";
    }
    out << "n\n";
  }
  return true;
}

template <typename GeomTraits>
inline void write_output_polylines_file(const std::string &file_name,
                                        const internal::Pslg<GeomTraits>& output_polylines)
{
  std::ofstream out(file_name);
  if (check_extension(file_name, ".dat")) {
    write_outputput_dat_file(out, output_polylines);
  }
}

} // namespace IO
} // namespace Alpha_wraps_2
} // namespace CGAL

#endif  // CGAL_ALPHA_WRAP_2_IO_H
