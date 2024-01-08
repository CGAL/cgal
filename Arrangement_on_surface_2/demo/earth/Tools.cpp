// Copyright(c) 2023, 2024 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#include "Tools.h"

#include <fstream>
#include <iostream>
#include <vector>

//! \brief
std::string read_file(const std::string& file_name) {
  const auto flags = std::ios::in | std::ios::binary | std::ios::ate;
  std::ifstream ifs(file_name.c_str(), flags);

  if (! ifs.is_open()) {
    std::cout << "could not open file: " << file_name << std::endl;
    return "";
  }

  std::ifstream::pos_type file_size = ifs.tellg();
  ifs.seekg(0, std::ios::beg);

  std::vector<char> bytes(file_size);
  ifs.read(&bytes[0], file_size);

  return std::string(&bytes[0], file_size);
}

//! \brief
std::ostream& operator << (std::ostream& os, const QVector2D& v) {
  os << v.x() << ", " << v.y();
  return os;
}

//! \brief
std::ostream& operator << (std::ostream& os, const QVector3D& v) {
  os << v.x() << ", " << v.y() << ", " << v.z();
  return os;
}

//! \brief
std::ostream& operator << (std::ostream& os, const QVector4D& v) {
  os << v.x() << ", " << v.y() << ", " << v.z() << ", " << v.w();
  return os;
}
