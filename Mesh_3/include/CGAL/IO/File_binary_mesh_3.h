// Copyright (c) 2012  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_IO_FILE_BINARY_MESH_3_H
#define CGAL_IO_FILE_BINARY_MESH_3_H

#include <CGAL/license/Triangulation_3.h>


#include <iostream>
#include <string>
#include <limits>
#include <CGAL/Mesh_3/io_signature.h>

namespace CGAL {

namespace Mesh_3 {

template <class C3T3>
bool
save_binary_file(std::ostream& os,
                 const C3T3& c3t3,
                 bool binary = true)
{
  typedef typename C3T3::Triangulation::Geom_traits::FT FT;
  if(binary) os << "binary ";
  os << "CGAL c3t3 " << CGAL::Get_io_signature<C3T3>()() << "\n";
  if(binary) {
    CGAL::set_binary_mode(os);
  } else {
    CGAL::set_ascii_mode(os);
    os.precision(std::numeric_limits<FT>::digits10+2);
  }
  return !!(os << c3t3);
  // call operator!() twice, because operator bool() is C++11
}

template <class C3T3>
bool load_binary_file(std::istream& is, C3T3& c3t3)
{
  std::string s;
  if(!(is >> s)) return false;
  bool binary = (s == "binary");
  if(binary) {
    if(!(is >> s)) return false;
  }
  if (s != "CGAL" ||
      !(is >> s) ||
      s != "c3t3")
  {
    return false;
  }
  std::getline(is, s);
  if(!s.empty()) {
    if(s[s.size()-1] == '\r') { // deal with Windows EOL
      s.resize(s.size() - 1);
    }
    if(s != std::string(" ") + CGAL::Get_io_signature<C3T3>()()) {
      std::cerr << "load_binary_file:"
                << "\n  expected format: " << CGAL::Get_io_signature<C3T3>()()
                << "\n       got format:" << s << std::endl;
      return false;
    }
  }
  if(binary) CGAL::set_binary_mode(is);
  is >> c3t3;
  return !!is;
  // call operator!() twice, because operator bool() is C++11
}

} // end namespace Mesh_3
} // end namespace CGAL

#endif // CGAL_IO_FILE_BINARY_MESH_3_H
