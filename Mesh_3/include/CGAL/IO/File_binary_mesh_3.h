// Copyright (c) 2012  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_IO_FILE_BINARY_MESH_3_H
#define CGAL_IO_FILE_BINARY_MESH_3_H

#include <iostream>
#include <string>
#include <CGAL/Mesh_3/io_signature.h>

namespace CGAL {

namespace Mesh_3 {

template <class C3T3>
bool
save_binary_file(std::ostream& os,
                 const C3T3& c3t3)
{
  os << "binary CGAL c3t3 " << CGAL::Get_io_signature<C3T3>()() << "\n";
  CGAL::set_binary_mode(os);
  return !!(os << c3t3);
  // call operator!() twice, because operator bool() is C++11
}

template <class C3T3>
bool load_binary_file(std::istream& is, C3T3& c3t3)
{
  std::string s;
  is >> s;
  if (s != "binary" ||
      !(is >> s) ||
      s != "CGAL" ||
      !(is >> s) ||
      s != "c3t3") 
  {
    return false;
  }
  std::getline(is, s);
  if(s != "") {
    if(s != std::string(" ") + CGAL::Get_io_signature<C3T3>()()) {
      std::cerr << "load_binary_file:"
                << "\n  expected format: " << CGAL::Get_io_signature<C3T3>()()
                << "\n       got format:" << s << std::endl;
      return false;
    }
  }
  CGAL::set_binary_mode(is);
  is >> c3t3;
  return !!is;
  // call operator!() twice, because operator bool() is C++11
}

} // end namespace Mesh_3
} // end namespace CGAL

#endif // CGAL_IO_FILE_BINARY_MESH_3_H
