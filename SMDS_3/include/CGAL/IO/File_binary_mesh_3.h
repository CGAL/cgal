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

#include <CGAL/license/SMDS_3.h>


#include <iostream>
#include <string>
#include <limits>
#include <CGAL/SMDS_3/io_signature.h>

namespace CGAL {

namespace IO {

  /**
   * @ingroup PkgSMDS3IOFunctions
   * @brief outputs a mesh complex to the CGAL binary file format (`.binary.cgal`).
   *
   * @tparam C3T3 Type of mesh complex, model of `MeshComplex_3InTriangulation_3`
   *
   * @param os the output stream, opened in binary mode
   * @param c3t3 the mesh complex
   *
   * @sa `CGAL::IO::load_binary_file()`
   */
template <class C3T3>
bool
save_binary_file(std::ostream& os,
                 const C3T3& c3t3
#ifdef DOXYGEN_RUNNING
                 )
#else
               , bool binary = true)
#endif
{
  typedef typename C3T3::Triangulation::Geom_traits::FT FT;
  if(binary) os << "binary ";
  os << "CGAL c3t3 " << CGAL::Get_io_signature<C3T3>()() << "\n";
  if(binary) {
    CGAL::IO::set_binary_mode(os);
  } else {
    CGAL::IO::set_ascii_mode(os);
    os.precision(std::numeric_limits<FT>::digits10+2);
  }
  return !!(os << c3t3);
  // call operator!() twice, because operator bool() is C++11
}

/**
 * @ingroup PkgSMDS3IOFunctions
 * @brief loads a mesh complex from a file written in CGAL binary file format (`.binary.cgal`).
 *
 * @tparam C3T3 Type of mesh complex, model of `MeshComplex_3InTriangulation_3`
 *
 * @param is the input stream, opened in binary mode
 * @param c3t3 the mesh complex
 *
 * @sa `CGAL::IO::save_binary_file()`
 */
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
  if(binary) CGAL::IO::set_binary_mode(is);
  is >> c3t3;
  return !!is;
  // call operator!() twice, because operator bool() is C++11
}

} // end namespace IO

#ifndef CGAL_NO_DEPRECATED_CODE
namespace Mesh_3 {
using IO::save_binary_file;
using IO::load_binary_file;
}
#endif

} // end namespace CGAL

#endif // CGAL_IO_FILE_BINARY_MESH_3_H
