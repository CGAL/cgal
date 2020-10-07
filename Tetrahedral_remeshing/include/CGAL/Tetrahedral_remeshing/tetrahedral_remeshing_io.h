// Copyright (c) 2020 GeometryFactory (France) and Telecom Paris (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois, Noura Faraj, Jean-Marc Thiery, Tamy Boubekeur

#ifndef CGAL_TETRAHEDRAL_REMESHING_IO_H
#define CGAL_TETRAHEDRAL_REMESHING_IO_H

#include <CGAL/license/Tetrahedral_remeshing.h>

#include <CGAL/IO/io.h>

#include <iostream>
#include <fstream>

namespace CGAL
{
template<typename T3>
bool load_triangulation(std::istream& is, T3& t3)
{
  std::string s;
  if (!(is >> s)) return false;
  bool binary = (s == "binary");
  if (binary) {
    if (!(is >> s)) return false;
  }
  if (s != "CGAL" || !(is >> s) || s != "c3t3")
    return false;

  std::getline(is, s);
  if (binary) CGAL::set_binary_mode(is);
  else        CGAL::set_ascii_mode(is);
  is >> t3;
  return bool(is);
}

template<typename T3>
bool save_binary_triangulation(std::ostream& os, const T3& t3)
{
  os << "binary CGAL c3t3\n";
  CGAL::set_binary_mode(os);
  return !!(os << t3);
}

template<typename T3>
bool save_ascii_triangulation(std::ostream& os, const T3& t3)
{
  os << "CGAL c3t3\n";
  CGAL::set_ascii_mode(os);
  return !!(os << t3);
}

}

#endif // CGAL_TETRAHEDRAL_REMESHING_IO_H
