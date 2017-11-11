// Copyright (c) 2016 GeometryFactory
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Andreas Fabri and Maxime Gimeno

#ifndef CGAL_IO_OBJ_READER_H
#define CGAL_IO_OBJ_READER_H

#include <istream>
#include <vector>


namespace CGAL {

template <class Point_3>
bool
read_OBJ( std::istream& input,
          std::vector<Point_3> &points,
          std::vector<std::vector<std::size_t> > &faces)
{
  Point_3 p;
  std::string line;
  while(getline(input, line)) {
    if(line[0] == 'v' && line[1] == ' ') {
      std::istringstream iss(line.substr(1));
      iss >> p;
      if(!iss)
        return false;
      points.push_back(p);
    }
    else if(line[0] == 'f') {
      std::istringstream iss(line.substr(1));
      int i;
      faces.push_back( std::vector<std::size_t>() );
      while(iss >> i)
      {
        faces.back().push_back(i-1);
        iss.ignore(256, ' ');
      }
    }
  }
  return true;
}

} // namespace CGAL

#endif // CGAL_IO_OBJ_READER_H
