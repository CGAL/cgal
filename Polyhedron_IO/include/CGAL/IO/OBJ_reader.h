// Copyright (c) 2016 GeometryFactory
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri and Maxime Gimeno

#ifndef CGAL_IO_OBJ_READER_H
#define CGAL_IO_OBJ_READER_H

#include <istream>
#include <sstream>
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
        if(i < 1)
        {
          faces.back().push_back(points.size()+i);//negative indices are relative references
        }
        else {
          faces.back().push_back(i-1);
        }
        iss.ignore(256, ' ');
      }
    }
    else
    {
      //std::cerr<<"ERROR : Cannnot read line beginning with "<<line[0]<<std::endl;
     continue;
    }
  }
  return true;
}

} // namespace CGAL

#endif // CGAL_IO_OBJ_READER_H
