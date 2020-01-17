// Copyright (c) 2016 GeometryFactory
//
// This file is part of CGAL (www.cgal.org);
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri and Maxime Gimeno

#ifndef CGAL_IO_OBJ_OBJ_READER_H
#define CGAL_IO_OBJ_OBJ_READER_H

#include <istream>
#include <sstream>
#include <vector>
#include <limits>


namespace CGAL {

//Points_3 a RandomAccessContainer of Point_3,
//Faces a RandomAccessContainer of RandomAccessContainer of std::size_t
template <class Points_3, class Faces>
bool
read_OBJ( std::istream& input,
          Points_3 &points,
          Faces &faces)
{
  typedef typename Points_3::value_type Point_3;
  int mini(1),
      maxi(-INT_MAX);
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
          if(i<mini)
            mini=i;
        }
        else {
          faces.back().push_back(i-1);
          if(i-1 > maxi)
            maxi = i-1;
        }
        iss.ignore(256, ' ');
      }
      if(!iss.good() && !iss.eof())
        return false;
    }
    else
    {
      //std::cerr<<"ERROR : Cannnot read line beginning with "<<line[0]<<std::endl;
     continue;
    }
  }
  if(maxi > points.size() || mini < -static_cast<int>(points.size())){
    std::cerr<<"a face index is invalid "<<std::endl;
    return false;
  }
  return true;
}

} // namespace CGAL

#endif // CGAL_IO_OBJ_OBJ_READER_H
