// Copyright (c) 2015  Geometry Factory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Lutz Kettner
//             Andreas Fabri
//             Maxime Gimeno

#ifndef CGAL_IO_OBJ_H
#define CGAL_IO_OBJ_H

#include <CGAL/IO/OBJ/File_writer_wavefront.h>

namespace CGAL {

//! \ingroup IOstreamFunctions
//!
/// reads the content of `input` into `points` and `faces`, using the `OBJ` format.
///
/// \tparam Points_3 a `RandomAccessContainer` of `Point_3,
/// \tparam Faces a `RandomAccessContainer` of `RandomAccessContainer` of `std::size_t`
///
/// \see \ref IOStreamOBJ
template <class Points_3, class Faces>
bool read_OBJ(std::istream& input,
              Points_3& points,
              Faces& faces)
{
  typedef typename Points_3::value_type                       Point_3;

  int mini(1),
      maxi(-INT_MAX);
  Point_3 p;
  std::string line;

  while(getline(input, line))
  {
    if(line[0] == 'v' && line[1] == ' ')
    {
      std::istringstream iss(line.substr(1));
      iss >> p;
      if(!iss)
        return false;
      points.push_back(p);
    }
    else if(line[0] == 'f')
    {
      std::istringstream iss(line.substr(1));
      int i;
      faces.push_back(std::vector<std::size_t>());
      while(iss >> i)
      {
        if(i < 1)
        {
          faces.back().push_back(points.size()+i);//negative indices are relative references
          if(i<mini)
            mini=i;
        }
        else
        {
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
      //std::cerr << "ERROR : Cannnot read line beginning with " << line[0] << std::endl;
     continue;
    }
  }

  if(maxi > points.size() || mini < -static_cast<int>(points.size()))
  {
    std::cerr << "a face index is invalid " << std::endl;
    return false;
  }

  return true;
}

} // namespace CGAL

#endif // CGAL_IO_OBJ_H
