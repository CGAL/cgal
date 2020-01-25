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

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read

//! \ingroup IOstreamFunctions
//!
/// reads the content of `input` into `points` and `faces`, using the `OBJ` format.
///
/// \tparam Points a `RandomAccessContainer` of `Point_3,
/// \tparam Faces a `RandomAccessContainer` of `RandomAccessContainer` of `std::size_t`
///
/// \see \ref IOStreamOBJ
template <class Points, class Faces>
bool read_OBJ(std::istream& input,
              Points& points,
              Faces& faces)
{
  typedef typename Points::value_type                          Point;

  int mini(1),
      maxi(-INT_MAX);
  Point p;
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
          faces.back().push_back(points.size()+i); // negative indices are relative references
          if(i < mini)
            mini = i;
        }
        else
        {
          faces.back().push_back(i-1);
          if(i-1 > maxi)
            maxi = i-1;
        }
        iss.ignore(256, ' ');
      }

      if(iss.fail())
        return false;
    }
    else
    {
      //std::cerr << "ERROR : Cannnot read line beginning with " << line[0] << std::endl;
     continue;
    }
  }

  if(maxi > static_cast<int>(points.size()) || mini < -static_cast<int>(points.size()))
  {
    std::cerr << "a face index is invalid " << std::endl;
    return false;
  }

  return !input.fail();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write


/*!
 * \ingroup IOstreamFunctions
 *
 * writes the content of `points` and `polygons` in `out`, in the OBJ format.
 *
 * \see \ref IOStreamOBJ
 */
template <class Point_3, class Polygon_3>
bool write_OBJ(std::ostream& out,
               std::vector<Point_3>& points,
               std::vector<Polygon_3>& polygons)
{
  CGAL::File_writer_wavefront writer;
  writer.write_header(out, points.size(), 0, polygons.size());

  for(std::size_t i = 0, end = points.size(); i<end; ++i)
  {
    const Point_3& p = points[i];
    writer.write_vertex(p.x(), p.y(), p.z());
  }

  writer.write_facet_header();
  for(std::size_t i=0, end=polygons.size(); i<end; ++i)
  {
    Polygon_3& polygon = polygons[i];
    const std::size_t size = polygon.size();

    writer.write_facet_begin(size);
    for(std::size_t j=0; j<size; ++j)
      writer.write_facet_vertex_index(polygon[j]);
    writer.write_facet_end();
  }
  writer.write_footer();

  return out.good();
}

} // namespace CGAL

#endif // CGAL_IO_OBJ_H
