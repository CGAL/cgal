// Copyright (c) 2015 GeometryFactory
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
//
// Author(s)     : Laurent Rineau and Sebastien Loriot

#ifndef CGAL_IO_OFF_READER_H
#define CGAL_IO_OFF_READER_H

#include <CGAL/IO/File_scanner_OFF.h>

#include <vector>
#include <iostream>

namespace CGAL{

  template <class Point_3, class Polygon_3>
  bool
  read_OFF( std::istream& in,
            std::vector< Point_3 >& points,
            std::vector< Polygon_3 >& polygons,
            bool /* verbose */ = false)
  {
    CGAL::File_scanner_OFF scanner(in);

    points.resize(scanner.size_of_vertices());
    polygons.resize(scanner.size_of_facets());
    for (std::size_t i = 0; i < scanner.size_of_vertices(); ++i) {
      double x, y, z, w;
      scanner.scan_vertex( x, y, z, w);
      CGAL_assertion(w!=0);
      points[i] = Point_3(x/w, y/w, z/w);
      scanner.skip_to_next_vertex( i);
    }
    if(!in)
      return false;

    for (std::size_t i = 0; i < scanner.size_of_facets(); ++i) {
      std::size_t no;

      scanner.scan_facet( no, i);
      polygons[i].resize(no);
      for(std::size_t j = 0; j < no; ++j) {
        std::size_t id;
        scanner.scan_facet_vertex_index(id, i);
        if(id < scanner.size_of_vertices())
        {
          polygons[i][j] = id;
        }
        else
          return false;
      }
    }
    return (bool) in;
  }

} // namespace CGAL

#endif // CGAL_IO_OFF_READER_H
