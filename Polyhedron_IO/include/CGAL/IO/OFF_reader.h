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
#include <CGAL/array.h>
#include <CGAL/assertions.h>
#include <CGAL/use.h>

namespace CGAL{

  namespace read_OFF_internal{
    template <class Point_3>
    void fill_point(double x, double y, double z, Point_3& pt)
    {
      pt = Point_3(x, y, z);
    }

    void fill_point(double x, double y, double z, CGAL::cpp11::array<double,3>& p)
    {
      p = CGAL::make_array(x,y,z);
    }

    template <class Polygon_3>
    void resize(Polygon_3& p, std::size_t size)
    {
      p.resize(size);
    }

    template <std::size_t N, class INT>
    void resize(CGAL::cpp11::array<INT, N>&, std::size_t size)
    {
      CGAL_USE(size);
      CGAL_assertion( size>=N );
    }
  }

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
      read_OFF_internal::fill_point( x/w, y/w, z/w, points[i] );
      scanner.skip_to_next_vertex( i);
    }
    if(!in)
      return false;

    for (std::size_t i = 0; i < scanner.size_of_facets(); ++i) {
      std::size_t no;

      scanner.scan_facet( no, i);
      read_OFF_internal::resize(polygons[i], no);
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
    return in.good();
  }

} // namespace CGAL

#endif // CGAL_IO_OFF_READER_H
