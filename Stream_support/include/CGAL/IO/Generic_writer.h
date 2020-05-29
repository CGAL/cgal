// Copyright (c) 1997
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>

#ifndef CGAL_IO_GENERIC_WRITER_H
#define CGAL_IO_GENERIC_WRITER_H

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <iostream>
#include <iterator>
#include <vector>

namespace CGAL {

template <typename Stream, typename FileWriter>
class Generic_writer
{
public:
  Generic_writer(Stream& out) : m_out(out) { }
  Generic_writer(Stream& out, FileWriter writer) : m_out(out), m_writer(writer) { }

  template <typename PointRange_3, typename PolygonRange_3>
  bool operator()(const PointRange_3& points,
                  const PolygonRange_3& polygons)
  {
    typedef typename PointRange_3::value_type Point_3;
    typedef typename PolygonRange_3::value_type Polygon_3;
    if(m_out.fail())
      return false;

    m_writer.write_header(m_out, points.size(), 0, polygons.size());
    for(std::size_t i=0, end=points.size(); i<end; ++i)
    {
      const Point_3& p = points[i];
      m_writer.write_vertex(p.x(), p.y(), p.z());
    }

    m_writer.write_facet_header();
    for(std::size_t i=0, end=polygons.size(); i<end; ++i)
    {
      const Polygon_3& polygon = polygons[i];
      const std::size_t size = polygon.size();

      m_writer.write_facet_begin(size);
      for(std::size_t j=0; j<size; ++j)
        m_writer.write_facet_vertex_index(polygon[j]);
      m_writer.write_facet_end();
    }
    m_writer.write_footer();

    return m_out.good();
  }

protected:
  Stream& m_out;
  FileWriter m_writer;
};

} // namespace CGAL

#endif // CGAL_IO_GENERIC_WRITER_H
