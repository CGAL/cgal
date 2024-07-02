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

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <iostream>
#include <iterator>
#include <vector>

namespace CGAL {

template <typename Stream, typename FileWriter>
class Generic_writer
{
public:
  Generic_writer(Stream& os) : m_os(os) { }
  Generic_writer(Stream& os, FileWriter writer) : m_os(os), m_writer(writer) { }

  template <typename PointRange, typename PolygonRange, typename NamedParameters>
  bool operator()(const PointRange& points,
                  const PolygonRange& polygons,
                  const NamedParameters& np = parameters::default_values())
  {
    typedef typename boost::range_value<PolygonRange>::type                   Poly;

    using parameters::choose_parameter;
    using parameters::get_parameter;

    typedef typename CGAL::GetPointMap<PointRange, NamedParameters>::type     PointMap;
    PointMap point_map = choose_parameter<PointMap>(get_parameter(np, internal_np::point_map));

    if(!m_os.good())
      return false;

    set_stream_precision_from_NP(m_os, np);

    m_writer.write_header(m_os, points.size(), 0, polygons.size());
    for(std::size_t i=0, end=points.size(); i<end; ++i)
    {
      const typename boost::property_traits<PointMap>::value_type& p = get(point_map, points[i]);
      m_writer.write_vertex(to_double(p.x()), to_double(p.y()), to_double(p.z()));
    }

    m_writer.write_facet_header();
    for(std::size_t i=0, end=polygons.size(); i<end; ++i)
    {
      const Poly& polygon = polygons[i];
      const std::size_t size = polygon.size();

      m_writer.write_facet_begin(size);
      for(std::size_t j=0; j<size; ++j)
        m_writer.write_facet_vertex_index(polygon[j]);
      m_writer.write_facet_end();
    }
    m_writer.write_footer();

    return m_os.good();
  }

protected:
  Stream& m_os;
  FileWriter m_writer;
};

} // namespace CGAL

#endif // CGAL_IO_GENERIC_WRITER_H
