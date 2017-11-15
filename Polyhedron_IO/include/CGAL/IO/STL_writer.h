// Copyright (c) 2017 GeometryFactory
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_IO_STL_WRITER_H
#define CGAL_IO_STL_WRITER_H

#include <CGAL/license/Polyhedron.h>

#include <CGAL/Kernel_traits.h>
#include <CGAL/boost/graph/properties.h>

#include <boost/cstdint.hpp>
#include <boost/foreach.hpp>
#include <boost/graph/graph_traits.hpp>



namespace CGAL{

template <class TriangleMesh>
std::ostream&
write_STL(const TriangleMesh& tm, std::ostream& out)
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::property_map<TriangleMesh, boost::vertex_point_t>::const_type Vpm;
  typedef typename boost::property_traits<Vpm>::reference Point_3_ref;
  typedef typename boost::property_traits<Vpm>::value_type Point_3;
  typedef typename Kernel_traits<Point_3>::Kernel::Vector_3 Vector_3;

  Vpm vpm = get(boost::vertex_point, tm);

  if (get_mode(out) == IO::BINARY)
  {
    out << "FileType: Binary                                                                ";
    const boost::uint32_t N32 = static_cast<boost::uint32_t>(faces(tm).size());
    out.write(reinterpret_cast<const char *>(&N32), sizeof(N32));

    BOOST_FOREACH(face_descriptor f, faces(tm))
    {
      halfedge_descriptor h = halfedge(f, tm);
      Point_3_ref p = get(vpm, target(h, tm));
      Point_3_ref q = get(vpm, target(next(h, tm), tm));
      Point_3_ref r = get(vpm, source(h, tm));

      Vector_3 n = collinear(p,q,r) ? Vector_3(1,0,0):
                                      unit_normal(p,q,r);

      const float coords[12]={
        static_cast<float>(n.x()), static_cast<float>(n.y()), static_cast<float>(n.z()),
        static_cast<float>(p.x()), static_cast<float>(p.y()), static_cast<float>(p.z()),
        static_cast<float>(q.x()), static_cast<float>(q.y()), static_cast<float>(q.z()),
        static_cast<float>(r.x()), static_cast<float>(r.y()), static_cast<float>(r.z()) };

      for (int i=0; i<12; ++i)
        out.write(reinterpret_cast<const char *>(&coords[i]), sizeof(coords[i]));
      out << "  ";
    }
  }
  else
  {
    out << "solid\n";
    BOOST_FOREACH(face_descriptor f, faces(tm))
    {
      halfedge_descriptor h = halfedge(f, tm);
      Point_3_ref p = get(vpm, target(h, tm));
      Point_3_ref q = get(vpm, target(next(h, tm), tm));
      Point_3_ref r = get(vpm, source(h, tm));

      Vector_3 n = collinear(p,q,r) ? Vector_3(1,0,0):
                                      unit_normal(p,q,r);
      out << "facet normal " << n << "\nouter loop\n";
      out << "vertex " << p << "\n";
      out << "vertex " << q << "\n";
      out << "vertex " << r << "\n";
      out << "endloop\nendfacet\n";
    }
    out << "endsolid\n";
  }
  return out;
}

} // end of namespace CGAL

#endif // CGAL_IO_STL_WRITER_H
