// Copyright (c) 2017 GeometryFactory
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
// Author(s)     : Sebastien Loriot

#ifndef CGAL_IO_STL_STL_WRITER_H
#define CGAL_IO_STL_STL_WRITER_H

#include <CGAL/Kernel_traits.h>
#include <CGAL/Kernel/global_functions.h>

#include <boost/cstdint.hpp>



namespace CGAL{

template <class PointRange, class TriangleRange>
std::ostream&
write_STL(const PointRange& points,
          const TriangleRange& facets,
          std::ostream& out)
{
  typedef typename PointRange::value_type Point;
  typedef typename CGAL::Kernel_traits<Point>::Kernel K;
  typedef typename K::Vector_3 Vector_3;


  if (get_mode(out) == IO::BINARY)
  {
    out << "FileType: Binary                                                                ";
    const boost::uint32_t N32 = static_cast<boost::uint32_t>(facets.size());
    out.write(reinterpret_cast<const char *>(&N32), sizeof(N32));

    for(auto face : facets)
    {
      const Point& p = points[face[0]];
      const Point& q = points[face[1]];
      const Point& r = points[face[2]];

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
    for(auto face : facets)
    {
      const Point& p = points[face[0]];
      const Point& q = points[face[1]];
      const Point& r = points[face[2]];

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

#endif // CGAL_IO_STL_STL_WRITER_H
