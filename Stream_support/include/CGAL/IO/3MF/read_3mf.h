// Copyright (c) 2019  Geometry Factory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Maxime Gimeno

#ifndef CGAL_IO_READ_3MF_H
#define CGAL_IO_READ_3MF_H

#ifdef CGAL_LINKED_WITH_LIB3MF

#include <CGAL/IO/Color.h>
#include <CGAL/Kernel_traits.h>

#include <boost/range/value_type.hpp>

#include <algorithm>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include <lib3mf_implicit.hpp>

namespace CGAL {

using namespace Lib3MF;

template<typename PointRange,
         typename TriangleRange,
         typename ColorRange>
bool extract_soups (PMeshObject pMeshObject,
                    const sTransform& transform,
                    PointRange& points,
                    TriangleRange& triangles,
                    ColorRange& colors,
                    std::string& name)
{
  typedef typename boost::range_value<PointRange>::type      Point_3;
  typedef typename boost::range_value<TriangleRange>::type   Triangle;
  typedef typename Kernel_traits<Point_3>::Kernel            Kernel;

  std::vector<char> pBuffer;

  name = pMeshObject->GetName();

  typename Kernel::Aff_transformation_3 t(
        transform.m_Fields[0][0], transform.m_Fields[0][1], transform.m_Fields[0][2], transform.m_Fields[0][3],
        transform.m_Fields[1][0], transform.m_Fields[1][1], transform.m_Fields[1][2], transform.m_Fields[1][3],
        transform.m_Fields[2][0], transform.m_Fields[2][1], transform.m_Fields[2][2], transform.m_Fields[2][3]
                                          );

  for(Lib3MF_uint32 vid = 0; vid < points.size(); ++vid)
  {
    sPosition pVertex = pMeshObject->GetVertex(vid);
    Point_3 p(pVertex.m_Coordinates[0],
              pVertex.m_Coordinates[1],
              pVertex.m_Coordinates[2]);
    points[vid] = t.transform(p);
  }

  // AF:  How to get from the property to the color ??
  // PColorGroupIterator pColorGroupIterator = model->GetColorGroups();
  // PColorGroup g;    g->GetResourceID(); // test for this one in sTriangleColor

  for(Lib3MF_uint32 pid = 0; pid < triangles.size(); ++pid)
  {
    sTriangle pTriangle = pMeshObject->GetTriangle(pid);
    Triangle triangle(3);
    for(Lib3MF_uint32 i = 0; i< 3; ++i)
      triangle[i] = pTriangle.m_Indices[i];

    triangles[pid] = triangle;
    sColor pColor;
    sTriangleProperties sTriangleColor;

    pMeshObject->GetTriangleProperties(pid, sTriangleColor);

    /*
    NMR::lib3mf_propertyhandler_getcolor(pPropertyHandler, pid, &pColor);
    NMR::MODELMESHCOLOR_SRGB mColor = pColor.m_Colors[0];
    colors[pid]=CGAL::IO::Color(mColor.m_Red, mColor.m_Green,
                            mColor.m_Blue, mColor.m_Alpha);
    */
  }

  return true;
}

} // namespace CGAL

#endif // CGAL_LINKED_WITH_LIB3MF

#endif // CGAL_IO_READ_3MF_H
