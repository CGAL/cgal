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


#include <CGAL/IO/Color.h>

#include <CGAL/Kernel_traits.h>

#include <boost/range/value_type.hpp>

#ifdef CGAL_LINKED_WITH_3MF
#include <Model/COM/NMR_DLLInterfaces.h>
#endif

#include <algorithm>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#ifdef CGAL_LINKED_WITH_3MF

namespace CGAL {
namespace transform_nmr_internal {

NMR::MODELTRANSFORM initMatrix()
{
  NMR::MODELTRANSFORM mMatrix;
  int i, j;
  for(i = 0; i < 4; i++) {
    for(j = 0; j < 3; j++) {
      mMatrix.m_fFields[j][i] = (i == j) ? 1.0f : 0.0f;
    }
  }

  return mMatrix;
}

} // namespace transform_nmr_internal

template<typename PointRange,
         typename TriangleRange,
         typename ColorRange>
bool extract_soups (NMR::PLib3MFModelMeshObject *pMeshObject,
                    const NMR::MODELTRANSFORM& transform,
                    PointRange& points,
                    TriangleRange& triangles,
                    ColorRange& colors,
                    std::string& name)
{
  typedef typename boost::range_value<PointRange>::type      Point_3;
  typedef typename boost::range_value<TriangleRange>::type   Triangle;
  typedef typename Kernel_traits<Point_3>::Kernel            Kernel;

  HRESULT hResult;
  DWORD nNeededChars;
  std::vector<char> pBuffer;

  // Retrieve Mesh Name Length
  hResult = NMR::lib3mf_object_getnameutf8(pMeshObject, NULL, 0, &nNeededChars);
  if(hResult != LIB3MF_OK)
  {
    std::cerr<<"Error during name extraction.";
    return false;
  }

  // Retrieve Mesh Name
  if(nNeededChars > 0)
  {
    pBuffer.resize(nNeededChars + 1);
    hResult = NMR::lib3mf_object_getnameutf8(pMeshObject, &pBuffer[0], nNeededChars + 1, NULL);
    pBuffer[nNeededChars] = 0;
    name = std::string(&pBuffer[0]);
  }
  else
  {
    name = std::string("Unknown Mesh");
  }

  typename Kernel::Aff_transformation_3 t(
        transform.m_fFields[0][0], transform.m_fFields[0][1], transform.m_fFields[0][2], transform.m_fFields[0][3],
      transform.m_fFields[1][0], transform.m_fFields[1][1], transform.m_fFields[1][2], transform.m_fFields[1][3],
      transform.m_fFields[2][0], transform.m_fFields[2][1], transform.m_fFields[2][2], transform.m_fFields[2][3]
      );

  NMR::PLib3MFPropertyHandler * pPropertyHandler;
  hResult = NMR::lib3mf_meshobject_createpropertyhandler(pMeshObject, &pPropertyHandler);
  if(hResult != LIB3MF_OK)
  {
    DWORD nErrorMessage;
    LPCSTR pszErrorMessage;
    std::cerr << "could not create property handler: " << std::hex << hResult << std::endl;
    NMR::lib3mf_getlasterror(pMeshObject, &nErrorMessage, &pszErrorMessage);
    std::cerr << "error #" << std::hex << nErrorMessage << ": " << pszErrorMessage << std::endl;
    NMR::lib3mf_release(pMeshObject);
    return false;
  }

  for(DWORD vid = 0; vid < points.size(); ++vid)
  {
    NMR::MODELMESHVERTEX pVertex;
    NMR::lib3mf_meshobject_getvertex(pMeshObject, vid, &pVertex);
    Point_3 p(pVertex.m_fPosition[0],
        pVertex.m_fPosition[1],
        pVertex.m_fPosition[2]);
    points[vid] = t.transform(p);
  }

  for(DWORD pid = 0; pid < triangles.size(); ++pid)
  {
    NMR::MODELMESHTRIANGLE pTriangle;
    NMR::lib3mf_meshobject_gettriangle(pMeshObject, pid, &pTriangle);
    Triangle triangle(3);
    for(DWORD i = 0; i< 3; ++i)
      triangle[i] = pTriangle.m_nIndices[i];
    triangles[pid] = triangle;
    NMR::MODELMESH_TRIANGLECOLOR_SRGB pColor;
    NMR::lib3mf_propertyhandler_getcolor(pPropertyHandler, pid, &pColor);
    NMR::MODELMESHCOLOR_SRGB mColor = pColor.m_Colors[0];
    colors[pid]=CGAL::IO::Color(mColor.m_Red, mColor.m_Green,
                            mColor.m_Blue, mColor.m_Alpha);
  }

  return true;
}

} // namespace CGAL

#endif // CGAL_LINKED_WITH_3MF

#endif // CGAL_IO_READ_3MF_H
