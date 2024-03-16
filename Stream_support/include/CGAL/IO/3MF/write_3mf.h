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

#ifndef CGAL_IO_WRITE_3MF_H
#define CGAL_IO_WRITE_3MF_H

#ifdef CGAL_LINKED_WITH_LIB3MF

#include <CGAL/IO/Color.h>

#include <CGAL/boost/graph/iterator.h>

#include <lib3mf_implicit.hpp>

#include <iostream>
#include <vector>
#include <string>

/*
 * \attention Only versions inferior to 2.0 of lib3mf are supported.
 * */
namespace CGAL {
namespace tmf_internal {

using namespace Lib3MF;

// Utility functions to create vertices and triangles
sPosition fnCreateVertex(float x, float y, float z)
{
  sPosition result;
  result.m_Coordinates[0] = x;
  result.m_Coordinates[1] = y;
  result.m_Coordinates[2] = z;
  return result;
}

sTriangle fnCreateTriangle(int v0, int v1, int v2)
{
  sTriangle result;
  result.m_Indices[0] = v0;
  result.m_Indices[1] = v1;
  result.m_Indices[2] = v2;
  return result;
}

sColor fnCreateColor(unsigned char red, unsigned char green,
                     unsigned char blue, unsigned char alpha=255)
{
  sColor result;
  result.m_Red = red;
  result.m_Green = green;
  result.m_Blue = blue;
  result.m_Alpha = alpha;
  return result;
}

sTriangleProperties fnCreateTriangleColor(PColorGroup colorGroup, Lib3MF_uint32 colorID1, Lib3MF_uint32 colorID2, Lib3MF_uint32 colorID3)
{
        sTriangleProperties sTriangleProperty;
        sTriangleProperty.m_ResourceID = colorGroup->GetResourceID();
        sTriangleProperty.m_PropertyIDs[0] = colorID1;
        sTriangleProperty.m_PropertyIDs[1] = colorID2;
        sTriangleProperty.m_PropertyIDs[2] = colorID3;
        return sTriangleProperty;
}

} // namespace tmf_internal

namespace IO {

using namespace Lib3MF;

bool add_build_item(PWrapper wrapper,
                    PModel pModel,
                    PMeshObject pMeshObject)
{
  /*
  HRESULT hResult;
  DWORD nErrorMessage;
  LPCSTR pszErrorMessage;
  */
  PBuildItem pBuildItem = pModel->AddBuildItem(pMeshObject.get(), wrapper->GetIdentityTransform());
  /*
  if(hResult != LIB3MF_OK)
  {
    std::cerr << "could not create build item: " << std::hex << hResult << std::endl;
    NMR::lib3mf_getlasterror(pModel, &nErrorMessage, &pszErrorMessage);
    std::cerr << "error #" << std::hex << nErrorMessage << ": " << pszErrorMessage << std::endl;
    NMR::lib3mf_release(pMeshObject);
    NMR::lib3mf_release(pModel);
    return false;
  }
  */
  return true;
}

bool export_model_to_file(const std::string& file_name,
                          PModel pModel)
{
  /*
  HRESULT hResult;
  DWORD nErrorMessage;
  LPCSTR pszErrorMessage;
  */
  // Output mesh as 3MF
  // Create Model Writer for 3MF
  PWriter writer = pModel->QueryWriter("3mf");
  writer->WriteToFile(file_name.c_str());
  /*
  if(hResult != LIB3MF_OK)
  {
    std::cerr << "could not create model reader: " << std::hex << hResult << std::endl;
    NMR::lib3mf_getlasterror(pModel, &nErrorMessage, &pszErrorMessage);
    std::cerr << "error #" << std::hex << nErrorMessage << ": " << pszErrorMessage << std::endl;
    NMR::lib3mf_release(pModel);
    return false;
  }

  // Export Model into File
  hResult = NMR::lib3mf_writer_writetofileutf8(p3MFWriter, file_name.c_str());
  if(hResult != LIB3MF_OK) {
    std::cerr << "could not write file: " << std::hex << hResult << std::endl;
    NMR::lib3mf_getlasterror(p3MFWriter, &nErrorMessage, &pszErrorMessage);
    std::cerr << "error #" << std::hex << nErrorMessage << ": " << pszErrorMessage << std::endl;
    NMR::lib3mf_release(pModel);
    NMR::lib3mf_release(p3MFWriter);
    return false;
  }
  */
  return true;
}

template<typename PointRange, typename TriangleRange, typename ColorRange>
bool write_mesh_to_model(const PointRange& points,
                         const TriangleRange& triangles,
                         const ColorRange& colors,
                         const std::string& name,
                         PMeshObject pMeshObject,
                         PModel  pModel,
                         PWrapper wrapper)
{
  /*
  DWORD nErrorMessage;
  LPCSTR pszErrorMessage;
  HRESULT hResult;
  */
  // Create mesh structure
  std::vector<sPosition> vertices;
  std::vector<sTriangle> striangles;

  for(const auto& point : points)
    vertices.push_back(tmf_internal::fnCreateVertex(point.x(), point.y(), point.z()));

  for(const auto& triangle : triangles)
    striangles.push_back(tmf_internal::fnCreateTriangle(triangle[0], triangle[1], triangle[2]));

  pMeshObject->SetGeometry(vertices, striangles);
  /*
  if(hResult != LIB3MF_OK)
  {
    std::cerr << "could not set mesh geometry: " << std::hex << hResult << std::endl;
    NMR::lib3mf_getlasterror(*pMeshObject, &nErrorMessage, &pszErrorMessage);
    std::cerr << "error #" << std::hex << nErrorMessage
              << ": " << pszErrorMessage << std::endl;
    NMR::lib3mf_release(*pMeshObject);
    NMR::lib3mf_release(pModel);
    return false;
  }
  */


  // define colors
  PColorGroup colorGroup = pModel->AddColorGroup();
  // AF: We should only create as many properties as we have faces
  for(std::size_t pid = 0; pid<colors.size(); ++pid)
  {
    Lib3MF_uint32 col = colorGroup->AddColor(wrapper->RGBAToColor(colors[pid].red(),
                                                                  colors[pid].green(),
                                                                  colors[pid].blue(),
                                                                  colors[pid].alpha()));

    sTriangleProperties sTriangleColor = tmf_internal::fnCreateTriangleColor(colorGroup, col, col, col);
    pMeshObject->SetTriangleProperties(pid, sTriangleColor);
  }

  // make sure to define a default property
  Lib3MF_uint32 lightgrey = colorGroup->AddColor(wrapper->RGBAToColor(211, 211, 212, 255));
  sTriangleProperties sTriangleColorLightgrey = tmf_internal::fnCreateTriangleColor(colorGroup, lightgrey, lightgrey, lightgrey);
  pMeshObject->SetObjectLevelProperty(sTriangleColorLightgrey.m_ResourceID, sTriangleColorLightgrey.m_PropertyIDs[0]);

  /*
  if(hResult != LIB3MF_OK)
  {
    std::cerr<< "could not create default property handler: " << std::hex << hResult << std::endl;
    NMR::lib3mf_getlasterror(*pMeshObject, &nErrorMessage, &pszErrorMessage);
    std::cerr<< "error #" << std::hex << nErrorMessage << ": " << pszErrorMessage << std::endl;
    NMR::lib3mf_release(*pMeshObject);
    NMR::lib3mf_release(pModel);
    return false;
  }
  */

  // Set name
  pMeshObject->SetName(name.c_str());
  /*
  if(hResult != LIB3MF_OK)
  {
    std::cerr << "could not set object name: " << std::hex << hResult << std::endl;
    NMR::lib3mf_getlasterror(*pMeshObject, &nErrorMessage, &pszErrorMessage);
    std::cerr << "error #" << std::hex << nErrorMessage << ": " << pszErrorMessage << std::endl;
    NMR::lib3mf_release(*pMeshObject);
    NMR::lib3mf_release(pModel);
    return false;
  }
  */
  //add a builditem to finish
  return add_build_item(wrapper, pModel, pMeshObject);
}

//remember that it adds 3 demmy vertices in the beginning, and a dummy triangle to be ignored.
template<typename PointRange, typename Color>
bool write_points(const PointRange& points,
                  const Color& color,
                  const std::string& name,
                  PMeshObject pMeshObject,
                  PModel pModel,
                  PWrapper wrapper)
{
  /*
  DWORD nErrorMessage;
  LPCSTR pszErrorMessage;
  HRESULT hResult;
  */
  // Create mesh structure
  std::vector<sPosition> vertices;
  std::vector<sTriangle> triangles(1);
  // Create Mesh Object
   pMeshObject = pModel->AddMeshObject();
   /*
  if(hResult != LIB3MF_OK)
  {
    std::cerr << "could not add mesh object: " << std::hex << hResult << std::endl;
    NMR::lib3mf_getlasterror(pModel, &nErrorMessage, &pszErrorMessage);
    std::cerr << "error #" << std::hex << nErrorMessage << ": " << pszErrorMessage << std::endl;
    NMR::lib3mf_release(pModel);
    return false;
  }
   */
  //add 3 dummy vertices to be sure to have a valid triangle and accept point sets with fewer than 3 vertices.
  for(int i = 0; i< 3; ++i)
    vertices.push_back(tmf_internal::fnCreateVertex(0,0,0));

  for(const auto& point : points)
    vertices.push_back(tmf_internal::fnCreateVertex(point.x(), point.y(), point.z()));

  triangles[0] =  tmf_internal::fnCreateTriangle(0,1,2); //add a triangle to avoid lib error.
  pMeshObject->SetGeometry(vertices, triangles);
  /*
  if(hResult != LIB3MF_OK)
  {
    std::cerr << "could not set mesh geometry: " << std::hex << hResult << std::endl;
    NMR::lib3mf_getlasterror(*pMeshObject, &nErrorMessage, &pszErrorMessage);
    std::cerr << "error #" << std::hex << nErrorMessage << ": " << pszErrorMessage << std::endl;
    NMR::lib3mf_release(*pMeshObject);
    NMR::lib3mf_release(pModel);
    return false;
  }
  */
  /*
  // Create color entries
  NMR::PLib3MFPropertyHandler * pPropertyHandler;
  hResult = NMR::lib3mf_meshobject_createpropertyhandler(*pMeshObject, &pPropertyHandler);
  if(hResult != LIB3MF_OK)
  {
    std::cerr << "could not create property handler: " << std::hex << hResult << std::endl;
    NMR::lib3mf_getlasterror(*pMeshObject, &nErrorMessage, &pszErrorMessage);
    std::cerr << "error #" << std::hex << nErrorMessage << ": " << pszErrorMessage << std::endl;
    NMR::lib3mf_release(*pMeshObject);
    NMR::lib3mf_release(pModel);
    return false;
  }

  // define colors
  sColor scolor = tmf_internal::fnCreateColor (color.red(),
                                               color.green(),
                                               color.blue(),
                                               color.alpha());
  // One-colored Triangles
  NMR::lib3mf_propertyhandler_setsinglecolor(pPropertyHandler, 0, &scolor);

  // make sure to define a default property
  NMR::PLib3MFDefaultPropertyHandler * pDefaultPropertyHandler;
  hResult = NMR::lib3mf_object_createdefaultpropertyhandler(*pMeshObject, &pDefaultPropertyHandler);
  if(hResult != LIB3MF_OK)
  {
    std::cerr<< "could not create default property handler: " << std::hex << hResult << std::endl;
    NMR::lib3mf_getlasterror(*pMeshObject, &nErrorMessage, &pszErrorMessage);
    std::cerr<< "error #" << std::hex << nErrorMessage << ": " << pszErrorMessage << std::endl;
    NMR::lib3mf_release(*pMeshObject);
    NMR::lib3mf_release(pModel);
    return false;
  }

  sColor default_color = tmf_internal::fnCreateColor(0,0,0,0);
  NMR::lib3mf_defaultpropertyhandler_setcolor(pDefaultPropertyHandler,
                                              &default_color);

  // release default property handler
  NMR::lib3mf_release(pDefaultPropertyHandler);
  */

  // Set name
  pMeshObject->SetName(name.c_str());
  /*
  if(hResult != LIB3MF_OK)
  {
    std::cerr << "could not set object name: " << std::hex << hResult << std::endl;
    NMR::lib3mf_getlasterror(*pMeshObject, &nErrorMessage, &pszErrorMessage);
    std::cerr << "error #" << std::hex << nErrorMessage << ": " << pszErrorMessage << std::endl;
    NMR::lib3mf_release(*pMeshObject);
    NMR::lib3mf_release(pModel);
    return false;
  }
  */

  return add_build_item(wrapper, pModel, pMeshObject);
}

template<typename PointRange, typename Color>
bool write_point_cloud_to_model(const PointRange& points,
                                const Color& color,
                                const std::string& name,
                                PMeshObject pMeshObject,
                                PModel pModel,
                                PWrapper wrapper)
{
  std::string pc_name = name;
  pc_name.append("_cgal_pc");
  return write_points(points, color, pc_name, pMeshObject, pModel, wrapper);
}

template<typename PointRange, typename Color>
bool write_polyline_to_model(const PointRange& points,
                             const Color& color,
                             const std::string& name,
                             PMeshObject pMeshObject,
                             PModel  pModel,
                             PWrapper wrapper)
{
  std::string pc_name = name;
  pc_name.append("_cgal_pl");
  return write_points(points, color, pc_name, pMeshObject, pModel, wrapper);
}

} // namespace IO
} // namespace CGAL

#endif // CGAL_LINKED_WITH_LIB3MF

#endif // CGAL_IO_WRITE_3MF_H
