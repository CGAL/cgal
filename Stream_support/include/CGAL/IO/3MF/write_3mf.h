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

#ifdef CGAL_LINKED_WITH_3MF

#include <CGAL/IO/Color.h>

#include <CGAL/boost/graph/iterator.h>

#include <Model/COM/NMR_DLLInterfaces.h>

#include <iostream>
#include <vector>
#include <string>

/*
 * \attention Only versions inferior to 2.0 of lib3mf are supported.
 * */
namespace CGAL {
namespace tmf_internal {

// Utility functions to create vertices and triangles
NMR::MODELMESHVERTEX fnCreateVertex(float x, float y, float z)
{
  NMR::MODELMESHVERTEX result;
  result.m_fPosition[0] = x;
  result.m_fPosition[1] = y;
  result.m_fPosition[2] = z;
  return result;
}

NMR::MODELMESHTRIANGLE fnCreateTriangle(int v0, int v1, int v2)
{
  NMR::MODELMESHTRIANGLE result;
  result.m_nIndices[0] = v0;
  result.m_nIndices[1] = v1;
  result.m_nIndices[2] = v2;
  return result;
}

NMR::MODELMESHCOLOR_SRGB fnCreateColor(unsigned char red, unsigned char green,
                                       unsigned char blue, unsigned char alpha=255)
{
  NMR::MODELMESHCOLOR_SRGB result;
  result.m_Red = red;
  result.m_Green = green;
  result.m_Blue = blue;
  result.m_Alpha = alpha;
  return result;
}

} // namespace tmf_internal

namespace IO {

bool add_build_item(NMR::PLib3MFModel * pModel,
                    NMR::PLib3MFModelMeshObject* pMeshObject)
{
  HRESULT hResult;
  DWORD nErrorMessage;
  LPCSTR pszErrorMessage;

  // Add Build Item for Mesh
  NMR::PLib3MFModelBuildItem * pBuildItem;
  hResult = NMR::lib3mf_model_addbuilditem(pModel, pMeshObject, NULL, &pBuildItem);

  if(hResult != LIB3MF_OK)
  {
    std::cerr << "could not create build item: " << std::hex << hResult << std::endl;
    NMR::lib3mf_getlasterror(pModel, &nErrorMessage, &pszErrorMessage);
    std::cerr << "error #" << std::hex << nErrorMessage << ": " << pszErrorMessage << std::endl;
    NMR::lib3mf_release(pMeshObject);
    NMR::lib3mf_release(pModel);
    return false;
  }

  // Release BuildItem and Mesh
  NMR::lib3mf_release(pMeshObject);
  NMR::lib3mf_release(pBuildItem);
  return true;
}

bool export_model_to_file(const std::string& file_name,
                          NMR::PLib3MFModel * pModel)
{
  HRESULT hResult;
  DWORD nErrorMessage;
  LPCSTR pszErrorMessage;

  // Output mesh as 3MF
  // Create Model Writer for 3MF
  NMR::PLib3MFModelWriter * p3MFWriter;
  hResult = NMR::lib3mf_model_querywriter(pModel, "3mf", &p3MFWriter);

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

  // Release Model Writer
  NMR::lib3mf_release(p3MFWriter);
  return true;
}

template<typename PointRange, typename TriangleRange, typename ColorRange>
bool write_mesh_to_model(const PointRange& points,
                         const TriangleRange& triangles,
                         const ColorRange& colors,
                         const std::string& name,
                         NMR::PLib3MFModelMeshObject** pMeshObject,
                         NMR::PLib3MFModel * pModel)
{
  DWORD nErrorMessage;
  LPCSTR pszErrorMessage;
  HRESULT hResult;

  // Create mesh structure
  std::vector<NMR::MODELMESHVERTEX> pVertices;
  std::vector<NMR::MODELMESHTRIANGLE> pTriangles;

  // Create Mesh Object
  hResult = NMR::lib3mf_model_addmeshobject(pModel, pMeshObject);
  if(hResult != LIB3MF_OK)
  {
    std::cerr << "could not add mesh object: " << std::hex << hResult << std::endl;
    NMR::lib3mf_getlasterror(pModel, &nErrorMessage, &pszErrorMessage);
    std::cerr << "error #" << std::hex << nErrorMessage << ": " << pszErrorMessage << std::endl;
    NMR::lib3mf_release(pModel);
    return false;
  }

  for(const auto& point : points)
    pVertices.push_back(tmf_internal::fnCreateVertex(point.x(), point.y(), point.z()));

  for(const auto& triangle : triangles)
    pTriangles.push_back(tmf_internal::fnCreateTriangle(triangle[0], triangle[1], triangle[2]));

  hResult = NMR::lib3mf_meshobject_setgeometry(*pMeshObject, pVertices.data(),
                                               pVertices.size(), pTriangles.data(),
                                               pTriangles.size());
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
  for(std::size_t pid = 0; pid<colors.size(); ++pid)
  {
    NMR::MODELMESHCOLOR_SRGB sColor = tmf_internal::fnCreateColor (colors[pid].red(),
                                                                   colors[pid].green(),
                                                                   colors[pid].blue(),
                                                                   colors[pid].alpha());
    // One-colored Triangles
    NMR::lib3mf_propertyhandler_setsinglecolor(pPropertyHandler, pid, &sColor);
  }

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

  NMR::MODELMESHCOLOR_SRGB default_color = tmf_internal::fnCreateColor(0,0,0,0);
  NMR::lib3mf_defaultpropertyhandler_setcolor(pDefaultPropertyHandler, &default_color);

  // release default property handler
  NMR::lib3mf_release(pDefaultPropertyHandler);

  // Set name
  hResult = NMR::lib3mf_object_setnameutf8(*pMeshObject, name.c_str());
  if(hResult != LIB3MF_OK)
  {
    std::cerr << "could not set object name: " << std::hex << hResult << std::endl;
    NMR::lib3mf_getlasterror(*pMeshObject, &nErrorMessage, &pszErrorMessage);
    std::cerr << "error #" << std::hex << nErrorMessage << ": " << pszErrorMessage << std::endl;
    NMR::lib3mf_release(*pMeshObject);
    NMR::lib3mf_release(pModel);
    return false;
  }

  //add a builditem to finish
  return add_build_item(pModel, *pMeshObject);
}

//remember that it adds 3 demmy vertices in the beginning, and a dummy triangle to be ignored.
template<typename PointRange, typename Color>
bool write_points(const PointRange& points,
                  const Color& color,
                  const std::string& name,
                  NMR::PLib3MFModelMeshObject** pMeshObject,
                  NMR::PLib3MFModel * pModel)
{
  DWORD nErrorMessage;
  LPCSTR pszErrorMessage;
  HRESULT hResult;

  // Create mesh structure
  std::vector<NMR::MODELMESHVERTEX> pVertices;

  // Create Mesh Object
  hResult = NMR::lib3mf_model_addmeshobject(pModel, pMeshObject);
  if(hResult != LIB3MF_OK)
  {
    std::cerr << "could not add mesh object: " << std::hex << hResult << std::endl;
    NMR::lib3mf_getlasterror(pModel, &nErrorMessage, &pszErrorMessage);
    std::cerr << "error #" << std::hex << nErrorMessage << ": " << pszErrorMessage << std::endl;
    NMR::lib3mf_release(pModel);
    return false;
  }

  //add 3 dummy vertices to be sure to have a valid triangle and accept point sets with fewer than 3 vertices.
  for(int i = 0; i< 3; ++i)
    pVertices.push_back(tmf_internal::fnCreateVertex(0,0,0));

  for(const auto& point : points)
    pVertices.push_back(tmf_internal::fnCreateVertex(point.x(), point.y(), point.z()));

  NMR::MODELMESHTRIANGLE dummy_triangle = tmf_internal::fnCreateTriangle(0,1,2); //add a triangle to avoid lib error.
  hResult = NMR::lib3mf_meshobject_setgeometry(*pMeshObject, pVertices.data(),
                                               pVertices.size(), &dummy_triangle, 1);
  if(hResult != LIB3MF_OK)
  {
    std::cerr << "could not set mesh geometry: " << std::hex << hResult << std::endl;
    NMR::lib3mf_getlasterror(*pMeshObject, &nErrorMessage, &pszErrorMessage);
    std::cerr << "error #" << std::hex << nErrorMessage << ": " << pszErrorMessage << std::endl;
    NMR::lib3mf_release(*pMeshObject);
    NMR::lib3mf_release(pModel);
    return false;
  }

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
  NMR::MODELMESHCOLOR_SRGB sColor = tmf_internal::fnCreateColor (color.red(),
                                                                 color.green(),
                                                                 color.blue(),
                                                                 color.alpha());
  // One-colored Triangles
  NMR::lib3mf_propertyhandler_setsinglecolor(pPropertyHandler, 0, &sColor);

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

  NMR::MODELMESHCOLOR_SRGB default_color = tmf_internal::fnCreateColor(0,0,0,0);
  NMR::lib3mf_defaultpropertyhandler_setcolor(pDefaultPropertyHandler,
                                              &default_color);

  // release default property handler
  NMR::lib3mf_release(pDefaultPropertyHandler);

  // Set name
  hResult = NMR::lib3mf_object_setnameutf8(*pMeshObject, name.c_str());
  if(hResult != LIB3MF_OK)
  {
    std::cerr << "could not set object name: " << std::hex << hResult << std::endl;
    NMR::lib3mf_getlasterror(*pMeshObject, &nErrorMessage, &pszErrorMessage);
    std::cerr << "error #" << std::hex << nErrorMessage << ": " << pszErrorMessage << std::endl;
    NMR::lib3mf_release(*pMeshObject);
    NMR::lib3mf_release(pModel);
    return false;
  }

  return add_build_item(pModel, *pMeshObject);
}

template<typename PointRange, typename Color>
bool write_point_cloud_to_model(const PointRange& points,
                                const Color& color,
                                const std::string& name,
                                NMR::PLib3MFModelMeshObject** pMeshObject,
                                NMR::PLib3MFModel * pModel)
{
  std::string pc_name = name;
  pc_name.append("_cgal_pc");
  return write_points(points, color, pc_name, pMeshObject, pModel);
}

template<typename PointRange, typename Color>
bool write_polyline_to_model(const PointRange& points,
                             const Color& color,
                             const std::string& name,
                             NMR::PLib3MFModelMeshObject** pMeshObject,
                             NMR::PLib3MFModel * pModel)
{
  std::string pc_name = name;
  pc_name.append("_cgal_pl");
  return write_points(points, color, pc_name, pMeshObject, pModel);
}

} // namespace IO
} // namespace CGAL

#endif // CGAL_LINKED_WITH_3MF

#endif // CGAL_IO_WRITE_3MF_H
