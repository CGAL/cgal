// Copyright (c) 2019  Geometry Factory
// All rights reserved.
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
// Author(s) : Maxime Gimeno

#ifndef READ_3MF_H
#define READ_3MF_H
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <functional>
#include <CGAL/IO/Color.h>
#include <Model/COM/NMR_DLLInterfaces.h>
namespace CGAL{

template<typename PointRange,
         typename PolygonRange,
         typename ColorRange>
bool extract_soups (NMR::PLib3MFModelMeshObject *pMeshObject,
           PointRange& points,
           PolygonRange& triangles,
           ColorRange& colors,
           std::string& name) {
  typedef typename PointRange::value_type Point_3;
  typedef typename PolygonRange::value_type Polygon;
  typedef typename ColorRange::value_type _Color;

  HRESULT hResult;
  DWORD nNeededChars;
  std::vector<char> pBuffer;
  // Retrieve Mesh Name Length
  hResult = NMR::lib3mf_object_getnameutf8(pMeshObject, NULL, 0, &nNeededChars);
  if (hResult != LIB3MF_OK)
  {
    std::cerr<<"Error during name extraction.";
    return false;
  }

  // Retrieve Mesh Name
  if (nNeededChars > 0) {
    pBuffer.resize(nNeededChars + 1);
    hResult = NMR::lib3mf_object_getnameutf8(pMeshObject, &pBuffer[0], nNeededChars + 1, NULL);
    pBuffer[nNeededChars] = 0;
    std::string temp(&pBuffer[0]);
    if(temp.find("_cgal_pc") != std::string::npos
       || temp.find("_cgal_pl")!= std::string::npos) //ignore point clouds and polylines
    {
      return false;
    }
    name = std::string(&pBuffer[0]);
  }
  else
    name = std::string("Unknown Mesh");

  NMR::PLib3MFPropertyHandler * pPropertyHandler;
  hResult = NMR::lib3mf_meshobject_createpropertyhandler(pMeshObject, &pPropertyHandler);
  if (hResult != LIB3MF_OK) {
    DWORD nErrorMessage;
    LPCSTR pszErrorMessage;
    std::cerr << "could not create property handler: " << std::hex << hResult << std::endl;
    NMR::lib3mf_getlasterror(pMeshObject, &nErrorMessage, &pszErrorMessage);
    std::cerr << "error #" << std::hex << nErrorMessage << ": " << pszErrorMessage << std::endl;
    NMR::lib3mf_release(pMeshObject);
    return -1;
  }

  for(DWORD vid = 0; vid < points.size(); ++vid)
  {
    NMR::MODELMESHVERTEX pVertex;
    NMR::lib3mf_meshobject_getvertex(pMeshObject, vid, &pVertex);
    points[vid] =
        Point_3(pVertex.m_fPosition[0],
        pVertex.m_fPosition[1],
        pVertex.m_fPosition[2]);
  }
  for(DWORD pid = 0; pid < triangles.size(); ++pid)
  {
    NMR::MODELMESHTRIANGLE pTriangle;
    NMR::lib3mf_meshobject_gettriangle(pMeshObject, pid, &pTriangle);
    Polygon triangle(3);
    for(DWORD i = 0; i< 3; ++i)
      triangle[i] = pTriangle.m_nIndices[i];
    triangles[pid] = triangle;
    NMR::MODELMESH_TRIANGLECOLOR_SRGB pColor;
    NMR::lib3mf_propertyhandler_getcolor(pPropertyHandler, pid, &pColor);
    NMR::MODELMESHCOLOR_SRGB mColor = pColor.m_Colors[0];
    colors[pid]=CGAL::Color(mColor.m_Red, mColor.m_Green,
                            mColor.m_Blue, mColor.m_Alpha);
  }
  return true;
}

template<typename PointRange,
         typename PolygonRange,
         typename ColorRange>
bool extract_polylines (NMR::PLib3MFModelMeshObject *pMeshObject,
           PointRange& points,
           PolygonRange&,
           ColorRange& colors,
           std::string& name) {
  typedef typename PointRange::value_type Point_3;
  typedef typename PolygonRange::value_type Polygon;

  HRESULT hResult;
  DWORD nNeededChars;
  std::vector<char> pBuffer;
  // Retrieve Mesh Name Length
  hResult = NMR::lib3mf_object_getnameutf8(pMeshObject, NULL, 0, &nNeededChars);
  if (hResult != LIB3MF_OK)
  {
    points.resize(0);
    std::cerr<<"Error during name extraction.";
    return false;
  }

  // Retrieve Mesh Name
  if (nNeededChars > 0) {
    pBuffer.resize(nNeededChars + 1);
    hResult = NMR::lib3mf_object_getnameutf8(pMeshObject, &pBuffer[0], nNeededChars + 1, NULL);
    pBuffer[nNeededChars] = 0;
    std::string temp(&pBuffer[0]);
    if(temp.find("_cgal_pl")== std::string::npos) //ignore not polylines
    {
      points.resize(0);
      return false;
    }
    name = std::string(&pBuffer[0]);
  }
  else
  {
    points.resize(0);
    return false;
  }

  NMR::PLib3MFPropertyHandler * pPropertyHandler;
  hResult = NMR::lib3mf_meshobject_createpropertyhandler(pMeshObject, &pPropertyHandler);
  if (hResult != LIB3MF_OK) {
    DWORD nErrorMessage;
    LPCSTR pszErrorMessage;
    std::cerr << "could not create property handler: " << std::hex << hResult << std::endl;
    NMR::lib3mf_getlasterror(pMeshObject, &nErrorMessage, &pszErrorMessage);
    std::cerr << "error #" << std::hex << nErrorMessage << ": " << pszErrorMessage << std::endl;
    NMR::lib3mf_release(pMeshObject);
    return -1;
  }

  points.resize(points.size()-3);//ignore dummy_vertices
  for(DWORD vid = 0; vid < points.size(); ++vid)
  {
    NMR::MODELMESHVERTEX pVertex;
    NMR::lib3mf_meshobject_getvertex(pMeshObject, vid+3, &pVertex);
    points[vid] =
        Point_3(pVertex.m_fPosition[0],
        pVertex.m_fPosition[1],
        pVertex.m_fPosition[2]);
  }

  NMR::MODELMESH_TRIANGLECOLOR_SRGB pColor;
  NMR::lib3mf_propertyhandler_getcolor(pPropertyHandler, 0, &pColor);//get color of the dummy triangle
  NMR::MODELMESHCOLOR_SRGB mColor = pColor.m_Colors[0];
  colors[0]=CGAL::Color(mColor.m_Red, mColor.m_Green,
                          mColor.m_Blue, mColor.m_Alpha);
  return true;
}

template<typename PointRange,
         typename PolygonRange,
         typename ColorRange>
bool extract_point_clouds (NMR::PLib3MFModelMeshObject *pMeshObject,
           PointRange& points,
           PolygonRange&,
           ColorRange& colors,
           std::string& name) {
  typedef typename PointRange::value_type Point_3;

  HRESULT hResult;
  DWORD nNeededChars;
  std::vector<char> pBuffer;
  // Retrieve Mesh Name Length
  hResult = NMR::lib3mf_object_getnameutf8(pMeshObject, NULL, 0, &nNeededChars);
  if (hResult != LIB3MF_OK)
  {
    std::cerr<<"Error during name extraction.";
    points.resize(0);
    return false;
  }

  // Retrieve Mesh Name
  if (nNeededChars > 0) {
    pBuffer.resize(nNeededChars + 1);
    hResult = NMR::lib3mf_object_getnameutf8(pMeshObject, &pBuffer[0], nNeededChars + 1, NULL);
    pBuffer[nNeededChars] = 0;
    std::string temp(&pBuffer[0]);
    if(temp.find("_cgal_pc")== std::string::npos) //ignore not point_cloud
    {
      points.resize(0);
      return false;
    }
    name = std::string(&pBuffer[0]);
  }
  else{
    points.resize(0);
    return false;
  }

  NMR::PLib3MFPropertyHandler * pPropertyHandler;
  hResult = NMR::lib3mf_meshobject_createpropertyhandler(pMeshObject, &pPropertyHandler);
  if (hResult != LIB3MF_OK) {
    DWORD nErrorMessage;
    LPCSTR pszErrorMessage;
    std::cerr << "could not create property handler: " << std::hex << hResult << std::endl;
    NMR::lib3mf_getlasterror(pMeshObject, &nErrorMessage, &pszErrorMessage);
    std::cerr << "error #" << std::hex << nErrorMessage << ": " << pszErrorMessage << std::endl;
    NMR::lib3mf_release(pMeshObject);
    return -1;
  }

  points.resize(points.size()-3);//ignore dummy_vertices
  for(DWORD vid = 0; vid < points.size(); ++vid)
  {
    NMR::MODELMESHVERTEX pVertex;
    NMR::lib3mf_meshobject_getvertex(pMeshObject, vid+3, &pVertex);
    points[vid] =
        Point_3(pVertex.m_fPosition[0],
        pVertex.m_fPosition[1],
        pVertex.m_fPosition[2]);
  }

  NMR::MODELMESH_TRIANGLECOLOR_SRGB pColor;
  NMR::lib3mf_propertyhandler_getcolor(pPropertyHandler, 0, &pColor);//get color of the dummy triangle
  NMR::MODELMESHCOLOR_SRGB mColor = pColor.m_Colors[0];
  colors[0]=CGAL::Color(mColor.m_Red, mColor.m_Green,
                          mColor.m_Blue, mColor.m_Alpha);
  return true;
}

template<typename PointRanges, typename PolygonRanges, typename ColorRanges,
         typename PointRange, typename PolygonRange, typename ColorRange>
int read_from_3mf(const std::string& file_name, PointRanges& all_points,
                        PolygonRanges& all_polygons, ColorRanges& all_colors,
                  std::vector<std::string>& names,
                  std::function<bool(NMR::PLib3MFModelMeshObject*,
                                     PointRange&,
                                     PolygonRange&,
                                     ColorRange&,
                                     std::string&)> func
                  )
{
  typedef typename PointRange::value_type Point_3;
  typedef typename PolygonRange::value_type Polygon;
  DWORD nInterfaceVersionMajor, nInterfaceVersionMinor, nInterfaceVersionMicro, nbVertices, nbPolygons;
  HRESULT hResult;
  NMR::PLib3MFModel * pModel;
  NMR::PLib3MFModelReader * pReader;
  // Extract Extension of filename
  std::string sReaderName("3mf");

  hResult = NMR::lib3mf_getinterfaceversion(&nInterfaceVersionMajor, &nInterfaceVersionMinor, &nInterfaceVersionMicro);
  if (hResult != LIB3MF_OK) {
    std::cerr << "could not get 3MF Library version: " << std::hex << hResult << std::endl;
    return -1;
  }

  // Create Model Instance
  hResult = NMR::lib3mf_createmodel(&pModel);
  if (hResult != LIB3MF_OK) {
    std::cerr << "could not create model: " << std::hex << hResult << std::endl;
    return -1;
  }

  // Create Model Reader
  hResult = NMR::lib3mf_model_queryreader(pModel, sReaderName.c_str(), &pReader);
  if (hResult != LIB3MF_OK) {
    std::cerr << "could not create model reader: " << std::hex << hResult << std::endl;
    NMR::lib3mf_release(pModel);
    return -1;
  }

  // Import Model from File
  hResult = NMR::lib3mf_reader_readfromfileutf8(pReader, file_name.c_str());
  if (hResult != LIB3MF_OK) {
    std::cerr << "could not parse file: " << std::hex << hResult << std::endl;
    NMR::lib3mf_release(pReader);
    NMR::lib3mf_release(pModel);
    return -1;
  }

  // Release Model Reader
  NMR::lib3mf_release(pReader);

  //Iterate Model

  BOOL pbHasNext;
  NMR::PLib3MFModelResourceIterator * pResourceIterator;

  hResult = NMR::lib3mf_model_getobjects(pModel, &pResourceIterator);
  if (hResult != LIB3MF_OK) {
    std::cerr << "could not get object: " << std::hex << hResult << std::endl;
    NMR::lib3mf_release(pModel);
    return -1;
  }
  hResult = NMR::lib3mf_resourceiterator_movenext(pResourceIterator, &pbHasNext);
  if (hResult != LIB3MF_OK) {
    std::cerr << "could not get next object: " << std::hex << hResult << std::endl;
    NMR::lib3mf_release(pResourceIterator);
    NMR::lib3mf_release(pModel);
    return -1;
  }
  while (pbHasNext) {
    NMR::PLib3MFModelResource * pResource;
    NMR::PLib3MFModelMeshObject * pMeshObject;
    NMR::PLib3MFModelComponentsObject * pComponentsObject;
    NMR::ModelResourceID ResourceID;

    // get current resource
    hResult = NMR::lib3mf_resourceiterator_getcurrent(pResourceIterator, &pResource);
    if (hResult != LIB3MF_OK) {
      std::cerr << "could not get resource: " << std::hex << hResult << std::endl;
      NMR::lib3mf_release(pResourceIterator);
      NMR::lib3mf_release(pModel);
      return -1;
    }

    // get resource ID
    hResult = NMR::lib3mf_resource_getresourceid(pResource, &ResourceID);
    if (hResult != LIB3MF_OK) {
      std::cerr << "could not get resource id: " << std::hex << hResult << std::endl;
      NMR::lib3mf_release(pResource);
      NMR::lib3mf_release(pResourceIterator);
      NMR::lib3mf_release(pModel);
      return -1;
    }
    // Query mesh interface
    BOOL bIsMeshObject;
    hResult = NMR::lib3mf_object_ismeshobject(pResource, &bIsMeshObject);
    if ((hResult == LIB3MF_OK) && (bIsMeshObject)) {
      pMeshObject = pResource;
      NMR::lib3mf_meshobject_getvertexcount(pMeshObject, &nbVertices);
      NMR::lib3mf_meshobject_gettrianglecount(pMeshObject, &nbPolygons);
      PointRange points (nbVertices);
      PolygonRange triangles(nbPolygons);
      ColorRange colors(nbPolygons);
      std::string name;

      if(func(pMeshObject, points, triangles, colors, name)){
        all_points.push_back(points);
        all_polygons.push_back(triangles);
        all_colors.push_back(colors);
        names.push_back(name);
      }
    }
    // free instances
    NMR::lib3mf_release(pResource);
    hResult = NMR::lib3mf_resourceiterator_movenext(pResourceIterator, &pbHasNext);
    if (hResult != LIB3MF_OK) {
      std::cerr << "could not get next object: " << std::hex << hResult << std::endl;
      return -1;
    }
  }
  return all_points.size();
}

/*!
 * \brief read_soups_from_3mf extracts ranges of points and triangles from the
 * `MeshObject`s contained in `file_name`
 * \tparam PointRanges a model of the concepts `RandomAccessContainer` and
 *  `BackInsertionSequence` whose `value type` is
 * a model of the concepts `RandomAccessContainer` and `BackInsertionSequence`
 *  whose `value type` is the point type.
 * \tparam PolygonRanges a model of the concept `RandomAccessContainer` whose
 *  `value_type` is a model of the concept `RandomAccessContainer`
 * whose `value_type` is a model of the concept `RandomAccessContainer` whose
 *  `value_type` is std::size_t.
 * \tparam ColorRanges a model of the concepts `RandomAccessContainer` and
 *  `BackInsertionSequence` whose `value type` is
 * a model of the concepts `RandomAccessContainer` and `BackInsertionSequence`
 *  whose `value type` is `CGAL::Color`.
 * \param file_name the name of the 3mf file to read.
 * \param all_points a `PointRanges` that will contain the points of the meshes
 *  in `file_name`.
 * Each of these meshes will add a range of its points.
 * \param all_polygons a `PolygonRanges` that will contain the triangles of the
 *  meshes in `file_name`.
 * Each of these meshes will add a range of its triangles. A `triangle` of
 *  all_polygons[i] contains the indices of its points in all_points[i].
 * \param all_colors will contain the color of each triangle for each soup.
 * \param names will contain the name of each mesh in `file_name` if any.
 *  If the i'th mesh has no name, it will be called "Unknown Mesh" in names.
 * \return the number of meshes processed in `file_name`.
 */
template<typename PointRanges, typename PolygonRanges, typename ColorRanges>
int read_soups_from_3mf(const std::string& file_name, PointRanges& all_points,
                        PolygonRanges& all_polygons, ColorRanges& all_colors,
                        std::vector<std::string>& names
                  )
{
  typedef typename PointRanges::value_type PointRange;
  typedef typename PolygonRanges::value_type PolygonRange;
  typedef typename ColorRanges::value_type ColorRange;
  return read_from_3mf<PointRanges,PolygonRanges,ColorRanges,
      PointRange, PolygonRange, ColorRange>
      (file_name, all_points, all_polygons,
       all_colors, names, extract_soups<PointRange, PolygonRange, ColorRange>);
}


template<typename PointRanges, typename ColorRanges>
int read_polylines_from_3mf(const std::string& file_name,
                            PointRanges& all_points,
                            ColorRanges& all_colors,
                            std::vector<std::string>& names
                  )
{
  typedef typename PointRanges::value_type PointRange;
  typedef std::vector<std::size_t> Polygon;
  typedef std::vector<Polygon> PolygonRange;
  typedef std::vector<CGAL::Color> ColorRange;
  std::vector<PolygonRange> all_polygons;
  return read_from_3mf<PointRanges,std::vector<PolygonRange>,
      std::vector<ColorRange>, PointRange, PolygonRange, ColorRange>
      (file_name, all_points, all_polygons, all_colors, names,
       extract_polylines<PointRange, PolygonRange, ColorRange>);
}


template<typename PointRanges, typename ColorRanges>
int read_point_clouds_from_3mf(const std::string& file_name,
                               PointRanges& all_points,
                               ColorRanges& all_colors,
                               std::vector<std::string>& names
                  )
{
  typedef typename PointRanges::value_type PointRange;
  typedef std::vector<std::size_t> Polygon;
  typedef std::vector<Polygon> PolygonRange;
  typedef std::vector<CGAL::Color> ColorRange;
  std::vector<PolygonRange> all_polygons;
  return read_from_3mf<PointRanges,std::vector<PolygonRange>,
      std::vector<ColorRange>, PointRange, PolygonRange, ColorRange>
      (file_name, all_points, all_polygons, all_colors, names,
       extract_point_clouds<PointRange, PolygonRange, ColorRange>);
}
}//end CGAL

#endif // READ_3MF_H

