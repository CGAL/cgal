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
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <functional>
#include <CGAL/IO/Color.h>
#include <CGAL/Kernel_traits.h>
#include <Model/COM/NMR_DLLInterfaces.h>
namespace CGAL{

namespace transform_nmr_internal{
NMR::MODELTRANSFORM initMatrix()
{
  NMR::MODELTRANSFORM mMatrix;
  int i, j;
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 3; j++) {
      mMatrix.m_fFields[j][i] = (i == j) ? 1.0f : 0.0f;
    }
  }

  return mMatrix;
}
}//end transform_nmr_internal

template<typename PointRange,
         typename PolygonRange,
         typename ColorRange>
bool extract_soups (NMR::PLib3MFModelMeshObject *pMeshObject,
                    const NMR::MODELTRANSFORM& transform,
                    PointRange& points,
                    PolygonRange& triangles,
                    ColorRange& colors,
                    std::string& name) {
  typedef typename PointRange::value_type Point_3;
  typedef typename PolygonRange::value_type Polygon;
  typedef typename Kernel_traits<Point_3>::Kernel Kernel;
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
    name = std::string(&pBuffer[0]);
  }
  else
    name = std::string("Unknown Mesh");

  typename Kernel::Aff_transformation_3 t(
        transform.m_fFields[0][0], transform.m_fFields[0][1], transform.m_fFields[0][2], transform.m_fFields[0][3],
        transform.m_fFields[1][0], transform.m_fFields[1][1], transform.m_fFields[1][2], transform.m_fFields[1][3],
        transform.m_fFields[2][0], transform.m_fFields[2][1], transform.m_fFields[2][2], transform.m_fFields[2][3]
      );

  NMR::PLib3MFPropertyHandler * pPropertyHandler;
  hResult = NMR::lib3mf_meshobject_createpropertyhandler(pMeshObject, &pPropertyHandler);
  if (hResult != LIB3MF_OK) {
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


template<typename PointRanges, typename PolygonRanges, typename ColorRanges,
         typename PointRange, typename PolygonRange, typename ColorRange>
int read_from_3mf(const std::string& file_name, PointRanges& all_points,
                  PolygonRanges& all_polygons, ColorRanges& all_colors,
                  std::vector<std::string>& names,
                  std::function<bool(NMR::PLib3MFModelMeshObject*,
                                     const NMR::MODELTRANSFORM&,
                                     PointRange&,
                                     PolygonRange&,
                                     ColorRange&,
                                     std::string&)> func
                  )
{
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

  /**************************************************
   **** Iterate Resources To Find Meshes ************
   **************************************************/

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
    if (hResult == LIB3MF_OK)
    {
      if(bIsMeshObject) {
        //skip it. Only get it through the components and buildItems.
      }
      else {
        BOOL bIsComponentsObject;
        hResult = NMR::lib3mf_object_iscomponentsobject(pResource, &bIsComponentsObject);
        if ((hResult == LIB3MF_OK) && (bIsComponentsObject)) {
          pComponentsObject = (NMR::PLib3MFModelComponentsObject*)pResource;
          DWORD nComponentCount;
          hResult = NMR::lib3mf_componentsobject_getcomponentcount(pComponentsObject, &nComponentCount);
          if (hResult != LIB3MF_OK)
            return -1;
          //for each component
          DWORD nIndex;
          for (nIndex = 0; nIndex < nComponentCount; ++nIndex) {
            NMR::PLib3MFModelResource * compResource;
            NMR::PLib3MFModelComponent * pComponent;
            hResult = NMR::lib3mf_componentsobject_getcomponent(pComponentsObject, nIndex, &pComponent);
            if (hResult != LIB3MF_OK) {
              return -1;
            }

            hResult = NMR::lib3mf_component_getobjectresource(pComponent, &compResource);
            if (hResult != LIB3MF_OK) {
              NMR::lib3mf_release(pComponent);
              return -1;
            }
            hResult = NMR::lib3mf_object_ismeshobject(compResource, &bIsMeshObject);
            if (hResult == LIB3MF_OK)
            {
              if(bIsMeshObject) {
                BOOL bHasTransform;
                NMR::MODELTRANSFORM Transform;
                hResult = NMR::lib3mf_component_hastransform(pComponent, &bHasTransform);
                if (hResult != LIB3MF_OK) {
                  NMR::lib3mf_release(pComponent);
                  return -1;
                }

                if (bHasTransform) {
                  // Retrieve Transform
                  hResult = NMR::lib3mf_component_gettransform(pComponent, &Transform);
                  if (hResult != LIB3MF_OK) {
                    NMR::lib3mf_release(pComponent);
                    return -1;
                  }
                }
                else {
                  Transform = transform_nmr_internal::initMatrix();
                }
                pMeshObject = compResource;
                NMR::lib3mf_meshobject_getvertexcount(pMeshObject, &nbVertices);
                NMR::lib3mf_meshobject_gettrianglecount(pMeshObject, &nbPolygons);
                PointRange points (nbVertices);
                PolygonRange triangles(nbPolygons);
                ColorRange colors(nbPolygons);
                std::string name;

                if(func(pMeshObject, Transform, points, triangles, colors, name)){
                  all_points.push_back(points);
                  all_polygons.push_back(triangles);
                  all_colors.push_back(colors);
                  names.push_back(name);
                }
              }
            }
          }
          //end component
        }
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


  /********************************************************
   **** Iterate Build items To Find Transformed Meshes ****
   ********************************************************/

  // Iterate through all the Build items
  NMR::PLib3MFModelBuildItemIterator * pBuildItemIterator;
  hResult = NMR::lib3mf_model_getbuilditems(pModel, &pBuildItemIterator);
  if (hResult != LIB3MF_OK) {
    std::cout << "could not get build items: " << std::hex << hResult << std::endl;
    NMR::lib3mf_release(pBuildItemIterator);
    NMR::lib3mf_release(pModel);
    return -1;
  }

  hResult = NMR::lib3mf_builditemiterator_movenext(pBuildItemIterator, &pbHasNext);
  if (hResult != LIB3MF_OK) {
    std::cout << "could not get next build item: " << std::hex << hResult << std::endl;
    NMR::lib3mf_release(pBuildItemIterator);
    NMR::lib3mf_release(pModel);
    return -1;
  }

  while (pbHasNext) {
    NMR::PLib3MFModelMeshObject * pMeshObject;
    NMR::MODELTRANSFORM Transform;
    NMR::PLib3MFModelBuildItem * pBuildItem;
    // Retrieve Build Item
    hResult = NMR::lib3mf_builditemiterator_getcurrent(pBuildItemIterator, &pBuildItem);
    if (hResult != LIB3MF_OK) {
      std::cout << "could not get build item: " << std::hex << hResult << std::endl;
      NMR::lib3mf_release(pBuildItemIterator);
      NMR::lib3mf_release(pModel);
      return -1;
    }

    // Retrieve Resource
    NMR::PLib3MFModelObjectResource * pObjectResource;
    hResult = NMR::lib3mf_builditem_getobjectresource(pBuildItem, &pObjectResource);
    if (hResult != LIB3MF_OK) {
      std::cout << "could not get build item resource: " << std::hex << hResult << std::endl;
      NMR::lib3mf_release(pBuildItem);
      NMR::lib3mf_release(pBuildItemIterator);
      NMR::lib3mf_release(pModel);
      return -1;
    }

    BOOL bIsMeshObject;
    hResult = NMR::lib3mf_object_ismeshobject(pObjectResource, &bIsMeshObject);
    if (hResult == LIB3MF_OK)
    {
      if(bIsMeshObject) {
        pMeshObject = pObjectResource;
        NMR::lib3mf_meshobject_getvertexcount(pMeshObject, &nbVertices);
        NMR::lib3mf_meshobject_gettrianglecount(pMeshObject, &nbPolygons);
        PointRange points (nbVertices);
        PolygonRange triangles(nbPolygons);
        ColorRange colors(nbPolygons);
        std::string name;


        // Check Object Transform
        BOOL bHasTransform;
        hResult = NMR::lib3mf_builditem_hasobjecttransform(pBuildItem, &bHasTransform);
        if (hResult != LIB3MF_OK) {
          NMR::lib3mf_release(pBuildItem);
          NMR::lib3mf_release(pBuildItemIterator);
          NMR::lib3mf_release(pModel);
          std::cerr << "could not check object transform: " << std::hex << hResult << std::endl;
          return -1;
        }

        if (bHasTransform) {
          // Retrieve Transform
          hResult = NMR::lib3mf_builditem_getobjecttransform(pBuildItem, &Transform);
          if (hResult != LIB3MF_OK) {
            NMR::lib3mf_release(pBuildItem);
            NMR::lib3mf_release(pBuildItemIterator);
            NMR::lib3mf_release(pModel);
            std::cerr << "could not get object transform: " << std::hex << hResult << std::endl;
            return -1;
          }
        }
        else {
          Transform = transform_nmr_internal::initMatrix();
        }

        if(func(pMeshObject, Transform, points, triangles, colors, name)){
          all_points.push_back(points);
          all_polygons.push_back(triangles);
          all_colors.push_back(colors);
          names.push_back(name);
        }
      }
    }
    // Release Object Resource ID
    NMR::lib3mf_release(pObjectResource);
    // Release Build Item
    NMR::lib3mf_release(pBuildItem);

    // Move to next Item
    hResult = NMR::lib3mf_builditemiterator_movenext(pBuildItemIterator, &pbHasNext);
    if (hResult != LIB3MF_OK) {
      std::cerr << "could not get next build item: " << std::hex << hResult << std::endl;
      return -1;
    }
  }
  // Release Build Item Iterator
  NMR::lib3mf_release(pBuildItemIterator);
  return all_points.size();
}

/*!
 * \brief extracts ranges of points and triangles from a 3mf file.
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
 *  `all_polygons[i]` contains the indices of its points in `all_points[i]`.
 * \param all_colors will contain the color of each triangle for each soup.
 * \param names will contain the name of each mesh in `file_name` if any.
 *  If the i'th mesh has no name, it will be called "Unknown Mesh" in names.
 * \return the number of soups read.
 *
 * \attention Only versions inferior to 2.0 of lib3mf are supported.
 */
template<typename PointRanges, typename PolygonRanges, typename ColorRanges>
int read_triangle_soups_from_3mf(const std::string& file_name, PointRanges& all_points,
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

}//end CGAL

#endif // CGAL_IO_READ_3MF_H

