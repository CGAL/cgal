// Copyright (c) 2015-2020 GeometryFactory
//
// This file is part of CGAL (www.cgal.org);
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Maxime Gimeno

#ifndef CGAL_IO_3MF_H
#define CGAL_IO_3MF_H

#include <CGAL/IO/3MF/read_3mf.h>
#include <CGAL/IO/3MF/write_3mf.h>
#include <CGAL/IO/Color.h>
#include <CGAL/IO/helpers.h>

#include <CGAL/boost/graph/iterator.h>

#include <boost/range/value_type.hpp>

#include <functional>
#include <iostream>
#include <string>
#include <vector>
#include <type_traits>

#if defined(CGAL_LINKED_WITH_LIB3MF) || defined(DOXYGEN_RUNNING)

namespace CGAL {

namespace IO {

using namespace Lib3MF;

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Read

/// \cond SKIP_IN_MANUAL

template<typename PointRanges, typename TriangleRanges, typename ColorRanges,
         typename PointRange, typename TriangleRange, typename ColorRange>
bool read_3MF(const std::string& fname,
              PointRanges& all_points,
              TriangleRanges& all_triangles,
              ColorRanges& all_colors,
              std::vector<std::string>& names,
              std::function<bool(PMeshObject,
                                 const sTransform&,
                                 PointRange&,
                                 TriangleRange&,
                                 ColorRange&,
                                 std::string&)> func)
{
  Lib3MF_uint32 nInterfaceVersionMajor, nInterfaceVersionMinor, nInterfaceVersionMicro, nbVertices, nbTriangles;


  PReader pReader;

  PWrapper wrapper = CWrapper::loadLibrary();
  PModel pModel = wrapper->CreateModel();

  // Extract Extension of filename
  std::string sReaderName("3mf");
  /*
  hResult = NMR::lib3mf_getinterfaceversion(&nInterfaceVersionMajor, &nInterfaceVersionMinor, &nInterfaceVersionMicro);
  if(hResult != LIB3MF_OK)
  {
    std::cerr << "could not get 3MF Library version: " << std::hex << hResult << std::endl;
    return false;
  }
  */

  // Create Model Reader
  pReader = pModel->QueryReader(sReaderName.c_str());

  // Import Model from File
  pReader->ReadFromFile(fname.c_str());

  //Iterate Model
  PObjectIterator pObjectIterator = pModel->GetObjects();


  /**************************************************
   **** Iterate Resources To Find Meshes ************
   **************************************************/

  while(pObjectIterator->MoveNext())
  {
    PMeshObject pMeshObject;
    PComponentsObject pComponentsObject;
    Lib3MF_uint32 ResourceID;

    // get current resource
    PObject pObject = pObjectIterator->GetCurrentObject();

    // Query mesh interface
    bool bIsMeshObject = pObject->IsMeshObject();

      if(bIsMeshObject)
      {
        //skip it. Only get it through the components and buildItems.
      }
      else
      {
        bool bIsComponentsObject = pObject->IsComponentsObject();
        if(bIsComponentsObject)
        {
          pComponentsObject = std::static_pointer_cast<CComponentsObject>(pObject);
          Lib3MF_uint32 nComponentCount = pComponentsObject->GetComponentCount();

          //for each component
          Lib3MF_uint32  nIndex;
          for(nIndex = 0; nIndex < nComponentCount; ++nIndex)
          {
            PComponent pComponent = pComponentsObject->GetComponent(nIndex);
            PObject pObject = pComponent->GetObjectResource();

            bIsMeshObject = pObject->IsMeshObject();
            if(bIsMeshObject)
              {
                bool bHasTransform = pComponent->HasTransform();
                sTransform Transform;
                if(bHasTransform)
                {
                  // Retrieve Transform
                  Transform = pComponent->GetTransform();

                }
                else
                {
                  Transform = wrapper->GetIdentityTransform();
                }

                pMeshObject = std::static_pointer_cast<CMeshObject>(pObject);
                nbVertices = pMeshObject->GetVertexCount();
                nbTriangles = pMeshObject->GetTriangleCount();
                PointRange points (nbVertices);
                TriangleRange triangles(nbTriangles);
                ColorRange colors(nbTriangles);
                std::string name;

                if(func(pMeshObject, Transform, points, triangles, colors, name))
                {
                  all_points.push_back(points);
                  all_triangles.push_back(triangles);
                  all_colors.push_back(colors);
                  names.push_back(name);
                }
              }

          }
          //end component
        }
      }
  }

  /********************************************************
   **** Iterate Build items To Find Transformed Meshes ****
   ********************************************************/

  // Iterate through all the Build items
  PBuildItemIterator pBuildItemIterator = pModel->GetBuildItems();


  while(pBuildItemIterator->MoveNext())
  {
    PMeshObject pMeshObject;
    sTransform Transform;
    PBuildItem pBuildItem = pBuildItemIterator->GetCurrent();


    // Retrieve Resource
    PObject pObject = pBuildItem->GetObjectResource();


      if(pObject->IsMeshObject()){
        pMeshObject = std::static_pointer_cast<CMeshObject>(pObject);
        nbVertices = pMeshObject->GetVertexCount();
        nbTriangles = pMeshObject->GetTriangleCount();
        PointRange points (nbVertices);
        TriangleRange triangles(nbTriangles);
        ColorRange colors(nbTriangles);
        std::string name;

        if(pBuildItem->HasObjectTransform())
        {
          Transform = pBuildItem->GetObjectTransform();
        }
        else
        {
          Transform = wrapper->GetIdentityTransform();
        }

        if(func(pMeshObject, Transform, points, triangles, colors, name)){
          all_points.push_back(points);
          all_triangles.push_back(triangles);
          all_colors.push_back(colors);
          names.push_back(name);
        }
      }
  }
  return true;
}

/// \endcond

/*!
 * \ingroup PkgStreamSupportIoFuncs3MF
 *
 * \brief reads ranges of points and triangles from an input file, using the \ref IOStream3MF.
 *
 * \attention The ranges not cleared, and the data from the file are appended.
 *
 * \tparam PointRanges a model of the concepts `RandomAccessContainer` and
 *                     `BackInsertionSequence` whose `value_type` is
 *                     a model of the concepts `RandomAccessContainer` and `BackInsertionSequence`
 *                     whose `value_type` is the point type
 * \tparam TriangleRanges a model of the concepts `RandomAccessContainer` and `BackInsertionSequence`
 *                        whose `value_type` is a model of the concepts `RandomAccessContainer`
 *                        and `BackInsertionSequence` whose `value_type` is a model of the concepts
 *                       `RandomAccessContainer` and `BackInsertionSequence` whose
 *                        `value_type` is an unsigned integer type convertible to `std::size_t`
 * \tparam ColorRanges a model of the concepts `RandomAccessContainer` and
 *                     `BackInsertionSequence` whose `value_type` is
 *                     a model of the concepts `RandomAccessContainer` and `BackInsertionSequence`
 *                     whose `value_type` is `CGAL::IO::Color`
 *
 * \param fname the name of the 3mf file to read
 * \param all_points a `PointRanges` that will contain the points of the meshes in `fname`.
 *                   Each of these meshes will add a range of its points.
 * \param all_triangles a `TriangleRanges` that will contain the triangles of the meshes in `fname`.
 *                      Each of these meshes will add a range of its triangles. A `triangle` of
 *                      `all_triangles[i]` contains the indices of its points in `all_points[i]`.
 * \param all_colors will contain the color of each triangle for each soup.
 * \param names will contain the name of each mesh in `fname` if any.
 *             If the i-th mesh has no name, it will be called "Unknown Mesh" in `names`.
 *
 * \returns `true` if reading was successful, `false` otherwise.
 */
template<typename PointRanges, typename TriangleRanges, typename ColorRanges>
bool read_3MF(const std::string& fname,
             PointRanges& all_points,
             TriangleRanges& all_triangles,
             ColorRanges& all_colors,
             std::vector<std::string>& names)
{
  typedef typename boost::range_value<PointRanges>::type    PointRange;
  typedef typename boost::range_value<TriangleRanges>::type TriangleRange;
  typedef typename boost::range_value<ColorRanges>::type    ColorRange;

  return read_3MF<PointRanges, TriangleRanges,ColorRanges,
                  PointRange, TriangleRange, ColorRange>(fname, all_points, all_triangles, all_colors, names,
                                                         extract_soups<PointRange, TriangleRange, ColorRange>);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Write

/*!
 * \ingroup PkgStreamSupportIoFuncs3MF
 *
 * \brief writes the triangle soups contained in `all_points` and `all_triangles` into the file `fname`, using the \ref IOStream3MF.
 *
 * \tparam PointRanges a model of the concept `RandomAccessContainer` whose `value_type` is
 *                     a model of the concept `RandomAccessContainer` whose `value_type` is the point type
 * \tparam TriangleRanges a model of the concept `RandomAccessContainer` whose
 *                        `value_type` is a model of the concept `RandomAccessContainer`
 *                        whose `value_type` is a model of the concept `RandomAccessContainer` whose
 *                        `value_type` is an unsigned integer type convertible to `std::size_t`
 *
 * \param fname the name of the 3mf file to write
 * \param all_points a `PointRanges` that contains the points of the soups to write
 * \param all_triangles a `TriangleRanges` that contains the triangles of the soups in `fname`
 * \param names a range of `std::string` associating a name to each soup, which will appear in the output
 *
 * \return `true` if the writing is successful, `false` otherwise.
 */
template<typename PointRanges, typename TriangleRanges>
bool write_3MF(const std::string& fname,
               const PointRanges& all_points,
               const TriangleRanges& all_triangles,
               const std::vector<std::string>& names)
{
  PWrapper wrapper = CWrapper::loadLibrary();
  // Create Model Instance
  PModel pModel =  wrapper->CreateModel();;
  /*
  if(hResult != LIB3MF_OK)
  {
    std::cerr << "could not create model: " << std::hex << hResult << std::endl;
    return false;
  }
  */
  for(std::size_t id = 0; id < all_points.size(); ++id)
  {
    PMeshObject pMeshObject = pModel->AddMeshObject();
    std::string name;
    if(names.size() > id && ! names[id].empty())
    {
      name=names[id];
    }
    else
    {
      name = std::string("");
    }

    // AF: FIX this all triangles will have the default color
    std::vector<CGAL::IO::Color> colors(all_triangles[id].size());
    colors[0] = CGAL::IO::Color(255,0,0,255);
    IO::write_mesh_to_model(all_points[id], all_triangles[id], colors, name, pMeshObject, pModel, wrapper);
  }

  return IO::export_model_to_file(fname, pModel);
}

/// \cond SKIP_IN_MANUAL

// convenience
template<typename PointRange, typename TriangleRange>
bool write_3MF(const std::string& fname,
               const PointRange& points,
               const TriangleRange& triangles,
               const std::string& name)
{
  std::vector<PointRange> all_points(1, points);
  std::vector<PointRange> all_triangles(1, triangles);
  std::vector<std::string> names(1, name);

  return write_3MF(fname, all_points, all_triangles, names);
}

template<typename PointRange, typename TriangleRange>
bool write_3MF(const std::string& fname, const PointRange& points, const TriangleRange& triangles,
               std::enable_if_t<internal::is_Range<PointRange>::value>* = nullptr)
{
  return write_3MF(fname, points, triangles, "anonymous");
}

/// \endcond

} // namespace IO

} // namespace CGAL

#endif // defined(CGAL_LINKED_WITH_LIB3MF) || defined(DOXYGEN_RUNNING)

#endif // CGAL_IO_3MF_H
