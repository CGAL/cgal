// Copyright (c) 2015 GeometryFactory
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

#include <CGAL/boost/graph/iterator.h>
#include <CGAL/IO/Color.h>
#include <CGAL/IO/3MF/read_3mf.h>
#include <CGAL/IO/3MF/write_3mf.h>

#include <Model/COM/NMR_DLLInterfaces.h>

#include <functional>
#include <iostream>
#include <string>
#include <vector>

namespace CGAL {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read

/*!
 * \brief extracts ranges of points and triangles from a 3mf file.
 *
 * \tparam PointRanges a model of the concepts `RandomAccessContainer` and
 *  `BackInsertionSequence` whose `value type` is
 *  a model of the concepts `RandomAccessContainer` and `BackInsertionSequence`
 *  whose `value type` is the point type.
 * \tparam PolygonRanges a model of the concept `RandomAccessContainer` whose
 *  `value_type` is a model of the concept `RandomAccessContainer`
 *  whose `value_type` is a model of the concept `RandomAccessContainer` whose
 *  `value_type` is std::size_t.
 * \tparam ColorRanges a model of the concepts `RandomAccessContainer` and
 *  `BackInsertionSequence` whose `value type` is
 *  a model of the concepts `RandomAccessContainer` and `BackInsertionSequence`
 *  whose `value type` is `CGAL::Color`.
 *
 * \param file_name the name of the 3mf file to read.
 * \param all_points a `PointRanges` that will contain the points of the meshes in `file_name`.
 *  Each of these meshes will add a range of its points.
 * \param all_polygons a `PolygonRanges` that will contain the triangles of the meshes in `file_name`.
 *  Each of these meshes will add a range of its triangles. A `triangle` of
 *  `all_polygons[i]` contains the indices of its points in `all_points[i]`.
 * \param all_colors will contain the color of each triangle for each soup.
 * \param names will contain the name of each mesh in `file_name` if any.
 *  If the i'th mesh has no name, it will be called "Unknown Mesh" in names.
 *
 * \return the number of soups read.
 */
template<typename PointRanges, typename PolygonRanges, typename ColorRanges>
int read_triangle_soups_from_3mf(const std::string& file_name,
                                 PointRanges& all_points,
                                 PolygonRanges& all_polygons,
                                 ColorRanges& all_colors,
                                 std::vector<std::string>& names)
{
  typedef typename PointRanges::value_type PointRange;
  typedef typename PolygonRanges::value_type PolygonRange;
  typedef typename ColorRanges::value_type ColorRange;
  return read_from_3mf<PointRanges,PolygonRanges,ColorRanges,
                       PointRange, PolygonRange, ColorRange>
          (file_name, all_points, all_polygons, all_colors, names,
           extract_soups<PointRange, PolygonRange, ColorRange>);
}

template<typename PointRanges, typename ColorRanges>
int read_polylines_from_3mf(const std::string& file_name,
                            PointRanges& all_points,
                            ColorRanges& all_colors,
                            std::vector<std::string>& names)
{
  typedef typename PointRanges::value_type PointRange;
  typedef std::vector<std::size_t> Polygon;
  typedef std::vector<Polygon> PolygonRange;
  typedef std::vector<CGAL::Color> ColorRange;

  std::vector<PolygonRange> all_polygons;
  return read_from_3mf<PointRanges, std::vector<PolygonRange>,
                       std::vector<ColorRange>, PointRange, PolygonRange, ColorRange>
          (file_name, all_points, all_polygons, all_colors, names,
           extract_polylines<PointRange, PolygonRange, ColorRange>);
}

template<typename PointRanges, typename ColorRanges>
int read_point_clouds_from_3mf(const std::string& file_name,
                               PointRanges& all_points,
                               ColorRanges& all_colors,
                               std::vector<std::string>& names)
{
  typedef typename PointRanges::value_type PointRange;
  typedef std::vector<std::size_t> Polygon;
  typedef std::vector<Polygon> PolygonRange;
  typedef std::vector<CGAL::Color> ColorRange;

  std::vector<PolygonRange> all_polygons;
  return read_from_3mf<PointRanges, std::vector<PolygonRange>,
                       std::vector<ColorRange>, PointRange, PolygonRange, ColorRange>
          (file_name, all_points, all_polygons, all_colors, names,
           extract_point_clouds<PointRange, PolygonRange, ColorRange>);
}

template<typename PointRanges, typename PolygonRanges, typename ColorRanges,
         typename PointRange, typename PolygonRange, typename ColorRange>
int read_from_3mf(const std::string& file_name,
                  PointRanges& all_points,
                  PolygonRanges& all_polygons,
                  ColorRanges& all_colors,
                  std::vector<std::string>& names,
                  std::function<bool(NMR::PLib3MFModelMeshObject*,
                                     const NMR::MODELTRANSFORM&,
                                     PointRange&,
                                     PolygonRange&,
                                     ColorRange&,
                                     std::string&)> func)
{
  DWORD nInterfaceVersionMajor, nInterfaceVersionMinor, nInterfaceVersionMicro, nbVertices, nbPolygons;
  HRESULT hResult;
  NMR::PLib3MFModel * pModel;
  NMR::PLib3MFModelReader * pReader;

  // Extract Extension of filename
  std::string sReaderName("3mf");

  hResult = NMR::lib3mf_getinterfaceversion(&nInterfaceVersionMajor, &nInterfaceVersionMinor, &nInterfaceVersionMicro);
  if(hResult != LIB3MF_OK)
  {
    std::cerr << "could not get 3MF Library version: " << std::hex << hResult << std::endl;
    return -1;
  }

  // Create Model Instance
  hResult = NMR::lib3mf_createmodel(&pModel);
  if(hResult != LIB3MF_OK)
  {
    std::cerr << "could not create model: " << std::hex << hResult << std::endl;
    return -1;
  }

  // Create Model Reader
  hResult = NMR::lib3mf_model_queryreader(pModel, sReaderName.c_str(), &pReader);
  if(hResult != LIB3MF_OK)
  {
    std::cerr << "could not create model reader: " << std::hex << hResult << std::endl;
    NMR::lib3mf_release(pModel);
    return -1;
  }

  // Import Model from File
  hResult = NMR::lib3mf_reader_readfromfileutf8(pReader, file_name.c_str());
  if(hResult != LIB3MF_OK)
  {
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
  if(hResult != LIB3MF_OK)
  {
    std::cerr << "could not get object: " << std::hex << hResult << std::endl;
    NMR::lib3mf_release(pModel);
    return -1;
  }

  hResult = NMR::lib3mf_resourceiterator_movenext(pResourceIterator, &pbHasNext);
  if(hResult != LIB3MF_OK)
  {
    std::cerr << "could not get next object: " << std::hex << hResult << std::endl;
    NMR::lib3mf_release(pResourceIterator);
    NMR::lib3mf_release(pModel);
    return -1;
  }

  /**************************************************
   **** Iterate Resources To Find Meshes ************
   **************************************************/

  while (pbHasNext)
  {
    NMR::PLib3MFModelResource * pResource;
    NMR::PLib3MFModelMeshObject * pMeshObject;
    NMR::PLib3MFModelComponentsObject * pComponentsObject;
    NMR::ModelResourceID ResourceID;

    // get current resource
    hResult = NMR::lib3mf_resourceiterator_getcurrent(pResourceIterator, &pResource);
    if(hResult != LIB3MF_OK)
    {
      std::cerr << "could not get resource: " << std::hex << hResult << std::endl;
      NMR::lib3mf_release(pResourceIterator);
      NMR::lib3mf_release(pModel);
      return -1;
    }

    // get resource ID
    hResult = NMR::lib3mf_resource_getresourceid(pResource, &ResourceID);
    if(hResult != LIB3MF_OK)
    {
      std::cerr << "could not get resource id: " << std::hex << hResult << std::endl;
      NMR::lib3mf_release(pResource);
      NMR::lib3mf_release(pResourceIterator);
      NMR::lib3mf_release(pModel);
      return -1;
    }

    // Query mesh interface
    BOOL bIsMeshObject;
    hResult = NMR::lib3mf_object_ismeshobject(pResource, &bIsMeshObject);
    if(hResult == LIB3MF_OK)
    {
      if(bIsMeshObject)
      {
        //skip it. Only get it through the components and buildItems.
      }
      else
      {
        BOOL bIsComponentsObject;
        hResult = NMR::lib3mf_object_iscomponentsobject(pResource, &bIsComponentsObject);
        if((hResult == LIB3MF_OK) && (bIsComponentsObject))
        {
          pComponentsObject = (NMR::PLib3MFModelComponentsObject*)pResource;
          DWORD nComponentCount;
          hResult = NMR::lib3mf_componentsobject_getcomponentcount(pComponentsObject, &nComponentCount);
          if(hResult != LIB3MF_OK)
            return -1;

          //for each component
          DWORD nIndex;
          for(nIndex = 0; nIndex < nComponentCount; ++nIndex)
          {
            NMR::PLib3MFModelResource * compResource;
            NMR::PLib3MFModelComponent * pComponent;
            hResult = NMR::lib3mf_componentsobject_getcomponent(pComponentsObject, nIndex, &pComponent);
            if(hResult != LIB3MF_OK)
              return -1;

            hResult = NMR::lib3mf_component_getobjectresource(pComponent, &compResource);
            if(hResult != LIB3MF_OK)
            {
              NMR::lib3mf_release(pComponent);
              return -1;
            }

            hResult = NMR::lib3mf_object_ismeshobject(compResource, &bIsMeshObject);
            if(hResult == LIB3MF_OK)
            {
              if(bIsMeshObject)
              {
                BOOL bHasTransform;
                NMR::MODELTRANSFORM Transform;
                hResult = NMR::lib3mf_component_hastransform(pComponent, &bHasTransform);
                if(hResult != LIB3MF_OK)
                {
                  NMR::lib3mf_release(pComponent);
                  return -1;
                }

                if(bHasTransform)
                {
                  // Retrieve Transform
                  hResult = NMR::lib3mf_component_gettransform(pComponent, &Transform);
                  if(hResult != LIB3MF_OK)
                  {
                    NMR::lib3mf_release(pComponent);
                    return -1;
                  }
                }
                else
                {
                  Transform = transform_nmr_internal::initMatrix();
                }

                pMeshObject = compResource;
                NMR::lib3mf_meshobject_getvertexcount(pMeshObject, &nbVertices);
                NMR::lib3mf_meshobject_gettrianglecount(pMeshObject, &nbPolygons);
                PointRange points (nbVertices);
                PolygonRange triangles(nbPolygons);
                ColorRange colors(nbPolygons);
                std::string name;

                if(func(pMeshObject, Transform, points, triangles, colors, name))
                {
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
    if(hResult != LIB3MF_OK)
    {
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
  if(hResult != LIB3MF_OK)
  {
    std::cout << "could not get build items: " << std::hex << hResult << std::endl;
    NMR::lib3mf_release(pBuildItemIterator);
    NMR::lib3mf_release(pModel);
    return -1;
  }

  hResult = NMR::lib3mf_builditemiterator_movenext(pBuildItemIterator, &pbHasNext);
  if(hResult != LIB3MF_OK)
  {
    std::cout << "could not get next build item: " << std::hex << hResult << std::endl;
    NMR::lib3mf_release(pBuildItemIterator);
    NMR::lib3mf_release(pModel);
    return -1;
  }

  while (pbHasNext)
  {
    NMR::PLib3MFModelMeshObject * pMeshObject;
    NMR::MODELTRANSFORM Transform;
    NMR::PLib3MFModelBuildItem * pBuildItem;
    // Retrieve Build Item
    hResult = NMR::lib3mf_builditemiterator_getcurrent(pBuildItemIterator, &pBuildItem);
    if(hResult != LIB3MF_OK)
    {
      std::cout << "could not get build item: " << std::hex << hResult << std::endl;
      NMR::lib3mf_release(pBuildItemIterator);
      NMR::lib3mf_release(pModel);
      return -1;
    }

    // Retrieve Resource
    NMR::PLib3MFModelObjectResource * pObjectResource;
    hResult = NMR::lib3mf_builditem_getobjectresource(pBuildItem, &pObjectResource);
    if(hResult != LIB3MF_OK)
    {
      std::cout << "could not get build item resource: " << std::hex << hResult << std::endl;
      NMR::lib3mf_release(pBuildItem);
      NMR::lib3mf_release(pBuildItemIterator);
      NMR::lib3mf_release(pModel);
      return -1;
    }

    BOOL bIsMeshObject;
    hResult = NMR::lib3mf_object_ismeshobject(pObjectResource, &bIsMeshObject);
    if(hResult == LIB3MF_OK)
    {
      if(bIsMeshObject)
      {
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
        if(hResult != LIB3MF_OK)
        {
          NMR::lib3mf_release(pBuildItem);
          NMR::lib3mf_release(pBuildItemIterator);
          NMR::lib3mf_release(pModel);
          std::cerr << "could not check object transform: " << std::hex << hResult << std::endl;
          return -1;
        }

        if(bHasTransform)
        {
          // Retrieve Transform
          hResult = NMR::lib3mf_builditem_getobjecttransform(pBuildItem, &Transform);
          if(hResult != LIB3MF_OK)
          {
            NMR::lib3mf_release(pBuildItem);
            NMR::lib3mf_release(pBuildItemIterator);
            NMR::lib3mf_release(pModel);
            std::cerr << "could not get object transform: " << std::hex << hResult << std::endl;
            return -1;
          }
        }
        else
        {
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
    if(hResult != LIB3MF_OK)
    {
      std::cerr << "could not get next build item: " << std::hex << hResult << std::endl;
      return -1;
    }
  }

  // Release Build Item Iterator
  NMR::lib3mf_release(pBuildItemIterator);
  return all_points.size();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write

/*!
 * \brief writes the triangle soups contained in `all_points` and
 *  `all_polygons` into the 3mf file `file_name`.
 *
 * \tparam PointRanges a model of the concepts `RandomAccessContainer` and
 *  `BackInsertionSequence` whose `value type` is
 * a model of the concepts `RandomAccessContainer` and `BackInsertionSequence`
 *  whose `value type` is the point type.
 * \tparam PolygonRanges a model of the concept `RandomAccessContainer` whose
 *  `value_type` is a model of the concept `RandomAccessContainer`
 * whose `value_type` is a model of the concept `RandomAccessContainer` whose
 *  `value_type` is std::size_t.
 *
 * \param file_name the name of the 3mf file to write.
 * \param all_points a `PointRanges` that contains the points of the soups to write.
 * \param all_polygons a `PolygonRanges` that contains the triangles of the soups in `file_name`.
 * \param names will contains the name of each mesh in `file_name`.
 *
 * \return `true` if the writing is successful, `false` otherwise.
 */
template<typename PointRanges, typename PolygonRanges>
bool write_triangle_soups_to_3mf(const std::string& file_name,
                                 const PointRanges& all_points,
                                 const PolygonRanges& all_polygons,
                                 const std::vector<std::string>& names)
{
  // Create Model Instance
  NMR::PLib3MFModel * pModel;
  HRESULT hResult = NMR::lib3mf_createmodel(&pModel);
  if(hResult != LIB3MF_OK)
  {
    std::cerr << "could not create model: " << std::hex << hResult << std::endl;
    return false;
  }

  for(std::size_t id = 0; id < all_points.size(); ++id)
  {
    NMR::PLib3MFModelMeshObject* pMeshObject;
    std::string name;
    if(names.size() > id && ! names[id].empty())
    {
      name=names[id];
    }
    else
    {
      name = std::string("");
    }

    std::vector<CGAL::Color> colors(all_polygons[id].size());
    write_mesh_to_model(all_points[id], all_polygons[id], colors, name, &pMeshObject, pModel);
  }

  return export_model_to_file(file_name, pModel);
}


/*!
 * \brief writes the triangle meshes contained in `tms` into the 3mf file `file_name`.
 *
 * \tparam TriangleMeshRange a model of the concepts `RandomAccessContainer`
 * and `BackInsertionSequence` whose `value type` is
 * a model of the concepts `FaceListGraph` and `HalfedgeListGraph`
 * that has only triangle faces.
 *
 * \param file_name the name of the 3mf file to write.
 * \param tms a `TriangleMeshRange` that contains the meshes
 *  to write. An internal property map for `CGAL::vertex_point_t`
 * must be available for each mesh.
 * \param names will contains the name of each mesh in `file_name`.
 *
 * \return `true` if the writing is successful, `false` otherwise.
 */
template<typename TriangleMeshRange>
bool write_triangle_meshes_to_3mf(const std::string& file_name,
                                  const TriangleMeshRange& tms,
                                  const std::vector<std::string>& names)
{
  typedef typename TriangleMeshRange::value_type                          Mesh;
  typedef typename boost::property_map<Mesh, boost::vertex_point_t>::type VPMap;
  typedef typename boost::property_traits<VPMap>::value_type              Point_3;

  typedef boost::graph_traits<Mesh>::vertex_descriptor                    vertex_descriptor;
  typedef boost::graph_traits<Mesh>::face_descriptor                      face_descriptor;

  typedef std::vector<std::size_t>                                        Polygon;
  typedef std::vector<Polygon>                                            PolygonRange;
  typedef std::vector<Point_3>                                            PointRange;

  std::vector<PointRange> all_points;
  std::vector<PolygonRange> all_polygons;

  for(const auto& tm : tms)
  {
    PointRange points;
    points.reserve(num_vertices(tm));
    PolygonRange triangles;
    triangles.reserve(num_faces(tm));

    VPMap vpm = get(boost::vertex_point, tm);
    std::unordered_map<typename boost::graph_traits<Mesh>::vertex_descriptor, std::size_t> vertex_id_map;

    std::size_t i = 0;
    for(vertex_descriptor v : vertices(tm))
    {
      points.push_back(get(vpm, v));
      vertex_id_map[v] = i++;
    }

    all_points.push_back(points);
    for(face_descriptor f : faces(tm))
    {
      Polygon triangle;
      for(vertex_descriptor vert : CGAL::vertices_around_face(halfedge(f, tm), tm))
        triangle.push_back(vertex_id_map[vert]);

      triangles.push_back(triangle);
    }

    all_polygons.push_back(triangles);
  }

  return write_triangle_soups_to_3mf(file_name, all_points, all_polygons, names);
}

} // namespace CGAL

#endif // CGAL_IO_3MF_H
