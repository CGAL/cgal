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
#include <iostream>
#include <vector>
#include <string>

#include <CGAL/IO/Color.h>

#include <Model/COM/NMR_DLLInterfaces.h>

namespace CGAL{
namespace tmf_internal{
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
} //end internal

bool add_build_item(NMR::PLib3MFModel * pModel,
                    NMR::PLib3MFModelMeshObject* pMeshObject)
{
  HRESULT hResult;
  DWORD nErrorMessage;
  LPCSTR pszErrorMessage;
  // Add Build Item for Mesh
  NMR::PLib3MFModelBuildItem * pBuildItem;
  hResult = NMR::lib3mf_model_addbuilditem(pModel, pMeshObject, NULL, &pBuildItem);
  if (hResult != LIB3MF_OK) {
    std::cerr << "could not create build item: "
              << std::hex << hResult << std::endl;
    NMR::lib3mf_getlasterror(pModel, &nErrorMessage, &pszErrorMessage);
    std::cerr << "error #" << std::hex << nErrorMessage
              << ": " << pszErrorMessage << std::endl;
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
  if (hResult != LIB3MF_OK) {
    std::cerr << "could not create model reader: "
              << std::hex << hResult << std::endl;
    NMR::lib3mf_getlasterror(pModel, &nErrorMessage, &pszErrorMessage);
    std::cerr << "error #" << std::hex << nErrorMessage
              << ": " << pszErrorMessage << std::endl;
    NMR::lib3mf_release(pModel);
    return false;
  }
  // Export Model into File
  hResult = NMR::lib3mf_writer_writetofileutf8(p3MFWriter, file_name.c_str());
  if (hResult != LIB3MF_OK) {
    std::cerr << "could not write file: " << std::hex << hResult << std::endl;
    NMR::lib3mf_getlasterror(p3MFWriter, &nErrorMessage, &pszErrorMessage);
    std::cerr << "error #" << std::hex << nErrorMessage << ": "
              << pszErrorMessage << std::endl;
    NMR::lib3mf_release(pModel);
    NMR::lib3mf_release(p3MFWriter);
    return false;
  }

  // Release Model Writer
  NMR::lib3mf_release(p3MFWriter);
  return true;
}

template<typename PointRange, typename PolygonRange, typename ColorRange>
bool write_mesh_to_model( const PointRange& points,
                          const PolygonRange& polygons,
                          const ColorRange& colors,
                          const std::string& name,
                          NMR::PLib3MFModelMeshObject** pMeshObject,
                          NMR::PLib3MFModel * pModel
                          )
{
  DWORD nErrorMessage;
  LPCSTR pszErrorMessage;
  HRESULT hResult;
  // Create mesh structure
  std::vector<NMR::MODELMESHVERTEX> pVertices;
  std::vector<NMR::MODELMESHTRIANGLE> pTriangles;
  // Create Mesh Object
  hResult = NMR::lib3mf_model_addmeshobject(pModel, pMeshObject);
  if (hResult != LIB3MF_OK) {
    std::cerr << "could not add mesh object: " << std::hex
              << hResult << std::endl;
    NMR::lib3mf_getlasterror(pModel, &nErrorMessage, &pszErrorMessage);
    std::cerr << "error #" << std::hex << nErrorMessage << ": "
              << pszErrorMessage << std::endl;
    NMR::lib3mf_release(pModel);
    return false;
  }
  for( auto point : points)
  {
    pVertices.push_back(tmf_internal::fnCreateVertex(point.x(), point.y(), point.z()));
  }

  for( auto triangle : polygons)
  {
    pTriangles.push_back(tmf_internal::fnCreateTriangle(triangle[0], triangle[1], triangle[2]));
  }

  hResult = NMR::lib3mf_meshobject_setgeometry(*pMeshObject, pVertices.data(),
                                          pVertices.size(), pTriangles.data(),
                                          pTriangles.size());
  if (hResult != LIB3MF_OK) {
    std::cerr << "could not set mesh geometry: "
              << std::hex << hResult << std::endl;
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
  if (hResult != LIB3MF_OK) {
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
  if (hResult != LIB3MF_OK) {
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
  if (hResult != LIB3MF_OK) {
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
                  NMR::PLib3MFModel * pModel
                  )
{
  DWORD nErrorMessage;
  LPCSTR pszErrorMessage;
  HRESULT hResult;
  // Create mesh structure
  std::vector<NMR::MODELMESHVERTEX> pVertices;
  // Create Mesh Object
  hResult = NMR::lib3mf_model_addmeshobject(pModel, pMeshObject);
  if (hResult != LIB3MF_OK) {
    std::cerr << "could not add mesh object: " << std::hex
              << hResult << std::endl;
    NMR::lib3mf_getlasterror(pModel, &nErrorMessage, &pszErrorMessage);
    std::cerr << "error #" << std::hex << nErrorMessage << ": "
              << pszErrorMessage << std::endl;
    NMR::lib3mf_release(pModel);
    return false;
  }
  //add 3 demmy vertices to be sure to have a valid triangle and accept point sets with less than 3 vertices.
  for(int i = 0; i< 3; ++i)
    pVertices.push_back(tmf_internal::fnCreateVertex(0,0,0));
  for( auto point : points)
  {
    pVertices.push_back(tmf_internal::fnCreateVertex(point.x(), point.y(), point.z()));
  }
  NMR::MODELMESHTRIANGLE dummy_triangle = tmf_internal::fnCreateTriangle(0,1,2); //add a triangle to avoid lib error.
  hResult = NMR::lib3mf_meshobject_setgeometry(*pMeshObject, pVertices.data(),
                                        pVertices.size(), &dummy_triangle, 1);
  if (hResult != LIB3MF_OK) {
    std::cerr << "could not set mesh geometry: "
              << std::hex << hResult << std::endl;
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
  if (hResult != LIB3MF_OK) {
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
  if (hResult != LIB3MF_OK) {
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
  if (hResult != LIB3MF_OK) {
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
                                NMR::PLib3MFModel * pModel
                                )
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
                             NMR::PLib3MFModel * pModel
                             )
{
  std::string pc_name = name;
  pc_name.append("_cgal_pl");
  return write_points(points, color, pc_name, pMeshObject, pModel);
}

/*!
 * \brief writes the triangle soups contained in `all_points` and
 *  `all_polygons` into the 3mf file `file_name`.
 * \tparam PointRanges a model of the concepts `RandomAccessContainer` and
 *  `BackInsertionSequence` whose `value type` is
 * a model of the concepts `RandomAccessContainer` and `BackInsertionSequence`
 *  whose `value type` is the point type.
 * \tparam PolygonRanges a model of the concept `RandomAccessContainer` whose
 *  `value_type` is a model of the concept `RandomAccessContainer`
 * whose `value_type` is a model of the concept `RandomAccessContainer` whose
 *  `value_type` is std::size_t.
 * \param file_name the name of the 3mf file to write.
 * \param all_points a `PointRanges` that contains the points of the soups
 *  to write.
 * \param all_polygons a `PolygonRanges` that contains the triangles of the
 *  soups in `file_name`.
 * \param names will contains the name of each mesh in `file_name`.
 * \return `true` if the writing is successful, `false` otherwise.
 *
 * \attention Only versions inferior to 2.0 of lib3mf are supported.
 */
template<typename PointRanges, typename PolygonRanges>
bool write_triangle_soups_to_3mf(const std::string& file_name,
                        const PointRanges& all_points,
                        const PolygonRanges& all_polygons,
                        const std::vector<std::string>& names)
{
  HRESULT hResult;

    // Create Model Instance
  NMR::PLib3MFModel * pModel;
  hResult = NMR::lib3mf_createmodel(&pModel);
  if (hResult != LIB3MF_OK) {
    std::cerr << "could not create model: " << std::hex << hResult << std::endl;
    return false;
  }
  for(std::size_t id = 0; id < all_points.size();  ++id)
  {
    NMR::PLib3MFModelMeshObject* pMeshObject;
    std::string name;
    if(names.size() > id
       && ! names[id].empty())
    {
      name=names[id];
    }
    else
      name = std::string("");
    std::vector<CGAL::Color> colors(all_polygons[id].size());
    write_mesh_to_model(all_points[id], all_polygons[id], colors,  name, &pMeshObject, pModel);
  }
  return export_model_to_file(file_name, pModel);
}


/*!
 * \brief writes the triangle meshes contained in `tms`
 * into the 3mf file `file_name`.
 * \tparam TriangleMeshRange a model of the concepts `RandomAccessContainer`
 * and `BackInsertionSequence` whose `value type` is
 * a model of the concepts `FaceListGraph` and `HalfedgeListGraph`
 * that has only triangle faces.
 * \param file_name the name of the 3mf file to write.
 * \param tms a `TriangleMeshRange` that contains the meshes
 *  to write. An internal property map for `CGAL::vertex_point_t`
 * must be available for each mesh.
 * \param names will contains the name of each mesh in `file_name`.
 * \return `true` if the writing is successful, `false` otherwise.
 *
 * \attention Only versions inferior to 2.0 of lib3mf are supported.
 */
template<typename TriangleMeshRange>
bool write_triangle_meshes_to_3mf(const std::string& file_name,
                                  const TriangleMeshRange& tms,
                                  const std::vector<std::string>& names)
{
  typedef typename TriangleMeshRange::value_type Mesh;
  typedef typename boost::property_map<Mesh, boost::vertex_point_t>::type VPMap;
  typedef typename boost::property_traits<VPMap>::value_type Point_3;

  typedef std::vector<std::size_t> Polygon;
  typedef std::vector<Polygon> PolygonRange;

  typedef std::vector<Point_3> PointRange;

  std::vector<PointRange> all_points;
  std::vector<PolygonRange> all_polygons;


  for(auto tm : tms)
  {
    PointRange points;
    PolygonRange triangles;
    VPMap vpm = get(boost::vertex_point, tm);
    std::unordered_map<typename boost::graph_traits<Mesh>::vertex_descriptor,
        std::size_t> vertex_id_map;
    std::size_t i = 0;
    for(auto v : vertices(tm))
    {
      points.push_back(get(vpm, v));
      vertex_id_map[v] = i++;
    }
    all_points.push_back(points);
    for(auto f : faces(tm))
    {
      Polygon triangle;
      for(auto vert : CGAL::vertices_around_face(halfedge(f, tm), tm))
      {
        triangle.push_back(vertex_id_map[vert]);
      }
      triangles.push_back(triangle);
    }
    all_polygons.push_back(triangles);
  }

  return write_triangle_soups_to_3mf(file_name, all_points, all_polygons, names);
}
}//end CGAL
#endif // CGAL_IO_WRITE_3MF_H
