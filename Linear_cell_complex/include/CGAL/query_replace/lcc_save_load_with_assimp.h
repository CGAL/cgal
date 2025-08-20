// Copyright (c) 2022 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of 3d-query-replace.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
////////////////////////////////////////////////////////////////////////////////
#ifndef LCC_SAVE_LOAD_WITH_ASSIMP_H
#define LCC_SAVE_LOAD_WITH_ASSIMP_H

#ifdef WITH_ASSIMP

#include <assimp/Importer.hpp>
#include <assimp/Exporter.hpp>
#include <assimp/scene.h>       // Output data structure
#include <assimp/postprocess.h> // Post processing flags
#include <unordered_map>
#include "My_linear_cell_complex_incremental_builder.h"

///////////////////////////////////////////////////////////////////////////////
template<class LCC>
bool load_object_3D_with_assimp(const std::string& filename, LCC& lcc)
{
  Assimp::Importer importer;
  importer.SetPropertyInteger(AI_CONFIG_PP_RVC_FLAGS,
                              aiComponent_NORMALS|
                              aiComponent_TANGENTS_AND_BITANGENTS|
                              aiComponent_COLORS|
                              aiComponent_TEXCOORDS|
                              aiComponent_BONEWEIGHTS|
                              aiComponent_ANIMATIONS|
                              aiComponent_TEXTURES|
                              aiComponent_LIGHTS|
                              aiComponent_CAMERAS|
                              aiComponent_MATERIALS);
  const aiScene* scene=importer.ReadFile(filename,
                                         aiProcess_RemoveComponent |
                                         aiProcess_JoinIdenticalVertices);
  if( scene==nullptr) return false; // the import failed

  // Retrieve the geometry node of a model from an input file
  aiMesh** allmeshes=scene->mMeshes;
  unsigned int n_meshes=scene->mNumMeshes;

  My_linear_cell_complex_incremental_builder_3<LCC> IB(lcc);
  for(unsigned int m=0; m<n_meshes; ++m)
  {
    IB.begin_surface();
    aiMesh* mesh=allmeshes[m];
    aiVector3D* vertices=mesh->mVertices;
    aiFace* faces=mesh->mFaces;

    unsigned int n_vertices=mesh->mNumVertices;
    unsigned int n_faces=mesh->mNumFaces;

    // Add all vertices
    for (unsigned int i=0; i<n_vertices; ++i)
    {
      IB.add_vertex(typename LCC::Point(vertices[i].x,
                                        vertices[i].y,
                                        vertices[i].z));
    }
    // Add all faces
    for (unsigned int i=0; i<n_faces; ++i)
    {
      IB.begin_facet();
      const unsigned int nbPts=faces[i].mNumIndices;
      for(unsigned int j=0; j<nbPts; ++j)
      { IB.add_vertex_to_facet(faces[i].mIndices[j]); }
      IB.end_facet();
    }
    IB.end_surface();
  }

  return true;
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
void save_object_3D_with_assimp(const std::string& filename, LCC& lcc)
{
  std::unordered_map<typename LCC::Vertex_attribute_handle, std::size_t> index;
  
  aiScene scene;
  scene.mRootNode=new aiNode();
  scene.mRootNode->mMeshes=new unsigned int[1];
  scene.mRootNode->mMeshes[0]=0;
  scene.mRootNode->mNumMeshes=1;

  scene.mNumMeshes=lcc.template one_dart_per_cell<3>().size();
  scene.mMeshes=new aiMesh*[scene.mNumMeshes];

  std::size_t nbvol=0, nb=0, nb2=0;
  typename LCC::Dart_handle dhcur;

  for(auto itvol=lcc.template one_dart_per_cell<3>().begin(),
      itvolend=lcc.template one_dart_per_cell<3>().end();
      itvol!=itvolend; ++itvol, ++nbvol)
  {
    scene.mMeshes[nbvol]=new aiMesh();
    scene.mMeshes[nbvol]->mMaterialIndex=0;

    auto pMesh=scene.mMeshes[nbvol];
    pMesh->mNumVertices=lcc.template one_dart_per_incident_cell<0,3,2>(itvol).size();
    pMesh->mVertices=new aiVector3D[pMesh->mNumVertices];

    // save each vertex
    nb=0;
    for(auto itv=lcc.template one_dart_per_incident_cell<0,3,2>(itvol).begin(),
        itvend=lcc.template one_dart_per_incident_cell<0,3,2>(itvol).end();
        itv!=itvend; ++itv, ++nb)
    {
      const auto& pt=lcc.point(itv);
      pMesh->mVertices[nb]=aiVector3D(pt.x(), pt.y(), pt.z());
      index[lcc.vertex_attribute(itv)]=nb;
    }

    pMesh->mNumFaces=lcc.template one_dart_per_incident_cell<2,3,2>(itvol).size();
    pMesh->mFaces=new aiFace[pMesh->mNumFaces];

    // save each face
    nb=0;
    for(auto itface=lcc.template one_dart_per_incident_cell<2,3,2>(itvol).begin(),
        itfaceend=lcc.template one_dart_per_incident_cell<2,3,2>(itvol).end();
        itface!=itfaceend; ++itface, ++nb)
    {
      aiFace& face=pMesh->mFaces[nb];
      nb2=0;
      dhcur=itface;
      do
      { ++nb2; dhcur=lcc.template beta<1>(dhcur); }
      while(dhcur!=itface);
      face.mNumIndices=nb2;
      face.mIndices=new unsigned int[face.mNumIndices];
      dhcur=itface; nb2=0;
      do
      {
        face.mIndices[nb2]=index[lcc.vertex_attribute(dhcur)];
        ++nb2; dhcur=lcc.template beta<1>(dhcur);
      }
      while(dhcur!=itface);
    }
  }

/*
'assbin'    *.assbin
'assxml'   Assxml Document - *.assxml
'3ds'      Autodesk 3DS (legacy) - *.3ds
'fbxa'     Autodesk FBX (ascii) - *.fbx
'fbx'      Autodesk FBX (binary) - *.fbx
'collada'  COLLADA - Digital Asset Exchange Schema - *.dae
'x3d'      Extensible 3D - *.x3d
'gltf'     GL Transmission Format - *.gltf
'glb'      GL Transmission Format (binary) - *.glb
'gltf2'    GL Transmission Format v. 2 - *.gltf
'glb2'     GL Transmission Format v. 2 (binary) - *.glb
'ply'      Stanford Polygon Library - *.ply
'plyb'     Stanford Polygon Library (binary) - *.ply
'stp'      Step Files - *.stp
'stl'      Stereolithography - *.stl
'stlb'     Stereolithography (binary) - *.stl
'3mf'      The 3MF-File-Format - *.3mf
'obj'      Wavefront OBJ format - *.obj
'objnomtl' Wavefront OBJ format without material file - *.obj
'x'        X Files - *.x */
  // Export(&scene, const std::string &pFormatId, const std::string &pPath);
  Assimp::Exporter exporter;
  // TODO exporter.Export(&scene, format, filename);
}

#endif // WITH_ASSIMP
///////////////////////////////////////////////////////////////////////////////
#endif // LCC_SAVE_LOAD_WITH_ASSIMP_H
