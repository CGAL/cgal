// Copyright (c) 2025 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
////////////////////////////////////////////////////////////////////////////////
#ifndef LCC_SAVE_LOAD_MESH_H
#define LCC_SAVE_LOAD_MESH_H

#include <fstream>
#include <unordered_map>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/IO/OBJ.h>
#include <CGAL/Element_topo.h>
#include <CGAL/Linear_cell_complex_incremental_builder_3.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Surface_mesh.h>

namespace CGAL::IO
{
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
void load_object_2D(const std::string& filename, LCC& lcc,
                    std::vector<typename LCC::Dart_descriptor>* tabdarts=nullptr)
{
  std::ifstream fi(filename.c_str());
  if(!fi.good()) return; // File open error

  unsigned int nbparticles, nbelements, dim;
  fi >> nbparticles >> nbelements >> dim;

  if(nbparticles==0 || dim!=2)
  {
    fi.close();
    return;
  }

  typename LCC::FT x, y, z;
  int index;
  unsigned int nb_sommets;
  std::size_t i, j;

  /*! Use incremental builder to create LCC from a 2D mesh */
  typedef CGAL::Linear_cell_complex_incremental_builder_3<LCC>
    IncrementalBuilder;
  IncrementalBuilder IB(lcc);

  IB.begin_surface(nbparticles, nbelements, 0);

  // Initialization of the particles
  for(i=0; i<nbparticles; i++)
  {
    fi>>x>>y>>z;
    IB.add_vertex(typename LCC::Point(x, y, z));
  }

  // Initialization of the elements
  for(i=0; i<nbelements; i++)
  {
    fi>>nb_sommets;
    IB.begin_facet();

    for(j=0; j<nb_sommets; j++)
    {
      fi>>index;
      IB.add_vertex_to_facet(index, tabdarts);
    }
    IB.end_facet();
  }

  IB.end_surface();
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
bool load_object_3D(std::istream& fi, LCC& lcc,
                    std::vector<typename LCC::Dart_descriptor>* tabdarts=nullptr)
{
  std::size_t nbparticles, nbvols, dim;
  fi >> nbparticles >> nbvols >> dim;

  if(nbparticles == 0 || dim != 3)
  { return false; }

  typename LCC::FT x, y, z;
  std::size_t index_element[8];
  std::ptrdiff_t nb_vertices_signed;
  std::size_t nb_faces, nb_vertices_in_face;
  std::size_t index;

  // Retrieve geometrical coordinates of particle;
  // add vertices in an incremental builder
  CGAL::Linear_cell_complex_incremental_builder_3<LCC> IB(lcc);
  for(std::size_t i = 0; i < nbparticles; ++i)
  {
    fi >> x >> y >> z;
    IB.add_vertex(typename LCC::Point(x, y, z));
  }

  for(std::size_t  i = 0; i < nbvols; ++i)
  {
    fi >> nb_vertices_signed;

    /* Convention used in the file (the same than the one of gmesh, vtk...)
     *      7----6
     *     /|   /|
     *    4----5 |
     *    | 3--|-2
     *    |/   |/
     *    0----1
     *
     *      3
     *     /|\
     *    4---5
     *    | | |
     *    | 0 |
     *    |/ \|
     *    1---2
     *
     *      4
     *     /|\
     *    0-|-3
     *    | | |
     *    1---2
     *
     *      3
     *     /|\
     *    0-|-2
     *     \|/
     *      1
     */
    if(nb_vertices_signed == 4)//tetra
    {
      fi >> index_element[0] >> index_element[1] >> index_element[2]
         >> index_element[3];

      make_tetrahedron_with_builder(IB, index_element[0], index_element[1],
          index_element[2], index_element[3], tabdarts);
    }

    else if(nb_vertices_signed == 5)//pyramid
    {
      fi >> index_element[0] >> index_element[1] >> index_element[2]
         >> index_element[3] >> index_element[4];

      make_pyramid_with_builder(IB, index_element[0], index_element[1],
          index_element[2], index_element[3], index_element[4], tabdarts);
    }

    else if(nb_vertices_signed == 6)//prism
    {
      fi >> index_element[0] >> index_element[1] >> index_element[2]
         >> index_element[3] >> index_element[4] >> index_element[5];

      make_prism_with_builder(IB, index_element[0], index_element[1],
          index_element[2], index_element[3], index_element[4],
          index_element[5], tabdarts);
    }

    else if(nb_vertices_signed == 8)// hexa
    {
      fi >> index_element[0] >> index_element[1] >> index_element[2]
         >> index_element[3] >> index_element[4] >> index_element[5]
         >> index_element[6] >> index_element[7];

      make_hexahedron_with_builder(IB, index_element[0], index_element[1],
                          index_element[2], index_element[3], index_element[4],
          index_element[5], index_element[6], index_element[7], tabdarts);
    }

    else if(nb_vertices_signed<0) // "generic" cell
    {
      fi>>nb_faces;
      IB.begin_surface();
      for(std::size_t j=0; j<nb_faces; ++j)
      {
        fi>>nb_vertices_in_face;
        IB.begin_facet();
        for(std::size_t k=0; k<nb_vertices_in_face; ++k)
        {
          fi>>index;
          IB.add_vertex_to_facet(index, tabdarts);
        }
        IB.end_facet();
      }
      IB.end_surface();
    }
  }

  return true;
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
bool load_object_3D(const std::string& filename, LCC& lcc,
                    std::vector<typename LCC::Dart_descriptor>* tabdarts=nullptr)
{
  std::ifstream fi(filename.c_str());
  if(!fi.good()) return false; // File open error
  return load_object_3D(fi, lcc, tabdarts);
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
void save_one_generic_face(const LCC& lcc,
                           typename LCC::Dart_const_descriptor dh,
                           std::unordered_map
                           <typename LCC::Vertex_attribute_const_descriptor,
                                              std::size_t>& index,
                           std::ostream& fo,
                           std::unordered_map<typename LCC::Dart_const_descriptor,
                                              std::size_t>* dartsids=nullptr)
{
  std::size_t nb=0;
  typename LCC::Dart_const_descriptor cur=dh;
  do
  {
    ++nb;
    cur=lcc.next(cur);
  }
  while(cur!=dh);
  if(dartsids!=nullptr) { dartsids->reserve(dartsids->size()+nb); }
  fo<<nb<<" ";
  do
  {
    fo<<index[lcc.vertex_attribute(cur)]<<" ";
    if(dartsids!=nullptr) { (*dartsids)[cur]=dartsids->size(); }
    cur=lcc.next(cur);
  }
  while(cur!=dh);
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
void push_darts_of_face(const LCC& lcc, typename LCC::Dart_const_descriptor dh,
                        std::unordered_map<typename LCC::Dart_const_descriptor,
                                           std::size_t>& dartsids)
{
  typename LCC::Dart_const_descriptor cur=dh;
  do
  {
    dartsids[cur]=dartsids.size();
    cur=lcc.next(cur);
  }
  while(cur!=dh);
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
void save_one_volume_internal
(const LCC& lcc, typename LCC::Dart_const_descriptor dh,
 std::unordered_map<typename LCC::Vertex_attribute_const_descriptor,
                                                 std::size_t>& index,
 std::ostream& fo,
 std::unordered_map<typename LCC::Dart_const_descriptor, std::size_t>* dartsids=nullptr)
{
  /* Convention used in CGAL LCC (different than the one used in the file, cf above)
   *      3
   *     /|\
   *    0-|-2
   *     \|/
   *      1
   *  Dart incident to p0, to edge p0,p1 and to facet (p0,p1,p2).
   *
   *      4
   *     /|\
   *    0-|-3
   *    | | |
   *    1---2
   *  Dart incident to p0 and to the facet (p0,p1,p2,p3).
   *
   *      3
   *     /|\
   *    4---5
   *    | | |
   *    | 0 |
   *    |/ \|
   *    1---2
   *  Dart incident to p0 and to the facet (p0,p1,p2).
   *
   *      7----6
   *     /|   /|
   *    4----5 |
   *    | 3--|-2
   *    |/   |/
   *    0----1
   *  Dart incident to p0 and to the facet (p0,p1,p2,p3).
   */
  using namespace CGAL::CMap::Element_topo;
  typename LCC::Dart_const_descriptor sd;
  cell_topo vol_type=get_cell_topo<3>(lcc, dh, sd);
  if(vol_type==TETRAHEDRON)
  {
    fo<<"4 "
      <<index[lcc.vertex_attribute(sd)]<<" "
      <<index[lcc.vertex_attribute(lcc.template beta<1>(sd))]<<" "
      <<index[lcc.vertex_attribute(lcc.template beta<0>(sd))]<<" "
      <<index[lcc.vertex_attribute(lcc.template beta<2, 0>(sd))]<<std::endl;

    if(dartsids!=nullptr)
    {
      dartsids->reserve(dartsids->size()+12);
      push_darts_of_face(lcc, sd, *dartsids);
      push_darts_of_face(lcc, lcc.template beta<2>(sd), *dartsids);
      push_darts_of_face(lcc, lcc.template beta<1,2>(sd), *dartsids);
      push_darts_of_face(lcc, lcc.template beta<0,2>(sd), *dartsids);
    }
  }
  else if(vol_type==PYRAMID)
  {
    fo<<"5 "
      <<index[lcc.vertex_attribute(sd)]<<" "
      <<index[lcc.vertex_attribute(lcc.template beta<1>(sd))]<<" "
      <<index[lcc.vertex_attribute(lcc.template beta<1,1>(sd))]<<" "
      <<index[lcc.vertex_attribute(lcc.template beta<0>(sd))]<<" "
      <<index[lcc.vertex_attribute(lcc.template beta<2,0>(sd))]<<std::endl;

    if(dartsids!=nullptr)
    {
      dartsids->reserve(dartsids->size()+16);
      push_darts_of_face(lcc, sd, *dartsids);
      push_darts_of_face(lcc, lcc.template beta<2>(sd), *dartsids);
      push_darts_of_face(lcc, lcc.template beta<1,2>(sd), *dartsids);
      push_darts_of_face(lcc, lcc.template beta<1,1,2>(sd), *dartsids);
      push_darts_of_face(lcc, lcc.template beta<0,2>(sd), *dartsids);
    }
  }
  else if(vol_type==PRISM)
  {
    fo<<"6 "
      <<index[lcc.vertex_attribute(sd)]<<" "
      <<index[lcc.vertex_attribute(lcc.template beta<1>(sd))]<<" "
      <<index[lcc.vertex_attribute(lcc.template beta<0>(sd))]<<" ";

    // Move to the up face
    typename LCC::Dart_const_descriptor d2=lcc.template beta<2, 1, 1, 2>(sd);
    fo<<index[lcc.vertex_attribute(lcc.template beta<1>(d2))]<<" "
      <<index[lcc.vertex_attribute(d2)]<<" "
      <<index[lcc.vertex_attribute(lcc.template beta<0>(d2))]<<std::endl;

    if(dartsids!=nullptr)
    {
      dartsids->reserve(dartsids->size()+18);
      push_darts_of_face(lcc, sd, *dartsids);
      push_darts_of_face(lcc, lcc.template beta<2>(sd), *dartsids);
      push_darts_of_face(lcc, lcc.template beta<1,2>(sd), *dartsids);
      push_darts_of_face(lcc, lcc.template beta<0,2>(sd), *dartsids);
      push_darts_of_face(lcc, lcc.template beta<2,1,1,2,0>(sd), *dartsids);
    }
  }
  else if(vol_type==HEXAHEDRON)
  {
    fo<<"8 ";
    // Darts associated with particles 0, 1, 2, 3
    for(unsigned int i=0; i<4; ++i)
    {
      fo<<index[lcc.vertex_attribute(sd)]<<" ";
      sd=lcc.template beta<1>(sd);
    }
    typename LCC::Dart_const_descriptor d2=lcc.template beta<2, 1, 1, 2, 1>(sd);
    // Darts associated with particles 4, 7, 6, 5
    for(unsigned int i = 0; i < 4; i++)
    {
      fo<<index[lcc.vertex_attribute(d2)]<<" ";
      d2 = lcc.template beta<0>(d2);
    }
    fo<<std::endl;

    if(dartsids!=nullptr)
    {
      dartsids->reserve(dartsids->size()+18);
      push_darts_of_face(lcc, sd, *dartsids);
      push_darts_of_face(lcc, lcc.template beta<2>(sd), *dartsids);
      push_darts_of_face(lcc, lcc.template beta<1,2>(sd), *dartsids);
      push_darts_of_face(lcc, lcc.template beta<1,1,2>(sd), *dartsids);
      push_darts_of_face(lcc, lcc.template beta<0,2>(sd), *dartsids);
      push_darts_of_face(lcc, lcc.template beta<2,1,1,2,1,1>(sd), *dartsids);
    }
  }
  else
  {
    // 1) number of vertices of the volume (negative number for generic cell)
    // and number of faces of the volume
    fo<<-static_cast<ptrdiff_t>(lcc.template one_dart_per_incident_cell<0,3,2>(sd).size())<<" "
      <<(lcc.template one_dart_per_incident_cell<2,3,2>(sd).size())<<std::endl;
    // 2) save each face
    for(auto itface=lcc.template one_dart_per_incident_cell<2,3,2>(sd).begin(),
          itfaceend=lcc.template one_dart_per_incident_cell<2,3,2>(sd).end();
        itface!=itfaceend; ++itface)
    {
      fo<<"  ";
      save_one_generic_face(lcc, itface, index, fo, dartsids);
      fo<<std::endl;
    }
  }
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
void save_object_3D(std::ostream& fo, const LCC& lcc,
                    std::unordered_map<typename LCC::Dart_const_descriptor,
                                       std::size_t>* dartsids=nullptr)
{
  if(lcc.is_empty())
  { fo<<"0 0 3"<<std::endl; return; }

  // #vertices #volumes dimension
  fo<<lcc.vertex_attributes().size()<<" "
    <<lcc.template one_dart_per_cell<3>().size()<<" 3"<<std::endl<<std::endl;

  std::unordered_map<typename LCC::Vertex_attribute_const_descriptor,
                     std::size_t> index;
  std::size_t nb=0;

  // All vertices
  for(auto itv=lcc.vertex_attributes().begin(),
        itvend=lcc.vertex_attributes().end(); itv!=itvend; ++itv)
  {
    fo<<itv->point()<<std::endl;
    index[itv]=nb++;
  }
  fo<<std::endl;

  // All volumes; using indices of vertices
  for(auto itvol=lcc.template one_dart_per_cell<3>().begin(),
        itvolend=lcc.template one_dart_per_cell<3>().end();
      itvol!=itvolend; ++itvol)
  { save_one_volume_internal(lcc, itvol, index, fo, dartsids); }
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
bool save_object_3D(const std::string& filename, const LCC& lcc,
                    std::unordered_map<typename LCC::Dart_const_descriptor,
                                       std::size_t>* dartsids=nullptr)
{
  std::ofstream fo(filename.c_str());
  if(!fo.good()) return false; // File open error
  save_object_3D(fo, lcc, dartsids);
  fo.close();
  return true;
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
void save_one_volume(std::ostream& fo, const LCC& lcc,
                     typename LCC::Dart_const_descriptor dh,
                     std::unordered_map<typename LCC::Dart_const_descriptor,
                                        std::size_t>* dartsids=nullptr)
{
  std::unordered_map<typename LCC::Vertex_attribute_const_descriptor,
                     std::size_t> index;
  std::size_t nb=0;
  typename LCC::Dart_const_descriptor sd;

  // #vertices #volumes dimension
  fo<<lcc.template one_dart_per_incident_cell<0,3>(dh).size()
    <<" 1 3"<<std::endl<<std::endl;

  // all vertices of volume(dh)
  auto itv=lcc.template one_dart_per_incident_cell<0,3>(dh).begin();
  for(auto itvend=lcc.template one_dart_per_incident_cell<0,3>(dh).end();
      itv!=itvend; ++itv)
  {
    fo<<lcc.point(itv)<<std::endl;
    index[lcc.vertex_attribute(itv)]=nb++;
  }
  fo<<std::endl;

  // The volume, using indices of vertices
  save_one_volume_internal(lcc, dh, index, fo, dartsids);
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
bool save_one_volume(const std::string& filename, const LCC& lcc,
                     typename LCC::Dart_const_descriptor dh,
                     std::unordered_map<typename LCC::Dart_const_descriptor,
                                        std::size_t>* dartsids=nullptr)
{
  std::ofstream fo(filename.c_str());
  if(!fo.good()) return false; // File open error
  save_one_volume(fo, lcc, dh, dartsids);
  fo.close();
  return true;
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
bool save_object_3D_gmsh(const std::string& filename, const LCC& lcc)
{
  using namespace CGAL::CMap::Element_topo;
  std::ofstream fo(filename.c_str());
  if(!fo.good()) return false; // File open error

  fo<<"$MeshFormat"<<std::endl;
  fo<<"2.2 0 8"<<std::endl;
  fo<<"$EndMeshFormat"<<std::endl;

  fo<<"$Nodes"<<std::endl;
  fo<<lcc.vertex_attributes().size()<<std::endl;
  std::unordered_map<typename LCC::Vertex_attribute_const_descriptor,
                     std::size_t> index;
  std::size_t nb=0;
  typename LCC::Dart_const_descriptor sd;

  for(auto itv=lcc.vertex_attributes().begin(),
        itvend=lcc.vertex_attributes().end(); itv!=itvend; ++itv)
  {
    fo<<nb<<" "<<itv->point()<<std::endl;
    index[itv]=nb++;
  }
  fo<<"$EndNodes"<<std::endl;

  nb=0;
  fo<<"$Elements"<<std::endl;

  // Count the number of regular elements, and write it in the file.
  for(auto itvol=lcc.template one_dart_per_cell<3>().begin(),
        itvolend=lcc.template one_dart_per_cell<3>().end();
      itvol!=itvolend; ++itvol)
  {
    cell_topo vol_type=get_cell_topo<3>(lcc, itvol, sd);
    if(vol_type==TETRAHEDRON || vol_type==PYRAMID ||
       vol_type==PRISM || vol_type==HEXAHEDRON)
    { ++nb; }
  }
  fo<<nb<<std::endl;

  // Write the different elements
  nb=0;
  for(auto itvol=lcc.template one_dart_per_cell<3>().begin(),
        itvolend=lcc.template one_dart_per_cell<3>().end();
      itvol!=itvolend; ++itvol)
  {
    cell_topo vol_type=get_cell_topo<3>(lcc, itvol, sd);
    if(vol_type==TETRAHEDRON)
    {
       fo<<nb++<<" 4 0 "
         <<index[lcc.vertex_attribute(sd)]<<" "
         <<index[lcc.vertex_attribute(lcc.template beta<1>(sd))]<<" "
         <<index[lcc.vertex_attribute(lcc.template beta<0>(sd))]<<" "
         <<index[lcc.vertex_attribute(lcc.template beta<2, 0>(sd))]<<std::endl;
    }
    else if(vol_type==PYRAMID)
    {
      fo<<nb++<<" 7 0 "
         <<index[lcc.vertex_attribute(sd)]<<" "
         <<index[lcc.vertex_attribute(lcc.template beta<1>(sd))]<<" "
         <<index[lcc.vertex_attribute(lcc.template beta<1,1>(sd))]<<" "
         <<index[lcc.vertex_attribute(lcc.template beta<0>(sd))]<<" "
         <<index[lcc.vertex_attribute(lcc.template beta<2,0>(sd))]<<std::endl;
    }
    else if(vol_type==PRISM)
    {
      fo<<nb++<<" 6 0 "
         <<index[lcc.vertex_attribute(sd)]<<" "
         <<index[lcc.vertex_attribute(lcc.template beta<1>(sd))]<<" "
         <<index[lcc.vertex_attribute(lcc.template beta<0>(sd))]<<" ";
      // Move to the up face
      typename LCC::Dart_const_descriptor d2=lcc.template beta<2, 1, 1, 2>(sd);
      fo<<index[lcc.vertex_attribute(lcc.template beta<1>(d2))]<<" "
        <<index[lcc.vertex_attribute(d2)]<<" "
        <<index[lcc.vertex_attribute(lcc.template beta<0>(d2))]<<std::endl;
    }
    else if(vol_type==HEXAHEDRON)
    {
      fo<<nb++<<" 5 0 ";
      for(unsigned int i=0; i<4; ++i)
      {
        fo<<index[lcc.vertex_attribute(sd)]<<" ";
        sd=lcc.template beta<1>(sd);
      }
      typename LCC::Dart_const_descriptor d2=lcc.template beta<2, 1, 1, 2, 1>(sd);
      // Darts associated with particles 4, 5, 6, 7
      for(unsigned int i = 0; i < 4; i++)
        {
          fo<<index[lcc.vertex_attribute(d2)]<<" ";
          d2 = lcc.template beta<0>(d2);
        }
      fo<<std::endl;
    }
    else
    { // TODO Generic case, not posible with gmsh format. is it true?
    }
  }
  fo<<"$EndElements"<<std::endl;
  return true;
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
bool read_object_3D_gmsh(const std::string& filename, LCC& lcc)
{
  std::ifstream fi(filename.c_str());
  if(!fi.good()) return false; // File open error

  std::string txt;
  double d, d2, d3;
  unsigned int i, i2, i3, i4;

  fi>>txt; // gmsh header
  if(txt!="$MeshFormat")
  {
    std::cerr<<"Error reading "<<filename<<"; not a gmsh format."<<std::endl;
    return false;
  }

  fi>>d; // gmsh version
  assert(d>=2.1);

  fi>>i; // 0 MSH file format version 2 (Legacy)
  if(i!=0)
  {
    std::cerr<<"Error reading "<<filename
             <<"; for now only MSH file format version 2 is supported."
             <<std::endl;
    return false;
  }
  fi>>i; // size of size_t
  assert(i==8);
  fi>>txt; // end of gmsh header
  assert(txt=="$EndMeshFormat");

  CGAL::Linear_cell_complex_incremental_builder_3<LCC> IB(lcc);

  // Read nodes
  do
  {
    if(fi.eof())
    {
      std::cerr<<"Error reading "<<filename<<"; no nodes."<<std::endl;
      return false;
    }
    fi>>txt; // begin nodes
  }
  while(txt!="$Nodes");

  // Map to give for each vertex id, its number in the incremental builder
  std::unordered_map<unsigned int, unsigned int> vertex_indices;
  unsigned int index=0;

  fi>>i; // Number of nodes
  for(unsigned int nb=0; nb<i; ++nb)
  {
    if(fi.eof())
    {
      std::cerr<<"Error reading "<<filename<<"; not enough nodes."<<std::endl;
      return false;
    }
    fi>>i2>>d>>d2>>d3; // node id; position of node x y z
    IB.add_vertex(typename LCC::Point(d, d2, d3));

    vertex_indices[i2]=index++;
  }

  fi>>txt; // end nodes
  assert(txt=="$EndNodes");

  // Read elements
  do
  {
    if(fi.eof())
    {
      std::cerr<<"Error reading "<<filename<<"; no elements."<<std::endl;
      return false;
    }
    fi>>txt; // begin nodes
  }
  while(txt!="$Elements");

  std::size_t index_element[8];
  fi>>i; // Number of elements
  for(unsigned int nb=0; nb<i; ++nb)
  {
    if(fi.eof())
    {
      std::cerr<<"Error reading "<<filename<<"; not enough elements."<<std::endl;
      return false;
    }
    fi>>i2>>i3>>i4; // element id; type of elements; number of tags
    // std::cout<<"Element "<<i2<<" "<<i3<<" "<<i4<<std::endl;
    // assert(i2==nb);
    for(unsigned int nb2=0; nb2<i4; ++nb2)
    { fi>>i2; } // Ignore all the tags

    switch(i3)
    {
    case 4: // Tetrahedron
       fi >> index_element[0] >> index_element[1] >> index_element[2]
         >> index_element[3];

      make_tetrahedron_with_builder(IB,
                                    vertex_indices[index_element[0]],
                                    vertex_indices[index_element[1]],
                                    vertex_indices[index_element[2]],
                                    vertex_indices[index_element[3]]);
     break;
    case 5: // Hexahedron
      fi >> index_element[0] >> index_element[1] >> index_element[2]
         >> index_element[3] >> index_element[4] >> index_element[5]
         >> index_element[6] >> index_element[7];

      make_hexahedron_with_builder(IB,
                                   vertex_indices[index_element[0]],
                                   vertex_indices[index_element[1]],
                                   vertex_indices[index_element[2]],
                                   vertex_indices[index_element[3]],
                                   vertex_indices[index_element[4]],
                                   vertex_indices[index_element[5]],
                                   vertex_indices[index_element[6]],
                                   vertex_indices[index_element[7]]);
      break;
    case 6: // Prism
      fi >> index_element[0] >> index_element[1] >> index_element[2]
         >> index_element[3] >> index_element[4] >> index_element[5];

      make_prism_with_builder(IB,
                              vertex_indices[index_element[0]],
                              vertex_indices[index_element[1]],
                              vertex_indices[index_element[2]],
                              vertex_indices[index_element[3]],
                              vertex_indices[index_element[4]],
                              vertex_indices[index_element[5]]);
      break;
    case 7: // Pyramid
      fi >> index_element[0] >> index_element[1] >> index_element[2]
         >> index_element[3] >> index_element[4];

      make_pyramid_with_builder(IB,
                              vertex_indices[index_element[0]],
                              vertex_indices[index_element[1]],
                              vertex_indices[index_element[2]],
                              vertex_indices[index_element[3]],
                              vertex_indices[index_element[4]]);
      break;
    default:
      std::cout<<"Element type "<<i2<<" not considered."<<std::endl;
      fi.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
  }

  fi>>txt; // end nodes
  assert(txt=="$EndElements");

  // It is possible that the mesh is non manifold => correct attributes.
  lcc.correct_invalid_attributes();

  return true;
}
///////////////////////////////////////////////////////////////////////////////
}
#endif // LCC_SAVE_LOAD_MESH_H
