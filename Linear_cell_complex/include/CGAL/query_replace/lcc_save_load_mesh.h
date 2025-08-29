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
#ifndef LCC_SAVE_LOAD_MESH_H
#define LCC_SAVE_LOAD_MESH_H

#include "Prism_and_pyramid_creation.h"
#include "Element_topo.h"
#include <fstream>
#include <unordered_map>
#include "My_linear_cell_complex_incremental_builder.h"

///////////////////////////////////////////////////////////////////////////////
template<class LCC>
void load_object_2D(const std::string& filename, LCC& lcc)
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
  typedef My_linear_cell_complex_incremental_builder_3<LCC>
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
      IB.add_vertex_to_facet(index);
    }
    IB.end_facet();
  }

  IB.end_surface();
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
bool load_object_3D(const std::string& filename, LCC& lcc)
{
  std::ifstream fi(filename.c_str());
  if(!fi.good()) return false; // File open error

  std::size_t nbparticles, nbvols, dim;
  fi >> nbparticles >> nbvols >> dim;

  if(nbparticles == 0 || dim != 3)
  {
    fi.close();
    return false;
  }

  typename LCC::FT x, y, z;
  std::size_t index_element[8];
  ptrdiff_t nb_vertices_signed;
  std::size_t nb_faces, nb_vertices_in_face;
  std::size_t index;

  // Retrieve geometrical coordinates of particle; add vertices in an incremental builder
  My_linear_cell_complex_incremental_builder_3<LCC> IB(lcc);
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
                       index_element[2], index_element[3]);
    }

    else if(nb_vertices_signed == 5)//pyramid
    {
      fi >> index_element[0] >> index_element[1] >> index_element[2]
         >> index_element[3] >> index_element[4];

      make_pyramid_with_builder(IB, index_element[0], index_element[1],
          index_element[2], index_element[3], index_element[4]);
    }

    else if(nb_vertices_signed == 6)//prism
    {
      fi >> index_element[0] >> index_element[1] >> index_element[2]
         >> index_element[3] >> index_element[4] >> index_element[5];

      make_prism_with_builder(IB, index_element[0], index_element[1],
          index_element[2], index_element[3], index_element[4], index_element[5]);
    }

    else if(nb_vertices_signed == 8)// hexa
    {
      fi >> index_element[0] >> index_element[1] >> index_element[2]
         >> index_element[3] >> index_element[4] >> index_element[5]
         >> index_element[6] >> index_element[7];

      make_hexahedron_with_builder(IB, index_element[0], index_element[1],
                          index_element[2], index_element[3], index_element[4],
                          index_element[5], index_element[6], index_element[7]);
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
          IB.add_vertex_to_facet(index);
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
void save_one_generic_face(LCC& lcc, typename LCC::Dart_handle dh,
                           std::unordered_map<typename LCC::Vertex_attribute_handle, std::size_t>& index,
                           std::ofstream& of)
{
  std::size_t nb=0;
  typename LCC::Dart_handle cur=dh;
  do
  {
    ++nb;
    cur=lcc.next(cur);
  }
  while(cur!=dh);
  of<<nb<<" ";
  do
  {
    of<<index[lcc.vertex_attribute(cur)]<<" ";
    cur=lcc.next(cur);
  }
  while(cur!=dh);
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
void save_object_3D(const std::string& filename, LCC& lcc)
{
  std::ofstream fo(filename.c_str());
  if(!fo.good()) return; // File open error

  // #vertices #volumes dimension
  fo<<lcc.vertex_attributes().size()<<" "
    <<lcc.template one_dart_per_cell<3>().size()<<" 3"<<std::endl<<std::endl;

  std::unordered_map<typename LCC::Vertex_attribute_handle, std::size_t> index;
  std::size_t nb=0;
  typename LCC::Dart_handle sd;
  
  // all vertices
  for(auto itv=lcc.vertex_attributes().begin(),
        itvend=lcc.vertex_attributes().end(); itv!=itvend; ++itv)
  {    
    fo<<itv->point()<<std::endl;
    index[itv]=nb++;
  }
  fo<<std::endl;
  
  // all volumes; using indices of vertices
  for(auto itvol=lcc.template one_dart_per_cell<3>().begin(),
        itvolend=lcc.template one_dart_per_cell<3>().end(); 
      itvol!=itvolend; ++itvol)
  {
    /* Convention used in CGAL LCC (different than the one used in the file, cf above)
     *      3
     *     /|\
     *    0-|-2
     *     \|/
     *      1
     *  Dart incident to p0, to edge p0,p1 and to facet p0,p1,p2.
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
     *  Dart incident to p0, to edge p0,p5 and to the facet (p0,p5,p6,p1).
     */
    cell_topo vol_type=get_cell_topo(lcc, itvol, sd);
    if(vol_type==TETRAHEDRON)
    {
       fo<<"4 "
         <<index[lcc.vertex_attribute(sd)]<<" "
         <<index[lcc.vertex_attribute(lcc.template beta<1>(sd))]<<" "
         <<index[lcc.vertex_attribute(lcc.template beta<0>(sd))]<<" "
         <<index[lcc.vertex_attribute(lcc.template beta<2, 0>(sd))]<<std::endl;
    }
    else if(vol_type==PYRAMID)
    {
      fo<<"5 "
        <<index[lcc.vertex_attribute(sd)]<<" "
        <<index[lcc.vertex_attribute(lcc.template beta<1>(sd))]<<" "
        <<index[lcc.vertex_attribute(lcc.template beta<1,1>(sd))]<<" "
        <<index[lcc.vertex_attribute(lcc.template beta<0>(sd))]<<" "
        <<index[lcc.vertex_attribute(lcc.template beta<2,0>(sd))]<<std::endl;
    }
    else if(vol_type==PRISM)
    {
      fo<<"6 "
        <<index[lcc.vertex_attribute(sd)]<<" "
        <<index[lcc.vertex_attribute(lcc.template beta<1>(sd))]<<" "
        <<index[lcc.vertex_attribute(lcc.template beta<0>(sd))]<<" ";
      
      // Move to the up face
      typename LCC::Dart_handle d2=lcc.template beta<2, 1, 1, 2>(sd);
      fo<<index[lcc.vertex_attribute(lcc.template beta<1>(d2))]<<" "
        <<index[lcc.vertex_attribute(d2)]<<" "
        <<index[lcc.vertex_attribute(lcc.template beta<0>(d2))]<<std::endl;
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
      typename LCC::Dart_handle d2=lcc.template beta<2, 1, 1, 2, 1>(sd);
      // Darts associated with particles 4, 5, 6, 7
      for(unsigned int i = 0; i < 4; i++)
        {
          fo<<index[lcc.vertex_attribute(d2)]<<" ";
          d2 = lcc.template beta<0>(d2);
        }
      fo<<std::endl;
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
        save_one_generic_face(lcc, itface, index, fo);
        fo<<std::endl;
      }
    }
  }

  fo.close();
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
void save_object_3D_gmsh(const std::string& filename, LCC& lcc)
{
  std::ofstream fo(filename.c_str());
  if(!fo.good()) return; // File open error

  fo<<"$MeshFormat"<<std::endl;
  fo<<"2.2 0 8"<<std::endl;
  fo<<"$EndMeshFormat"<<std::endl;

  fo<<"$Nodes"<<std::endl;
  fo<<lcc.vertex_attributes().size()<<std::endl;
  std::unordered_map<typename LCC::Vertex_attribute_handle, std::size_t> index;
  std::size_t nb=0;
  typename LCC::Dart_handle sd;

  for(auto itv=lcc.vertex_attributes().begin(),
        itvend=lcc.vertex_attributes().end(); itv!=itvend; ++itv)
  {
    fo<<nb<<" "<<itv->point()<<std::endl;
    index[itv]=nb++;
  }
  fo<<"$EndNodes"<<std::endl;

  nb=0;
  fo<<"$Elements"<<std::endl;
  fo<<lcc.template one_dart_per_cell<3>().size()<<std::endl;
  for(auto itvol=lcc.template one_dart_per_cell<3>().begin(),
        itvolend=lcc.template one_dart_per_cell<3>().end();
      itvol!=itvolend; ++itvol, ++nb)
  {
    cell_topo vol_type=get_cell_topo(lcc, itvol, sd);
    if(vol_type==TETRAHEDRON)
    {
       fo<<nb<<" 4 2 0 1 "
         <<index[lcc.vertex_attribute(sd)]<<" "
         <<index[lcc.vertex_attribute(lcc.template beta<1>(sd))]<<" "
         <<index[lcc.vertex_attribute(lcc.template beta<0>(sd))]<<" "
         <<index[lcc.vertex_attribute(lcc.template beta<2, 0>(sd))]<<std::endl;
    }
    else if(vol_type==PYRAMID)
    {
      fo<<nb<<" 7 2 0 1 "
         <<index[lcc.vertex_attribute(sd)]<<" "
         <<index[lcc.vertex_attribute(lcc.template beta<1>(sd))]<<" "
         <<index[lcc.vertex_attribute(lcc.template beta<1,1>(sd))]<<" "
         <<index[lcc.vertex_attribute(lcc.template beta<0>(sd))]<<" "
         <<index[lcc.vertex_attribute(lcc.template beta<2,0>(sd))]<<std::endl;
    }
    else if(vol_type==PRISM)
    {
      fo<<nb<<" 6 2 0 1 "
         <<index[lcc.vertex_attribute(sd)]<<" "
         <<index[lcc.vertex_attribute(lcc.template beta<1>(sd))]<<" "
         <<index[lcc.vertex_attribute(lcc.template beta<0>(sd))]<<" ";
      // Move to the up face
      typename LCC::Dart_handle d2=lcc.template beta<2, 1, 1, 2>(sd);
      fo<<index[lcc.vertex_attribute(lcc.template beta<1>(d2))]<<" "
        <<index[lcc.vertex_attribute(d2)]<<" "
        <<index[lcc.vertex_attribute(lcc.template beta<0>(d2))]<<std::endl;
    }
    else if(vol_type==HEXAHEDRON)
    {
      fo<<nb<<" 5 2 0 1 ";
      for(unsigned int i=0; i<4; ++i)
      {
        fo<<index[lcc.vertex_attribute(sd)]<<" ";
        sd=lcc.template beta<1>(sd);
      }
      typename LCC::Dart_handle d2=lcc.template beta<2, 1, 1, 2, 1>(sd);
      // Darts associated with particles 4, 5, 6, 7
      for(unsigned int i = 0; i < 4; i++)
        {
          fo<<index[lcc.vertex_attribute(d2)]<<" ";
          d2 = lcc.template beta<0>(d2);
        }
      fo<<std::endl;
    }
    else
    { // TODO Generic case, not posible with gmsh format
    }
  }
  fo<<"$EndElements"<<std::endl;
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
void save_lcc_surface_into_off(const std::string& filename, LCC& lcc)
{
  std::ofstream fo(filename.c_str());
  if(!fo.good()) return; // File open error

  auto vertex_3free=lcc.get_new_mark(), face_3free=lcc.get_new_mark();

  // count and mark all 3 free vertices and faces.
  std::size_t nbvertices=0, nbfaces=0;
  for(auto it=lcc.darts().begin(), itend=lcc.darts().end(); it!=itend; ++it)
  {
    if(lcc.template is_free<3>(it))
    {
      if(!lcc.is_marked(it, vertex_3free))
      { lcc.template mark_cell<0>(it, vertex_3free); ++nbvertices; }
      if(!lcc.is_marked(it, face_3free))
      { lcc.template mark_cell<2,2>(it, face_3free); ++nbfaces; }
    }
  }

  // #vertices #faces 0 (0 to ignore number of edges)
  fo<<"OFF"<<std::endl<<nbvertices<<" "<<nbfaces<<" 0"<<std::endl<<std::endl;

  std::unordered_map<typename LCC::Vertex_attribute_handle, std::size_t> index;

  // all vertices
  std::size_t nb=0;
  // For this loop, we cannot iterate through vertex attributes since they may
  // not contain one of its dart. In such a case, it is not possible to ummark
  // the 3free vertices.
  for(auto it=lcc.darts().begin(), itend=lcc.darts().end(); it!=itend; ++it)
  {
    if (lcc.is_marked(it, vertex_3free))
    {
      fo<<lcc.point(it)<<std::endl;
      index[lcc.vertex_attribute(it)]=nb++;
      lcc.template unmark_cell<0>(it, vertex_3free);
    }
  }
  fo<<std::endl;

  // all 3free faces
  for(auto it=lcc.darts().begin(), itend=lcc.darts().end(); it!=itend; ++it)
  {
    if (lcc.is_marked(it, face_3free))
    {
      save_one_generic_face(lcc, it, index, fo);
      fo<<std::endl;
      lcc.template unmark_cell<2,2>(it, face_3free);
    }
  }

  lcc.free_mark(vertex_3free);
  lcc.free_mark(face_3free);
}
///////////////////////////////////////////////////////////////////////////////
#endif // LCC_SAVE_LOAD_MESH_H
