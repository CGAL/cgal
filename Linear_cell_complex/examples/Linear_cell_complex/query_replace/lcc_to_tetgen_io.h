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
#ifndef LCC_TO_TETGEN_IO_H
#define LCC_TO_TETGEN_IO_H

#include <tetgen.h>
#include <unordered_map>
#include <vector>
#include <iostream>
////////////////////////////////////////////////////////////////////////////////
/** Internal function that process one vertex of the lcc
 */
template <typename LCC>
void lcc_to_tetgen_process_vertex(const LCC& lcc,
                                  typename LCC::Dart_const_handle dh,
                                  typename LCC::size_type markv,
                                  tetgenio& io, std::size_t nbv,
                                  std::unordered_map
                                  <typename LCC::Vertex_attribute_const_handle,
                                  std::size_t>& vertices_map)
{
  io.pointlist[nbv*3]=lcc.point(dh).x();
  io.pointlist[(nbv*3)+1]=lcc.point(dh).y();
  io.pointlist[(nbv*3)+2]=lcc.point(dh).z();
  vertices_map[lcc.vertex_attribute(dh)]=nbv;
  lcc.template unmark_cell<0>(dh, markv);
}
////////////////////////////////////////////////////////////////////////////////
/** Internal function that process one face of the lcc
 */
template <typename LCC>
void lcc_to_tetgen_process_face(const LCC& lcc,
                                typename LCC::Dart_const_handle dh,
                                typename LCC::size_type markf,
                                tetgenio& io, std::size_t nbf,
                                std::unordered_map
                                <typename LCC::Vertex_attribute_const_handle,
                                std::size_t>& vertices_map)
{
  std::size_t nb=0;
  typename LCC::Dart_const_handle cur=dh;
  do
  {
    lcc.unmark(cur, markf);
    ++nb;
    cur=lcc.next(cur);
  }
  while(cur!=dh);

  tetgenio::facet *f=&io.facetlist[nbf];
  // Initialize the fields of this facet.
  //   There is one polygon, no hole.
  f->numberofpolygons=1;
  f->polygonlist=new tetgenio::polygon[1]; // Allocate memory for polygons.
  f->numberofholes=0;
  f->holelist=nullptr;
  tetgenio::polygon *p=&f->polygonlist[0];
  p->numberofvertices=nb;
  p->vertexlist=new int[nb]; // Allocate memory for vertices.
  nb=0;
  do
  {
    p->vertexlist[nb++]=vertices_map[lcc.vertex_attribute(cur)];
    cur=lcc.next(cur);
  }
  while(cur!=dh);
}
////////////////////////////////////////////////////////////////////////////////
/** Build a tetgenio frwom a given volume of an LCC.
   * @pre all the faces of the volume must be triangles
   *       (we can call constrained_delaunay_triangulation as pre-process).
   */
template <typename LCC>
void lcc_to_tetgen_one_volume_(const LCC& lcc,
                               typename LCC::Dart_const_handle dh,
                               tetgenio& io)
{
  // All indices start from 0.
  io.firstnumber=0;

  // First copy vertices
  std::unordered_map<typename LCC::Vertex_attribute_const_handle, std::size_t>
      vertices_map;
  auto markv=lcc.get_new_mark();
  auto markf=lcc.get_new_mark();
  auto treatedv=lcc.get_new_mark();

  io.numberofpoints=0;
  io.numberoffacets=0;
  for(auto itvol=lcc.template darts_of_cell_basic<3>(dh, treatedv).begin(),
      itvolend=lcc.template darts_of_cell_basic<3>(dh, treatedv).end();
      itvol!=itvolend; ++itvol)
  {
    lcc.mark(itvol, treatedv);
    if(!lcc.is_marked(itvol, markv))
    {
      ++(io.numberofpoints);
      lcc.template mark_cell<0>(itvol, markv);
    }
    if(!lcc.is_marked(itvol, markf))
    {
      ++(io.numberoffacets);
      lcc.template mark_cell<2>(itvol, markf);
    }
  }

  lcc.negate_mark(treatedv);
  io.pointlist=new REAL[io.numberofpoints*3];
  std::size_t nb=0;
  std::vector<typename LCC::Dart_const_handle> faces;
  faces.reserve(io.numberoffacets);
  for(auto itvol=lcc.template darts_of_cell_basic<3>(dh, treatedv).begin(),
      itvolend=lcc.template darts_of_cell_basic<3>(dh, treatedv).end();
      itvol!=itvolend; ++itvol)
  {
    lcc.mark(itvol, treatedv);
    if(lcc.is_marked(itvol, markv))
    { lcc_to_tetgen_process_vertex(lcc, itvol, markv, io, nb++, vertices_map); }
    if(lcc.is_marked(itvol, markf))
    {
      faces.push_back(itvol);
      lcc.template unmark_cell<2>(itvol, markf);
    }
  }
  lcc.negate_mark(treatedv);

  io.facetlist=new tetgenio::facet[io.numberoffacets];
  nb=0;
  for(auto dh: faces)
 { lcc_to_tetgen_process_face(lcc, dh, markf, io, nb++, vertices_map); }

  assert(lcc.is_whole_map_unmarked(markv));
  assert(lcc.is_whole_map_unmarked(markf));
  assert(lcc.is_whole_map_unmarked(treatedv));
  lcc.free_mark(markv);
  lcc.free_mark(markf);
  lcc.free_mark(treatedv);
}
////////////////////////////////////////////////////////////////////////////////
/** Build a tetgenio from an LCC, using only marked volumes
 *   (a volume is considered marked if one of its dart is marked)
 *  Face between 2 marked volumes are ignored (to satisty tetgen constraint)
 *  and thus only preserve the external boundary of the marked volumes.
 * @pre all the faces of marked volumes must be triangles
 *       (we can call constrained_delaunay_triangulation as pre-process).
 */
template <typename LCC>
void lcc_to_tetgen(const LCC& lcc, typename LCC::size_type amark, tetgenio& io)
{
  typedef typename LCC::Dart_const_handle DH;
  typedef typename LCC::Vertex_attribute_const_handle VH;

  // All indices start from 0.
  io.firstnumber=0;

  // First copy vertices
  std::unordered_map<VH, std::size_t> vertices_map;
  auto markv=lcc.get_new_mark();
  auto markf=lcc.get_new_mark();
  auto treatedv=lcc.get_new_mark();

  io.numberoffacets=0;
  // 1) We mark all faces
  for(auto it=lcc.darts().begin(), itend=lcc.darts().end(); it!=itend; ++it)
  {
    if(lcc.is_marked(it, amark) && !lcc.is_marked(it, treatedv))
    {
      for(auto itvol=lcc.template darts_of_cell_basic<3>(it, treatedv).begin(),
          itvolend=lcc.template darts_of_cell_basic<3>(it, treatedv).end();
          itvol!=itvolend; ++itvol)
      {
        lcc.mark(itvol, treatedv);
        if(!lcc.is_marked(itvol, markf))
        {
          ++(io.numberoffacets);
          lcc.template mark_cell<2,2>(itvol, markf);
        }
      }
    }
  }

  // 2) We unmark faces that are incident to two marked volumes,
  //    and mark and count vertices incident fo marked faces.
  io.numberofpoints=0;
  for(auto it=lcc.darts().begin(), itend=lcc.darts().end(); it!=itend; ++it)
  {
    if(lcc.is_marked(it, markf))
    {
      if(lcc.is_marked(lcc.template beta<3>(it), markf))
      {
        lcc.template unmark_cell<2>(it, markf);
        (io.numberoffacets)-=2;
      }
      else
      {
        DH cur=it;
        do
        {
          if(!lcc.is_marked(cur, markv))
          {
            ++(io.numberofpoints);
            lcc.template mark_cell<0>(cur, markv);
          }
          cur=lcc.next(cur);
        }
        while(cur!=it);
      }
    }
  }

  // 3) We create the array of points and fill it.
  io.pointlist=new REAL[io.numberofpoints*3];
  std::size_t nb=0;
  for(auto it=lcc.darts().begin(), itend=lcc.darts().end(); it!=itend; ++it)
  {
    lcc.unmark(it, treatedv);
    if(lcc.is_marked(it, markv))
    { lcc_to_tetgen_process_vertex(lcc, it, markv, io, nb++, vertices_map); }
  }

  // 4) We create the array of faces and fill it.
  io.facetlist=new tetgenio::facet[io.numberoffacets];
  nb=0;
  for(auto it=lcc.darts().begin(), itend=lcc.darts().end(); it!=itend; ++it)
   {
    if(lcc.is_marked(it, markf))
    { lcc_to_tetgen_process_face(lcc, it, markf, io, nb++, vertices_map); }
  }

  assert(lcc.is_whole_map_unmarked(markv));
  assert(lcc.is_whole_map_unmarked(markf));
  assert(lcc.is_whole_map_unmarked(treatedv));
  lcc.free_mark(markv);
  lcc.free_mark(markf);
  lcc.free_mark(treatedv);
}
////////////////////////////////////////////////////////////////////////////////
/** Build a tetgenio from the entire given LCC. Preserve only the external
 *  surface of the volumic mesh (and not the internal faces).
 * @pre all the faces of marked volumes must be triangles
 *       (we can call constrained_delaunay_triangulation as pre-process).
 */
template <typename LCC>
void lcc_to_tetgen(const LCC& lcc, tetgenio& io)
{
  auto amark=lcc.get_new_mark();
  lcc.negate_mark(amark);
  lcc_to_tetgen(lcc, amark, io);
  lcc.negate_mark(amark);
}
////////////////////////////////////////////////////////////////////////////////
template <typename LCC>
typename LCC::Dart_handle dart_of_ith_face_of_tetra(LCC& lcc,
                                                    typename LCC::Dart_handle dh,
                                                    int i)
{
  switch(i)
  {
    case 3: return dh;
    case 1: return lcc.template beta<2>(dh);
    case 0: return lcc.template beta<1,2>(dh);
    case 2: return lcc.template beta<0,2>(dh);
    default: break;
  }
  std::cerr<<"ERROR dart_of_ith_face_of_tetra"<<std::endl;
  return nullptr;
}
////////////////////////////////////////////////////////////////////////////////
// Return the dart of the volume containing dh2 that matches dh1
template <typename LCC>
typename LCC::Dart_handle find_dart_that_match(LCC& lcc,
                                               typename LCC::Dart_handle dh1,
                                               typename LCC::Dart_handle dh2)
{
  typename LCC::Vertex_attribute_handle vh1=lcc.vertex_attribute(dh1);
  typename LCC::Vertex_attribute_handle vh2=lcc.vertex_attribute(lcc.next(dh1));
  typename LCC::Vertex_attribute_handle vh3=lcc.vertex_attribute(lcc.previous(dh1));

  /*if (lcc.template attributes<0>().index(vh1)==0 ||
      lcc.template attributes<0>().index(vh2)==0 ||
      lcc.template attributes<0>().index(vh3)==0)
    std::cout<<"Face to match: "<<lcc.template attributes<0>().index(vh1)<<" "
             <<lcc.template attributes<0>().index(vh2)<<" "
             <<lcc.template attributes<0>().index(vh3)<<std::endl;*/

  for (auto it=lcc.template darts_of_cell<3>(dh2).begin(),
       itend=lcc.template darts_of_cell<3>(dh2).end(); it!=itend; ++it)
  {
    /*if (lcc.template attributes<0>().index(vh1)==0 ||
        lcc.template attributes<0>().index(vh2)==0 ||
        lcc.template attributes<0>().index(vh3)==0)
      std::cout<<"     tested: "<<lcc.template attributes<0>().index(lcc.vertex_attribute(lcc.next(it)))<<" "
               <<lcc.template attributes<0>().index(lcc.vertex_attribute(it))<<" "
               <<lcc.template attributes<0>().index(lcc.vertex_attribute(lcc.previous(it)))<<std::endl;*/

    if (vh1==lcc.vertex_attribute(lcc.next(it)) &&
        vh2==lcc.vertex_attribute(it) &&
        vh3==lcc.vertex_attribute(lcc.previous(it)))
    { return it; }
  }
  std::cerr<<"ERROR in find_dart_that_match."<<std::endl;
  return nullptr;
}
////////////////////////////////////////////////////////////////////////////////
template <typename LCC>
void tetgen_to_lcc(const tetgenio& io, LCC& lcc)
{
  if (io.numberofcorners!=4)
  {
    std::cout<<"Conversion from tetgen to lcc for "<<io.numberofcorners
             <<" not implemented."<<std::endl;
    return;
  }

  std::vector<typename LCC::Vertex_attribute_handle> TV;
  TV.resize(io.numberofpoints);
  for (int i=0; i<io.numberofpoints; ++i)
  {
    TV[i]=lcc.create_vertex_attribute(typename LCC::Point(io.pointlist[3*i],
                                                          io.pointlist[(3*i)+1],
                                                          io.pointlist[(3*i)+2]));
  }

  typename LCC::Dart_handle dh, nh;
  std::vector<typename LCC::Dart_handle> TH; // One dart per tetrahedra
  TH.resize(io.numberoftetrahedra);
  for (int i=0; i<io.numberoftetrahedra; ++i)
  {
    /*if (lcc.template attributes<0>().index(TV[io.tetrahedronlist[i*4]])==0 ||
        lcc.template attributes<0>().index(TV[io.tetrahedronlist[(i*4)+1]])==0 ||
        lcc.template attributes<0>().index(TV[io.tetrahedronlist[(i*4)+2]])==0 ||
        lcc.template attributes<0>().index(TV[io.tetrahedronlist[(i*4)+3]])==0)
      std::cout<<"make_tetrahedron "
               <<lcc.template attributes<0>().index(TV[io.tetrahedronlist[i*4]])<<" "
               <<lcc.template attributes<0>().index(TV[io.tetrahedronlist[(i*4)+1]])<<" "
               <<lcc.template attributes<0>().index(TV[io.tetrahedronlist[(i*4)+2]])<<" "
               <<lcc.template attributes<0>().index(TV[io.tetrahedronlist[(i*4)+3]])<<std::endl;*/

    TH[i]=lcc.make_tetrahedron(TV[io.tetrahedronlist[i*4]],
        TV[io.tetrahedronlist[(i*4)+2]],
        TV[io.tetrahedronlist[(i*4)+1]],
        TV[io.tetrahedronlist[(i*4)+3]]);
    for (int j=0; j<4; ++j)
    {
      if (io.neighborlist[(4*i)+j]!=-1 && io.neighborlist[(4*i)+j]<i) // jth neighboor of ith tetra
      { // We have one dart of this tetra
        dh=dart_of_ith_face_of_tetra(lcc, TH[i], j);
        nh=find_dart_that_match(lcc, dh, TH[io.neighborlist[(4*i)+j]]);
        assert(nh!=nullptr);
        lcc.template topo_sew<3>(nh, dh);
      }
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
#endif //  CGAL_LCC_TO_TETGEN_IO_H
////////////////////////////////////////////////////////////////////////////////
