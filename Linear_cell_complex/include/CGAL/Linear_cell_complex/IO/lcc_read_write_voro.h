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
#ifndef LCC_READ_WRITE_VORO_H
#define LCC_READ_WRITE_VORO_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <CGAL/Linear_cell_complex_incremental_builder_3.h>
#include <CGAL/Linear_cell_complex/lcc_sew3_similar_facets.h>

namespace CGAL::IO
{
// Mesh format for VORO++ (or a variant for VEM++)
// * vem++
//  -> VoronoiCell X
//     nbv // number of vertices
//     x0 y0 z0
//     x1 y1 z1
//     ...
//     nbf // number of faces
//     iv0_0 iv1_0 iv2_0 ... #
//     iv0_1 iv1_1 iv2_1 ... #
//     ...
//     id0 id1 ... id(nbf-1) // tags for faces, negative if border face
//     Then another VoronoiCell X ...
//
// * voro++
// -> exemple ./voro++ -c "%w %P %s %t" -10 10 -10 10 -10 10 pack_irregular
// One cell per line:
// 22 (-7.38612,-4.79217,-10) (7.96003,-10,-9.55694) (0.769875,3.19235,-10) (0.285365,-0.021303,-2.823) (-3.37854,-0.288531,-4.77004) (-6.07213,-3.52433,-8.1543) (-1.80528,-0.422809,-3.3515) (-8.44931,-1.53507,-10) (-6.61129,0.977114,-10) (2.80333,-4.22525,-4.51965) (8.44242,-10,-10) (-7.18355,-1.34615,-8.58203) (-2.96211,-7.04195,-10) (6.74654,-10,-10) (-0.140831,-0.269506,-2.54162) (1.03397,-1.47014,-2.75092) (-1.58578,-0.608579,-3.19964) (6.01626,-1.06559,-10) (-0.11516,-0.581078,-2.46046) (1.53373,-0.856098,-3.17792) (2.18414,-1.07794,-3.78686) (8.50369,-8.92064,-10) 13 (1,9,15,19,20,21,10) (1,13,12,9) (1,10,13) (2,17,20,19,3) (2,8,7,0,12,13,10,21,17) (2,3,14,6,4,8) (3,19,15,18,14) (4,11,7,8) (4,6,16,5,11) (5,0,7,11) (5,16,18,15,9,12,0) (6,14,18,16) (17,21,20)
// nbv (x0,y0,z0) (x1,y1,z1) ... nbf (iv0_0,iv1_0,iv2_0,...)(iv0_1,iv1_1,iv2_1,...)
//
////////////////////////////////////////////////////////////////////////////////
template<typename LCC>
typename LCC::Dart_descriptor read_voro(std::ifstream& file, LCC& lcc)
{
  typename LCC::Dart_descriptor res=lcc.null_descriptor;
  std::string s;
  std::size_t n, index;
  double x, y, z;
  while(file.good())
  {
    file>>s;
    if(s!="VoronoiCell")
    {
      std::cerr<<"[ERROR] read_voro"<<std::endl;
      break;
    }
    file>>n; // id of voronoi cell; ignored
    CGAL::Linear_cell_complex_incremental_builder_3<LCC> ib(lcc);
    ib.begin_surface();
    file>>n; // number of points
    for(std::size_t i=0; i<n; ++i)
    {
      file>>x>>y>>z;
      ib.add_vertex(typename LCC::Point(x, y, z));
    }

    file>>n; // number of faces
    for(std::size_t i=0; i<n; ++i)
    {
      do
      {
        std::getline(file, s); // indices of faces, end by #
      }
      while((s=="" || s.find_first_not_of(' ')==std::string::npos) &&
             file.good());
      if(file.good())
      {
        std::stringstream ss(s);
        ib.begin_facet();
        while(ss.good())
        {
          ss>>index;
          if(ss.good())
          { // it was an index and not #
            ib.add_vertex_to_facet(index);
          }
        }
        ib.end_facet();
      }
    }

    for(std::size_t i=0; i<n; ++i)
    { file>>x; } // Ignore these values (serve to tag some faces, negative for border faces)

    typename LCC::Dart_descriptor d=ib.end_surface();
    if(res==lcc.null_descriptor) { res=d; }
  }

  sew3_similar_facets(lcc, 1e-9);
  return res;
}////////////////////////////////////////////////////////////////////////////////
template<typename LCC>
typename LCC::Dart_descriptor read_vol(std::ifstream& file, LCC& lcc)
{
  typename LCC::Dart_descriptor res=lcc.null_descriptor;
  std::string s;
  char c;
  std::size_t n, index;
  double x, y, z;
  file>>n; // number of points (of the first cell)
  while(file.good())
  {
    CGAL::Linear_cell_complex_incremental_builder_3<LCC> ib(lcc);
    ib.begin_surface();
    for(std::size_t i=0; i<n; ++i)
    {
      file>>c>>x>>c>>y>>c>>z>>c; // (x,y,z)
      ib.add_vertex(typename LCC::Point(x, y, z));
    }

    file>>n; // number of faces
    for(std::size_t i=0; i<n; ++i)
    {
      file>>c; // (
      ib.begin_facet();
      do
      {
        file>>index>>c; // index then , or )
        ib.add_vertex_to_facet(index);
      }
      while(c!=')' && file.good());
      ib.end_facet();
    }

    typename LCC::Dart_descriptor d=ib.end_surface();
    if(res==lcc.null_descriptor) { res=d; }

    file>>n; // number of points of the next cell (or bad bit on the file)
  }

  sew3_similar_facets(lcc, 1e-9);
  return res;
}
////////////////////////////////////////////////////////////////////////////////
template<typename LCC>
typename LCC::Dart_descriptor read_voro(const std::string& filename, LCC& lcc)
{
  std::ifstream file(filename);
  if(!file.is_open())
  {
    std::cerr<<"[ERROR] read_voro: cannot open file "<<filename<<std::endl;
    return lcc.null_descriptor;
  }

  return read_voro(file, lcc);
}
////////////////////////////////////////////////////////////////////////////////
template<typename LCC>
typename LCC::Dart_descriptor read_vol(const std::string& filename, LCC& lcc)
{
  std::ifstream file(filename);
  if(!file.is_open())
  {
    std::cerr<<"[ERROR] read_vol: cannot open file "<<filename<<std::endl;
    return lcc.null_descriptor;
  }

  return read_vol(file, lcc);
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
void write_voro_one_cell(std::ofstream& fo, const LCC& lcc,
                         typename LCC::Dart_const_descriptor dh, std::size_t nb)
{
  fo<<"VoronoiCell "<<nb<<std::endl;
  fo<<lcc.template one_dart_per_incident_cell<0,3>(dh).size()<<std::endl;

  std::unordered_map<typename LCC::Vertex_attribute_const_descriptor,
                     std::size_t> pmap;
  nb=0;
  for(auto it=lcc.template one_dart_per_incident_cell<0,3>(dh).begin(),
        itend=lcc.template one_dart_per_incident_cell<0,3>(dh).end();
      it!=itend; ++it, ++nb)
  {
    pmap[lcc.vertex_attribute(it)]=nb;
    fo<<lcc.point(it)<<std::endl;
  }
  fo<<std::endl;

  fo<<lcc.template one_dart_per_incident_cell<2,3>(dh).size()<<std::endl;
  for(auto it=lcc.template one_dart_per_incident_cell<2,3>(dh).begin(),
        itend=lcc.template one_dart_per_incident_cell<2,3>(dh).end();
      it!=itend; ++it, ++nb)
  {
    typename LCC::Dart_const_descriptor cur=it;
    do
    {
      fo<<pmap[lcc.vertex_attribute(cur)]<<" ";
      cur=lcc.template beta<1>(cur);
    }
    while(cur!=it);
    fo<<"#"<<std::endl;
  }
  fo<<std::endl;

  nb=0;
  for(auto it=lcc.template one_dart_per_incident_cell<2,3>(dh).begin(),
       itend=lcc.template one_dart_per_incident_cell<2,3>(dh).end();
       it!=itend; ++it, ++nb)
  {
    // Not sure this is correct !! (but for now I don't use these tags...)
    if(lcc.template is_free<3>(it)) { fo<<-nb<<" "; }
    else { fo<<nb<<" "; }
  }
  fo<<std::endl<<std::endl;
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
bool write_voro(std::ofstream& fo, const LCC& lcc)
{
  std::size_t nbcells=0;
  for(auto itvol=lcc.template one_dart_per_cell<3>().begin(),
       itvolend=lcc.template one_dart_per_cell<3>().end();
       itvol!=itvolend; ++itvol, ++nbcells)
  { write_voro_one_cell(fo, lcc, itvol, nbcells); }

  return true;
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
bool write_voro(const std::string& filename, const LCC& lcc)
{
  std::ofstream fo(filename.c_str());
  if(!fo.good()) return false; // File open error

  return write_voro(fo, lcc);
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
void write_vol_one_cell(std::ofstream& fo, const LCC& lcc,
                         typename LCC::Dart_const_descriptor dh, std::size_t nb)
{
  fo<<lcc.template one_dart_per_incident_cell<0,3>(dh).size()<<" ";

  std::unordered_map<typename LCC::Vertex_attribute_const_descriptor,
                     std::size_t> pmap;
  nb=0;
  for(auto it=lcc.template one_dart_per_incident_cell<0,3>(dh).begin(),
       itend=lcc.template one_dart_per_incident_cell<0,3>(dh).end();
       it!=itend; ++it, ++nb)
  {
    pmap[lcc.vertex_attribute(it)]=nb;
    fo<<"("<<lcc.point(it).x()<<","<<lcc.point(it).y()<<","
      <<lcc.point(it).z()<<") ";
  }

  fo<<lcc.template one_dart_per_incident_cell<2,3>(dh).size()<<" ";
  for(auto it=lcc.template one_dart_per_incident_cell<2,3>(dh).begin(),
       itend=lcc.template one_dart_per_incident_cell<2,3>(dh).end();
       it!=itend; ++it, ++nb)
  {
    fo<<"(";
    typename LCC::Dart_const_descriptor cur=it;
    do
    {
      fo<<pmap[lcc.vertex_attribute(cur)];
      cur=lcc.template beta<1>(cur);
      if(cur==it) { fo<<") "; }
      else { fo<<","; }
    }
    while(cur!=it);
  }
  fo<<std::endl;
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
bool write_vol(std::ofstream& fo, const LCC& lcc)
{
  std::size_t nbcells=0;
  for(auto itvol=lcc.template one_dart_per_cell<3>().begin(),
       itvolend=lcc.template one_dart_per_cell<3>().end();
       itvol!=itvolend; ++itvol, ++nbcells)
  { write_vol_one_cell(fo, lcc, itvol, nbcells); }

  return true;
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
bool write_vol(const std::string& filename, const LCC& lcc)
{
  std::ofstream fo(filename.c_str());
  if(!fo.good()) return false; // File open error

  return write_vol(fo, lcc);
}
///////////////////////////////////////////////////////////////////////////////
}
#endif // LCC_READ_WRITE_VORO_H
