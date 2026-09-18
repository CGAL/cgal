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
#ifndef LCC_SAVE_DIRECT_H
#define LCC_SAVE_DIRECT_H

#include <fstream>
#include <unordered_map>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>

namespace CGAL::IO
{
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
void build_faces_from_lcc(const LCC& lcc,
                          std::vector<typename LCC::Point>& points,
                          std::vector<std::vector<std::size_t>>& polygons)
{
  points.clear();
  polygons.clear();
  auto markface=lcc.get_new_mark();

  // All faces
  std::unordered_map<typename LCC::Vertex_attribute_const_descriptor,
                     std::size_t> index;
  for(auto it=lcc.darts().begin(), itend=lcc.darts().end(); it!=itend; ++it)
  {
    if (!lcc.is_marked(it, markface))
    {
      polygons.push_back(std::vector<std::size_t>());
      std::vector<std::size_t>& poly=polygons.back();
      std::size_t nb=0;
      typename LCC::Dart_const_descriptor cur=it;
      do
      {
        auto it=index.find(lcc.vertex_attribute(cur));
        if(it==index.end())
        {
          index[lcc.vertex_attribute(cur)]=points.size();
          nb=points.size();
          points.push_back(lcc.point(cur));
        }
        else
        { nb=it->second; }
        poly.push_back(nb);

        lcc.mark(cur, markface);
        if(!lcc.template is_free<3>(cur))
        { lcc.mark(lcc.template beta<3>(cur), markface); }
        cur=lcc.next(cur);
      }
      while(cur!=it);
    }
  }
  lcc.free_mark(markface);
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
bool save_lcc_into_gocad(const std::string& filename, const LCC& lcc)
{
  std::ofstream fo(filename.c_str());
  if(!fo.good()) { return false; } // File open error

  std::vector<typename LCC::Point> points;
  std::vector<std::vector<std::size_t>> polygons;
  build_faces_from_lcc(lcc, points, polygons);

  bool res=CGAL::IO::write_GOCAD(fo, points, polygons);
  fo.close();
  return res;
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
bool save_lcc_into_obj(const std::string& filename, const LCC& lcc)
{
  std::ofstream fo(filename.c_str());
  if(!fo.good()) { return false; } // File open error

  std::vector<typename LCC::Point> points;
  std::vector<std::vector<std::size_t>> polygons;
  build_faces_from_lcc(lcc, points, polygons);

  bool res=CGAL::IO::write_OBJ(fo, points, polygons);
  fo.close();
  return res;
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
bool save_lcc_into_off(const std::string& filename, const LCC& lcc)
{
  std::ofstream fo(filename.c_str());
  if(!fo.good()) { return false; } // File open error

  std::vector<typename LCC::Point> points;
  std::vector<std::vector<std::size_t>> polygons;
  build_faces_from_lcc(lcc, points, polygons);

  bool res=CGAL::IO::write_OFF(fo, points, polygons);
  fo.close();
  return res;
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
bool save_lcc_into_ply(const std::string& filename, const LCC& lcc)
{
  std::ofstream fo(filename.c_str());
  if(!fo.good()) { return false; } // File open error

  std::vector<typename LCC::Point> points;
  std::vector<std::vector<std::size_t>> polygons;
  build_faces_from_lcc(lcc, points, polygons);

  bool res=CGAL::IO::write_PLY(fo, points, polygons);
  fo.close();
  return res;
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
bool save_lcc_into_stl(const std::string& filename, const LCC& lcc)
{
  std::ofstream fo(filename.c_str());
  if(!fo.good()) { return false; } // File open error

  std::vector<typename LCC::Point> points;
  std::vector<std::vector<std::size_t>> polygons;
  build_faces_from_lcc(lcc, points, polygons);

  bool res=CGAL::IO::write_STL(fo, points, polygons);
  fo.close();
  return res;
}
///////////////////////////////////////////////////////////////////////////////
}
#endif // LCC_SAVE_DIRECT_H
