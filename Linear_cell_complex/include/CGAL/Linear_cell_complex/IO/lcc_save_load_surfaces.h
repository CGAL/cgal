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
#ifndef LCC_SAVE_LOAD_SURFACES_H
#define LCC_SAVE_LOAD_SURFACES_H

#include <fstream>
#include <unordered_map>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/Linear_cell_complex_incremental_builder_3.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Surface_mesh.h>
#include "lcc_save_load_mesh.h"

namespace CGAL::IO
{
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
bool save_lcc_surface_into_off(const std::string& filename, const LCC& lcc)
{
  std::ofstream fo(filename.c_str());
  if(!fo.good()) return false; // File open error

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

  std::unordered_map<typename LCC::Vertex_attribute_const_descriptor,
                     std::size_t> index;

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
  return true;
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
void build_manifold_surface_from_lcc(const LCC& lcc,
                                     CGAL::Surface_mesh<typename LCC::Point>& mesh,
                                     bool orient_ccs=true)
{
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

  std::unordered_map<typename LCC::Vertex_attribute_const_descriptor,
                     std::size_t> index;
  std::vector<typename LCC::Point> points;
  std::vector<std::vector<std::size_t>> polygons;

  // all vertices
  std::size_t nb=0;
  // For this loop, we cannot iterate through vertex attributes since they may
  // not contain one of its dart. In such a case, it is not possible to ummark
  // the 3free vertices.
  points.reserve(nbvertices);
  for(auto it=lcc.darts().begin(), itend=lcc.darts().end(); it!=itend; ++it)
  {
    if (lcc.is_marked(it, vertex_3free))
    {
      points.push_back(lcc.point(it));
      index[lcc.vertex_attribute(it)]=nb++;
      lcc.template unmark_cell<0>(it, vertex_3free);
    }
  }

  // all 3free faces
  polygons.resize(nbfaces);
  nbfaces=0;
  for(auto it=lcc.darts().begin(), itend=lcc.darts().end(); it!=itend; ++it)
  {
    if (lcc.is_marked(it, face_3free))
    {
      nb=0;
      typename LCC::Dart_const_descriptor cur=it;
      do
      { ++nb; cur=lcc.next(cur); }
      while(cur!=it);
      polygons[nbfaces].reserve(nb);
      do
      {
        polygons[nbfaces].push_back(index[lcc.vertex_attribute(cur)]);
        lcc.unmark(cur, face_3free);
        cur=lcc.next(cur);
      }
      while(cur!=it);
      ++nbfaces;
    }
  }
  lcc.free_mark(vertex_3free);
  lcc.free_mark(face_3free);

  mesh.clear();
  std::vector<CGAL::Surface_mesh<typename LCC::Point>> tabmeshes;
  CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons);
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, mesh);

  if(orient_ccs)
  {
    CGAL::Polygon_mesh_processing::split_connected_components(mesh, tabmeshes);
    mesh.clear();
    for(auto& m: tabmeshes)
    {
      if (CGAL::is_closed(m) && (!CGAL::Polygon_mesh_processing::is_outward_oriented(m)))
      { CGAL::Polygon_mesh_processing::reverse_face_orientations(m); }
      CGAL::copy_face_graph(m, mesh);
    }
  }
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
bool save_lcc_surface_into_manifold_gocad(const std::string& filename,
                                          const LCC& lcc)
{
  std::ofstream fo(filename.c_str());
  if(!fo.good()) { return false; } // File open error

  CGAL::Surface_mesh<typename LCC::Point> mesh;
  build_manifold_surface_from_lcc(lcc, mesh);

  bool res=CGAL::IO::write_GOCAD(fo, mesh);
  fo.close();
  return res;
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
bool save_lcc_surface_into_manifold_obj(const std::string& filename,
                                        const LCC& lcc)
{
  std::ofstream fo(filename.c_str());
  if(!fo.good()) { return false; } // File open error

  CGAL::Surface_mesh<typename LCC::Point> mesh;
  build_manifold_surface_from_lcc(lcc, mesh);

  bool res=CGAL::IO::write_OBJ(fo, mesh);
  fo.close();
  return res;
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
bool save_lcc_surface_into_manifold_off(const std::string& filename,
                                        const LCC& lcc)
{
  std::ofstream fo(filename.c_str());
  if(!fo.good()) { return false; } // File open error

  CGAL::Surface_mesh<typename LCC::Point> mesh;
  build_manifold_surface_from_lcc(lcc, mesh);

  fo<<mesh;
  fo.close();
  return true;
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
bool save_lcc_surface_into_manifold_ply(const std::string& filename,
                                        const LCC& lcc)
{
  std::ofstream fo(filename.c_str());
  if(!fo.good()) { return false; } // File open error

  CGAL::Surface_mesh<typename LCC::Point> mesh;
  build_manifold_surface_from_lcc(lcc, mesh);

  bool res=CGAL::IO::write_PLY(fo, mesh);
  fo.close();
  return res;
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
bool save_lcc_surface_into_manifold_stl(const std::string& filename,
                                        const LCC& lcc)
{
  std::ofstream fo(filename.c_str());
  if(!fo.good()) { return false; } // File open error

  CGAL::Surface_mesh<typename LCC::Point> mesh;
  build_manifold_surface_from_lcc(lcc, mesh);

  bool res=CGAL::IO::write_STL(fo, mesh);
  fo.close();
  return res;
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
bool save_lcc_surface_into_manifold_vtk(const std::string& filename,
                                        const LCC& lcc)
{
  std::ofstream fo(filename.c_str());
  if (!fo.good())
  {
     // File open error
    return false;
  }

  fo << "# vtk DataFile Version 3.0\n";
  fo << "LCC Surface Mesh\n";
  fo << "ASCII\n";
  fo << "DATASET POLYDATA\n";
  std::vector<unsigned int> count = lcc.count_all_cells();
  fo << "POINTS " << count[0] << " float\n";
  fo << std::setprecision(10);

  // Dump point coordinates and build index map
  std::unordered_map<typename LCC::Vertex_attribute_handle, std::size_t> index;
  std::size_t nb = 0;
  typename LCC::Dart_descriptor sd;
  for (auto itv = lcc.vertex_attributes().begin(), itvend = lcc.vertex_attributes().end();
       itv != itvend; ++itv)
  {
    fo << itv->point() << "\n";
    index[itv] = nb++;
  }

  // First pass on faces to count the total size of the cell indices
  std::size_t data_size = 0;
  for (auto it = lcc.template one_dart_per_cell<2>().begin();
       it != lcc.template one_dart_per_cell<2>().end(); ++it)
  {
    data_size += 1 + lcc.template darts_of_orbit<1>(it).size();
  }

  fo << "POLYGONS " << count[2] << " " << data_size << "\n";

  // Second pass on faces to write the cell indices
  for (auto it = lcc.template one_dart_per_cell<2>().begin();
       it != lcc.template one_dart_per_cell<2>().end(); ++it)
  {
    auto darts = lcc.template darts_of_orbit<1>(it);
    fo << darts.size() << " ";
    for (auto iti = darts.begin(); iti != darts.end(); ++iti)
    {
      fo << index[lcc.vertex_attribute(lcc.template beta<0>(iti))] << " ";
    }
    fo << "\n";
  }

  fo.close();
  return true;
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
void build_lcc_from_polygons(LCC& lcc,
                             const std::vector<typename LCC::Point>& points,
                             const std::vector<std::vector<std::size_t>>& polygons)
{
  CGAL::Linear_cell_complex_incremental_builder_3<LCC> B(lcc);
  B.begin_surface();
  for(auto& p: points)
  { B.add_vertex(p); }
  for(auto& p: polygons)
  {
    B.begin_facet();
    for(auto i: p)
    { B.add_vertex_to_facet(static_cast<typename LCC::size_type>(i)); }
    B.end_facet();
  }
  B.end_surface();
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
bool lcc_read_obj(LCC& lcc, const std::string& name)
{
  std::vector<typename LCC::Point> points;
  std::vector<std::vector<std::size_t>> polygons;
  if(!CGAL::IO::read_OBJ(name, points, polygons))
  { return false; }
  build_lcc_from_polygons(lcc, points, polygons);
  return true;
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
bool lcc_read_off(LCC& lcc, const std::string& name)
{
  std::vector<typename LCC::Point> points;
  std::vector<std::vector<std::size_t>> polygons;
  if(!CGAL::IO::read_OFF(name, points, polygons))
  { return false; }
  build_lcc_from_polygons(lcc, points, polygons);
  return true;
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
bool lcc_read_ply(LCC& lcc, const std::string& name)
{
  std::vector<typename LCC::Point> points;
  std::vector<std::vector<std::size_t>> polygons;
  if(!CGAL::IO::read_PLY(name, points, polygons))
  { return false; }
  build_lcc_from_polygons(lcc, points, polygons);
  return true;
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
bool lcc_read_gocad(LCC& lcc, const std::string& name)
{
  std::vector<typename LCC::Point> points;
  std::vector<std::vector<std::size_t>> polygons;
  if(!CGAL::IO::read_GOCAD(name, points, polygons))
  { return false; }
  build_lcc_from_polygons(lcc, points, polygons);
  return true;
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
bool lcc_read_stl(LCC& lcc, const std::string& name)
{
  std::vector<typename LCC::Point> points;
  std::vector<std::vector<std::size_t>> polygons;
  if(!CGAL::IO::read_STL(name, points, polygons))
  { return false; }
  build_lcc_from_polygons(lcc, points, polygons);
  return true;
}
///////////////////////////////////////////////////////////////////////////////
}
#endif // LCC_SAVE_LOAD_SURFACES_H
