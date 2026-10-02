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
#ifndef LCC_READ_WRITE_DEPENDING_EXTENSION_H
#define LCC_READ_WRITE_DEPENDING_EXTENSION_H

#include <filesystem>
#include <string>

#include <CGAL/Combinatorial_map_save_load.h>
#include <CGAL/Linear_cell_complex_constructors.h>

#include <CGAL/Linear_cell_complex/IO/MOKA.h>
#include <CGAL/Linear_cell_complex/IO/lcc_read_write_voro.h>
#include <CGAL/Linear_cell_complex/IO/lcc_save_direct.h>
#include <CGAL/Linear_cell_complex/IO/lcc_save_load_mesh.h>
#include <CGAL/Linear_cell_complex/IO/lcc_save_load_surfaces.h>
#include <CGAL/Linear_cell_complex/IO/VTK.h>

namespace CGAL::IO
{
///////////////////////////////////////////////////////////////////////////////
inline
bool is_an_lcc_known_extension(const std::string& filename)
{
  std::filesystem::path p(filename);
  std::string ext=p.extension().string();
  return
    ext==".3map" ||    // My own format *
    ext==".gmsh" ||    // Gmsh format *
    ext==".gocad" ||   // cgal reader *
#ifdef USE_IFC
    ext==".ifc" ||     // Ifc format *
#endif // USE_IFC
    //ext==".medit" ||   // experimental reader medit files usually are .mesh but is conflicting here
    ext==".mesh" ||    // My own format *
    ext==".msh" ||     // Gmsh format *
    ext==".moka" ||    // Moka format *
    ext==".obj" ||     // cgal reader *
    ext==".off" ||     // cgal reader *
    ext==".ply" ||     // cgal reader *
    ext==".stl" ||     // cgal reader *
    ext==".toposim" || // "Toposim format"=mesh format *
    ext==".vol" ||     // voro++ format *
    ext==".voro" ||    // vem++ format *
    ext==".vtk";       // vtk my own reader *
}
///////////////////////////////////////////////////////////////////////////////
template<typename LCC>
typename LCC::Dart_descriptor read_depending_extension(const std::string& filename,
                                                       LCC& lcc)
{
  if(!is_an_lcc_known_extension(filename)) { return lcc.null_descriptor; }

  bool res=false;
  typename LCC::size_type amark=LCC::INVALID_MARK;
  if(!lcc.is_empty())
  {
    amark=lcc.get_new_mark();
    lcc.negate_mark(amark); // All old darts are marked
  }
  std::filesystem::path p(filename);
  std::string ext=p.extension().string();
  if(ext==".3map")
  { res=CGAL::load_combinatorial_map(filename.c_str(), lcc); }
  else if(ext==".gmsh" || ext==".msh")
  { res=read_object_3D_gmsh(filename, lcc); }
  else if(ext==".gocad")
  { res=lcc_read_gocad(lcc, filename); }
#ifdef USE_IFC
  else if(ext==".ifc")
  { res=(lcc_read_ifc(lcc, filename)!=lcc.null_descriptor); }
#endif // USE_IFC
  //else if(ext==".medit")// experimental reader medit files usually are .mesh but is conflicting here
  //{res = read_medit_as_lcc(filename, lcc);}
  else if(ext==".mesh" || ext==".toposim")
  { res=load_object_3D(filename, lcc); }
  else if(ext==".moka")
  { res=read_MOKA(lcc, filename.c_str()); }
  else if(ext==".obj")
  { res=lcc_read_obj(lcc, filename); }
  else if(ext==".off")
  { res=lcc_read_off(lcc, filename); }
  else if(ext==".ply")
  { res=lcc_read_ply(lcc, filename); }
  else if(ext==".stl")
  { res=lcc_read_stl(lcc, filename); }
  else if(ext==".vol")
  { res=(read_vol(filename, lcc)!=lcc.null_descriptor); }
  else if(ext==".voro")
  { res=(read_voro(filename, lcc)!=lcc.null_descriptor); }
  else if(ext==".vtk")
  { res=CGAL::IO::read_VTK(filename.c_str(), lcc);}

  if(!res) { return lcc.null_descriptor; }
  if(amark==LCC::INVALID_MARK) { return lcc.darts().begin(); }

  typename LCC::Dart_descriptor resdh=lcc.null_descriptor;
  for(auto it=lcc.darts().begin(), itend=lcc.darts().end(); it!=itend; ++it)
  {
    if(lcc.is_marked(it, amark))
    {
      lcc.unmark(it, amark);
      if(resdh==lcc.null_descriptor) { resdh=it; }
    }
  }
  lcc.free_mark(amark);
  return resdh;
}
///////////////////////////////////////////////////////////////////////////////
template<typename LCC>
bool write_depending_extension(const std::string& filename,
                               const LCC& lcc,
                               bool transform_lcc_into_manifold=false)
{
  if(!is_an_lcc_known_extension(filename)) { return false; }

  bool res=false;
  std::filesystem::path p(filename);
  std::string ext=p.extension().string();
  if(ext==".3map")
  { res=CGAL::save_combinatorial_map(lcc, filename.c_str()); }
  else if(ext==".gmsh" || ext==".msh")
  { res=save_object_3D_gmsh(filename, lcc); }
  else if(ext==".gocad")
  { res=save_lcc_surface_into_manifold_gocad(filename, lcc); }
#ifdef USE_IFC
  else if(ext==".ifc")
  { std::cerr<<"[write_depending_extension] this is not (yet) possible to "
             <<"write an lcc into an ifc file."<<std::endl; }
#endif // USE_IFC
  else if(ext==".mesh" || ext==".toposim")
  { res=save_object_3D(filename, lcc); }
  else if(ext==".moka")
  { res=write_MOKA(lcc, filename.c_str()); }
  else if(ext==".obj")
  {
    if(transform_lcc_into_manifold)
    { res=save_lcc_surface_into_manifold_obj(filename, lcc); }
    else
    { res=save_lcc_into_obj(filename, lcc); }
  }
  else if(ext==".off")
  {
    if(transform_lcc_into_manifold)
    { res=save_lcc_surface_into_manifold_off(filename, lcc); }
    else
    { res=save_lcc_into_off(filename, lcc); }
  }
  else if(ext==".ply")
  {
    if(transform_lcc_into_manifold)
    { res=save_lcc_surface_into_manifold_ply(filename, lcc); }
    else
    { res=save_lcc_into_ply(filename, lcc); }
  }
  else if(ext==".stl")
  {
    if(transform_lcc_into_manifold)
    { res=save_lcc_surface_into_manifold_stl(filename, lcc); }
    else
    { res=save_lcc_into_stl(filename, lcc); }
  }
  else if(ext==".vol")
  { res=write_vol(filename, lcc); }
  else if(ext==".voro")
  { res=write_voro(filename, lcc); }
  else if(ext==".vtk")
  { res=CGAL::IO::write_VTK(filename.c_str(), lcc); }

  return res;
}
///////////////////////////////////////////////////////////////////////////////
}
#endif // LCC_READ_WRITE_DEPENDING_EXTENSION_H
