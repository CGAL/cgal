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
#ifndef LCC_READ_DEPENDING_EXTENSION_H
#define LCC_READ_DEPENDING_EXTENSION_H

#include <filesystem>
#include <string>

#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Combinatorial_map_save_load.h>

#include "import_moka.h"
#include "lcc_read_from_vtk.h"
#include "lcc_save_load_mesh.h"
#include "lcc_save_load_with_assimp.h"

///////////////////////////////////////////////////////////////////////////////
bool is_lcc_readable_file(const std::string& filename)
{
  std::filesystem::path p(filename);
  std::string ext=p.extension().string();
  return ext==".3map" ||
      ext==".moka" ||
      ext==".off" ||
      ext==".mesh" ||
      ext==".toposim" ||
      ext==".vtk" ||
      ext==".obj" ||
      ext==".ply"; // More possibilities from collada ?
}
///////////////////////////////////////////////////////////////////////////////
template<typename LCC>
typename LCC::Dart_handle read_depending_extension(const std::string& filename,
                                                   LCC& lcc)
{
  if(!is_lcc_readable_file(filename)) { return nullptr; }

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
  else if(ext==".moka")
  { res=import_from_moka(lcc, filename.c_str()); }
  else if(ext==".off")
  { res=CGAL::load_off(lcc, filename.c_str()); }
  else if(ext==".mesh" || ext==".toposim")
  { res=load_object_3D(filename, lcc); }
  else if(ext==".vtk")
  { res=(read_vtk(filename, lcc)!=nullptr); }
#ifdef WITH_ASSIMP
  else // By default use load with collada
  { res=load_object_3D_with_assimp(filename, lcc); }
#endif // WITH_ASSIMP

  if(!res) { return nullptr; }
  if(amark==LCC::INVALID_MARK) { return lcc.darts().begin(); }

  typename LCC::Dart_handle resdh=nullptr;
  for(auto it=lcc.darts().begin(), itend=lcc.darts().end(); it!=itend; ++it)
  {
    if(lcc.is_marked(it, amark))
    {
      lcc.unmark(it, amark);
      if(resdh==nullptr) { resdh=it; }
    }
  }
  lcc.free_mark(amark);
  return resdh;
}
///////////////////////////////////////////////////////////////////////////////
#endif // LCC_READ_DEPENDING_EXTENSION_H
