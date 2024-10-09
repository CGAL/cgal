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
// Prog that uses the query/replace method to create an hexahedral mesh from
// a tetrahedral one.

#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_operations.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/draw_linear_cell_complex.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/aff_transformation_tags.h>

#include <map>
#include <sys/types.h>
#include <unistd.h>

#include "cmap_copy.h"
#include "cmap_query_replace.h"
#include "Compute_stats.h"
#include "init_to_preserve_for_query_replace.h"
#include "lcc_read_depending_extension.h"
#include "lcc_save_load_mesh.h"
#include "Print_txt.h"
#include "Tetrahedral_tools.h"

////////////////////////////////////////////////////////////////////////////////
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT       FT;
typedef Kernel::Point_3  Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::Linear_cell_complex_for_combinatorial_map<3,3> LCC3;
typedef typename LCC3::Dart_handle Dart_handle;
typedef typename LCC3::Vertex_attribute_handle Vertex_handle;
typedef typename LCC3::size_type   size_type;

////////////////////////////////////////////////////////////////////////////////
[[ noreturn ]] void usage(int /*argc*/, char** argv)
{
  // Name
  std::cout<<"Name"<<std::endl;
  std::cout<<"        "<<argv[0]<<" - subdivides the given off file in hexahedra.";
  std::cout<<std::endl<<std::endl;
  // Synopsis
  std::cout<<"SYNOPSIS"<<std::endl;
  std::cout<<"        "<<argv[0]<<" [--help|-h|-?] "
           <<"[-draw] [-save] filename"
           <<std::endl<<std::endl;
  // Description
  std::cout<<"DESCRIPTION"<<std::endl;
  std::cout<<"        "<<" subdivides the given off file in hexahedra, by first computing a tetrahedral mesh, then transform all tetrahedra into hexahedra."
           <<std::endl<<std::endl;
  // Options
  std::cout<<"        --help, -h, -?"<<std::endl
           <<"                display this help and exit."
           <<std::endl<<std::endl;
  std::cout<<"        -draw"<<std::endl
           <<"                draw the final lcc."
           <<std::endl<<std::endl;
  std::cout<<"        -save"<<std::endl
           <<"                save the lcc resulting of the subdvision (format mesh)."
           <<std::endl;
  exit(EXIT_FAILURE);
}
////////////////////////////////////////////////////////////////////////////////
[[ noreturn ]] void error_command_line(int argc, char** argv, const char* msg)
{
  std::cout<<"ERROR: "<<msg<<std::endl;
  usage(argc, argv);
}
////////////////////////////////////////////////////////////////////////////////
void process_command_line(int argc, char** argv,
                          std::string& filename,
                          bool& draw,
                          bool& save
                          )
{
  filename="";
  draw=false;
  save=false;

  bool helprequired=false;
  std::string arg;
  for (int i=1; i<argc; ++i)
  {
    arg=std::string(argv[i]);
    if(arg==std::string("-h") || arg==std::string("--help") || arg==std::string("-?"))
    { helprequired=true; }
    else if(arg=="-draw")
    { draw=true; }
    else if(arg=="-save")
    { save=true; }
    else if(arg[0]=='-')
    { std::cout<<"Unknown option "<<arg<<", ignored."<<std::endl; }
    else { filename=arg; }
  }
  if (helprequired || filename.empty()) { usage(argc, argv); }
}
////////////////////////////////////////////////////////////////////////////////
void subdivide_edges(LCC3& lcc)
{
  std::vector<Dart_handle> edges_to_subdivide;
  edges_to_subdivide.reserve(lcc.number_of_darts()/3); // a la louche ;)

  for(auto itd=lcc.one_dart_per_cell<1>().begin(),
        itdend=lcc.one_dart_per_cell<1>().end(); itd!=itdend; ++itd)
  { edges_to_subdivide.push_back(itd); }

  for(Dart_handle itd: edges_to_subdivide)
  { lcc.insert_barycenter_in_cell<1>(itd); }
}
////////////////////////////////////////////////////////////////////////////////
void subdivide_faces(LCC3& lcc,
                     size_type corner_mark,
                     Pattern_substituer<LCC3>& ps)
{
  std::vector<Dart_handle> faces_to_subdivide;
  faces_to_subdivide.reserve(lcc.number_of_darts()/9); // a la louche ;)
  for(auto itd=lcc.one_dart_per_cell<2>().begin(),
        itdend=lcc.one_dart_per_cell<2>().end(); itd!=itdend; ++itd)
  {
    if(lcc.is_marked(itd, corner_mark))
    { faces_to_subdivide.push_back(itd); }
    else
    {
      assert(lcc.is_marked(lcc.beta<0>(itd), corner_mark));
      faces_to_subdivide.push_back(lcc.beta<0>(itd));
    }
  }

  for(Dart_handle itd: faces_to_subdivide)
  {
    ps.replace_one_face_from_dart(lcc, itd, ps.m_fpatterns[0],
                                  ps.m_fsignatures.begin()->second.first);
  }
}
////////////////////////////////////////////////////////////////////////////////
void subdivide_volumes(LCC3& lcc,
                       size_type corner_mark,
                       Pattern_substituer<LCC3>& ps)
{
  std::vector<Dart_handle> vols_to_subdivide;
  vols_to_subdivide.reserve(lcc.number_of_darts()/24);
  for(auto itd=lcc.one_dart_per_cell<3>().begin(),
        itdend=lcc.one_dart_per_cell<3>().end(); itd!=itdend; ++itd)
  {
    if(lcc.is_marked(itd, corner_mark))
    { vols_to_subdivide.push_back(itd); }
    else
    {
      assert(lcc.is_marked(lcc.beta<0>(itd), corner_mark));
      vols_to_subdivide.push_back(lcc.beta<0>(itd));
    }
  }

  for(Dart_handle itd: vols_to_subdivide)
  {
    ps.replace_one_volume_from_dart(lcc, itd, ps.m_vpatterns[0],
                                    ps.m_vsignatures.begin()->second.first);
  }
}
////////////////////////////////////////////////////////////////////////////////
bool create_mesh_from_off(const std::string& filename,
                          bool draw,
                          bool save
                          )
{
  std::chrono::system_clock::time_point start=std::chrono::system_clock::now();
  std::chrono::system_clock::time_point start2=std::chrono::system_clock::now();

  LCC3 lcc;
  size_type corner_mark=lcc.get_new_mark();
  if(read_depending_extension(filename, lcc)==nullptr)
  {
    std::cout<<"[ERROR] problem when reading file "<<filename<<std::endl;
    return false;
  }
  tetrahedralize_with_tetgen(lcc);
  lcc.negate_mark(corner_mark); // At the beginning, all darts are corners
  std::chrono::duration<double> diff=std::chrono::system_clock::now()-start2;
  print_txt_with_endl("    Tetrahedral mesh computation: ", diff, "s.");

  start2=std::chrono::system_clock::now();

  Pattern_substituer<LCC3> ps1;
  ps1.load_fpatterns(CGAL::data_file_path("query_replace/tetra-to-hexa/fpattern/"),
                     mark_fpattern_corners<LCC3>);
  ps1.load_vpatterns(CGAL::data_file_path("query_replace/tetra-to-hexa/vpattern/"),
                     mark_vpattern_corners<LCC3>);

  subdivide_edges(lcc);
  subdivide_faces(lcc, corner_mark, ps1);
  subdivide_volumes(lcc, corner_mark, ps1);

  diff=std::chrono::system_clock::now()-start2;
  print_txt_with_endl("    Subdivision: ", diff, "s.");

  diff=std::chrono::system_clock::now()-start;
  print_txt_with_endl("Create mesh TOTAL: ", diff, "s.");
  std::cout<<"Final map: "; display_stats(lcc);
  assert(lcc.is_valid());

  if(draw)
  { CGAL::draw(lcc); }

  if(save)
  {
    std::filesystem::path p(filename);
    save_object_3D(p.stem().string()+"-hexa.mesh", lcc);
  }

  lcc.free_mark(corner_mark);
  return true;
}
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  if (argc<2) // We need at least one filename
  { usage(argc, argv); }

  std::string filename;
  bool draw;
  bool save;

  process_command_line(argc, argv,
                       filename,
                       draw,
                       save);

  if (!create_mesh_from_off(filename,
                            draw,
                            save))
  { return EXIT_FAILURE; }

  return EXIT_SUCCESS;
}
////////////////////////////////////////////////////////////////////////////////
