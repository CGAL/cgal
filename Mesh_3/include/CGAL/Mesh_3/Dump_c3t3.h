// Copyright (c) 2012  GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_MESH_3_DUMP_C3T3_H
#define CGAL_MESH_3_DUMP_C3T3_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Mesh_3/io_signature.h>
#include <CGAL/Mesh_3/Mesh_complex_3_in_triangulation_3_base.h>
#include <CGAL/IO/File_medit.h>

#include <CGAL/is_streamable.h>

#include <fstream>
#include <string>

namespace CGAL {

template <typename C3t3,
          bool is_streamable =
            is_streamable<typename C3t3::Triangulation::Vertex>::value &&
            is_streamable<typename C3t3::Triangulation::Cell>::value
            &&
            (is_streamable<typename C3t3::Surface_patch_index>::value ||
             Output_rep<typename C3t3::Surface_patch_index>::is_specialized)
            &&
            (is_streamable<typename C3t3::Subdomain_index>::value ||
             Output_rep<typename C3t3::Subdomain_index>::is_specialized)
          >
struct Dump_c3t3 {
  void dump_c3t3(const C3t3& c3t3, std::string prefix) const
  {
    std::clog<<"======dump c3t3===== to: " << prefix << std::endl;
    std::ofstream medit_file((prefix+".mesh").c_str());
    medit_file.precision(17);
    CGAL::output_to_medit(medit_file, c3t3, false /*rebind*/, true /*show_patches*/);
    medit_file.close();

    std::string bin_filename = prefix;
    bin_filename += ".binary.cgal";
    std::ofstream bin_file(bin_filename.c_str(),
                           std::ios_base::binary | std::ios_base::out);
    std::string signature = CGAL::Get_io_signature<C3t3>()();
    CGAL_assertion(signature != std::string());
    bin_file << "binary CGAL c3t3 " << signature << "\n";
    CGAL::set_binary_mode(bin_file);
    bin_file << c3t3;
  }
}; // end struct template Dump_c3t3<C3t3, bool>

template <typename C3t3>
struct Dump_c3t3<C3t3, false>
{
  void dump_c3t3(const C3t3&, std::string) {
    std::cerr << "Warning " << __FILE__ << ":" << __LINE__ << "\n"
              << "  the c3t3 object of following type:\n"
              << typeid(C3t3).name() << std::endl
              << "  cannot be dumped because some types are not streamable:\n";
    if(!is_streamable<typename C3t3::Triangulation::Vertex>::value) {
      std::cerr << "     - C3t3::Triangulation::Vertex is not streamble\n";
      std::cerr << "       "
                << typeid(typename C3t3::Triangulation::Vertex).name()
                << "\n";
    }

    if(!is_streamable<typename C3t3::Triangulation::Cell>::value) {
      std::cerr << "     - C3t3::Triangulation::Cell is not streamble\n";
      std::cerr << "       "
                << typeid(typename C3t3::Triangulation::Cell).name()
                << "\n";
    }

    if(!is_streamable<typename C3t3::Surface_patch_index>::value &&
       !CGAL::Output_rep<typename C3t3::Surface_patch_index>::is_specialized)
    {
      std::cerr << "     - C3t3::Surface_patch_index is not streamable\n";
      std::cerr << "       "
                << typeid(typename C3t3::Surface_patch_index).name()
                << "\n";
    }
    if(!is_streamable<typename C3t3::Subdomain_index>::value &&
       !CGAL::Output_rep<typename C3t3::Subdomain_index>::is_specialized)
    {
      std::cerr << "     - C3t3::Subdomain_index is not streamable\n";
      std::cerr << "       "
                << typeid(typename C3t3::Subdomain_index).name()
                << "\n";
    }
  }
}; // end struct template specialization Dump_c3t3<C3t3, false>

template <typename C3t3>
void dump_c3t3_edges(const C3t3& c3t3, std::string prefix)
{
  typename C3t3::Triangulation::Geom_traits::Construct_point_3 cp =
    c3t3.triangulation().geom_traits().construct_point_3_object();

  std::ofstream file((prefix+".polylines.txt").c_str());
  file.precision(17);
  for(typename C3t3::Edges_in_complex_iterator
        edge_it = c3t3.edges_in_complex_begin(),
        end     = c3t3.edges_in_complex_end();
      edge_it != end; ++edge_it)
  {
    const typename C3t3::Triangulation::Cell_handle c = edge_it->first;
    const int i = edge_it->second;
    const int j = edge_it->third;
    const typename C3t3::Triangulation::Weighted_point& ei = c3t3.triangulation().point(c, i);
    const typename C3t3::Triangulation::Weighted_point& ej = c3t3.triangulation().point(c, j);
    file << "2 " << cp(ei) << " "  << cp(ej) << "\n";
  }
}
template <typename C3t3>
void dump_c3t3(const C3t3& c3t3, std::string prefix)
{
  if(!prefix.empty()) {
    Dump_c3t3<C3t3> dump;
    dump.dump_c3t3(c3t3, prefix);
  }
}

} // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_MESH_3_DUMP_C3T3_H
