// Copyright (c) 2012  GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_MESH_3_DUMP_C3T3_H
#define CGAL_MESH_3_DUMP_C3T3_H

#include <string>
#include <CGAL/Mesh_3/io_signature.h>
#include <CGAL/Mesh_3/Mesh_complex_3_in_triangulation_3_base.h>
#include <CGAL/is_streamable.h>
#include <fstream>

namespace CGAL {

template <typename C3t3, 
          bool is_streamable = 
            is_streamable<typename C3t3::Triangulation::Vertex>::value &&
            is_streamable<typename C3t3::Triangulation::Cell>::value &&
            is_streamable<typename C3t3::Surface_patch_index>::value &&
            is_streamable<typename C3t3::Subdomain_index>::value 
          >
struct Dump_c3t3 {
  void dump_c3t3(const C3t3& c3t3, std::string prefix) const {
    std::clog<<"======dump c3t3===== to: " << prefix << std::endl;
    std::ofstream medit_file((prefix+".mesh").c_str());
    medit_file.precision(17);
    CGAL::output_to_medit(medit_file, c3t3, false, true);
    medit_file.close();

    std::string bin_filename = prefix;
    bin_filename += ".binary.cgal";
    std::ofstream bin_file(bin_filename.c_str(),
                           std::ios_base::binary | std::ios_base::out);
    bin_file << "binary CGAL c3t3 " << CGAL::Get_io_signature<C3t3>()() << "\n";
    CGAL::set_binary_mode(bin_file);
    bin_file << c3t3;
  }
}; // end struct template Dump_c3t3<C3t3, bool>

template <typename C3t3>
struct Dump_c3t3<C3t3, false> {
  void dump_c3t3(const C3t3&, std::string) {
    std::cerr << "Warning " << __FILE__ << ":" << __LINE__ << "\n"
              << "  the c3t3 object cannot be dumped because some types are"
              << " not streamable:\n";
    if(!is_streamable<typename C3t3::Triangulation::Vertex>::value)
      std::cerr << "     - C3t3::Triangulation::Vertex is not streamble\n";

    if(!is_streamable<typename C3t3::Triangulation::Cell>::value)
      std::cerr << "     - C3t3::Triangulation::Cell is not streamble\n";

    if(!is_streamable<typename C3t3::Surface_patch_index>::value)
      std::cerr << "     - C3t3::Surface_patch_index is not streamable\n";
      
    if(!is_streamable<typename C3t3::Subdomain_index>::value)
      std::cerr << "     - C3t3::Subdomain_index is not streamable\n";      
  }
}; // end struct template specialization Dump_c3t3<C3t3, false>

template <typename C3t3>
void dump_c3t3(const C3t3& c3t3, std::string prefix) {
  if(!prefix.empty()) {
    Dump_c3t3<C3t3> dump;
    dump.dump_c3t3(c3t3, prefix);
  }
}

} // end namespace CGAL

#endif // CGAL_MESH_3_DUMP_C3T3_H
