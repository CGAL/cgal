// Copyright (C) 2020  GeometryFactory Sarl
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//

#ifndef CGAL_TMDS_3_MESH_COMPLEX_3_IN_TRIANGULATION_3_FWD_H
#define CGAL_TMDS_3_MESH_COMPLEX_3_IN_TRIANGULATION_3_FWD_H

/// \file Mesh_complex_3_in_triangulation_3_fwd.h
/// Forward declarations of the MDS_3 package.

#ifndef DOXYGEN_RUNNING
namespace CGAL {

// fwdS for the public interface
template <typename Tr,
          typename CornerIndex = int,
          typename CurveIndex = int>
class Mesh_complex_3_in_triangulation_3;

namespace IO {
  template <class C3T3>
  void output_to_medit(std::ostream& os,
                       const C3T3& c3t3,
                       bool rebind = false,
                       bool show_patches = false
#ifndef DOXYGEN_RUNNING
                     , bool all_vertices = true
                     , bool all_cells = false
#endif
  );
}

namespace TMDS_3 {

  template<class Tr>
  bool build_triangulation_from_file(std::istream& is,
                                     Tr& tr,
                                     bool verbose = false,
                                     bool replace_domain_0 = false);

} // MDS_3
} // CGAL
#endif

#endif /* CGAL_TMDS_3_MESH_COMPLEX_3_IN_TRIANGULATION_3_FWD_H */


