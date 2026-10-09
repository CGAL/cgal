// Copyright (C) 2020  GeometryFactory Sarl
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//

#ifndef CGAL_MESH_3_MESH_COMPLEX_3_IN_TRIANGULATION_3_FWD_H
#define CGAL_MESH_3_MESH_COMPLEX_3_IN_TRIANGULATION_3_FWD_H

/// \file Mesh_complex_3_in_triangulation_3_fwd.h
/// Forward declarations of the Mesh_3 package.

#include <istream>

#ifndef DOXYGEN_RUNNING
namespace CGAL {

// fwdS for the public interface
template <typename Tr,
          typename CornerIndex = int,
          typename CurveIndex = int>
class Mesh_complex_3_in_triangulation_3;

template <class C3T3, bool c3t3_loader_failed>
bool build_mesh_complex_from_file(std::istream& is,
                                  C3T3& c3t3,
                                  bool replace_domain_0);

template <class C3T3, bool c3t3_loader_failed>
bool build_mesh_complex_from_file(std::istream& is,
                                  C3T3& c3t3);
} // CGAL
#endif

#endif /* CGAL_MESH_3_MESH_COMPLEX_3_IN_TRIANGULATION_3_FWD_H */


