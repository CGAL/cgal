// Copyright (C) 2014-2018  GeometryFactory Sarl
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//

#ifndef CGAL_SURFACE_MESH_FWD_H
#define CGAL_SURFACE_MESH_FWD_H

/// \file Surface_mesh_fwd.h
/// Forward declarations of the Surface_mesh package.

#ifndef DOXYGEN_RUNNING
namespace CGAL {

// fwdS for the public interface
template<typename P>
class Surface_mesh;

template<typename T>
struct is_derived_from_Surface_mesh : std::false_type {};

template<typename U>
struct is_derived_from_Surface_mesh<CGAL::Surface_mesh<U>> : std::true_type {};

} // CGAL
#endif

#endif /* CGAL_SURFACE_MESH_FWD_H */


