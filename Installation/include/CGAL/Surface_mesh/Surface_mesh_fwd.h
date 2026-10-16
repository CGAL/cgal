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

#include <boost/mpl/has_xxx.hpp>

#ifndef DOXYGEN_RUNNING

namespace CGAL {
BOOST_MPL_HAS_XXX_TRAIT_DEF(Point)

// fwdS for the public interface
template<typename P>
class Surface_mesh;

template<typename T, typename = void>
struct is_surface_mesh : std::false_type {};

template <typename T>
struct is_surface_mesh<T, std::enable_if_t<CGAL::has_Point<T>::value>>: std::is_base_of< CGAL::Surface_mesh<typename T::Point>, T>
{};


} // CGAL
#endif

#endif /* CGAL_SURFACE_MESH_FWD_H */


