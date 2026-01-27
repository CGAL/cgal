// Copyright (c) 2025 Geometry Factory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Léo Valque
//

#ifndef CGAL_CONVEX_HULL_INTERNAL_HELPERS_H
#define CGAL_CONVEX_HULL_INTERNAL_HELPERS_H

#include <CGAL/license/Convex_hull_3.h>

#ifdef CGAL_PROFILE_CONVEX_HULL_DO_INTERSECT
      size_t nb_visited=0;
#endif

namespace CGAL::Convex_hull_3::internal{

template <class T, template <class...> class Template>
inline constexpr bool is_instance_of_v = false;

template <template <class...> class Template, class... Args>
inline constexpr bool is_instance_of_v<Template<Args...>, Template> = true;

// template class to deduce the GT from Convex and NamedParameters
template<class Convex,
         class NamedParameters,
         bool Is_range=CGAL::IO::internal::is_Range_v<Convex> >
struct GetGeomTraitsFromConvex{
  typedef typename GetGeomTraits<Convex, NamedParameters>::type type;
};

template<class Convex, class NamedParameters>
struct GetGeomTraitsFromConvex<Convex, NamedParameters, true>{
  typedef typename Point_set_processing_3_np_helper<Convex, NamedParameters>::Geom_traits type;
};

template<class Mesh, class NamedParameters>
struct GetGeomTraitsFromConvex<Convex_hull_hierarchy<Mesh>, NamedParameters, false>{
  typedef typename GetGeomTraits<Mesh, NamedParameters>::type type;
};

}

#endif // CGAL_CONVEX_HULL_INTERNAL_HELPERS_H