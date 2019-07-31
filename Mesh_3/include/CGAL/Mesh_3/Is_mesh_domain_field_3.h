// Copyright (c) 2019 GeometryFactory Sarl  (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Laurent Rineau
//
//******************************************************************************
// File Description :
// Mesh_facet_criteria_3 class.
//******************************************************************************

#ifndef CGAL_MESH_3_IS_MESH_DOMAIN_FIELD_3_H
#define CGAL_MESH_3_IS_MESH_DOMAIN_FIELD_3_H

#include <CGAL/license/Mesh_3.h>

#include <boost/config.hpp>
#if BOOST_VERSION >= 106600
#  include <boost/callable_traits/is_invocable.hpp>
#else
#  include <boost/mpl/has_xxx.hpp>
#endif

#include <CGAL/tags.h>

namespace CGAL {
  namespace Mesh_3 {
#if BOOST_VERSION >= 106600
    template <typename Tr, typename Type>
    struct Is_mesh_domain_field_3 :
      public CGAL::Boolean_tag
      <
        boost::callable_traits::is_invocable_r<
          typename Tr::FT,
          Type,
          typename Tr::Bare_point,
          int,
          typename Tr::Vertex::Index
        >::value
      >
    {};
#else // Boost before 1.66
    BOOST_MPL_HAS_XXX_TRAIT_DEF(FT)
    template <typename Tr, typename Type>
    struct Is_mesh_domain_field_3 : public Boolean_tag<has_FT<Type>::value>
    {};
#endif // Boost before 1.66
  } // end namespace Mesh_3
} // end namespace CGAL

#endif // CGAL_MESH_3_IS_MESH_DOMAIN_FIELD_3_H
