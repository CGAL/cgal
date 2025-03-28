// Copyright (c) 2019 GeometryFactory Sarl  (France).
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
//
//******************************************************************************
// File Description :
// Mesh_facet_criteria_3 class.
//******************************************************************************

#ifndef CGAL_MESH_3_IS_MESH_DOMAIN_FIELD_3_H
#define CGAL_MESH_3_IS_MESH_DOMAIN_FIELD_3_H

#include <CGAL/license/Mesh_3.h>

#include <boost/config.hpp>
#include <boost/callable_traits/is_invocable.hpp>

#include <CGAL/tags.h>
#include <CGAL/type_traits.h>

namespace CGAL {
  namespace Mesh_3 {
    template <typename Tr, typename Type>
    struct Is_mesh_domain_field_3 :
      public CGAL::Boolean_tag
      <
        boost::callable_traits::is_invocable_r<
          typename Tr::Geom_traits::FT,
          Type,
          Bare_point_type_t<Tr>,
          int,
          typename Tr::Vertex::Index
        >::value
      >
    {};
  } // end namespace Mesh_3
} // end namespace CGAL

#endif // CGAL_MESH_3_IS_MESH_DOMAIN_FIELD_3_H
