// Copyright (c) 2015 GeometryFactory (France).
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
// Author(s)     : Jane Tournois
//
//******************************************************************************
// File Description :
//
//
//******************************************************************************

#ifndef CGAL_INTERNAL_MESH_3_CHECK_WEIGHTS_H
#define CGAL_INTERNAL_MESH_3_CHECK_WEIGHTS_H

#include <CGAL/license/Mesh_3.h>


#include <CGAL/tags.h>
#include <CGAL/Mesh_3/Has_features.h>
#include <boost/type_traits/is_same.hpp>

namespace CGAL {
namespace internal {
namespace Mesh_3 {


template<typename Triangulation, typename MeshDomain>
bool has_non_protecting_weights(const Triangulation& tr,
                                const MeshDomain&)
{
  bool with_features =
    boost::is_same<Has_features<MeshDomain>, CGAL::Tag_true>::value;

  for (typename Triangulation::Finite_vertices_iterator
        vv = tr.finite_vertices_begin();
        vv != tr.finite_vertices_end();
        ++vv)
  {
    if (vv->point().weight() != 0.)
    {
      if (with_features && vv->in_dimension() > 1)
        return true;
    }
  }
  return false;
}


}//end namespace Mesh_3
}//end namespace internal
}//end namespace CGAL

#endif //CGAL_INTERNAL_MESH_3_CHECK_WEIGHTS_H

