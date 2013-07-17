// Copyright (c) 2013  GeometryFactory (France). All rights reserved.
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
// Author(s)     : Sebastien Loriot
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GEOMETRIC_TEST_SKIPPER_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GEOMETRIC_TEST_SKIPPER_H 1

namespace CGAL {

namespace Surface_mesh_simplification
{
  template <class GetPlacement>
  struct Geometric_test_skipper: public GetPlacement{
    struct Skip_geom_valid_test{};
  };
}

}

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GEOMETRIC_TEST_SKIPPER_H