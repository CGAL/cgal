// Copyright (c) 2023   GeometryFactory Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_TESTSUITE_DERIVED_SURFACE_MESH_H
#define CGAL_TESTSUITE_DERIVED_SURFACE_MESH_H

namespace CGAL { namespace Testsuite {

  template <typename P>
  struct DerivedSurfaceMesh: public CGAL::Surface_mesh<P> {
    typedef CGAL::Surface_mesh<P> Base;
    std::string name;
  };


  } // namespace Testsuite
} // namespace CGAL


#define CGAL_GRAPH_TRAITS_INHERITANCE_TEMPLATE_PARAMS typename P
#define CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME CGAL::Testsuite::DerivedSurfaceMesh<P>
#define CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME CGAL::Surface_mesh<P>
#include <CGAL/boost/graph/graph_traits_inheritance_macros.h>

#endif //CGAL_TESTSUITE_DERIVED_SURFACE_MESH_H
