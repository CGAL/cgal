// Copyright (c) 2015 GeometryFactory
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

#ifndef CGAL_MESH_3_NULL_SUBDOMAIN_INDEX
#define CGAL_MESH_3_NULL_SUBDOMAIN_INDEX

#include <CGAL/license/Mesh_3.h>

#ifndef DOXYGEN_RUNNING

namespace CGAL {
  struct Null_subdomain_index {
    template <typename T>
    bool operator()(const T& x) const { return 0 == x; }
  };
}

#endif

#endif //CGAL_MESH_3_NULL_SUBDOMAIN_INDEX
