// Copyright (c) 2011, 2015 GeometryFactory (France).
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
//
//
// Author(s)     : Sebastien Loriot and Andreas Fabri

// temporarily needed so that Mean_curvature_skeleton test works with/without PMP

#ifndef CGAL_PMP_POLYHEDRON_SUBSET_EXTRACTION_H
#define CGAL_PMP_POLYHEDRON_SUBSET_EXTRACTION_H

#include <CGAL/Polygon_mesh_processing/connected_components.h>

namespace CGAL {

  namespace internal {
    
    template <typename Polyhedron,typename Output_iterator>
    void extract_connected_components(const Polyhedron& P,Output_iterator out)
    {
      corefinement::extract_connected_components(P,out);
    }
  }
}


#endif // CGAL_PMP_POLYHEDRON_SUBSET_EXTRACTION_H
