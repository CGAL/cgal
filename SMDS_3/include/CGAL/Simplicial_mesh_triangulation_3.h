// Copyright (c) 2006-2007  INRIA Sophia-Antipolis (France).
// Copyright (c) 2008,2011 GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau, Stephane Tayeb, Andreas Fabri, Jane Tournois

#ifndef CGAL_SIMPLICIAL_MESH_TRIANGULATION_3_H
#define CGAL_SIMPLICIAL_MESH_TRIANGULATION_3_H

#include <CGAL/license/SMDS_3.h>

#include <CGAL/Simplicial_mesh_vertex_base_3.h>
#include <CGAL/Simplicial_mesh_cell_base_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_3.h>

namespace CGAL
{
  /**
  *\ingroup PkgSMDS3Classes
  * `Simplicial_mesh_triangulation_3`
  * @todo
  */
  template<typename K,
           typename SubdomainIndex = int,
           typename SurfacePatchIndex = int,
           typename CurveIndex = int,
           typename CornerIndex = int>
  using Simplicial_mesh_triangulation_3 =
      CGAL::Triangulation_3<K,
        CGAL::Triangulation_data_structure_3<
          CGAL::Simplicial_mesh_vertex_base_3<K,
                                              SubdomainIndex,
                                              SurfacePatchIndex,
                                              CurveIndex,
                                              CornerIndex>,
          CGAL::Simplicial_mesh_cell_base_3<K, SubdomainIndex, SurfacePatchIndex>,
          CGAL::Sequential_tag
        >
      >;
};

#endif // CGAL_SIMPLICIAL_MESH_TRIANGULATION_3_H


