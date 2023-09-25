// Copyright (c) 2023 University of Genova (Italy).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Claudio Mancinelli

#ifndef CGAL_POLYGON_MESH_PROCESSING_BSURF_LOCALLY_SHORTEST_PATH_H
#define CGAL_POLYGON_MESH_PROCESSING_BSURF_LOCALLY_SHORTEST_PATH_H

// #include <CGAL/license/Polygon_mesh_processing/bsurf.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/graph/graph_traits.hpp>

#include <CGAL/Polygon_mesh_processing/locate.h>

namespace CGAL{
namespace Polygon_mesh_processing {

template <class TriangleMesh, class FT>
struct Edge_location
{
  typename boost::graph_traits<TriangleMesh>::edge_descriptor ed;
  FT alpha;
};

template <class FT, class TriangleMesh, class EdgeLocationRange>
void
locally_shortest_path(const Face_location<TriangleMesh, FT>& src,
                      const Face_location<TriangleMesh, FT>& tgt,
                      const TriangleMesh& tmesh,
                      EdgeLocationRange& edge_locations)
{

}


} } // CGAL::Polygon_mesh_processing

#endif
