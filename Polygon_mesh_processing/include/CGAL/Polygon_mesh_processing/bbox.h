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
// 
//
// Author(s)     : Jane Tournois

#ifndef CGAL_POLYGON_MESH_PROCESSING_BOUNDING_BOX_H
#define CGAL_POLYGON_MESH_PROCESSING_BOUNDING_BOX_H

#include <CGAL/Bbox_3.h>

#include <boost/graph/graph_traits.hpp>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#include <boost/foreach.hpp>

namespace CGAL {

  namespace Polygon_mesh_processing {

    /*!
    * \ingroup PkgPolygonMeshProcessing
    *  computes a bounding box of a polygon mesh.
    *
    * @tparam PolygonMesh a model of `HalfedgeListGraph`
    *   that has an internal property map for `CGAL::vertex_point_t`
    * @tparam NamedParameters a sequence of \ref namedparameters
    *
    * @param pmesh a polygon mesh
    * @param np optional sequence of \ref namedparameters among the ones listed below
    *
    * \cgalNamedParamsBegin
    *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh` \cgalParamEnd
    * \cgalNamedParamsEnd
    *
    * @return a bounding box of `pmesh`
    */
    template<typename PolygonMesh, typename NamedParameters>
    CGAL::Bbox_3 bbox_3(const PolygonMesh& pmesh,
                        const NamedParameters& np)
    {
      using boost::choose_const_pmap;
      using boost::get_param;
      typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type
        vpm = choose_const_pmap(get_param(np, CGAL::vertex_point),
                                pmesh,
                                CGAL::vertex_point);

      typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

      halfedge_descriptor h0 = *(halfedges(pmesh).first);
      CGAL::Bbox_3 bb = get(vpm, target(h0, pmesh)).bbox();
      BOOST_FOREACH(halfedge_descriptor h, halfedges(pmesh))
      {
        bb += get(vpm, target(h, pmesh)).bbox();
      }
      return bb;
    }

    template<typename PolygonMesh>
    CGAL::Bbox_3 bbox_3(const PolygonMesh& pmesh)
    {
      return bbox_3(pmesh,
        CGAL::Polygon_mesh_processing::parameters::all_default());
    }


  }
}

#endif //CGAL_POLYGON_MESH_PROCESSING_BOUNDING_BOX_H

