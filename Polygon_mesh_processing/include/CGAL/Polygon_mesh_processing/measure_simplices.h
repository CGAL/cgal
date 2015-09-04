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
// Author(s)     : Andreas Fabri

#ifndef CGAL_POLYGON_MESH_PROCESSING_MEASURE_SIMPLICES_H
#define CGAL_POLYGON_MESH_PROCESSING_MEASURE_SIMPLICES_H

#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/properties.h>
#include <boost/graph/graph_traits.hpp>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

namespace CGAL {

namespace Polygon_mesh_processing {

  template<typename PolygonMesh,
           typename NamedParameters>
  double length(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h
              , const PolygonMesh& pmesh
              , const NamedParameters& np)
  {
    using boost::choose_const_pmap;
    using boost::get_param;

    typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type
    vpm = choose_const_pmap(get_param(np, CGAL::vertex_point),
                            pmesh,
                            CGAL::vertex_point);

    return CGAL::sqrt(CGAL::squared_distance(get(vpm, source(h, pmesh)),
                                             get(vpm, target(h, pmesh))));
  }

  template<typename PolygonMesh>
  double length(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h
    , const PolygonMesh& pmesh)
  {
    return length(h, pmesh,
      CGAL::Polygon_mesh_processing::parameters::all_default());
  }

  template<typename PolygonMesh,
           typename NamedParameters>
  double
  border_length(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h
              , const PolygonMesh& pmesh
              , const NamedParameters& np)
  {
    double result = 0.;
    BOOST_FOREACH(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor haf,
                  halfedges_around_face(h, pmesh))
    {
      result += length(haf, pmesh, np);
    }
    return result;
  }

  template<typename PolygonMesh>
  double
  border_length(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h
              , const PolygonMesh& pmesh)
  {
    return border_length(h, pmesh,
      CGAL::Polygon_mesh_processing::parameters::all_default());
  }


  template<typename PolygonMesh,
           typename NamedParameters>
  double area(typename boost::graph_traits<PolygonMesh>::face_descriptor f
            , const PolygonMesh& pmesh
            , const NamedParameters& np)
  {
    using boost::choose_const_pmap;
    using boost::get_param;

    typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type
    vpm = choose_const_pmap(get_param(np, CGAL::vertex_point),
                            pmesh,
                            CGAL::vertex_point);

    typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
    halfedge_descriptor hd = halfedge(f, pmesh);
    halfedge_descriptor nhd = next(hd, pmesh);

    return CGAL::sqrt(CGAL::squared_area(get(vpm, source(hd, pmesh)),
                                         get(vpm, target(hd, pmesh)),
                                         get(vpm, target(nhd, pmesh))));
  }

  template<typename PolygonMesh>
  double area(typename boost::graph_traits<PolygonMesh>::face_descriptor f
            , const PolygonMesh& pmesh)
  {
    return area(f, pmesh,
      CGAL::Polygon_mesh_processing::parameters::all_default());
  }

  template<typename FaceRange,
           typename PolygonMesh,
           typename NamedParameters>
  double area(FaceRange face_range
            , const PolygonMesh& pmesh
            , const NamedParameters& np)
  {
    typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
    double result = 0.;
    BOOST_FOREACH(face_descriptor f, face_range)
    {
      result += area(f, pmesh, np);
    }
    return result;
  }

  template<typename PolygonMesh, typename FaceRange>
  double area(FaceRange face_range, const PolygonMesh& pmesh)
  {
    return area(face_range, pmesh,
      CGAL::Polygon_mesh_processing::parameters::all_default());
  }

}
}

#endif // CGAL_POLYGON_MESH_PROCESSING_MEASURE_SIMPLICES_H
