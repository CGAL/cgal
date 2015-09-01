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

#ifndef CGAL_POLYGON_MESH_PROCESSING_DIMENSIONS_H
#define CGAL_POLYGON_MESH_PROCESSING_DIMENSIONS_H

#include <CGAL/boost/graph/iterator.h>
#include <CGAL/property_map.h>

namespace CGAL {

namespace Polygon_mesh_processing {

  template<typename PolygonMesh>
  double
  border_length(PolygonMesh& pmesh,
                typename boost::graph_traits<PolygonMesh>::halfedge_descriptor hd)
  {
    double result = 0;
    typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
    typedef typename boost::property_map<PolygonMesh, boost::vertex_point_t>::const_type VPM;
    VPM vpm = get (CGAL::vertex_point, pmesh);
    BOOST_FOREACH(halfedge_descriptor haf, halfedges_around_face(hd,pmesh)){
      result += sqrt(CGAL::squared_distance(get(vpm, source(haf,pmesh)),
                                            get(vpm, target(haf,pmesh))));
    }
    return result;
  }



  template<typename PolygonMesh, typename FaceRange>
  double
  area(PolygonMesh& pmesh,
       FaceRange fr)
  {
    double result = 0;
    typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;

    typedef typename boost::property_map<PolygonMesh, boost::vertex_point_t>::const_type VPM;
    VPM vpm = get (CGAL::vertex_point, pmesh);
    
    BOOST_FOREACH(face_descriptor f, fr){
      halfedge_descriptor hd = halfedge(f, pmesh);
      halfedge_descriptor nhd = next(hd,pmesh);
      
      result += sqrt(CGAL::squared_area(get(vpm, source(hd,pmesh)),
                                        get(vpm, target(hd,pmesh)),
                                        get(vpm, target(nhd,pmesh))));
    }
    return result;
  }

}
}

#endif // CGAL_POLYGON_MESH_PROCESSING_DIMENSIONS_H
