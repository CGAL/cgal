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

#ifndef CGAL_POLYGON_MESH_PROCESSING_GET_BORDER_H
#define CGAL_POLYGON_MESH_PROCESSING_GET_BORDER_H

#include <boost/graph/graph_traits.hpp>

#include <boost/foreach.hpp>
#include <set>

namespace CGAL{
namespace Polygon_mesh_processing {

  /**
  * collects the border of a face range
  * @param faces the range of face descriptors around which the
  *              border is computed
  * @param out the output iterator that collects edges that form the border
  *            of `faces`, seen from inside the surface patch
  *
  * @todo code : what shall we do for more than one connected components
  * @todo code : make sure the halfedges returned actually belong to a face
  * from `faces`. It's not the case yet
  */
  template<typename PolygonMesh
         , typename FaceRange
         , typename HalfedgeOutputIterator>
  void get_border(const PolygonMesh& pmesh
                , const FaceRange& faces
                , HalfedgeOutputIterator out)
  {
    typedef PolygonMesh PM;
    typedef typename boost::graph_traits<PM>::edge_descriptor     edge_descriptor;
    typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;
    
    typedef typename boost::graph_traits<PM>::face_descriptor     face_descriptor;

    //collect halfedges that appear only once
    std::set<halfedge_descriptor> border;
    BOOST_FOREACH(face_descriptor f, faces)
    {
      BOOST_FOREACH(halfedge_descriptor h,
                    halfedges_around_face(halfedge(f, pmesh), pmesh))
      {
        //halfedge_descriptor is model of `LessThanComparable`
        halfedge_descriptor he = (h < opposite(h, pmesh))
          ? h
          : opposite(h, pmesh);
        if (border.find(he) != border.end())
          border.erase(he); //even number of appearances
        else
          border.insert(he);//odd number of appearances
      }
    }
    //copy them in out
    BOOST_FOREACH(halfedge_descriptor h, border)
    {
      *out++ = h;
    }
  }

} } // end of namespace CGAL::Polygon_mesh_processing


#endif //CGAL_POLYGON_MESH_PROCESSING_GET_BORDER_H
