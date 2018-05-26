// Copyright (c) 2018  Carnegie Mellon University (USA), GeometryFactory (France)
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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s) : Christina Vaz, Keenan Crane, Andreas Fabri

#ifndef CGAL_HEAT_METHOD_3_HEAT_METHOD_3_H
#define CGAL_HEAT_METHOD_3_HEAT_METHOD_3_H

#include <CGAL/license/Heat_method_3.h>

#include <CGAL/disable_warnings.h>
#include <set>

#include <Eigen/Cholesky>

#include <Eigen/Dense>

#include <boost/foreach.hpp>

namespace CGAL {
namespace Heat_method_3 {

  // This class will later go into another file
  // It encapsulates what we use from Eigen so that one potentially can use another LA library
  struct Heat_method_Eigen_traits_3 {
    typedef Eigen::Matrix3d Matrix;
  };


  /**
   * Class `Heat_method_3` is a ...
   * \tparam TriangleMesh a triangulated surface mesh, model of `FaceGraph` and `HalfedgeListGraph`
   * \tparam Traits a model of HeatMethodTraits_3
   * \tparam VertexPointMap a model of `ReadablePropertyMap` with
   *        `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key and
   *        `Traits::Point_3` as value type.
   *        The default is `typename boost::property_map< TriangleMesh, vertex_point_t>::%type`.
   *
   */
  template <typename TriangleMesh,
            typename Traits,
            typename VertexPointMap = typename boost::property_map< TriangleMesh, vertex_point_t>::type,
            typename LA = Heat_method_Eigen_traits_3>
  class Heat_method_3
  {
    /// Polygon_mesh typedefs
    typedef typename boost::graph_traits<TriangleMesh>               graph_traits;
    typedef typename graph_traits::vertex_descriptor            vertex_descriptor;
    typedef typename graph_traits::edge_descriptor                edge_descriptor;
    typedef typename graph_traits::halfedge_descriptor        halfedge_descriptor;
    typedef typename graph_traits::face_descriptor                face_descriptor;

    /// Geometric typedefs
    typedef typename Traits::Point_3                                      Point_3;
    typedef typename Traits::FT                                                FT;

    typedef typename LA::Matrix Matrix;
  public:

    Heat_method_3(const TriangleMesh& tm)
      : tm(tm), vpm(get(vertex_point,tm))
    {
      build();
    }

    Heat_method_3(const TriangleMesh& tm, VertexPointMap vpm)
      : tm(tm), vpm(vpm)
    {
      build();
    }

    /**
     * add `vd` to the source set, returning `false` if `vd` is already in the set.
     */
    bool add_source(vertex_descriptor vd)
    {
      return sources.insert(vd).second;
    }

    /**
     * remove 'vd' from the source set, returning 'true' if 'vd' was in the set
     */
    bool remove_source(vertex_descriptor vd)
    {
      if(sources.find(vd))
      {
        sources.erase(vd);
        return true;
      }
      else
      {
        return false;
      }
    }

    /**
     *return current source set
    */
    vertex_descriptor getSources()
    {
      return sources;
    }

    /**
     * clear the current source set
     */
    void clearSources()
    {
      sources.clear();
      return;
    }

    /**
     * return vertex_descriptor to first vertex in the source set
     */
    vertex_descriptor sourcesBegin()
    {
      return sources.begin();
    }
    /**
     * return vertex_descriptor to last vertex in the source set
     */
    vertex_descriptor sourcesEnd()
    {
      return sources.end();
    }

    /**
     * get distance from the current source set to a vertex ` vd`.
     */
    double distance(vertex_descriptor vd)
    {
      return 0;
    }

  private:

    void build()
    {
      BOOST_FOREACH(vertex_descriptor vd, vertices(tm)){
        BOOST_FOREACH(vertex_descriptor one_ring_vd, vertices_around_target(halfedge(vd,tm),tm)){
          Point_3 p = get(vpm,one_ring_vd);
          std::cout << p << std::endl;
        }
      }
    }

    const TriangleMesh& tm;
    VertexPointMap vpm;
    std::set<vertex_descriptor> sources;

    Matrix m;
  };

} // namespace Heat_method_3
} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif CGAL_HEAT_METHOD_3_HEAT_METHOD_3_H
