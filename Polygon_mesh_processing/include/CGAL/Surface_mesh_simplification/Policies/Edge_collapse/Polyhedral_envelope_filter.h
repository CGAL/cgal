// Copyright (c) 2020  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_POLYHEDRAL_ENVELOPE_FILTER_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_POLYHEDRAL_ENVELOPE_FILTER_H

#include <CGAL/license/Polygon_mesh_processing/Polyhedral_envelope.h>

#include <CGAL/assertions.h>
#include <CGAL/Default.h>
#include <CGAL/intersections.h>
#include <CGAL/boost/graph/named_params_helper.h>


#include <CGAL/Polyhedral_envelope.h>

#include <boost/optional.hpp>

#include <vector>
#include <type_traits>

namespace CGAL {
namespace Surface_mesh_simplification {

namespace internal {

  struct Dummy_filter2 {
  template <typename Profile>
  inline
  const boost::optional<typename Profile::Point>
  operator()(const Profile&, const boost::optional<typename Profile::Point>& op) const
  {
    return op;
  }

};

} // namesapce internal

template<typename GeomTraits,typename BaseFilter = internal::Dummy_filter2>
class Polyhedral_envelope_filter
{
  typedef GeomTraits                                                          Geom_traits;
  typedef typename Geom_traits::FT                                            FT;
  typedef typename Geom_traits::Point_3                                       Point_3;

  typedef CGAL::Polyhedral_envelope<GeomTraits>                               Envelope;
  typedef typename Envelope::Vector3i                                         Vector3i;

private:
  template <typename Profile>
  void initialize_envelope(const Profile& profile) const
  {
    CGAL_static_assertion((std::is_same<GeomTraits, typename Profile::Geom_traits>::value));

    typedef typename Profile::Triangle_mesh                                   Triangle_mesh;
    typedef typename boost::graph_traits<Triangle_mesh>::halfedge_descriptor  halfedge_descriptor;
    typedef typename boost::graph_traits<Triangle_mesh>::face_descriptor      face_descriptor;
    typedef typename boost::graph_traits<Triangle_mesh>::vertex_descriptor    vertex_descriptor;

    const Triangle_mesh& tm = profile.surface_mesh();

    m_vertices.reserve(num_vertices(tm));
    m_faces.reserve(num_faces(tm));

    for(vertex_descriptor v : vertices(tm)){
      m_vertices.emplace_back(get(profile.vertex_point_map(),v));
    }

    auto vim = get(vertex_index, tm);
    for(face_descriptor f : faces(tm))
    {
      halfedge_descriptor h = halfedge(f, tm);
      CGAL_assertion(!is_border(h, tm));
      int i = get(vim, source(h, tm));
      int j = get(vim, target(h, tm));
      int k = get(vim, target(next(h, tm), tm));

      Vector3i face = { i, j, k };
      m_faces.push_back(face);
    }

    m_envelope = new Envelope(m_vertices, m_faces, m_dist);
  }


public:
  Polyhedral_envelope_filter(const FT dist,
                             const BaseFilter& filter = BaseFilter())
    :
      m_dist(dist),
      m_envelope(nullptr),
      m_base_filter(filter)
  {}

  ~Polyhedral_envelope_filter()
  {
    if(m_envelope != nullptr){
      delete m_envelope;
    }
  }


  template <typename Profile>
  boost::optional<typename Profile::Point>
  operator()(const Profile& profile, boost::optional<typename Profile::Point> op) const
  {
    typedef typename Profile::Point Point;
    typedef typename Profile::vertex_descriptor_vector Link;
    typedef typename Profile::Triangle_mesh Triangle_mesh;
    typedef typename boost::graph_traits<Triangle_mesh>::vertex_descriptor vertex_descriptor;

    op = m_base_filter(profile, op);
    if(op)
    {
      if(m_envelope == nullptr){
        initialize_envelope(profile);
      }
      const Point& p = *op;

      if(! (*m_envelope)(p)){
        // the new placement is outside envelope
        return boost::none;
      }

      const Link link = profile.link();

      vertex_descriptor v = link.back();
      Point pv = get(profile.vertex_point_map(),v);

      for(std::size_t i=0; i!= link.size();i++){
        vertex_descriptor w = link[i];
        Point pw = get(profile.vertex_point_map(),w);

        if(! (*m_envelope)(p, pv, pw)){
          // the triange intersects the envelope
          return boost::none;
        }
        pv = pw;

      }
    }
    return op;
  }


private:
  const FT m_dist;
  mutable Envelope* m_envelope;
  mutable std::vector<Point_3> m_vertices;
  mutable std::vector<Vector3i> m_faces;

  const BaseFilter m_base_filter;
};



} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_POLYHEDRAL_ENVELOPE_FILTER_H
