// Copyright (c) 2020 GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois

#ifndef CGAL_PMP_REMESHING_ADAPTIVE_SIZING_FIELD_H
#define CGAL_PMP_REMESHING_ADAPTIVE_SIZING_FIELD_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>

#include <CGAL/Polygon_mesh_processing/internal/Isotropic_remeshing/Sizing_field.h>

#include <CGAL/Polygon_mesh_processing/interpolated_corrected_curvatures.h>

#include <CGAL/number_utils.h>

namespace CGAL
{
namespace Polygon_mesh_processing
{
template <class PolygonMesh>
class Adaptive_sizing_field : public CGAL::Sizing_field<PolygonMesh>
{
private:
  typedef CGAL::Sizing_field<PolygonMesh> Base;

public:
  typedef typename Base::K          K;
  typedef typename Base::FT         FT;
  typedef typename Base::Point_3    Point_3;
  typedef typename Base::face_descriptor     face_descriptor;
  typedef typename Base::halfedge_descriptor halfedge_descriptor;
  typedef typename Base::vertex_descriptor   vertex_descriptor;

  typedef typename CGAL::dynamic_vertex_property_t<FT> Vertex_property_tag;
  typedef typename boost::property_map<PolygonMesh,
                                       Vertex_property_tag>::type Vertex_sizing_map;

  template <typename FaceRange>
  Adaptive_sizing_field(const double tol
                      , const std::pair<FT, FT>& edge_len_min_max
                      , const FaceRange& face_range
                      , PolygonMesh& pmesh)
    : tol(tol)
    , m_short(edge_len_min_max.first)
    , m_long(edge_len_min_max.second)
    , m_pmesh(pmesh)
  {
    m_vertex_sizing_map = get(Vertex_property_tag(), m_pmesh);

    if (face_range.size() == faces(m_pmesh).size())
    {
      // calculate curvature from the whole mesh
      calc_sizing_map(m_pmesh);
    }
    else
    {
      // expand face selection and calculate curvature from it
      std::vector<face_descriptor> selection(face_range.begin(), face_range.end());
      auto is_selected = get(CGAL::dynamic_face_property_t<bool>(), m_pmesh);
      for (face_descriptor f : faces(m_pmesh)) put(is_selected, f, false);
      for (face_descriptor f : face_range)  put(is_selected, f, true);
      CGAL::expand_face_selection(selection, m_pmesh, 1
                                , is_selected, std::back_inserter(selection));
      CGAL::Face_filtered_graph<PolygonMesh> ffg(m_pmesh, selection);

      calc_sizing_map(ffg);
    }
  }

private:
  template <typename FaceGraph>
  void calc_sizing_map(FaceGraph& face_graph)
  {
    //todo ip: please check if this is good enough to store curvature
    typedef Principal_curvatures_and_directions<K> Principal_curvatures;
    typedef typename CGAL::dynamic_vertex_property_t<Principal_curvatures> Vertex_curvature_tag;
    typedef typename boost::property_map<FaceGraph,
      Vertex_curvature_tag>::type Vertex_curvature_map;

#ifdef CGAL_PMP_REMESHING_VERBOSE
    int oversize  = 0;
    int undersize = 0;
    int insize    = 0;
    std::cout << "Calculating sizing field..." << std::endl;
#endif

    Vertex_curvature_map vertex_curvature_map = get(Vertex_curvature_tag(), face_graph);
    interpolated_corrected_principal_curvatures_and_directions(face_graph
                                                             , vertex_curvature_map);
    // calculate vertex sizing field L(x_i) from the curvature field
    for(vertex_descriptor v : vertices(face_graph))
    {
      auto vertex_curv = get(vertex_curvature_map, v);
      const FT max_absolute_curv = CGAL::max(CGAL::abs(vertex_curv.max_curvature)
                                           , CGAL::abs(vertex_curv.min_curvature));
      const FT vertex_size_sq = 6 * tol / max_absolute_curv - 3 * CGAL::square(tol);
      if (vertex_size_sq > CGAL::square(m_long))
      {
        put(m_vertex_sizing_map, v, m_long);
#ifdef CGAL_PMP_REMESHING_VERBOSE
        ++oversize;
#endif
      }
      else if (vertex_size_sq < CGAL::square(m_short))
      {
        put(m_vertex_sizing_map, v, m_short);
#ifdef CGAL_PMP_REMESHING_VERBOSE
        ++undersize;
#endif
      }
      else
      {
        put(m_vertex_sizing_map, v, CGAL::approximate_sqrt(vertex_size_sq));
#ifdef CGAL_PMP_REMESHING_VERBOSE
        ++insize;
#endif
      }
    }
#ifdef CGAL_PMP_REMESHING_VERBOSE
    std::cout << " done (" << insize << " from curvature, "
              << oversize  << " set to max, "
              << undersize << " set to min)" << std::endl;
#endif
  }

  FT sqlength(const vertex_descriptor va,
              const vertex_descriptor vb) const
  {
    typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::const_type
      vpmap = get(CGAL::vertex_point, m_pmesh);
    return FT(CGAL::squared_distance(get(vpmap, va), get(vpmap, vb)));
  }

  FT sqlength(const halfedge_descriptor& h) const
  {
    return sqlength(target(h, m_pmesh), source(h, m_pmesh));
  }

public:
  FT get_sizing(const vertex_descriptor v) const {
      CGAL_assertion(get(m_vertex_sizing_map, v));
      return get(m_vertex_sizing_map, v);
    }

  boost::optional<FT> is_too_long(const halfedge_descriptor h) const
  {
    const FT sqlen = sqlength(h);
    FT sqtarg_len = CGAL::square(4./3. * CGAL::min(get(m_vertex_sizing_map, source(h, m_pmesh)),
                                                   get(m_vertex_sizing_map, target(h, m_pmesh))));
    CGAL_assertion(get(m_vertex_sizing_map, source(h, m_pmesh)));
    CGAL_assertion(get(m_vertex_sizing_map, target(h, m_pmesh)));
    if(sqlen > sqtarg_len)
      return sqlen;
    else
      return boost::none;
  }

  boost::optional<FT> is_too_long(const vertex_descriptor va,
                                  const vertex_descriptor vb) const
  {
    const FT sqlen = sqlength(va, vb);
    FT sqtarg_len = CGAL::square(4./3. * CGAL::min(get(m_vertex_sizing_map, va),
                                                   get(m_vertex_sizing_map, vb)));
    CGAL_assertion(get(m_vertex_sizing_map, va));
    CGAL_assertion(get(m_vertex_sizing_map, vb));
    if (sqlen > sqtarg_len)
      return sqlen;
    else
      return boost::none;
  }

  boost::optional<FT> is_too_short(const halfedge_descriptor h) const
  {
    const FT sqlen = sqlength(h);
    FT sqtarg_len = CGAL::square(4./5. * CGAL::min(get(m_vertex_sizing_map, source(h, m_pmesh)),
                                                   get(m_vertex_sizing_map, target(h, m_pmesh))));
    CGAL_assertion(get(m_vertex_sizing_map, source(h, m_pmesh)));
    CGAL_assertion(get(m_vertex_sizing_map, target(h, m_pmesh)));
    if (sqlen < sqtarg_len)
      return sqlen;
    else
      return boost::none;
  }

  virtual Point_3 split_placement(const halfedge_descriptor h) const
  {
    typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::const_type
      vpmap = get(CGAL::vertex_point, m_pmesh);
    return CGAL::midpoint(get(vpmap, target(h, m_pmesh)),
                          get(vpmap, source(h, m_pmesh)));
  }

  void update_sizing_map(const vertex_descriptor v)
  {
    // calculating it as the average of two vertices on other ends
    // of halfedges as updating is done during an edge split
    FT vertex_size = 0;
    CGAL_assertion(CGAL::halfedges_around_target(v, m_pmesh).size() == 2);
    for (halfedge_descriptor ha: CGAL::halfedges_around_target(v, m_pmesh))
    {
      vertex_size += get(m_vertex_sizing_map, source(ha, m_pmesh));
    }
    vertex_size /= CGAL::halfedges_around_target(v, m_pmesh).size();

    put(m_vertex_sizing_map, v, vertex_size);
  }

  //todo ip: is_protected_constraint_too_long() from PR

private:
  const FT tol;
  const FT m_short;
  const FT m_long;
  PolygonMesh& m_pmesh;
  Vertex_sizing_map m_vertex_sizing_map;
};

}//end namespace Polygon_mesh_processing
}//end namespace CGAL

#endif //CGAL_PMP_REMESHING_ADAPTIVE_SIZING_FIELD_H
