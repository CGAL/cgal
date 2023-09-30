// Copyright (c) 2023 GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Ivan Paden

#ifndef CGAL_PMP_REMESHING_ADAPTIVE_SIZING_FIELD_H
#define CGAL_PMP_REMESHING_ADAPTIVE_SIZING_FIELD_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>

#include <CGAL/Polygon_mesh_processing/internal/Sizing_field_base.h>
#include <CGAL/Polygon_mesh_processing/interpolated_corrected_curvatures.h>

#include <CGAL/boost/graph/selection.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>

#include <CGAL/number_utils.h>

namespace CGAL
{
namespace Polygon_mesh_processing
{
/*!
* \ingroup PMP_meshing_grp
* provides criteria for isotropic remeshing to achieve variable mesh edge lengths as a function
* of local discrete curvatures.
*
* The local discrete curvatures are calculated using the
* `CGAL::Polygon_mesh_processing::interpolated_corrected_principal_curvatures_and_directions()` function.
*
* Edges longer than the local target edge length are split in two, while
* edges shorter than the local target edge length are collapsed.
*
* \cgalModels{PMPSizingField}
*
* \sa `isotropic_remeshing()`
* \sa `Uniform_sizing_field`
*
* @tparam PolygonMesh model of `MutableFaceGraph` that
*         has an internal property map for `CGAL::vertex_point_t`.
*/
template <class PolygonMesh,
          class VPMap = typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::const_type>
class Adaptive_sizing_field
#ifndef DOXYGEN_RUNNING
: public internal::Sizing_field_base<PolygonMesh, VPMap>
#endif
{
private:
  typedef internal::Sizing_field_base<PolygonMesh, VPMap> Base;
  typedef typename CGAL::dynamic_vertex_property_t<typename Base::FT> Vertex_property_tag;
  typedef typename boost::property_map<PolygonMesh,
                                       Vertex_property_tag>::type VertexSizingMap;

public:
  typedef typename Base::K          K;
  typedef typename Base::FT         FT;
  typedef typename Base::Point_3    Point_3;
  typedef typename Base::face_descriptor     face_descriptor;
  typedef typename Base::halfedge_descriptor halfedge_descriptor;
  typedef typename Base::vertex_descriptor   vertex_descriptor;

  /// \name Creation
  /// @{
  /*!
  * Constructor
  *
  * @tparam FaceRange range of `boost::graph_traits<PolygonMesh>::%face_descriptor`,
  *         model of `Range`. Its iterator type is `ForwardIterator`.
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param tol the error tolerance, used together with curvature to derive target edge length.
  *        Lower tolerance values will result in shorter mesh edges.
  * @param edge_len_min_max contains the bounds for minimum and maximum
  *        edge lengths
  * @param face_range the range of triangular faces defining one or several surface patches
  *        to be remeshed. It should be the same as the range of faces passed to `isotropic_remeshing()`.
  * @param pmesh a polygon mesh with triangulated surface patches to be remeshed. It should be the
  *              same mesh as the one passed to `isotropic_remeshing()`.
  * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
  *     \cgalParamType{a class model of `ReadWritePropertyMap` with
  *                   `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
  *                    as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
  *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
  *                     must be available in `PolygonMesh`.}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{ball_radius}
  *     \cgalParamDescription{`ball_radius` parameter passed to `interpolated_corrected_curvatures()`}
  *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  */
  template <typename FaceRange
          , typename NamedParameters = parameters::Default_named_parameters>
  Adaptive_sizing_field(const FT tol
                      , const std::pair<FT, FT>& edge_len_min_max
                      , const FaceRange& face_range
                      , PolygonMesh& pmesh
                      , const NamedParameters& np = parameters::default_values())
    : tol(tol)
    , m_short(edge_len_min_max.first)
    , m_long(edge_len_min_max.second)
    , m_vpmap(parameters::choose_parameter(
                parameters::get_parameter(np, internal_np::vertex_point),
                get_property_map(vertex_point, pmesh)))
    , m_vertex_sizing_map(get(Vertex_property_tag(), pmesh))
  {
    if (face_range.size() == faces(pmesh).size())
    {
      // calculate curvature from the whole mesh
      calc_sizing_map(pmesh, np);
    }
    else
    {
      // expand face selection and calculate curvature from it
      std::vector<face_descriptor> selection(face_range.begin(), face_range.end());
      auto is_selected = get(CGAL::dynamic_face_property_t<bool>(), pmesh);
      for (face_descriptor f : faces(pmesh)) put(is_selected, f, false);
      for (face_descriptor f : face_range)  put(is_selected, f, true);
      expand_face_selection(selection, pmesh, 1,
                            is_selected, std::back_inserter(selection));
      Face_filtered_graph<PolygonMesh> ffg(pmesh, selection);

      calc_sizing_map(ffg, np);
    }
  }

  ///@}

private:
  template <typename FaceGraph
          , typename NamedParameters = parameters::Default_named_parameters>
  void calc_sizing_map(FaceGraph& face_graph
                     , const NamedParameters& np)
  {
    //todo ip: please check if this is good enough to store curvature
    typedef Principal_curvatures_and_directions<K> Principal_curvatures;
    typedef typename CGAL::dynamic_vertex_property_t<Principal_curvatures> Vertex_curvature_tag;
    typedef typename boost::property_map<FaceGraph,
      Vertex_curvature_tag>::type Vertex_curvature_map;

    using parameters::choose_parameter;
    using parameters::get_parameter;
    typename Base::FT radius = choose_parameter(get_parameter(np, internal_np::ball_radius), -1);

#ifdef CGAL_PMP_REMESHING_VERBOSE
    int oversize  = 0;
    int undersize = 0;
    int insize    = 0;
    std::cout << "Calculating sizing field..." << std::endl;
#endif

    Vertex_curvature_map vertex_curvature_map = get(Vertex_curvature_tag(), face_graph);
    interpolated_corrected_principal_curvatures_and_directions(face_graph,
                                                               vertex_curvature_map,
                                                               parameters::ball_radius(radius));
    // calculate vertex sizing field L(x_i) from the curvature field
    for(vertex_descriptor v : vertices(face_graph))
    {
      auto vertex_curv = get(vertex_curvature_map, v);
      const FT max_absolute_curv = (CGAL::max)(CGAL::abs(vertex_curv.max_curvature),
                                               CGAL::abs(vertex_curv.min_curvature));
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
    return FT(CGAL::squared_distance(get(m_vpmap, va), get(m_vpmap, vb)));
  }

  FT sqlength(const halfedge_descriptor& h, const PolygonMesh& pmesh) const
  {
    return sqlength(target(h, pmesh), source(h, pmesh));
  }

public:
  FT get_sizing(const vertex_descriptor v) const
  {
    CGAL_assertion(get(m_vertex_sizing_map, v));
    return get(m_vertex_sizing_map, v);
  }

  std::optional<FT> is_too_long(const vertex_descriptor va, const vertex_descriptor vb) const
  {
    const FT sqlen = sqlength(va, vb);
    FT sqtarg_len = CGAL::square(4./3. * (CGAL::min)(get(m_vertex_sizing_map, va),
                                                     get(m_vertex_sizing_map, vb)));
    CGAL_assertion(get(m_vertex_sizing_map, va));
    CGAL_assertion(get(m_vertex_sizing_map, vb));
    if (sqlen > sqtarg_len)
      return sqlen / sqtarg_len;
    else
      return std::nullopt;
  }

  std::optional<FT> is_too_short(const halfedge_descriptor h, const PolygonMesh& pmesh) const
  {
    const FT sqlen = sqlength(h, pmesh);
    FT sqtarg_len = CGAL::square(4./5. * (CGAL::min)(get(m_vertex_sizing_map, source(h, pmesh)),
                                                     get(m_vertex_sizing_map, target(h, pmesh))));
    CGAL_assertion(get(m_vertex_sizing_map, source(h, pmesh)));
    CGAL_assertion(get(m_vertex_sizing_map, target(h, pmesh)));
    if (sqlen < sqtarg_len)
      return sqlen / sqtarg_len;
    else
      return std::nullopt;
  }

  Point_3 split_placement(const halfedge_descriptor h, const PolygonMesh& pmesh) const
  {
    return midpoint(get(m_vpmap, target(h, pmesh)),
                    get(m_vpmap, source(h, pmesh)));
  }

  void update(const vertex_descriptor v, const PolygonMesh& pmesh)
  {
    // calculating it as the average of two vertices on other ends
    // of halfedges as updating is done during an edge split
    FT vertex_size = 0;
    CGAL_assertion(CGAL::halfedges_around_target(v, pmesh).size() == 2);
    for (halfedge_descriptor ha: CGAL::halfedges_around_target(v, pmesh))
    {
      vertex_size += get(m_vertex_sizing_map, source(ha, pmesh));
    }
    vertex_size /= CGAL::halfedges_around_target(v, pmesh).size();

    put(m_vertex_sizing_map, v, vertex_size);
  }

private:
  const FT tol;
  const FT m_short;
  const FT m_long;
  const VPMap m_vpmap;
  VertexSizingMap m_vertex_sizing_map;
};

}//end namespace Polygon_mesh_processing
}//end namespace CGAL

#endif //CGAL_PMP_REMESHING_ADAPTIVE_SIZING_FIELD_H
