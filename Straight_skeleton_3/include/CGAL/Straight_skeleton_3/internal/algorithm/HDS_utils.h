// Copyright (c) 2024-2025 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_HDS_UTILS_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_HDS_UTILS_H

#include <CGAL/Straight_skeleton_3/internal/debug.h>
#include <CGAL/Straight_skeleton_3/internal/kernel/Kernel_wrapper.h>
#include <CGAL/Straight_skeleton_3/internal/HDS/Polyhedron.h>
#include <CGAL/Straight_skeleton_3/Straight_skeleton_3.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/Polyhedron_transformation.h>

#include <CGAL/assertions.h>

#include <memory>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {
namespace algorithm {

template <typename GeomTraits>
class Polyhedron_transformation;

template <typename GeomTraits>
class Hds_utils
{
  using FT = typename GeomTraits::FT;
  using Point_3 = typename GeomTraits::Point_3;
  using Vector_3 = typename GeomTraits::Vector_3;
  using Line_3 = typename GeomTraits::Line_3;
  using Plane_3 = typename GeomTraits::Plane_3;

private:
  using Polyhedron = HDS::Polyhedron<GeomTraits>;
  using PolyhedronSPtr = typename Polyhedron::PolyhedronSPtr;

  using Vertex = typename Polyhedron::Vertex;
  using VertexSPtr = typename Polyhedron::VertexSPtr;
  using Skeleton_vertex_data = typename Polyhedron::Skeleton_vertex_data;
  using SkelVertexDataSPtr = typename Polyhedron::SkelVertexDataSPtr;
  using Edge = typename Polyhedron::Edge;
  using EdgeWPtr = typename Polyhedron::EdgeWPtr;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;
  using Skeleton_edge_data = typename Polyhedron::Skeleton_edge_data;
  using SkelEdgeDataSPtr = typename Polyhedron::SkelEdgeDataSPtr;
  using Facet = typename Polyhedron::Facet;
  using FacetWPtr = typename Polyhedron::FacetWPtr;
  using FacetSPtr = typename Polyhedron::FacetSPtr;
  using Skeleton_facet_data = typename Polyhedron::Skeleton_facet_data;
  using SkelFacetDataSPtr = typename Polyhedron::SkelFacetDataSPtr;

private:
  using Straight_skeleton_3 = Straight_skeleton_3<GeomTraits>;
  using StraightSkeletonWPtr = std::weak_ptr<Straight_skeleton_3>;
  using StraightSkeletonSPtr = std::shared_ptr<Straight_skeleton_3>;

  using NodeSPtr = typename Straight_skeleton_3::NodeSPtr;
  using ArcSPtr = typename Straight_skeleton_3::ArcSPtr;
  using SheetSPtr = typename Straight_skeleton_3::SheetSPtr;

private:
  using Kernel_wrapper = kernel::Kernel_wrapper<GeomTraits>;
  using Transformation = algorithm::Polyhedron_transformation<GeomTraits>;

public:
  static Line_3 line(const EdgeSPtr& edge)
  {
    CGAL_SS3_DEBUG_SPTR(edge);
    Line_3 result;
    VertexSPtr vertex_src = edge->source();
    VertexSPtr vertex_tgt = edge->target();
    if (vertex_src->point() == vertex_tgt->point()) {
      FacetSPtr facet_l = edge->get_facet_L();
      FacetSPtr facet_r = edge->get_facet_R();
      FacetSPtr facet_src = edge->get_facet_src();
      FacetSPtr facet_tgt = edge->get_facet_tgt();

      Plane_3 offset_plane_l = Transformation::shift_plane(facet_l, -1);
      Plane_3 offset_plane_r = Transformation::shift_plane(facet_r, -1);
      Plane_3 offset_plane_src = Transformation::shift_plane(facet_src, -1);
      Plane_3 offset_plane_tgt = Transformation::shift_plane(facet_tgt, -1);

      std::optional<Point_3> p_src = Kernel_wrapper::intersection(offset_plane_src, offset_plane_l, offset_plane_r);
      std::optional<Point_3> p_tgt = Kernel_wrapper::intersection(offset_plane_tgt, offset_plane_l, offset_plane_r);
      if (!p_src || !p_tgt) {
        return { };
      }
      Vector_3 v_dir { *p_tgt - *p_src };
      result = Line_3 { vertex_src->point(), v_dir };
    } else {
      result = edge->line();
    }
    return result;
  }

public:
  static ArcSPtr get_arc(const VertexSPtr& vertex)
  {
    CGAL_SS3_DEBUG_SPTR(vertex);
    CGAL_precondition(vertex->has_data());
    SkelVertexDataSPtr data = std::dynamic_pointer_cast<Skeleton_vertex_data>(vertex->get_data());
    CGAL_SS3_DEBUG_SPTR(data);
    const ArcSPtr& arc = data->get_arc();
    CGAL_SS3_DEBUG_SPTR(arc);
    return arc;
  }

  static void set_arc(const VertexSPtr& vertex,
                      const ArcSPtr& arc)
  {
    CGAL_SS3_DEBUG_SPTR(vertex);
    CGAL_SS3_DEBUG_SPTR(arc);
    CGAL_precondition(vertex->has_data());
    SkelVertexDataSPtr data = std::dynamic_pointer_cast<Skeleton_vertex_data>(vertex->get_data());
    CGAL_SS3_DEBUG_SPTR(data);
    data->set_arc(arc);
  }

  static NodeSPtr get_node(const VertexSPtr& vertex)
  {
    CGAL_SS3_DEBUG_SPTR(vertex);
    CGAL_precondition(vertex->has_data());
    SkelVertexDataSPtr data = std::dynamic_pointer_cast<Skeleton_vertex_data>(vertex->get_data());
    CGAL_SS3_DEBUG_SPTR(data);
    const NodeSPtr& node = data->get_node();
    CGAL_SS3_DEBUG_SPTR(node);
    return node;
  }

  static void set_node(const VertexSPtr& vertex,
                       const NodeSPtr& node)
  {
    CGAL_SS3_DEBUG_SPTR(vertex);
    CGAL_SS3_DEBUG_SPTR(node);
    CGAL_precondition(vertex->has_data());
    SkelVertexDataSPtr data = std::dynamic_pointer_cast<Skeleton_vertex_data>(vertex->get_data());
    CGAL_SS3_DEBUG_SPTR(data);
    data->set_node(node);
  }

  static bool has_final_point(const VertexSPtr& vertex)
  {
    CGAL_SS3_DEBUG_SPTR(vertex);
    CGAL_precondition(vertex->has_data());
    SkelVertexDataSPtr data = std::dynamic_pointer_cast<Skeleton_vertex_data>(vertex->get_data());
    CGAL_SS3_DEBUG_SPTR(data);
    return data->has_final_point();
  }

  static const Point_3& get_final_point(const VertexSPtr& vertex)
  {
    CGAL_SS3_DEBUG_SPTR(vertex);
    CGAL_precondition(vertex->has_data());
    SkelVertexDataSPtr data = std::dynamic_pointer_cast<Skeleton_vertex_data>(vertex->get_data());
    CGAL_SS3_DEBUG_SPTR(data);
    return data->get_final_point();
  }

  static void set_final_point(const VertexSPtr& vertex,
                              const std::optional<Point_3>& opt)
  {
    CGAL_SS3_DEBUG_SPTR(vertex);
    CGAL_precondition(vertex->has_data());
    SkelVertexDataSPtr data = std::dynamic_pointer_cast<Skeleton_vertex_data>(vertex->get_data());
    CGAL_SS3_DEBUG_SPTR(data);
    return data->set_final_point(opt);
  }

  /**
    * updates the final position of the vertex of a polyhedron,
    * computed from the planes of its incident faces.
    */
  static bool reset_final_point(const VertexSPtr& vertex,
                                const FT& offset_future_bound)
  {
    CGAL_SS3_TRANSF_TRACE_V(16, "reset_final_point() of " << vertex->to_string());
    CGAL_SS3_DEBUG_SPTR(vertex);

    std::array<const Plane_3*, 3> final_planes;
    unsigned int i = 0;
    typename std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
    while (i < 3 && it_f != vertex->facets().end()) {
      FacetWPtr facet_wptr = *it_f++;
      if (FacetSPtr facet = facet_wptr.lock()) {
        if (!has_final_plane(facet)) {
          reset_final_plane(facet, offset_future_bound);
        }

        final_planes[i++] = &(get_final_plane(facet));

        CGAL_SS3_TRANSF_TRACE_V(32,"  Facet " << facet->id() << " [" << get_final_plane(facet) << "]");
      }
    }
    CGAL_postcondition(i == 3);

    std::optional<Point_3> point = Kernel_wrapper::intersection(*(final_planes[0]),
                                                                *(final_planes[1]),
                                                                *(final_planes[2]));
    if (!point) {
      CGAL_SS3_TRANSF_TRACE_V(1,"Warning: triplet of planes does not define a point!");
      return false;
    }

    CGAL_SS3_TRANSF_TRACE_V(32,"  New final point = " << *point);

    CGAL_precondition(vertex->has_data());
    SkelVertexDataSPtr data = std::dynamic_pointer_cast<Skeleton_vertex_data>(vertex->get_data());
    CGAL_SS3_DEBUG_SPTR(data);
    data->set_final_point(*point);

    return true;
  }

  static const Point_3& get_final_point(const VertexSPtr& vertex,
                                        const FT& offset_future_bound)
  {
    CGAL_SS3_DEBUG_SPTR(vertex);
    CGAL_precondition(vertex->has_data());
    SkelVertexDataSPtr data = std::dynamic_pointer_cast<Skeleton_vertex_data>(vertex->get_data());
    CGAL_SS3_DEBUG_SPTR(data);

    // @todo technically we should check if it's the correct point
    if (!data->has_final_point()) {
      reset_final_point(vertex, offset_future_bound);
    }
    return data->get_final_point();
  }

public:
  static std::optional<FT> get_vanish_time(const EdgeSPtr& edge)
  {
    CGAL_SS3_DEBUG_SPTR(edge);
    CGAL_precondition(edge->has_data());
    SkelEdgeDataSPtr data = std::dynamic_pointer_cast<Skeleton_edge_data>(edge->get_data());
    CGAL_SS3_DEBUG_SPTR(data);
    return data->get_vanish_time();
  }

  static void set_vanish_time(const EdgeSPtr& edge,
                              const std::optional<FT>& vanish_time)
  {
    CGAL_SS3_DEBUG_SPTR(edge);
    CGAL_precondition(edge->has_data());
    SkelEdgeDataSPtr data = std::dynamic_pointer_cast<Skeleton_edge_data>(edge->get_data());
    CGAL_SS3_DEBUG_SPTR(data);
    data->set_vanish_time(vanish_time);
  }

  static SheetSPtr get_sheet(const EdgeSPtr& edge)
  {
    CGAL_SS3_DEBUG_SPTR(edge);
    CGAL_precondition(edge->has_data());
    SkelEdgeDataSPtr data = std::dynamic_pointer_cast<Skeleton_edge_data>(edge->get_data());
    CGAL_SS3_DEBUG_SPTR(data);
    const SheetSPtr& sheet = data->get_sheet();
    CGAL_SS3_DEBUG_SPTR(data);
    return sheet;
  }

  static void set_sheet(const EdgeSPtr& edge,
                        const SheetSPtr& sheet)
  {
    CGAL_SS3_DEBUG_SPTR(edge);
    CGAL_SS3_DEBUG_SPTR(sheet);
    CGAL_precondition(edge->has_data());
    SkelEdgeDataSPtr data = std::dynamic_pointer_cast<Skeleton_edge_data>(edge->get_data());
    CGAL_SS3_DEBUG_SPTR(data);
    data->set_sheet(sheet);
  }

  static void clear_sheet(const EdgeSPtr& edge)
  {
    CGAL_SS3_DEBUG_SPTR(edge);
    CGAL_precondition(edge->has_data());
    SkelEdgeDataSPtr data = std::dynamic_pointer_cast<Skeleton_edge_data>(edge->get_data());
    CGAL_SS3_DEBUG_SPTR(data);
    data->set_sheet(SheetSPtr());
  }

public:
  // just for convenience to avoid all the verbosity of extracting the data
  static const FT& get_speed(const FacetSPtr& facet)
  {
    CGAL_SS3_DEBUG_SPTR(facet);
    CGAL_precondition(facet->has_data());
    SkelFacetDataSPtr data = std::dynamic_pointer_cast<Skeleton_facet_data>(facet->get_data());
    CGAL_SS3_DEBUG_SPTR(data);
    return data->get_speed();
  }

  static void set_speed(const FacetSPtr& facet,
                        const FT& speed)
  {
    CGAL_SS3_DEBUG_SPTR(facet);
    CGAL_precondition(facet->has_data());
    SkelFacetDataSPtr data = std::dynamic_pointer_cast<Skeleton_facet_data>(facet->get_data());
    CGAL_SS3_DEBUG_SPTR(data);
    data->set_speed(speed);
  }

  static FacetSPtr get_facet_origin(const FacetSPtr& facet)
  {
    CGAL_SS3_DEBUG_SPTR(facet);
    CGAL_precondition(facet->has_data());
    SkelFacetDataSPtr data = std::dynamic_pointer_cast<Skeleton_facet_data>(facet->get_data());
    CGAL_SS3_DEBUG_SPTR(data);
    return data->get_facet_origin();
  }

  static const Plane_3& get_base_plane(const FacetSPtr& facet)
  {
    CGAL_SS3_DEBUG_SPTR(facet);
    CGAL_precondition(facet->has_data());
    SkelFacetDataSPtr data = std::dynamic_pointer_cast<Skeleton_facet_data>(facet->get_data());
    CGAL_SS3_DEBUG_SPTR(data);
    return data->get_base_plane();
  }

  static void set_base_plane(const FacetSPtr& facet, const Plane_3& plane)
  {
    CGAL_SS3_DEBUG_SPTR(facet);
    CGAL_SS3_DEBUG_SPTR(plane);
    CGAL_precondition(facet->has_data());
    SkelFacetDataSPtr data = std::dynamic_pointer_cast<Skeleton_facet_data>(facet->get_data());
    CGAL_SS3_DEBUG_SPTR(data);
    data->set_base_plane(plane);
    CGAL_SS3_DEBUG_SPTR(get_base_plane(facet));
  }

  static bool has_final_plane(const FacetSPtr& facet)
  {
    CGAL_SS3_DEBUG_SPTR(facet);
    CGAL_precondition(facet->has_data());
    SkelFacetDataSPtr data = std::dynamic_pointer_cast<Skeleton_facet_data>(facet->get_data());
    CGAL_SS3_DEBUG_SPTR(data);
    return data->has_final_plane();
  }

  static const Plane_3& get_final_plane(const FacetSPtr& facet)
  {
    CGAL_SS3_DEBUG_SPTR(facet);
    CGAL_precondition(facet->has_data());
    SkelFacetDataSPtr data = std::dynamic_pointer_cast<Skeleton_facet_data>(facet->get_data());
    CGAL_SS3_DEBUG_SPTR(data);
    return data->get_final_plane();
  }

  static void set_final_plane(const FacetSPtr& facet, const Plane_3& plane)
  {
    CGAL_SS3_DEBUG_SPTR(facet);
    CGAL_SS3_DEBUG_SPTR(plane);
    CGAL_precondition(facet->has_data());
    SkelFacetDataSPtr data = std::dynamic_pointer_cast<Skeleton_facet_data>(facet->get_data());
    CGAL_SS3_DEBUG_SPTR(data);
    return data->set_final_plane(plane);
  }

    /**
    * updates the final coefficients of the plane of a facet
    */
  // resets using the first 3 final planes, even if there are more
  static void reset_final_plane(const FacetSPtr& facet,
                                const FT& offset_future_bound)
  {
    CGAL_SS3_DEBUG_SPTR(facet);
    Plane_3 plane = Transformation::shift_plane(facet, offset_future_bound);

    CGAL_precondition(facet->has_data());
    SkelFacetDataSPtr data = std::dynamic_pointer_cast<Skeleton_facet_data>(facet->get_data());
    CGAL_SS3_DEBUG_SPTR(data);
    data->set_final_plane(plane);
  }

  static const Plane_3& get_final_plane(const FacetSPtr& facet,
                                        const FT& offset_future_bound)
  {
    CGAL_SS3_DEBUG_SPTR(facet);
    CGAL_precondition(facet->has_data());
    SkelFacetDataSPtr data = std::dynamic_pointer_cast<Skeleton_facet_data>(facet->get_data());
    CGAL_SS3_DEBUG_SPTR(data);

    // @todo technically we should check if it's the correct plane
    if (!data->has_final_plane()) {
      reset_final_plane(facet, offset_future_bound);
    }
    return data->get_final_plane();
  }

public:
  static bool is_reflex(const EdgeSPtr& edge,
                        const bool future_facing = true)
  {
    CGAL_SS3_DEBUG_SPTR(edge);

    if (edge->reflex_status()) {
      return *(edge->reflex_status());
    }

    bool result = false;

    VertexSPtr vertex_src = edge->source();
    VertexSPtr vertex_tgt = edge->target();

    // @todo is degenerate edge? pointer comparison or actual position comparison?
    if (vertex_src->point() == vertex_tgt->point()) {
      FacetSPtr facet_l = edge->get_facet_L();
      FacetSPtr facet_r = edge->get_facet_R();
      FacetSPtr facet_src = edge->get_facet_src();
      FacetSPtr facet_tgt = edge->get_facet_tgt();

      FT od = future_facing ? -1 : 1;
      const Plane_3 offset_plane_l = Transformation::shift_plane(facet_l, od);
      const Plane_3 offset_plane_r = Transformation::shift_plane(facet_r, od);
      const Plane_3 offset_plane_src = Transformation::shift_plane(facet_src, od);
      const Plane_3 offset_plane_tgt = Transformation::shift_plane(facet_tgt, od);

      std::optional<Point_3> p_src = Kernel_wrapper::intersection(offset_plane_src, offset_plane_l, offset_plane_r);
      std::optional<Point_3> p_tgt = Kernel_wrapper::intersection(offset_plane_tgt, offset_plane_l, offset_plane_r);
      CGAL_assertion(p_src && p_tgt);

      const Vector_3 v_dir { *p_src, *p_tgt };
      CGAL_assertion(v_dir != CGAL::NULL_VECTOR);

      const Vector_3 normal_l = offset_plane_l.orthogonal_vector();
      CGAL_assertion(normal_l != CGAL::NULL_VECTOR);

      const Vector_3 v_cross = CGAL::cross_product(normal_l, v_dir);
      const Point_3 p = *p_src + v_cross;
      if (Kernel_wrapper::side(offset_plane_r, p) > 0) {
        result = true;
      }
    } else {
      result = edge->is_reflex();
    }
    return result;
  }

  static bool is_reflex(const VertexSPtr& vertex)
  {
    CGAL_SS3_DEBUG_SPTR(vertex);
    if (vertex->degree() == 0) {
      return false;
    }
    bool result = true;
    for (EdgeWPtr edge_wptr : vertex->edges()) {
      if (EdgeSPtr edge = edge_wptr.lock()) {
        if (!is_reflex(edge)) {
          result = false;
          break;
        }
      }
    }
    return result;
  }

  static bool is_convex(const VertexSPtr& vertex)
  {
    CGAL_SS3_DEBUG_SPTR(vertex);
    if (vertex->degree() == 0) {
      return false;
    }
    bool result = true;
    for (EdgeWPtr edge_wptr : vertex->edges()) {
      if (EdgeSPtr edge = edge_wptr.lock()) {
        if (is_reflex(edge)) {
          result = false;
          break;
        }
      }
    }
    return result;
  }

  /**
    * Checks if the given edge is part of a triangle
    * on the surface of the polyhedron.
    * Used by: nextEdgeEvent, nextTriangleEvent
    */
  static bool is_triangle(const FacetSPtr& facet,
                          EdgeSPtr edge_begin)
  {
    CGAL_SS3_DEBUG_SPTR(facet);
    CGAL_SS3_DEBUG_SPTR(edge_begin);
    bool result = false;
    if (!facet->has_incident_edge(edge_begin)) {
      return false;
    }
    unsigned int i = 0;
    // VertexSPtr vertices[3];
    EdgeSPtr edge = edge_begin;
    for (i = 0; i < 3; ++i) {
      // vertices[i] = edge->src(facet);
      edge = edge->next(facet);
      if (!edge) {
        break;
      }
    }
    if (i == 3 && edge_begin == edge) {
      result = true;
      // if (vertices[0]->point() != vertices[1]->point() &&
      //     vertices[1]->point() != vertices[2]->point() &&
      //     vertices[2]->point() != vertices[0]->point()) {
      //   // check if triangle is a hole
      //   Plane_3 plane { vertices[0]->point(),
      //                   vertices[1]->point(),
      //                   vertices[2]->point() };
      //   Vector_3 v1 = plane.orthogonal_vector();
      //   Vector_3 v2 = facet->get_plane().orthogonal_vector();
      //   if (CGAL::scalar_product(v1, v2) > 0) {
      //     // angle between orthogonal vectors of planes < CGAL_PI/2.0
      //     // not a hole
      //     result = true;
      //   }
      // }
    }
    return result;
  }

  /**
    * Checks if the given edge is part of a tetrahedron.
    * Used by: nextTriangleEvent, nextTetrahedronEvent
    */
  static bool is_tetrahedron(const EdgeSPtr& edge_begin)
  {
    CGAL_SS3_DEBUG_SPTR(edge_begin);
    bool result = true;
    VertexSPtr vertices[4];
    for (unsigned int i = 0; i < 4; ++i) {
      vertices[i] = VertexSPtr();
    }
    vertices[0] = edge_begin->source();
    vertices[1] = edge_begin->target();
    EdgeSPtr edge_l = edge_begin->next(edge_begin->get_facet_L());
    vertices[2] = edge_l->tgt(edge_begin->get_facet_L());
    EdgeSPtr edge_r = edge_begin->next(edge_begin->get_facet_R());
    vertices[3] = edge_r->tgt(edge_begin->get_facet_R());
    for (unsigned int i = 0; i < 4; ++i) {
      if (!vertices[i]) {
        result = false;
        return result;
      }
    }
    for (unsigned int i = 0; i < 4; ++i) {
      for (unsigned int j = 1; j < 4; ++j) {
        EdgeSPtr edge = vertices[i]->find_edge(vertices[(i+j)%4]);
        if (!edge) {
          result = false;
        }
      }
    }
    return result;
  }
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_HDS_UTILS_H */
