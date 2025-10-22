// Copyright (c) 2024-2025 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

/**
 * file   algo/3d/PolyhedronTransformation.h
 * author Gernot Walzl
 * date   2012-09-01
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_POLYHEDRON_TRANSFORMATION_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_POLYHEDRON_TRANSFORMATION_H

#include <CGAL/Straight_skeleton_3/internal/debug.h>
#include <CGAL/Straight_skeleton_3/internal/kernel/Kernel_factory.h>
#include <CGAL/Straight_skeleton_3/internal/kernel/Kernel_wrapper.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/Geom_utils.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/HDS_utils.h>
#include <CGAL/Straight_skeleton_3/internal/HDS/Polyhedron.h>
#include <CGAL/Straight_skeleton_3/IO/String_factory.h>
#include <CGAL/Straight_skeleton_3/Configuration.h>

#ifdef CGAL_SS3_DUMP_FILES
# include <CGAL/Straight_skeleton_3/IO/OBJ.h>
#endif

#include <CGAL/Projection_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/mark_domain_in_triangulation.h>
#include <CGAL/simplest_rational_in_interval.h>
#include <CGAL/Polygon_mesh_processing/autorefinement.h> // only for double_ceil
#include <CGAL/Polygon_mesh_processing/internal/triangle_soup_snap_rounding.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <list>
#include <map>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <stack>
#include <unordered_map>
#include <vector>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {
namespace algorithm {

template <typename Traits>
class Hds_utils;

template <typename Traits>
class Polyhedron_transformation
{
  using FT = typename Traits::FT;
  using Point_3 = typename Traits::Point_3;
  using Segment_3 = typename Traits::Segment_3;
  using Vector_3 = typename Traits::Vector_3;
  using Line_3 = typename Traits::Line_3;
  using Plane_3 = typename Traits::Plane_3;

  using Point3SPtr = std::shared_ptr<Point_3>;
  using Segment3SPtr = std::shared_ptr<Segment_3>;
  using Vector3SPtr = std::shared_ptr<Vector_3>;
  using Line3SPtr = std::shared_ptr<Line_3>;
  using Plane3SPtr = std::shared_ptr<Plane_3>;

private:
  using Polyhedron = HDS::Polyhedron<Traits>;
  using PolyhedronSPtr = typename Polyhedron::PolyhedronSPtr;

  using Vertex = typename Polyhedron::Vertex;
  using VertexSPtr = typename Polyhedron::VertexSPtr;
  using Edge = typename Polyhedron::Edge;
  using EdgeWPtr = typename Polyhedron::EdgeWPtr;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;
  using Facet = typename Polyhedron::Facet;
  using FacetWPtr = typename Polyhedron::FacetWPtr;
  using FacetSPtr = typename Polyhedron::FacetSPtr;

  using Skeleton_vertex_data = typename Polyhedron::Skeleton_vertex_data;
  using SkelVertexDataSPtr = typename Polyhedron::SkelVertexDataSPtr;
  using Skeleton_edge_data = typename Polyhedron::Skeleton_edge_data;
  using SkelEdgeDataSPtr = typename Polyhedron::SkelEdgeDataSPtr;
  using Skeleton_facet_data = typename Polyhedron::Skeleton_facet_data;
  using SkelFacetDataSPtr = typename Polyhedron::SkelFacetDataSPtr;

private:
  using Kernel_factory = kernel::Kernel_factory<Traits>;
  using Kernel_wrapper = kernel::Kernel_wrapper<Traits>;
  using Geom_utils = algorithm::Geom_utils<Traits>;
  using Hds_utils = algorithm::Hds_utils<Traits>;

public:
  static Point3SPtr bounding_box_min(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    FT p_min[3];
    for (unsigned int i = 0; i < 3; ++i) {
      p_min[i] = std::numeric_limits<double>::max();
    }
    for (const VertexSPtr& vertex : polyhedron->vertices()) {
      Point3SPtr p = vertex->point();
      for (unsigned int i = 0; i < 3; ++i) {
        if ((*p)[i] < p_min[i]) {
          p_min[i] = (*p)[i];
        }
      }
    }
    Point3SPtr result = Kernel_factory::createPoint3(p_min[0], p_min[1], p_min[2]);
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  static Point3SPtr bounding_box_max(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    FT p_max[3];
    for (unsigned int i = 0; i < 3; ++i) {
      p_max[i] = -std::numeric_limits<double>::max();
    }
    for (const VertexSPtr& vertex : polyhedron->vertices()) {
      Point3SPtr p = vertex->point();
      for (unsigned int i = 0; i < 3; ++i) {
        if ((*p)[i] > p_max[i]) {
          p_max[i] = (*p)[i];
        }
      }
    }
    Point3SPtr result = Kernel_factory::createPoint3(p_max[0], p_max[1], p_max[2]);
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  static bool is_inside_bbox(const PolyhedronSPtr& polyhedron,
                             const Point3SPtr& p_box_min, const Point3SPtr& p_box_max)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    CGAL_SS3_DEBUG_SPTR(p_box_min);
    CGAL_SS3_DEBUG_SPTR(p_box_max);
    bool result = true;
    for (const VertexSPtr& vertex : polyhedron->vertices()) {
      Point3SPtr p = vertex->point();
      for (unsigned int i = 0; i < 3; ++i) {
        if (!((*p_box_min)[i] <= (*p)[i] && (*p)[i] <= (*p_box_max)[i])) {
          result = false;
          // CGAL_SS3_TRANSF_TRACE(*p << " is not in the box " << *p_box_min << " " << *p_box_max);
          break;
        }
      }
      if (!result) {
        break;
      }
    }
    return result;
  }

  static void translate(const PolyhedronSPtr& polyhedron,
                        const Vector3SPtr& v_t)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    CGAL_SS3_DEBUG_SPTR(v_t);

    for (const VertexSPtr& vertex : polyhedron->vertices()) {
      Point3SPtr p = vertex->point();
      Point3SPtr p_t = Kernel_factory::createPoint3(*p + *v_t);
      vertex->set_point(p_t);
    }

    polyhedron->init_planes();
  }

  static void scale(const PolyhedronSPtr& polyhedron,
                    const Vector3SPtr& v_s)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    CGAL_SS3_DEBUG_SPTR(v_s);

    for (const VertexSPtr& vertex : polyhedron->vertices()) {
      Point3SPtr p = vertex->point();
      Point3SPtr p_s = Kernel_factory::createPoint3((*p)[0] * (*v_s)[0],
                                                   (*p)[1] * (*v_s)[1],
                                                   (*p)[2] * (*v_s)[2]);
      vertex->set_point(p_s);
    }

    polyhedron->init_planes();
  }

  static void translate_and_scale(const PolyhedronSPtr& polyhedron,
                                  const Point3SPtr& p_box_min, const Point3SPtr& p_box_max)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    CGAL_SS3_DEBUG_SPTR(p_box_min);
    CGAL_SS3_DEBUG_SPTR(p_box_max);

    Vector3SPtr v_box_min = Kernel_factory::createVector3(p_box_min);
    Vector3SPtr v_box_max = Kernel_factory::createVector3(p_box_max);
    Vector3SPtr v_size = Kernel_factory::createVector3(*v_box_max - *v_box_min);
    Vector3SPtr v_center = Kernel_factory::createVector3((*v_box_min + *v_box_max) / 2.0);

    Point3SPtr p_box_min_curr = bounding_box_min(polyhedron);
    Point3SPtr p_box_max_curr = bounding_box_max(polyhedron);
    Vector3SPtr v_box_min_curr = Kernel_factory::createVector3(p_box_min_curr);
    Vector3SPtr v_box_max_curr = Kernel_factory::createVector3(p_box_max_curr);
    Vector3SPtr v_size_curr = Kernel_factory::createVector3(*v_box_max_curr - *v_box_min_curr);
    Vector3SPtr v_center_curr = Kernel_factory::createVector3((*v_box_min_curr + *v_box_max_curr) / 2.0);

    FT scale_factor = std::numeric_limits<double>::max(); // do not put FT
    for (unsigned int i = 0; i < 3; ++i) {
      FT s = (*v_size)[i]/(*v_size_curr)[i];
      if (scale_factor > s) {
        scale_factor = s;
      }
    }
    scale_factor = floor(CGAL::to_double(scale_factor)*1000.0)/1000.0;
    Vector3SPtr v_s = Kernel_factory::createVector3(scale_factor, scale_factor, scale_factor);

    Vector3SPtr v_t = Kernel_factory::createVector3((*v_center_curr) * -1.0);
    if (v_t->squared_length() > 0.0) {
      translate(polyhedron, v_t);
    }
    if (scale_factor != 1.0) {
      scale(polyhedron, v_s);
    }
    if (v_center->squared_length() > 0.0) {
      translate(polyhedron, v_center);
    }
  }

  static void truncate_precision(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);

    ConfigurationSPtr config = Configuration::get_instance();
    double range = 1e-10;
    if (config->is_loaded()) {
      range = config->get_double("Preprocessing", "truncate_precision");
    }

    if (range == 0.) {
      return;
    }

    CGAL_SS3_TRANSF_TRACE_V(4, "Lower precision of input polyhedron");
    CGAL_SS3_TRANSF_TRACE_V(8, " truncate precision: " << range);

    double exp = std::ceil(std::log2(1.0 / range));
    double scale = std::pow(2.0, exp);
    CGAL_SS3_TRANSF_TRACE_V(8,"  scale = " << scale);

    for (const VertexSPtr& vertex : polyhedron->vertices()) {
      Point3SPtr p = vertex->point();
      CGAL_SS3_TRANSF_TRACE_V(32,"    Truncated from: " << *(vertex->point()));

      double rx = CGAL::Polygon_mesh_processing::autorefine_impl::double_ceil(p->x() * scale) / scale;
      double ry = CGAL::Polygon_mesh_processing::autorefine_impl::double_ceil(p->y() * scale) / scale;
      double rz = CGAL::Polygon_mesh_processing::autorefine_impl::double_ceil(p->z() * scale) / scale;

      vertex->set_point(Kernel_factory::createPoint3(rx, ry, rz));
      CGAL_SS3_TRANSF_TRACE_V(32,"    Truncated to: " << *(vertex->point()));
    }
  }

  static bool has_coplanar_facets(const EdgeSPtr& edge,
                                  const double epsilon)
  {
    CGAL_SS3_DEBUG_SPTR(edge);

    bool result = false;
    FacetSPtr facet_l = edge->get_facet_L();
    FacetSPtr facet_r = edge->get_facet_R();
    if (facet_l && facet_r) {
      if (epsilon == 0.) {
        return (*(facet_l->plane()) == *(facet_r->plane())); // planes are normalized
      }

      Vector3SPtr normal_l = Kernel_factory::createVector3(facet_l->plane());
      Vector3SPtr normal_r = Kernel_factory::createVector3(facet_r->plane());
      FT length_l = 0.0;
      FT length_r = 0.0;
      for (unsigned int i = 0; i < 3; ++i) {
        length_l += (*normal_l)[i] * (*normal_l)[i];
        length_r += (*normal_r)[i] * (*normal_r)[i];
      }
      // tolerate this sqrt because it does not matter for robustness
      length_l = CGAL::approximate_sqrt(length_l);
      length_r = CGAL::approximate_sqrt(length_r);
      FT diff = 0.0;
      FT diff_sq_length = 0.0;
      for (unsigned int i = 0; i < 3; ++i) {
        diff = ((*normal_l)[i]/length_l) - ((*normal_r)[i]/length_r);
        diff_sq_length += diff*diff;
      }
      result = (diff_sq_length < CGAL::square(epsilon));
    }
    return result;
  }

  static void merge_facets(const EdgeSPtr& edge,
                           const FacetSPtr& facet_into,
                           const FacetSPtr& facet_from,
                           const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_DEBUG_SPTR(edge);
    CGAL_SS3_DEBUG_SPTR(facet_into);
    CGAL_SS3_DEBUG_SPTR(facet_from);
    CGAL_precondition(facet_into != facet_from);
    CGAL_precondition(edge->get_facet_L() == facet_from || edge->get_facet_R() == facet_from);
    CGAL_precondition(edge->get_facet_L() == facet_into || edge->get_facet_R() == facet_into);

    CGAL_SS3_TRANSF_TRACE_V(16, "Merging F" << facet_from->get_ID() << " into F" << facet_into->get_ID() <<
                                " Common edge E" << edge->get_ID() << " [V" << edge->get_vertex_src()->get_ID()
                                                                  << " - V" << edge->get_vertex_dst()->get_ID() << "]");
    CGAL_SS3_TRANSF_TRACE_V(16, "  FROM normal: " << *(Kernel_factory::createVector3(facet_from->get_plane())));
    CGAL_SS3_TRANSF_TRACE_V(16, "  INTO normal: " << *(Kernel_factory::createVector3(facet_into->get_plane())));

    // First remove the facet incidence info from the edge such that the facets are not deleted
    // when Polyhedron::remove_edge() is called
    facet_into->remove_edge(edge);
    facet_from->remove_edge(edge);
    polyhedron->remove_edge(edge);

    facet_into->merge(facet_from);
    polyhedron->remove_facet(facet_from);

    CGAL_postcondition(polyhedron->is_consistent());
  }

  static void merge_facets(const EdgeSPtr& edge,
                          const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_DEBUG_SPTR(edge);
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    return merge_facets(edge, edge->get_facet_L(), edge->get_facet_R(), polyhedron);
  }

  static int merge_coplanar_facets(const PolyhedronSPtr& polyhedron,
                                   const double epsilon)
  {
    CGAL_SS3_TRANSF_TRACE_V(4, "\nMerging coplanar faces with epsilon = " << epsilon);
    CGAL_SS3_TRANSF_TRACE_V(4, "  initial facet count: " << polyhedron->facets().size());

    CGAL_SS3_DEBUG_SPTR(polyhedron);

#ifdef CGAL_SS3_DUMP_FILES
    IO::OBJFile::save("results/coplanar_merge_before.obj", polyhedron, false /*do not triangulate*/);
#endif

    int result = 0;
    std::list<EdgeWPtr> edges_toremove;
    for (const EdgeSPtr& edge : polyhedron->edges()) {
      if (has_coplanar_facets(edge, epsilon)) {
        edges_toremove.push_back(edge);
      }
    }

    CGAL_SS3_TRANSF_TRACE(edges_toremove.size() << " edges to merge");
    CGAL_SS3_TRANSF_TRACE_CODE(if (edges_toremove.size() > 0))
    CGAL_SS3_TRANSF_TRACE("Adjacent facets of the following edges are detected to be coplanar and will be merged.");

    for (EdgeWPtr edge_w : edges_toremove) {
      if (EdgeSPtr edge = edge_w.lock()) {
        merge_facets(edge, polyhedron);
        ++result;
      }
    }

    IO::OBJFile::save("results/post_merge.obj", polyhedron, false /*do not triangulate*/);
    std::cout << polyhedron->to_string() << std::endl;

    // There can be multiple defects after merging, for example:
    // - when we merge two facets sharing multiple edges, this can leave dangling vertices
    // - we can end up with degree 2 vertices
    // - ...

    sanitize(polyhedron); // remove degenerate vertices and facets
    CGAL_postcondition(polyhedron->is_consistent());

    polyhedron->initialize_all_IDs();

    CGAL_SS3_TRANSF_TRACE_V(4, "  final facet count: " << polyhedron->facets().size());

#ifdef CGAL_SS3_DUMP_FILES
    IO::OBJFile::save("results/coplanar_merge_after.obj", polyhedron,
                      false /*do not triangulate*/);
#endif

    return result;
  }

  static int merge_coplanar_facets(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    double epsilon = 0.0;
    ConfigurationSPtr config = Configuration::get_instance();
    if (config->is_loaded()) {
      std::string section("Preprocessing");
      std::string key("coplanarity_epsilon");
      if (config->contains(section, key)) {
        epsilon = config->get_double(section, key);
      }
    }

    return merge_coplanar_facets(polyhedron, epsilon);
  }

  static int remove_vertices_deg_lt3(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_TRANSF_TRACE("Remove vertices with degree < 3");

    CGAL_SS3_DEBUG_SPTR(polyhedron);
    CGAL_SS3_TRANSF_TRACE("  initial vertex count: " << polyhedron->vertices().size());

    int result = 0;
    std::list<VertexSPtr> vertices_toremove;
    typename std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
      VertexSPtr vertex = *it_v++;
      if (vertex->degree() < 3) {
        vertices_toremove.push_back(vertex);
        CGAL_SS3_TRANSF_TRACE("Enlist " << vertex->to_string());
        for (FacetWPtr wf : vertex->facets()) {
          FacetSPtr facet = wf.lock();
          CGAL_SS3_TRANSF_TRACE("  Incident facet with: " << facet->vertices().size() << " vertices");
        }
      }
    }
    it_v = vertices_toremove.begin();
    while (it_v != vertices_toremove.end()) {
      VertexSPtr vertex = *it_v++;
      CGAL_SS3_TRANSF_TRACE("Removing " << vertex->to_string());

      if (vertex->degree() == 0) {
        // degree 0 so there are no incident edges, but it might still be a vertex incident to a face...
        typename std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
        while (it_f != vertex->facets().end()) { // no C++11
          FacetWPtr facet_wptr = *it_f++;
          if (FacetSPtr facet = facet_wptr.lock()) {
            facet->remove_vertex(vertex);
          }
        }
        polyhedron->remove_vertex(vertex);
      } else if (vertex->degree() == 1) {
        EdgeSPtr edge = vertex->edges().front().lock();
        CGAL_SS3_DEBUG_SPTR(edge);

        typename std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
        while (it_f != vertex->facets().end()) { // no C++11
          FacetWPtr facet_wptr = *it_f++;
          if (FacetSPtr facet = facet_wptr.lock()) {
            facet->remove_vertex(vertex);
          }
        }
        for (FacetWPtr facet_wptr : { edge->get_facet_L(), edge->get_facet_R() }) {
          if (FacetSPtr facet = facet_wptr.lock()) {
            facet->remove_edge(edge);
          }
        }

        polyhedron->remove_vertex(vertex); // removes the edge too
      } else if (vertex->degree() == 2) {
        EdgeSPtr edge_src = vertex->first_edge();
        EdgeSPtr edge_dst = edge_src->next(vertex);
        CGAL_assertion(edge_src->has_vertex(vertex));
        CGAL_assertion(edge_dst->has_vertex(vertex));

        VertexSPtr vertex_src = edge_src->other(vertex);
        VertexSPtr vertex_dst = edge_dst->other(vertex);

        FacetSPtr fL = edge_src->get_facet_L();
        FacetSPtr fR = edge_src->get_facet_R();
        CGAL_assertion(fL && fR);
        CGAL_assertion(fL != fR);

        if (fL->vertices().size() == 3) {
          CGAL_SS3_TRANSF_TRACE("Deg 2 vertex is the apex of a triangle facet (fL=" << fL->get_ID() << ")");
          EdgeSPtr third_edge;
          for (const EdgeSPtr& edge : fL->edges()) {
            if (edge != edge_src && edge != edge_dst) {
              third_edge = edge;
              break;
            }
          }
          CGAL_assertion(third_edge != nullptr);
          merge_facets(third_edge, polyhedron);
        }

        if (fR->vertices().size() == 3) {
          CGAL_SS3_TRANSF_TRACE("Deg 2 vertex is the apex of a triangle facet (fR=" << fR->get_ID() << ")");
          EdgeSPtr third_edge;
          for (const EdgeSPtr& edge : fR->edges()) {
            if (edge != edge_src && edge != edge_dst) {
              third_edge = edge;
              break;
            }
          }
          CGAL_assertion(third_edge != nullptr);
          merge_facets(third_edge, polyhedron);
        }

        typename std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
        while (it_f != vertex->facets().end()) { // no C++11
          FacetWPtr facet_wptr = *it_f++;
          if (FacetSPtr facet = facet_wptr.lock()) {
            facet->remove_vertex(vertex);
          }
        }

        edge_dst->get_facet_L()->remove_edge(edge_dst);
        edge_dst->get_facet_R()->remove_edge(edge_dst);
        polyhedron->remove_edge(edge_dst);

        if (edge_src->get_vertex_dst() == vertex) {
          edge_src->replace_vertex_dst(vertex_dst);
        } else if (edge_src->get_vertex_src() == vertex) {
          edge_src->replace_vertex_src(vertex_dst);
        } else {
          CGAL_assertion(false);
        }

        polyhedron->remove_vertex(vertex);
      }

      ++result;
    }

    CGAL_SS3_TRANSF_TRACE("  final vertex count: " << polyhedron->vertices().size());
    CGAL_postcondition(polyhedron->is_consistent());

    return result;
  }

  static int remove_facets_deg_lt3(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_TRANSF_TRACE("Remove facets with size < 3");

    CGAL_SS3_DEBUG_SPTR(polyhedron);
    CGAL_SS3_TRANSF_TRACE("  initial facet count: " << polyhedron->facets().size());

    int result = 0;
    std::list<FacetSPtr> facets_tomerge;
    for (const FacetSPtr& facet : polyhedron->facets()) {
      if (facet->vertices().size() < 3) {
        facets_tomerge.push_back(facet);
      }
    }
    for (const FacetSPtr& facet : facets_tomerge) {
      // Facet could have grown from another merge, so check again
      if (facet->vertices().size() >= 3) {
        continue;
      }

      // Out of the two edges of the facet, find the edge that is incident to the facet
      // that is the largest, and merge into that one
      EdgeSPtr best_edge = nullptr;
      FacetSPtr best_neighbor = nullptr;
      std::size_t best_size = 0;

      for (const EdgeSPtr& edge : facet->edges()) {
        FacetSPtr neighbor = nullptr;
        if (edge->get_facet_L() == facet && edge->get_facet_R() && edge->get_facet_R() != facet) {
          neighbor = edge->get_facet_R();
        } else if (edge->get_facet_R() == facet && edge->get_facet_L() && edge->get_facet_L() != facet) {
          neighbor = edge->get_facet_L();
        }
        if (neighbor && neighbor->vertices().size() > best_size) {
          best_edge = edge;
          best_neighbor = neighbor;
          best_size = neighbor->vertices().size();
        }
      }

      if (best_neighbor) {
        merge_facets(best_edge, best_neighbor, facet, polyhedron);
      } else {
        // No neighbor, just remove
        for (const EdgeSPtr& edge : facet->edges()) {
          polyhedron->remove_edge(edge);
        }
      }
      ++result;
    }

    CGAL_SS3_TRANSF_TRACE("  final facet count: " << polyhedron->vertices().size());
    CGAL_postcondition(polyhedron->is_consistent());

    return result;
  }

  static int sanitize(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_TRANSF_TRACE("Sanitizing polyhedron...");

    CGAL_SS3_DEBUG_SPTR(polyhedron);

    // - remove_vertices_deg_lt3 can create facets with fewer than 3 vertices
    // - remove_facets_deg_lt3 removes facets with fewer than 3 vertices
    // so loop till nothing is done anymore
    int result = 0;
    for (;;) {
      std::size_t vlt3 = remove_vertices_deg_lt3(polyhedron);
      std::size_t flt3 = remove_facets_deg_lt3(polyhedron);
      int partial = vlt3 + flt3;
      if (partial == 0) {
        break;
      }
      result += vlt3 + flt3;
    }
    return result;
  }

  static Point3SPtr offset_point_from_base(const VertexSPtr& vertex, const FT& time)
  {
    CGAL_SS3_TRANSF_TRACE_V(16, "absolute offset of " << vertex->to_string());
    CGAL_SS3_DEBUG_SPTR(vertex);

    std::array<Plane3SPtr, 3> planes;
    unsigned int i = 0;
    typename std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
    while (i < 3 && it_f != vertex->facets().end()) {
      FacetWPtr facet_wptr = *it_f++;
      if (FacetSPtr facet = facet_wptr.lock()) {
        CGAL_SS3_TRANSF_TRACE_V(16, "  Facet " << facet->get_ID());
        CGAL_assertion(facet->has_data());
        Plane3SPtr base_plane = Hds_utils::get_base_plane(facet);
        const FT& speed = Hds_utils::get_speed(facet);
        planes[i++] = Geom_utils::offset_plane(base_plane, speed*time);
      }
    }
    CGAL_postcondition(i == 3);

    Point3SPtr point = Kernel_wrapper::intersection(planes[0], planes[1], planes[2]);
    if (!point) {
      CGAL_SS3_TRANSF_TRACE_V(1, "Warning: triplet of offset planes does not define a point!");
      Point3SPtr result = Point3SPtr();
      CGAL_SS3_DEBUG_SPTR(result);
      return { };
    }

    CGAL_SS3_TRANSF_TRACE_V(16, "  New point = " << *point);

    return point;
  }

  static Segment3SPtr offset_edge_from_base(const EdgeSPtr& edge, const FT& time)
  {
    CGAL_SS3_DEBUG_SPTR(edge);
    return Kernel_factory::createSegment3(offset_point_from_base(edge->get_vertex_src(), time),
                                         offset_point_from_base(edge->get_vertex_dst(), time));
  }

  static Plane3SPtr offset_plane_from_base(const FacetSPtr& facet, const FT& time)
  {
    CGAL_SS3_DEBUG_SPTR(facet);
    CGAL_precondition(facet->has_data());
    Plane3SPtr base_plane = Hds_utils::get_base_plane(facet);
    const FT& speed = Hds_utils::get_speed(facet);
    return Geom_utils::offset_plane(base_plane, speed*time);
  }

  /**
    * updates the position of the vertex of a polyhedron, computed from the planes of
    * its incident faces.
    */
  static bool reset_point(const VertexSPtr& vertex,
                          const std::array<Plane3SPtr, 3>& planes)
  {
    CGAL_SS3_DEBUG_SPTR(vertex);
    CGAL_SS3_DEBUG_SPTR(planes[0]);
    CGAL_SS3_DEBUG_SPTR(planes[1]);
    CGAL_SS3_DEBUG_SPTR(planes[2]);

    Point3SPtr point = Kernel_wrapper::intersection(planes[0], planes[1], planes[2]);
    if (!point) {
      CGAL_SS3_TRANSF_TRACE_V(1, "Warning: triplet of planes does not define a point!");
      Point3SPtr result = Point3SPtr();
      CGAL_SS3_DEBUG_SPTR(result);
      return false;
    }

    vertex->set_point(point);
    CGAL_SS3_TRANSF_TRACE_V(16, "  New point = " << *point);

    CGAL_postcondition_code(for (FacetWPtr facet_wptr : vertex->facets()) {)
    CGAL_postcondition_code(    if (FacetSPtr facet = facet_wptr.lock()) {)
    CGAL_postcondition(             facet->get_plane()->has_on(*point));
    CGAL_postcondition_code(    })
    CGAL_postcondition_code(})

    return true;
  }

  // resets using the first 3 planes, even if there are more
  static bool reset_point(const VertexSPtr& vertex)
  {
    CGAL_SS3_TRANSF_TRACE_V(16, "reset_point() of " << vertex->to_string());
    CGAL_SS3_DEBUG_SPTR(vertex);

    std::array<Plane3SPtr, 3> planes;
    unsigned int i = 0;
    typename std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
    while (i < 3 && it_f != vertex->facets().end()) {
      FacetWPtr facet_wptr = *it_f++;
      if (FacetSPtr facet = facet_wptr.lock()) {
        planes[i++] = facet->get_plane();
        CGAL_SS3_TRANSF_TRACE_V(16, "  Facet " << facet->get_ID() << " [" << *(facet->get_plane()) << "]");
      }
    }
    CGAL_postcondition(i == 3);

    return reset_point(vertex, { planes[0], planes[1], planes[2] });
  }

  /**
    * updates the positions of the vertices of a polyhedron, computed from the planes of
    * their incident faces.
    */
  static bool reset_points(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_TRANSF_TRACE("Reset point positions");
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    for (VertexSPtr vertex : polyhedron->vertices()) {
      if (!reset_point(vertex)) {
        CGAL_SS3_TRANSF_TRACE_V(1, "Warning: failed to reset vertex " << vertex->to_string());
        return false;
      }
    }
    return true;
  }

  /**
    * returns the shifted position of the vertex of a polyhedron
    * \pre vertex->degree() == 3
    */
  static Point3SPtr shift_point(const VertexSPtr& vertex,
                                const FT& time)
  {
    CGAL_SS3_TRANSF_TRACE_V(32, "shift " << vertex->to_string() << "\ntime = " << time);
    CGAL_SS3_DEBUG_SPTR(vertex);
    CGAL_precondition(vertex->degree() >= 3);

    Plane3SPtr planes[3];
    unsigned int i = 0;
    typename std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
    while (i < 3 && it_f != vertex->facets().end()) {
      FacetWPtr facet_wptr = *it_f++;
      if (FacetSPtr facet = facet_wptr.lock()) {
        Plane3SPtr plane = facet->plane();

        if (vertex->degree() > 3) {
          // planes are _offset_ planes, but it doesn't matter for the tests
          bool independent = true;
          if (i == 1) {
            independent = !(CGAL::parallel(*(planes[0]), *plane));
          } else if (i == 2) {
            independent = !is_zero(CGAL::determinant(planes[0]->a(), planes[0]->b(), planes[0]->c(),
                                                     planes[1]->a(), planes[1]->b(), planes[1]->c(),
                                                     plane->a(), plane->b(), plane->c()));
          }

          if (!independent) {
            continue;
          }
        }

        planes[i] = shift_plane(facet, time);

        CGAL_SS3_TRANSF_TRACE_V(64, "facet[" << i << "] = " << facet->get_ID());
        CGAL_SS3_TRANSF_TRACE_V(64, "Offset Plane[" << i << "] = " << *(planes[i]));

        ++i;
      }
    }

    if (i < 3) {
      CGAL_SS3_TRANSF_TRACE_V(1, "Warning: Could not find three independent planes for high-degree vertex" << vertex->to_string());
      return { };
    }

    Point3SPtr point = Kernel_wrapper::intersection(planes[0], planes[1], planes[2]);
    if (!point) {
      CGAL_SS3_TRANSF_TRACE_V(1, "Error: triplet of planes doesn't define a point!");
      Point3SPtr result = Point3SPtr();
      CGAL_SS3_DEBUG_SPTR(result);
      return result;
    }

    CGAL_SS3_TRANSF_TRACE_V(32, "  New point = " << *point);

    CGAL_assertion_code(for (Plane3SPtr pi : planes))
    CGAL_assertion(pi->has_on(*point));

    return point;
  }

  /**
    * returns the shifted position of the edge of a polyhedron
    */
  static Segment3SPtr shift_edge(const EdgeSPtr& edge,
                                 const FT& time)
  {
    CGAL_SS3_DEBUG_SPTR(edge);

    FacetSPtr facet_l = edge->get_facet_L();
    FacetSPtr facet_r = edge->get_facet_R();
    FacetSPtr facet_src = edge->get_facet_L()->next(edge->get_vertex_src());
    FacetSPtr facet_dst = edge->get_facet_R()->next(edge->get_vertex_dst());

    const FT& speed_l = Hds_utils::get_speed(facet_l);
    const FT& speed_r = Hds_utils::get_speed(facet_r);
    const FT& speed_src = Hds_utils::get_speed(facet_src);
    const FT& speed_dst = Hds_utils::get_speed(facet_dst);

    // Offset the two common planes
    Plane3SPtr offset_plane_l = Geom_utils::offset_plane(facet_l->plane(), speed_l*time);
    Plane3SPtr offset_plane_r = Geom_utils::offset_plane(facet_r->plane(), speed_r*time);
    Plane3SPtr offset_plane_src = Geom_utils::offset_plane(facet_src->plane(), speed_src*time);
    Plane3SPtr offset_plane_dst = Geom_utils::offset_plane(facet_dst->plane(), speed_dst*time);

#if 0
    // leaving it here because it's not that intuitive: factoring the intersection of the
    // two common planes is much slower than computing two 3-plane intersections
    Line3SPtr common_line = Kernel_wrapper::intersection(offset_plane_l, offset_plane_r);
    Point3SPtr src_point = Kernel_wrapper::intersection(offset_plane_src, common_line);
    Point3SPtr dst_point = Kernel_wrapper::intersection(offset_plane_dst, common_line);
#else
    Point3SPtr src_point = Kernel_wrapper::intersection(offset_plane_src, offset_plane_l, offset_plane_r);
    Point3SPtr dst_point = Kernel_wrapper::intersection(offset_plane_dst, offset_plane_l, offset_plane_r);
#endif

    CGAL_SS3_DEBUG_SPTR(src_point);
    CGAL_SS3_DEBUG_SPTR(dst_point);
    return Kernel_factory::createSegment3(src_point, dst_point);
  }

  /**
    * returns the shifted position of the facet of a polyhedron
    */
  static Plane3SPtr shift_plane(const FacetSPtr& facet,
                                const FT& time)
  {
    CGAL_SS3_DEBUG_SPTR(facet);
    const Plane3SPtr& plane = facet->get_plane();
    const FT& speed = Hds_utils::get_speed(facet);
    return Geom_utils::offset_plane(plane, speed*time);
  }

  /**
    * Offsets the polyhedron `polyhedron`
    * Negative offset points to the interior of the polyhedron.
    * This function is for the main shift in the event loop.
    */
  static void shift_facets(const PolyhedronSPtr& polyhedron,
                           const FT& time,
                           const bool recompute_positions = true)
  {
    CGAL_SS3_TRANSF_TRACE_V(32, "~~~~ Shift polyhedron by " << time << " [in place]");
    CGAL_SS3_DEBUG_SPTR(polyhedron);

    typename std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
      FacetSPtr facet = *it_f++;
      Plane3SPtr offset_plane = shift_plane(facet, time);
      facet->set_plane(offset_plane);
    }

    typename std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
      VertexSPtr vertex = *it_v++;

      // See comment at the top of the function
      CGAL_assertion(vertex->degree() == 3);

      Point3SPtr old_point = vertex->point(), new_point;
      if (time != 0 || recompute_positions) {
        CGAL_assertion_code(bool ok =)
        reset_point(vertex);
        CGAL_assertion(ok);
      }
    }
  }

  /**
  * Offsets the polyhedron `polyhedron`, which may have degree 1 vertices.
  * Negative offset points to the interior of the polyhedron.
  */
  static bool shift_facets_deg1(const PolyhedronSPtr& polyhedron,
                                const FT& time)
  {
    CGAL_SS3_TRANSF_TRACE_V(32, "~~~~ Shift polyhedron by " << time << " [in place w/ degree 1]");
    CGAL_SS3_DEBUG_SPTR(polyhedron);

    // Maps to store shifted planes and points
    std::unordered_map<FacetSPtr, Plane3SPtr> facet_to_shifted_plane;
    std::unordered_map<VertexSPtr, Point3SPtr> vertex_to_shifted_point;

    for (const FacetSPtr& facet : polyhedron->facets()) {
      Plane3SPtr offset_plane = shift_plane(facet, time);
      facet_to_shifted_plane[facet] = offset_plane;
    }

    // Degree 1 vertices are shifted by translating the shifted adjacent degree 3+ vertex adjacent
    // to the degree 1 vertex by the same direction and distance as the unshifted vertices.
    // So, before treating degree 1 vertices, we must treat the degree 3 vertices.
    for (const VertexSPtr& vertex : polyhedron->vertices()) {
      if (vertex->degree() >= 3) {
        std::array<Plane3SPtr, 3> shifted_planes;
        unsigned int i = 0;
        typename std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
        while (i < 3 && it_f != vertex->facets().end()) {
          FacetWPtr facet_wptr = *it_f++;
          if (FacetSPtr facet = facet_wptr.lock()) {
            shifted_planes[i++] = facet_to_shifted_plane[facet];
          }
        }
        CGAL_postcondition(i == 3);

        Point3SPtr shifted_point = Kernel_wrapper::intersection(shifted_planes[0], shifted_planes[1], shifted_planes[2]);
        if (!shifted_point) {
          CGAL_SS3_TRANSF_TRACE_V(1, "Error: triplet of shifted planes doesn't define a point!");
          return false;
        }
        vertex_to_shifted_point[vertex] = shifted_point;
      }
    }

    // Now we can shift degree 1 vertices
    for (const VertexSPtr& vertex : polyhedron->vertices()) {
      if (vertex->degree() == 1) {
        EdgeSPtr edge = nullptr;
        CGAL_assertion_code(unsigned int i = 0;)
        for (EdgeWPtr edge_wptr : vertex->edges()) {
          if ((edge = edge_wptr.lock())) {
            CGAL_assertion_code(++i;)
          }
        }
        CGAL_assertion(i == 1);

        VertexSPtr vertex_other = edge->other(vertex);
        CGAL_assertion(vertex_other->degree() >= 3);

        FacetSPtr facet1 = edge->get_facet_L();
        FacetSPtr facet2 = edge->get_facet_R();
        Plane3SPtr plane1 = facet_to_shifted_plane[facet1];
        Plane3SPtr plane2 = facet_to_shifted_plane[facet2];
        Line3SPtr intersection_line = Kernel_wrapper::intersection(plane1, plane2);

        // Determine the direction using the unshifted positions
        Point3SPtr point_other = vertex_to_shifted_point[vertex_other];
        Point3SPtr point = vertex->point();
        CGAL_assertion(*point != *point_other);
        Vector_3 direction = intersection_line->to_vector();
        if (CGAL::scalar_product(direction, *point - *vertex_other->point()) < 0) {
          direction = -direction;
        }

        // Apply the shift to the vertex
        Point3SPtr offset_point = Kernel_factory::createPoint3(*point_other + direction);
        if (!offset_point) {
          return false;
        }
        vertex_to_shifted_point[vertex] = offset_point;
      }
    }

    // Apply all shifts
    for (const FacetSPtr& facet : polyhedron->facets()) {
      facet->set_plane(facet_to_shifted_plane[facet]);
    }
    for (const VertexSPtr& vertex : polyhedron->vertices()) {
      if (vertex_to_shifted_point.count(vertex)) {
        vertex->set_point(vertex_to_shifted_point[vertex]);
      }
    }

    return true;
  }

    // The tag is a template parameter because for debugging outputs,
  // we might want to tolerate intersections
  template <typename CDT2Tag>
  static auto triangulate_facet(const FacetSPtr& facet)
  {
    CGAL_SS3_DEBUG_SPTR(facet);

    using Itag = CDT2Tag;
    using PK = CGAL::Projection_traits_3<Traits>;
    using PVbb = CGAL::Triangulation_vertex_base_with_info_2<VertexSPtr, PK>;
    using PVb = CGAL::Triangulation_vertex_base_2<PK, PVbb>;
    using PFb = CGAL::Constrained_triangulation_face_base_2<PK>;
    using PTDS = CGAL::Triangulation_data_structure_2<PVb,PFb>;
    using PCDT = CGAL::Constrained_Delaunay_triangulation_2<PK, PTDS, Itag>;
    using PCDT_VH = typename PCDT::Vertex_handle;

    Vector3SPtr n = Kernel_factory::createVector3(facet->plane());
    CGAL_precondition(*n != CGAL::NULL_VECTOR);

    PK projection_traits(*n);
    PCDT pcdt(projection_traits);

    std::map<VertexSPtr, PCDT_VH> face_vhs; // might have multiple vertices at the same position

    typename std::list<VertexSPtr>::iterator it_v = facet->vertices().begin();
    while (it_v != facet->vertices().end()) {
      VertexSPtr vertex = *it_v++;
      auto res = face_vhs.emplace(vertex, PCDT_VH());
      if (res.second) // first time seeing this point
      {
        PCDT_VH vh = pcdt.insert(*(vertex->point()));
        res.first->second = vh;
        vh->info() = vertex;
      }
    }

    typename std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
    while (it_e != facet->edges().end()) {
      EdgeSPtr edge = *it_e++;
      VertexSPtr v0 = edge->src(facet);
      VertexSPtr v1 = edge->dst(facet);
      CGAL_assertion(*(v0->point()) != *(v1->point()));

      PCDT_VH vh0 = face_vhs.at(v0);
      PCDT_VH vh1 = face_vhs.at(v1);

      try {
        pcdt.insert_constraint(vh0, vh1);
      } catch(const typename PCDT::Intersection_of_constraints_exception&) {
        CGAL_SS3_TRANSF_TRACE("Error: Intersection of constraint w/ " << vh0->point() << " " << vh1->point());
        CGAL_SS3_TRANSF_TRACE(facet->to_string());
        CGAL_assertion_msg(false, "Intersections in CDT2 are not allowed");
        return PCDT(projection_traits);
      }
    }

    return pcdt;
  }

  /**
    * Triangulate the facet 'f' and returns vertices and all newly created triangle facets
    */
  static std::pair<std::list<VertexSPtr>, std::list<FacetSPtr> >
  triangulate_facet(const FacetSPtr& facet,
                    const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_TRANSF_TRACE("Triangulate facet " << facet->get_ID() << " of polyhedron "
                          << polyhedron->get_ID() << " with " << facet->vertices().size()
                          << " vertices");
    CGAL_SS3_DEBUG_SPTR(facet);
    CGAL_SS3_DEBUG_SPTR(polyhedron);

    std::list<VertexSPtr> facet_vertices = facet->vertices();

    CGAL_precondition(facet && polyhedron && facet->vertices().size() >= 3);

    if (facet_vertices.size() == 3) {
      return { facet_vertices, { facet } };
    }

    facet->sort_vertices();

    using CDT2_Tag = CGAL::No_constraint_intersection_tag;
    auto pcdt = triangulate_facet<CDT2_Tag>(facet);

    using PCDT = decltype(pcdt);
    using PCDT_FH = typename PCDT::Face_handle;

    std::unordered_map<PCDT_FH, bool> in_domain_map;
    boost::associative_property_map<std::unordered_map<PCDT_FH, bool>> in_domain(in_domain_map);
    CGAL::mark_domain_in_triangulation(pcdt, in_domain);

    // Speed of the subdivided facet is applied to all the subfacets
    const FT& parent_speed = Hds_utils::get_speed(facet);

    polyhedron->remove_facet(facet);

    // Create new facets and edges for each triangle
    std::list<FacetSPtr> created_facets;
    for (auto fh : pcdt.finite_face_handles()) {
      if (!get(in_domain, fh)) {
        continue;
      }

      VertexSPtr v0 = fh->vertex(0)->info();
      VertexSPtr v1 = fh->vertex(1)->info();
      VertexSPtr v2 = fh->vertex(2)->info();
      std::vector<VertexSPtr> verts = {v0, v1, v2};
      FacetSPtr new_facet = Facet::create(verts);
      Plane3SPtr plane = Kernel_factory::createPlane3(v0->point(),
                                                      v1->point(),
                                                      v2->point());
      new_facet->set_plane(plane);
      normalize_plane_coefficients(new_facet);
      SkelFacetDataSPtr new_data = Skeleton_facet_data::create(new_facet);
      new_data->set_speed(parent_speed);
      polyhedron->add_facet(new_facet);
      created_facets.push_back(new_facet);
    }

    return { facet_vertices , created_facets };
  }

  /**
    * checks if the plane is normalized
    */
  static bool has_normalized_plane(const FacetSPtr& facet)
  {
    CGAL_SS3_DEBUG_SPTR(facet);
    const FT& a = facet->plane()->a();
    const FT& b = facet->plane()->b();
    const FT& c = facet->plane()->c();
    return (a*a + b*b + c*c - 1) <= 1e-5;
  }

  /**
    * normalizes the plane coefficients to obtain a canonical plane representation
    */
  static bool normalize_plane_coefficients(const FacetSPtr& facet)
  {
    CGAL_SS3_DEBUG_SPTR(facet);
    const FT& a = facet->plane()->a();
    const FT& b = facet->plane()->b();
    const FT& c = facet->plane()->c();
    const FT& d = facet->plane()->d();

    // this should be the only place with unavoidable SQRTs
    const FT n = CGAL::approximate_sqrt(CGAL::square(a) + CGAL::square(b) + CGAL::square(c));

    if (!is_zero(n)) {
      facet->set_plane(Kernel_factory::createPlane3(a/n, b/n, c/n, d/n)); // @todo to_double() it here too?
      return true;
    } else {
      facet->set_plane(Kernel_factory::createPlane3(a, b, c, d));
      return false;
    }
  }

  /**
    * normalizes facet planes
    */
  static bool normalize_facet_planes(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    bool result = true;
    for (const FacetSPtr& facet : polyhedron->facets()) {
      result &= normalize_plane_coefficients(facet);
    }
    return result;
  }
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_POLYHEDRON_TRANSFORMATION_H */
