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
class PolyhedronTransformation
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

  using Vertex = typename Polyhedron::template Vertex<Traits>;
  using VertexSPtr = typename Polyhedron::VertexSPtr;
  using Edge = typename Polyhedron::template Edge<Traits>;
  using EdgeWPtr = typename Polyhedron::EdgeWPtr;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;
  using Facet = typename Polyhedron::template Facet<Traits>;
  using FacetWPtr = typename Polyhedron::FacetWPtr;
  using FacetSPtr = typename Polyhedron::FacetSPtr;

  using SkelVertexData = typename Polyhedron::SkelVertexData;
  using SkelVertexDataSPtr = typename Polyhedron::SkelVertexDataSPtr;
  using SkelEdgeData = typename Polyhedron::SkelEdgeData;
  using SkelEdgeDataSPtr = typename Polyhedron::SkelEdgeDataSPtr;
  using SkelFacetData = typename Polyhedron::SkelFacetData;
  using SkelFacetDataSPtr = typename Polyhedron::SkelFacetDataSPtr;

private:
  using KernelFactory = kernel::KernelFactory<Traits>;
  using KernelWrapper = kernel::KernelWrapper<Traits>;
  using GeomUtils = algorithm::GeomUtils<Traits>;

public:
  static Point3SPtr boundingBoxMin(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    FT p_min[3];
    for (unsigned int i = 0; i < 3; ++i) {
      p_min[i] = std::numeric_limits<double>::max();
    }
    for (const VertexSPtr& vertex : polyhedron->vertices()) {
      Point3SPtr p = vertex->getPoint();
      for (unsigned int i = 0; i < 3; ++i) {
        if ((*p)[i] < p_min[i]) {
          p_min[i] = (*p)[i];
        }
      }
    }
    Point3SPtr result = KernelFactory::createPoint3(p_min[0], p_min[1], p_min[2]);
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  static Point3SPtr boundingBoxMax(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    FT p_max[3];
    for (unsigned int i = 0; i < 3; ++i) {
      p_max[i] = -std::numeric_limits<double>::max();
    }
    for (const VertexSPtr& vertex : polyhedron->vertices()) {
      Point3SPtr p = vertex->getPoint();
      for (unsigned int i = 0; i < 3; ++i) {
        if ((*p)[i] > p_max[i]) {
          p_max[i] = (*p)[i];
        }
      }
    }
    Point3SPtr result = KernelFactory::createPoint3(p_max[0], p_max[1], p_max[2]);
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  static bool isInsideBox(const PolyhedronSPtr& polyhedron,
                          const Point3SPtr& p_box_min, const Point3SPtr& p_box_max)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    CGAL_SS3_DEBUG_SPTR(p_box_min);
    CGAL_SS3_DEBUG_SPTR(p_box_max);
    bool result = true;
    for (const VertexSPtr& vertex : polyhedron->vertices()) {
      Point3SPtr p = vertex->getPoint();
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
      Point3SPtr p = vertex->getPoint();
      Point3SPtr p_t = KernelFactory::createPoint3(*p + *v_t);
      vertex->setPoint(p_t);
    }

    polyhedron->initPlanes();
  }

  static void scale(const PolyhedronSPtr& polyhedron,
                    const Vector3SPtr& v_s)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    CGAL_SS3_DEBUG_SPTR(v_s);

    for (const VertexSPtr& vertex : polyhedron->vertices()) {
      Point3SPtr p = vertex->getPoint();
      Point3SPtr p_s = KernelFactory::createPoint3((*p)[0] * (*v_s)[0],
                                                   (*p)[1] * (*v_s)[1],
                                                   (*p)[2] * (*v_s)[2]);
      vertex->setPoint(p_s);
    }

    polyhedron->initPlanes();
  }

  static void translateNscale(const PolyhedronSPtr& polyhedron,
                              const Point3SPtr& p_box_min, const Point3SPtr& p_box_max)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    CGAL_SS3_DEBUG_SPTR(p_box_min);
    CGAL_SS3_DEBUG_SPTR(p_box_max);

    Vector3SPtr v_box_min = KernelFactory::createVector3(p_box_min);
    Vector3SPtr v_box_max = KernelFactory::createVector3(p_box_max);
    Vector3SPtr v_size = KernelFactory::createVector3(*v_box_max - *v_box_min);
    Vector3SPtr v_center = KernelFactory::createVector3((*v_box_min + *v_box_max) / 2.0);

    Point3SPtr p_box_min_curr = boundingBoxMin(polyhedron);
    Point3SPtr p_box_max_curr = boundingBoxMax(polyhedron);
    Vector3SPtr v_box_min_curr = KernelFactory::createVector3(p_box_min_curr);
    Vector3SPtr v_box_max_curr = KernelFactory::createVector3(p_box_max_curr);
    Vector3SPtr v_size_curr = KernelFactory::createVector3(*v_box_max_curr - *v_box_min_curr);
    Vector3SPtr v_center_curr = KernelFactory::createVector3((*v_box_min_curr + *v_box_max_curr) / 2.0);

    FT scale_factor = std::numeric_limits<double>::max(); // do not put FT
    for (unsigned int i = 0; i < 3; ++i) {
      FT s = (*v_size)[i]/(*v_size_curr)[i];
      if (scale_factor > s) {
        scale_factor = s;
      }
    }
    scale_factor = floor(CGAL::to_double(scale_factor)*1000.0)/1000.0;
    Vector3SPtr v_s = KernelFactory::createVector3(scale_factor, scale_factor, scale_factor);

    Vector3SPtr v_t = KernelFactory::createVector3((*v_center_curr) * -1.0);
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

  static void truncatePrecision(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);

    ConfigurationSPtr config = Configuration::getInstance();
    double range = 1e-10;
    if (config->isLoaded()) {
      range = config->getDouble("Preprocessing", "truncate_precision");
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
      Point3SPtr p = vertex->getPoint();
      CGAL_SS3_TRANSF_TRACE_V(32,"    Truncated from: " << *(vertex->getPoint()));

      double rx = CGAL::Polygon_mesh_processing::autorefine_impl::double_ceil(p->x() * scale) / scale;
      double ry = CGAL::Polygon_mesh_processing::autorefine_impl::double_ceil(p->y() * scale) / scale;
      double rz = CGAL::Polygon_mesh_processing::autorefine_impl::double_ceil(p->z() * scale) / scale;

      vertex->setPoint(KernelFactory::createPoint3(rx, ry, rz));
      CGAL_SS3_TRANSF_TRACE_V(32,"    Truncated to: " << *(vertex->getPoint()));
    }
  }

  static bool hasCoplanarFacets(const EdgeSPtr& edge,
                                const double epsilon)
  {
    CGAL_SS3_DEBUG_SPTR(edge);

    bool result = false;
    FacetSPtr facet_l = edge->getFacetL();
    FacetSPtr facet_r = edge->getFacetR();
    if (facet_l && facet_r) {
      if (epsilon == 0.) {
        return (*(facet_l->plane()) == *(facet_r->plane())); // planes are normalized
      }

      Vector3SPtr normal_l = KernelFactory::createVector3(facet_l->plane());
      Vector3SPtr normal_r = KernelFactory::createVector3(facet_r->plane());
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

  static void mergeFacets(const EdgeSPtr& edge,
                          const FacetSPtr& facet_into,
                          const FacetSPtr& facet_from,
                          const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_DEBUG_SPTR(edge);
    CGAL_SS3_DEBUG_SPTR(facet_into);
    CGAL_SS3_DEBUG_SPTR(facet_from);
    CGAL_precondition(facet_into != facet_from);
    CGAL_precondition(edge->getFacetL() == facet_from || edge->getFacetR() == facet_from);
    CGAL_precondition(edge->getFacetL() == facet_into || edge->getFacetR() == facet_into);

    CGAL_SS3_TRANSF_TRACE_V(16, "Merging F" << facet_from->getID() << " into F" << facet_into->getID() <<
                                " Common edge E" << edge->getID() << " [V" << edge->getVertexSrc()->getID()
                                                                  << " - V" << edge->getVertexDst()->getID() << "]");
    CGAL_SS3_TRANSF_TRACE_V(16, "  FROM normal: " << *(KernelFactory::createVector3(facet_from->getPlane())));
    CGAL_SS3_TRANSF_TRACE_V(16, "  INTO normal: " << *(KernelFactory::createVector3(facet_into->getPlane())));

    // remove the facet incidence info from the edge such that the facets are not deleted
    // when Polyhedron::removeEdge() is called
    facet_into->removeEdge(edge);
    facet_from->removeEdge(edge);

    polyhedron->removeEdge(edge);
    if (facet_into && facet_from && facet_into != facet_from) {
      facet_into->merge(facet_from);
      polyhedron->removeFacet(facet_from);
    }

    CGAL_postcondition(polyhedron->isConsistent());
  }

  static void mergeFacets(const EdgeSPtr& edge,
                          const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_DEBUG_SPTR(edge);
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    return mergeFacets(edge, edge->getFacetL(), edge->getFacetR(), polyhedron);
  }

  static int mergeCoplanarFacets(const PolyhedronSPtr& polyhedron,
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
      if (hasCoplanarFacets(edge, epsilon)) {
        edges_toremove.push_back(edge);
      }
    }

    CGAL_SS3_TRANSF_TRACE(edges_toremove.size() << " edges to merge");
    CGAL_SS3_TRANSF_TRACE_CODE(if (edges_toremove.size() > 0))
    CGAL_SS3_TRANSF_TRACE("Adjacent facets of the following edges are detected to be coplanar and will be merged.");

    for (EdgeWPtr edge_w : edges_toremove) {
      if (EdgeSPtr edge = edge_w.lock()) {
        mergeFacets(edge, polyhedron);
        ++result;
      }
    }

    sanitize(polyhedron);

    polyhedron->initializeAllIDs();

    CGAL_SS3_TRANSF_TRACE_V(4, "  final facet count: " << polyhedron->facets().size());

#ifdef CGAL_SS3_DUMP_FILES
    IO::OBJFile::save("results/coplanar_merge_after.obj", polyhedron,
                      false /*do not triangulate*/);
#endif

    return result;
  }

  static int mergeCoplanarFacets(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    double epsilon = 0.0;
    ConfigurationSPtr config = Configuration::getInstance();
    if (config->isLoaded()) {
      std::string section("Preprocessing");
      std::string key("coplanarity_epsilon");
      if (config->contains(section, key)) {
        epsilon = config->getDouble(section, key);
      }
    }

    return mergeCoplanarFacets(polyhedron, epsilon);
  }

  static int removeVerticesDegLt3(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_TRANSF_TRACE("Remove Vertices with degree < 3");
    CGAL_SS3_TRANSF_TRACE("  initial vertex count: " << polyhedron->vertices().size());

    CGAL_SS3_DEBUG_SPTR(polyhedron);

    int result = 0;
    std::list<VertexSPtr> vertices_toremove;
    typename std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
      VertexSPtr vertex = *it_v++;
      if (vertex->degree() < 3) {
        vertices_toremove.push_back(vertex);
        CGAL_SS3_TRANSF_TRACE("Enlist: V" << vertex->getID());
        for (FacetWPtr wf : vertex->facets()) {
          FacetSPtr facet = wf.lock();
          CGAL_SS3_TRANSF_TRACE("  Incident facet with: " << facet->vertices().size() << " vertices");
        }
      }
    }
    it_v = vertices_toremove.begin();
    while (it_v != vertices_toremove.end()) {
      VertexSPtr vertex = *it_v++;
      CGAL_SS3_TRANSF_TRACE("Removing " << vertex->toString());

      typename std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
      while (it_f != vertex->facets().end()) {
        FacetWPtr facet_wptr = *it_f++;
        if (FacetSPtr facet = facet_wptr.lock()) {
          facet->removeVertex(vertex);
        }
      }

      if (vertex->degree() == 2) {
        EdgeSPtr edge_src = vertex->firstEdge();
        EdgeSPtr edge_dst = edge_src->next(vertex);

        VertexSPtr vertex_src = edge_src->getVertexSrc();
        if (vertex_src == vertex) {
          vertex_src = edge_src->getVertexDst();
        }
        VertexSPtr vertex_dst = edge_dst->getVertexSrc();
        if (vertex_dst == vertex) {
          vertex_dst = edge_dst->getVertexDst();
        }

        FacetSPtr fL = edge_src->getFacetL();
        FacetSPtr fR = edge_src->getFacetR();
        CGAL_assertion(fL != fR);

        // if there is any facet with degree 3, put it in fL (both could be with degree 3,
        // but then the swap does not change anything)
        //
        // It's "== 2" because we have already removed the vertex from its incident facets
        CGAL_assertion(!fL->hasVertex(vertex));
        CGAL_assertion(!fR->hasVertex(vertex));
        if (fR->vertices().size() == 2) {
          std::swap(fL, fR);
        }

        if (fL->vertices().size() == 2) {
          if (fR->vertices().size() == 2) {
            CGAL_SS3_TRANSF_TRACE("Vertex is the apex of two facets of degree 3");
            // both facets have degree 3, so remove everything (both facets, both edges,
            // and one of the other edges + setting up incident facets properly)
            VertexSPtr other_vertex = edge_src->other(vertex);
            EdgeSPtr third_edge_1 = edge_src->next(other_vertex);
            EdgeSPtr third_edge_2 = edge_src->prev(other_vertex);
            CGAL_assertion(third_edge_1 != third_edge_2);
            mergeFacets(third_edge_1, polyhedron);
            mergeFacets(third_edge_2, polyhedron);

            edge_dst->getFacetL()->removeEdge(edge_dst);
            edge_dst->getFacetR()->removeEdge(edge_dst);
            polyhedron->removeEdge(edge_dst);

            if (edge_src->getVertexDst() == vertex) {
              edge_src->replaceVertexDst(vertex_dst);
            } else if (edge_src->getVertexSrc() == vertex) {
              edge_src->replaceVertexSrc(vertex_dst);
            }
          } else {
            CGAL_SS3_TRANSF_TRACE("Vertex is the apex of one facet of degree 3");
            // one facet has degree 3, so remove the vertex, the two incident edges, and the facet.
            EdgeSPtr third_edge;
            for (const EdgeSPtr& edge : fR->edges()) {
              if (edge != edge_src && edge != edge_dst) {
                third_edge = edge;
                break;
              }
            }
            CGAL_assertion(third_edge != nullptr);

            fR->removeVertex(vertex);
            fR->removeEdge(edge_src);
            fR->removeEdge(edge_dst);
            fL->removeEdge(third_edge);
            if (third_edge->getFacetL() == fL) {
              third_edge->setFacetL(fR);
            } else {
              third_edge->setFacetR(fR);
            }
          }
        } else {
          edge_dst->getFacetL()->removeEdge(edge_dst);
          edge_dst->getFacetR()->removeEdge(edge_dst);
          polyhedron->removeEdge(edge_dst);

          if (edge_src->getVertexDst() == vertex) {
            edge_src->replaceVertexDst(vertex_dst);
          } else if (edge_src->getVertexSrc() == vertex) {
            edge_src->replaceVertexSrc(vertex_dst);
          }
        }
      }

      polyhedron->removeVertex(vertex);
      CGAL_postcondition(polyhedron->isConsistent());
      ++result;
    }

    CGAL_SS3_TRANSF_TRACE("  final vertex count: " << polyhedron->vertices().size());

    return result;
  }

  static int removeFacetsDegLt3(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);
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
        if (edge->getFacetL() == facet && edge->getFacetR() && edge->getFacetR() != facet) {
          neighbor = edge->getFacetR();
        } else if (edge->getFacetR() == facet && edge->getFacetL() && edge->getFacetL() != facet) {
          neighbor = edge->getFacetL();
        }
        if (neighbor && neighbor->vertices().size() > best_size) {
          best_edge = edge;
          best_neighbor = neighbor;
          best_size = neighbor->vertices().size();
        }
      }

      if (best_neighbor) {
        mergeFacets(best_edge, best_neighbor, facet, polyhedron);
      } else {
        // No neighbor, just remove
        for (const EdgeSPtr& edge : facet->edges()) {
          polyhedron->removeEdge(edge);
        }
      }
      ++result;
    }
    return result;
  }

  static int sanitize(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_TRANSF_TRACE("Sanitizing polyhedron...");

    CGAL_SS3_DEBUG_SPTR(polyhedron);

    // - removeVerticesDegLt3 can create facets with fewer than 3 vertices
    // - removeFacetsDegLt3 removes facets with fewer than 3 vertices
    // so loop till nothing is done anymore
    int result = 0;
    for (;;) {
      std::size_t vlt3 = removeVerticesDegLt3(polyhedron);
      std::size_t flt3 = removeFacetsDegLt3(polyhedron);
      int partial = vlt3 + flt3;
      if (partial == 0) {
        break;
      }
      result += vlt3 + flt3;
    }
    return result;
  }

  /**
    * updates the position of the vertex of a polyhedron, computed from the planes of
    * its incident faces.
    */
  static bool resetPoint(const VertexSPtr& vertex,
                         const std::array<Plane3SPtr, 3>& planes)
  {
    CGAL_SS3_DEBUG_SPTR(vertex);
    CGAL_SS3_DEBUG_SPTR(planes[0]);
    CGAL_SS3_DEBUG_SPTR(planes[1]);
    CGAL_SS3_DEBUG_SPTR(planes[2]);

    Point3SPtr point = KernelWrapper::intersection(planes[0], planes[1], planes[2]);
    if (!point) {
      CGAL_SS3_TRANSF_TRACE_V(1, "Warning: triplet of planes does not define a point!");
      Point3SPtr result = Point3SPtr();
      CGAL_SS3_DEBUG_SPTR(result);
      return false;
    }

    vertex->setPoint(point);
    CGAL_SS3_TRANSF_TRACE_V(16, "  New point = " << *point);

    CGAL_postcondition_code(for (FacetWPtr facet_wptr : vertex->facets()) {)
    CGAL_postcondition_code(    if (FacetSPtr facet = facet_wptr.lock()) {)
    CGAL_postcondition(             facet->getPlane()->has_on(*point));
    CGAL_postcondition_code(    })
    CGAL_postcondition_code(})

    return true;
  }

  // resets using the first 3 planes, even if there are more
  static bool resetPoint(const VertexSPtr& vertex)
  {
    CGAL_SS3_TRANSF_TRACE_V(16, "resetPoint() of " << vertex->toString());
    CGAL_SS3_DEBUG_SPTR(vertex);

    std::array<Plane3SPtr, 3> planes;
    unsigned int i = 0;
    typename std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
    while (i < 3 && it_f != vertex->facets().end()) {
      FacetWPtr facet_wptr = *it_f++;
      if (FacetSPtr facet = facet_wptr.lock()) {
        planes[i++] = facet->getPlane();
        CGAL_SS3_TRANSF_TRACE_V(16, "  Facet " << facet->getID() << " [" << *(facet->getPlane()) << "]");
      }
    }
    CGAL_postcondition(i == 3);

    return resetPoint(vertex, { planes[0], planes[1], planes[2] });
  }

  /**
    * updates the positions of the vertices of a polyhedron, computed from the planes of
    * their incident faces.
    */
  static bool resetPoints(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_TRANSF_TRACE("Reset point positions");
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    for (VertexSPtr vertex : polyhedron->vertices()) {
      if (!resetPoint(vertex)) {
        CGAL_SS3_TRANSF_TRACE_V(1, "Warning: failed to reset vertex " << vertex->toString());
        return false;
      }
    }
    return true;
  }

  /**
    * returns the shifted position of the vertex of a polyhedron
    * \pre vertex->degree() == 3
    */
  static Point3SPtr shiftPoint(const VertexSPtr& vertex,
                               const FT& offset)
  {
    CGAL_SS3_TRANSF_TRACE_V(32, "shift " << vertex->toString() << "\noffset = " << offset);
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

        planes[i] = shiftPlane(facet, offset);

        CGAL_SS3_TRANSF_TRACE_V(64, "facet[" << i << "] = " << facet->getID());
        CGAL_SS3_TRANSF_TRACE_V(64, "Offset Plane[" << i << "] = " << *(planes[i]));

        ++i;
      }
    }

    if (i < 3) {
      CGAL_SS3_TRANSF_TRACE_V(1, "Warning: Could not find three independent planes for high-degree vertex" << vertex->toString());
      return { };
    }

    Point3SPtr point = KernelWrapper::intersection(planes[0], planes[1], planes[2]);
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
  static Segment3SPtr shiftEdge(const EdgeSPtr& edge,
                                const FT& offset)
  {
    CGAL_SS3_DEBUG_SPTR(edge);

    FacetSPtr facet_l = edge->getFacetL();
    FacetSPtr facet_r = edge->getFacetR();
    FacetSPtr facet_src = edge->getFacetL()->next(edge->getVertexSrc());
    FacetSPtr facet_dst = edge->getFacetR()->next(edge->getVertexDst());

    FT speed_l = 1.0;
    if (facet_l->hasData()) {
      speed_l = std::dynamic_pointer_cast<SkelFacetData>(facet_l->getData())->getSpeed();
    }

    FT speed_r = 1.0;
    if (facet_r->hasData()) {
      speed_r = std::dynamic_pointer_cast<SkelFacetData>(facet_r->getData())->getSpeed();
    }

    FT speed_src = 1.0;
    if (facet_src->hasData()) {
      speed_src = std::dynamic_pointer_cast<SkelFacetData>(facet_src->getData())->getSpeed();
    }

    FT speed_dst = 1.0;
    if (facet_dst->hasData()) {
      speed_dst = std::dynamic_pointer_cast<SkelFacetData>(facet_dst->getData())->getSpeed();
    }

    // Offset the two common planes
    Plane3SPtr offset_plane_l = GeomUtils::offsetPlane(facet_l->plane(), offset*speed_l);
    Plane3SPtr offset_plane_r = GeomUtils::offsetPlane(facet_r->plane(), offset*speed_r);
    Plane3SPtr offset_plane_src = GeomUtils::offsetPlane(facet_src->plane(), offset*speed_src);
    Plane3SPtr offset_plane_dst = GeomUtils::offsetPlane(facet_dst->plane(), offset*speed_dst);

#if 0
    // leaving it here because it's not that intuitive: factoring the intersection of the
    // two common planes is much slower than computing two 3-plane intersections
    Line3SPtr common_line = KernelWrapper::intersection(offset_plane_l, offset_plane_r);
    Point3SPtr src_point = KernelWrapper::intersection(offset_plane_src, common_line);
    Point3SPtr dst_point = KernelWrapper::intersection(offset_plane_dst, common_line);
#else
    Point3SPtr src_point = KernelWrapper::intersection(offset_plane_src, offset_plane_l, offset_plane_r);
    Point3SPtr dst_point = KernelWrapper::intersection(offset_plane_dst, offset_plane_l, offset_plane_r);
#endif

    CGAL_assertion(bool(src_point) && bool(dst_point));
    return KernelFactory::createSegment3(src_point, dst_point);
  }

  /**
    * returns the shifted position of the facet of a polyhedron
    */
  static Plane3SPtr shiftPlane(const FacetSPtr& facet,
                               const FT& offset)
  {
    CGAL_SS3_DEBUG_SPTR(facet);

    FT facet_speed = 1.0;
    if (facet->hasData()) {
      facet_speed = std::dynamic_pointer_cast<SkelFacetData>(facet->getData())->getSpeed();
    }

    return GeomUtils::offsetPlane(facet->plane(), facet_speed*offset);
  }

  /**
    * Offsets the polyhedron `polyhedron`
    * Negative offset points to the interior of the polyhedron.
    * This function is for the main shift in the event loop.
    */
  static void shiftFacets(const PolyhedronSPtr& polyhedron,
                          const FT& offset,
                          const bool recompute_positions = true)
  {
    CGAL_SS3_TRANSF_TRACE("~~~~ Shift polyhedron by " << offset << " [in place]");
    CGAL_SS3_DEBUG_SPTR(polyhedron);

    typename std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
      FacetSPtr facet = *it_f++;
      Plane3SPtr offset_plane = shiftPlane(facet, offset);
      facet->setPlane(offset_plane);
    }

    typename std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
      VertexSPtr vertex = *it_v++;

      // See comment at the top of the function
      CGAL_assertion(vertex->degree() == 3);

      Point3SPtr old_point = vertex->getPoint(), new_point;
      if (offset != 0 || recompute_positions) {
        CGAL_assertion_code(bool ok =)
        resetPoint(vertex);
        CGAL_assertion(ok);
      }
    }
  }

  /**
  * Offsets the polyhedron `polyhedron`, which may have degree 1 vertices.
  * Negative offset points to the interior of the polyhedron.
  */
  static bool shiftFacetsDegree1(const PolyhedronSPtr& polyhedron,
                                 const FT& offset)
  {
    CGAL_SS3_TRANSF_TRACE("~~~~ Shift polyhedron by " << offset << " [in place w/ degree 1]");
    CGAL_SS3_DEBUG_SPTR(polyhedron);

    // Maps to store shifted planes and points
    std::unordered_map<FacetSPtr, Plane3SPtr> facet_to_shifted_plane;
    std::unordered_map<VertexSPtr, Point3SPtr> vertex_to_shifted_point;

    for (const FacetSPtr& facet : polyhedron->facets()) {
      Plane3SPtr offset_plane = shiftPlane(facet, offset);
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

        Point3SPtr shifted_point = KernelWrapper::intersection(shifted_planes[0], shifted_planes[1], shifted_planes[2]);
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
        unsigned int i = 0;
        for (EdgeWPtr edge_wptr : vertex->edges()) {
          if ((edge = edge_wptr.lock())) {
            ++i;
          }
        }
        CGAL_assertion(i == 1);

        VertexSPtr vertex_other = edge->other(vertex);
        CGAL_assertion(vertex_other->degree() >= 3);

        FacetSPtr facet1 = edge->getFacetL();
        FacetSPtr facet2 = edge->getFacetR();
        Plane3SPtr plane1 = facet_to_shifted_plane[facet1];
        Plane3SPtr plane2 = facet_to_shifted_plane[facet2];
        Line3SPtr intersection_line = KernelWrapper::intersection(plane1, plane2);
        CGAL_assertion(bool(intersection_line));

        // Determine the direction using the unshifted positions
        Point3SPtr point_other = vertex_to_shifted_point[vertex_other];
        Point3SPtr point = vertex->getPoint();
        CGAL_assertion(*point != *point_other);
        Vector_3 direction = intersection_line->to_vector();
        if (CGAL::scalar_product(direction, *point - *vertex_other->getPoint()) < 0) {
          direction = -direction;
        }

        // Apply the shift to the vertex
        Point3SPtr offset_point = KernelFactory::createPoint3(*point_other + direction);
        if (!offset_point) {
          return false;
        }
        vertex_to_shifted_point[vertex] = offset_point;
      }
    }

    // Apply all shifts
    for (const FacetSPtr& facet : polyhedron->facets()) {
      facet->setPlane(facet_to_shifted_plane[facet]);
    }
    for (const VertexSPtr& vertex : polyhedron->vertices()) {
      if (vertex_to_shifted_point.count(vertex)) {
        vertex->setPoint(vertex_to_shifted_point[vertex]);
      }
    }

    return true;
  }

    // The tag is a template parameter because for debugging outputs,
  // we might want to tolerate intersections
  template <typename CDT2Tag>
  static auto constructFacetTriangulation(const FacetSPtr& facet)
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

    Vector3SPtr n = KernelFactory::createVector3(facet->plane());
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
        PCDT_VH vh = pcdt.insert(*(vertex->getPoint()));
        res.first->second = vh;
        vh->info() = vertex;
      }
    }

    typename std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
    while (it_e != facet->edges().end()) {
      EdgeSPtr edge = *it_e++;
      VertexSPtr v0 = edge->src(facet);
      VertexSPtr v1 = edge->dst(facet);
      CGAL_assertion(*(v0->getPoint()) != *(v1->getPoint()));

      PCDT_VH vh0 = face_vhs.at(v0);
      PCDT_VH vh1 = face_vhs.at(v1);

      try {
        pcdt.insert_constraint(vh0, vh1);
      } catch(const typename PCDT::Intersection_of_constraints_exception&) {
        CGAL_SS3_TRANSF_TRACE("Error: Intersection of constraint w/ " << vh0->point() << " " << vh1->point());
        CGAL_SS3_TRANSF_TRACE(facet->toString());
        CGAL_assertion_msg(false, "Intersections in CDT2 are not allowed");
        return PCDT(projection_traits);
      }
    }

    return pcdt;
  }

  /**
    * Triangulate the facet 'f' and returns vertices and all newly created triangle facets
    */
  static std::pair<std::list<VertexSPtr>,
                   std::list<FacetSPtr> > triangulateFacet(const FacetSPtr& facet,
                                                           const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_TRANSF_TRACE("Triangulate facet " << facet->getID() << " of polyhedron "
                          << polyhedron->getID() << " with " << facet->vertices().size()
                          << " vertices");
    CGAL_SS3_DEBUG_SPTR(facet);
    CGAL_SS3_DEBUG_SPTR(polyhedron);

    std::list<VertexSPtr> facet_vertices = facet->vertices();

    CGAL_precondition(facet && polyhedron && facet->vertices().size() >= 3);

    if (facet_vertices.size() == 3) {
      return { facet_vertices, { facet } };
    }

    facet->sortVertices();

    using CDT2_Tag = CGAL::No_constraint_intersection_tag;
    auto pcdt = constructFacetTriangulation<CDT2_Tag>(facet);

    using PCDT = decltype(pcdt);
    using PCDT_FH = typename PCDT::Face_handle;

    std::unordered_map<PCDT_FH, bool> in_domain_map;
    boost::associative_property_map<std::unordered_map<PCDT_FH, bool>> in_domain(in_domain_map);
    CGAL::mark_domain_in_triangulation(pcdt, in_domain);

    // Speed of the subdivided facet is applied to all the subfacets
    FT parent_speed = 1.0;
    if (facet->hasData()) {
      parent_speed = std::dynamic_pointer_cast<SkelFacetData>(facet->getData())->getSpeed();
    }

    polyhedron->removeFacet(facet);

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
      Plane3SPtr plane = KernelFactory::createPlane3(v0->getPoint(),
                                                      v1->getPoint(),
                                                      v2->getPoint());
      new_facet->setPlane(plane);
      normalizePlaneCoefficients(new_facet);
      SkelFacetDataSPtr new_data = SkelFacetData::create(new_facet);
      new_data->setSpeed(parent_speed);
      polyhedron->addFacet(new_facet);
      created_facets.push_back(new_facet);
    }

    return { facet_vertices , created_facets };
  }

  /**
  * Checks if the plane is normalized
  */
  static bool hasNormalizedPlane(const FacetSPtr& facet)
  {
    CGAL_SS3_DEBUG_SPTR(facet);
    const FT& a = facet->plane()->a();
    const FT& b = facet->plane()->b();
    const FT& c = facet->plane()->c();
    return (a*a + b*b + c*c - 1) <= 1e-5;
  }

  /**
  * Normalizes the plane coefficients to obtain a canonical plane representation
  */
  static bool normalizePlaneCoefficients(const FacetSPtr& facet)
  {
    CGAL_SS3_DEBUG_SPTR(facet);
    const FT& a = facet->plane()->a();
    const FT& b = facet->plane()->b();
    const FT& c = facet->plane()->c();
    const FT& d = facet->plane()->d();

    // this should be the only place with unavoidable SQRTs
    const FT n = CGAL::approximate_sqrt(CGAL::square(a) + CGAL::square(b) + CGAL::square(c));

    if (!is_zero(n)) {
      facet->setPlane(KernelFactory::createPlane3(a/n, b/n, c/n, d/n)); // @todo to_double() it here too?
      return true;
    } else {
      facet->setPlane(KernelFactory::createPlane3(a, b, c, d));
      return false;
    }
  }

  /**
    * Normalize facet planes
  */
  static bool normalizeFacetPlanes(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    bool result = true;
    for (const FacetSPtr& facet : polyhedron->facets()) {
      result &= normalizePlaneCoefficients(facet);
    }
    return result;
  }
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_POLYHEDRON_TRANSFORMATION_H */
